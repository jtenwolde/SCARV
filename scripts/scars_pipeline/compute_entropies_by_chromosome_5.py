from scars import *

import pybedtools
import pyranges as pr
from keras.models import load_model
import keras_genomics
import joblib
import numpy as np
import pandas as pd
import sys


chrom = sys.argv[1]
chrom_base = "chrX" if "X" in chrom else chrom

window_size = int(sys.argv[2])
ancestry = sys.argv[3]
pop_size = int(sys.argv[4])
outFile = sys.argv[5]

k = 13
genome = pybedtools.genome_registry.hg38
reference_fasta = "/rds/project/who1000-1/rds-who1000-cbrc/ref/UCSC/hg38/hg38.fa"

snvs = scars_queries.load_variants("/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/scars_pipeline_gnomad_hg38/" + ancestry +"/variants/snvs_by_chr/pass_snvs_" + chrom_base + ".bed")
deletions = scars_queries.load_variants("/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/scars_pipeline_gnomad_hg38/" + ancestry +"/variants/pass_deletions.bed")
reliable_sites = pr.read_bed("/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/scars_pipeline_gnomad_hg38/" + ancestry + "/quality_filtering/reliable_sites_by_chr/reliable_sites_" + chrom + ".bed")

snvs_df = snvs.as_df()
snvs_df = snvs_df.set_index('Start')

preds_by_kmer_file = "/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/scars_pipeline_gnomad_hg38/" + ancestry + "/preds_by_kmer.pkl"
preds_by_kmer = pd.read_pickle(preds_by_kmer_file)

chrom_gr = pr.from_dict({'Chromosome': [chrom_base], 'Start': [0], 'End': [genome[chrom_base][1]]})
chrom_seq = pr.get_fasta(chrom_gr, reference_fasta)[0].upper()
chrom_seq_ohe = np.array(pd.get_dummies(pd.Series(list(chrom_seq)))[['A','C','G','T']])


# Compute chromosome-wide predictions
k_mers = pd.Series([chrom_seq[(pos-k//2):(pos+k//2+1)] for pos in range(genome[chrom_base][1])])

preds = np.empty(shape=(genome[chrom_base][1], 4))
preds[:] = np.nan

reliable_indices = reliable_sites.tile(1).Start.to_numpy()
identified_nucleotides = np.where(["N" not in seq for seq in k_mers])[0]
query_indices = np.intersect1d(reliable_indices, identified_nucleotides)

preds[query_indices] = preds_by_kmer.loc[k_mers[query_indices]].to_numpy()


# Compute chromosome-wide allele frequencies
allele_counts = 2 * pop_size * chrom_seq_ohe.astype('int32')

deletion_allele_counts = deletions[chrom_base].tile(1).as_df().groupby('Start')['ac'].sum()
allele_counts[deletion_allele_counts.index] -=  deletion_allele_counts.to_numpy()[:, np.newaxis] * chrom_seq_ohe[deletion_allele_counts.index] 

alternate_alleles = pd.get_dummies(snvs_df.alt).mul(snvs_df.ac, axis=0)\
                                               .groupby('Start')\
                                               .sum()\
                                               .to_numpy()

reference_alleles = pd.get_dummies(snvs_df.ref).mul(snvs_df.an, axis=0)\
                                               .groupby('Start')\
                                               .first()\
                                               .to_numpy()

reference_alleles -= np.sum(alternate_alleles, axis=1)[:, np.newaxis]
reference_alleles[reference_alleles < 0] = 0

allele_counts[np.unique(snvs.Start)] = reference_alleles + alternate_alleles
allele_counts[np.setdiff1d(np.arange(allele_counts.shape[0]), query_indices)] = 0   # remove contribution from unreliable or unidentified sequences
allele_frequencies = allele_counts / allele_counts.sum(axis=1)[:, np.newaxis]


# Generate score
observed_entropy = scars_assess.get_entropy(allele_frequencies)
observed_entropy_sliding = np.convolve(observed_entropy, [1] * window_size, 'same')

predicted_entropy = scars_assess.get_entropy(preds)
predicted_entropy_sliding = np.convolve(predicted_entropy, [1] * window_size, 'same')

coverage_indicator = np.zeros(genome[chrom_base][1], dtype=np.uint8)
coverage_indicator[query_indices] = 1
coverage_sliding = np.convolve(coverage_indicator, [1] * window_size, 'same')


# Print to file
output = pd.DataFrame({'chrom': chrom_base, 'start': np.arange(genome[chrom_base][1]), 'end': np.arange(genome[chrom_base][1]) + 1, 
    'obs': observed_entropy_sliding, 'exp': predicted_entropy_sliding, 'coverage': coverage_sliding})

output.to_csv(outFile, sep='\t', header=False, index=False, compression='gzip')




