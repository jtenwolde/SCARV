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

k = 25
genome = pybedtools.genome_registry.hg38
reference_fasta = "/rds/project/who1000-1/rds-who1000-cbrc/ref/UCSC/hg38/hg38.fa"

snvs = scars_queries.load_variants("/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/scars_pipeline_gnomad_hg38/" + ancestry +"/variants/snvs_by_chr/pass_snvs_" + chrom_base + ".bed")

deletions = scars_queries.load_variants("/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/scars_pipeline_gnomad_hg38/" + ancestry +"/variants/pass_deletions.bed")[chrom_base]
deletions.ref = deletions.ref.str.slice(1)
deletions.alt = "X"
deletions = deletions[deletions.ref.str.len() == 1]
deletions.Start += 1
deletions = deletions[deletions.lengths() == 1]

snvs_and_deletions = pr.PyRanges(pd.concat([snvs.as_df(), deletions.as_df()]))

insertions = scars_queries.load_variants("/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/scars_pipeline_gnomad_hg38/" + ancestry +"/variants/pass_insertions.bed")[chrom_base]
insertions.ref = "X"
insertions.alt = insertions.alt.str.slice(1)
insertions = insertions[insertions.alt.str.len() == 1]

reliable_sites = pr.read_bed("/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/scars_pipeline_gnomad_hg38/" + ancestry + "/quality_filtering/reliable_sites_by_chr/reliable_sites_" + chrom + ".bed")


preds_by_kmer_file = "/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/scars_pipeline_gnomad_hg38/" + ancestry + "/preds_by_kmer.pkl"
preds_by_kmer = pd.read_pickle(preds_by_kmer_file)

chrom_gr = pr.from_dict({'Chromosome': [chrom_base], 'Start': [0], 'End': [genome[chrom_base][1]]})
chrom_seq = pr.get_fasta(chrom_gr, reference_fasta)[0].upper()
chrom_seq_padded = "X".join(chrom_seq)
chrom_seq_ohe = np.array(pd.get_dummies(pd.Series(list(chrom_seq_padded)))[['A','C','X','G','T']])


# Compute chromosome-wide predictions
k_mers = pd.Series([chrom_seq_padded[(pos-k//2):(pos+k//2+1)] for pos in range(len(chrom_seq_padded)-k//2)])


preds = np.empty(shape=(2*genome[chrom_base][1]-1, 5))
preds[:] = np.nan

reliable_indices = np.sort(np.concatenate([2*reliable_sites.tile(1).Start.to_numpy(), 2*reliable_sites.tile(1).Start.to_numpy()+1])) 
identified_nucleotides = np.where(["N" not in seq for seq in k_mers])[0]
query_indices = np.intersect1d(reliable_indices, identified_nucleotides)

preds[query_indices] = preds_by_kmer.loc[k_mers[query_indices]].to_numpy()


# Compute chromosome-wide allele frequencies
allele_counts = 2 * pop_size * chrom_seq_ohe.astype('int32')

allele_counts = scars_assess.process_variants (insertions, allele_counts, insertion=True)
allele_counts = scars_assess.process_variants (snvs_and_deletions, allele_counts, insertion=False)

allele_counts[np.setdiff1d(np.arange(allele_counts.shape[0]), query_indices)] = 0   # remove contribution from unreliable or unidentified sequences
allele_frequencies = allele_counts / allele_counts.sum(axis=1)[:, np.newaxis]


# Generate score
observed_entropy = scars_assess.get_entropy(allele_frequencies)
observed_entropy_sliding = np.convolve(observed_entropy, [1] * 2 * window_size, 'same')[::2]

predicted_entropy = scars_assess.get_entropy(preds)
predicted_entropy_sliding = np.convolve(predicted_entropy, [1] * 2 * window_size, 'same')[::2]

coverage_indicator = np.zeros(2*genome[chrom_base][1]-1, dtype=np.uint8)
coverage_indicator[query_indices] = 1
coverage_sliding = np.convolve(coverage_indicator, [1] * 2 * window_size, 'same')[::2]


# Print to file
output = pd.DataFrame({'chrom': chrom_base, 'start': np.arange(genome[chrom_base][1]), 'end': np.arange(genome[chrom_base][1]) + 1, 
    'obs': observed_entropy_sliding, 'exp': predicted_entropy_sliding, 'coverage': coverage_sliding})

output.to_csv(outFile, sep='\t', header=False, index=False, compression='gzip')


