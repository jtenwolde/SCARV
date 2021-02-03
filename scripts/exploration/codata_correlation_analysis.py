def load_snvs(fn):
    import pyranges as pr

    snvs = pr.read_bed(fn) 
    snvs.columns = ["Chromosome", "Start", "End", "ref", "alt", "ac", "an"]
    snvs.ac = snvs.ac.astype(int)

    return snvs


def query_bigWigs(gr, filename_list, impute_values):
    import pyranges as pr
    import pyBigWig
    import numpy as np

    N = gr.length
    out = np.empty(shape=(N, len(filename_list)))
    for ix, file in enumerate(filename_list):
        bw = pyBigWig.open(file)
        scores_nested = [bw.values(row.Chromosome, row.Start, row.End) for row in gr.as_df().itertuples()]
        scores = np.array([item for sublist in scores_nested for item in sublist])
        scores[scores==np.nan] = impute_values[ix]
        out[:,ix] = scores
        bw.close()

    return out


def get_correlations(df, colname, minBinSize, nBins):
	df_tmp = df[['scars', 'phyloP', colname]]
	df_tmp = df_tmp[~np.isnan(df_tmp[colname])]
	df_tmp['binned_midpoints'] = [(a.left + a.right)/2 for a in pd.qcut(df_tmp[colname], nBins, duplicates='drop')]
	grpd_means = df_tmp.groupby('binned_midpoints')[['scars', 'phyloP']].mean()
	grpd_means['count'] = df_tmp.groupby('binned_midpoints')[colname].count()
	reliable_grpd_means = grpd_means[grpd_means['count'] >= minBinSize]
	corr_scars = np.corrcoef(reliable_grpd_means.scars, reliable_grpd_means.index)[0,1]
	corr_phyloP = np.corrcoef(reliable_grpd_means.phyloP, reliable_grpd_means.index)[0,1]
	return corr_scars, corr_phyloP, reliable_grpd_means.shape[0]



import pyranges as pr
from scars import *
import pandas as pd
import numpy as np


flank = 6
pop_size = 32399
n_samples = int(1e5)

training_sites = pr.read_bed("/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/scars_pipeline_gnomad_hg38/nfe/training_sites/training_loci_phyloP_min3_to_0.bed")
sampled_training_loci = scars_queries.sample_loci_from_pr(training_sites, n_samples)


# annotate with observed

# load singletons
rareSnvsFn = "/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/scars_pipeline_gnomad_hg38/nfe/variants/rare_snvs_corrected.bed"
snvs = load_snvs(rareSnvsFn)
singletons = snvs[snvs.ac/snvs.an < 1.5/(2*pop_size)]
singletons = singletons[singletons.ac > 0]

index = pd.Series(range(len(sampled_training_loci)), name='ix')
sampled_training_loci = sampled_training_loci.insert(index)
hits = sampled_training_loci.join(singletons)
acByIndex = hits.as_df().groupby('ix')['ac'].sum()

obs = np.zeros(n_samples)
obs[acByIndex.index] = acByIndex


# query sequence
reference_fasta = "/rds/project/who1000-1/rds-who1000-cbrc/ref/UCSC/hg38/hg38.fa"
sequence = scars_queries.query_sequence(sampled_training_loci, reference_fasta, flank)
scars_fit.anchor_on_AC(sequence)

# load models
from keras.models import load_model
import keras_genomics
import joblib

cnn = load_model('/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/scars_pipeline_gnomad_hg38/nfe/cnn/cnn_singleton_input_N1e6_seqonly.h5', custom_objects={'RevCompConv1D'\
: keras_genomics.layers.RevCompConv1D, 'DenseAfterRevcompConv1D': keras_genomics.layers.DenseAfterRevcompConv1D})
calibration_model = joblib.load("/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/scars_pipeline_gnomad_hg38/nfe/calibration/calibration_singleton_input_N1e6_seqonly.h5")

# annotate with expected
exp = scars_assess.get_expected(sequence, np.empty(shape=(sequence.shape[0], 0)), cnn, calibration_model, pop_size)


# annotate with phyloP 
phyloP_bw = "/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/hg38.phyloP100way.bw"
phyloP_scores = query_bigWigs(sampled_training_loci, [phyloP_bw], [np.nan])


# annotate with codata 
import os
import pybedtools 

genome = pybedtools.genome_registry.hg38

chr_list = ['chr' + str(i) for i in np.arange(1,23,1)]
chr_list.append('chrX')
chr_list.append('chrY')

epidata_folder = "/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/epigenomic_tracks/"
epidata_files = os.listdir(epidata_folder)
non_methylation_bw_files = sorted([epidata_folder + f for f in epidata_files if "WGBS" not in f])
human_ovary_impute_vals_raw = list(map(lambda x: scars_queries.average_bw (x, chr_list, genome), non_methylation_bw_files))
human_ovary_impute_vals = [0 if "DNase" in fn else x for x,fn in zip(human_ovary_impute_vals_raw, non_methylation_bw_files)]
CpG_methylation_files = [epidata_folder + f for f in epidata_files if "WGBS" in f]

codata = query_bigWigs(sampled_training_loci, non_methylation_bw_files, human_ovary_impute_vals)
codata_meth = query_bigWigs(sampled_training_loci, CpG_methylation_files, [np.nan, np.nan])


M = np.c_[obs, exp, phyloP_scores, codata, codata_meth]

df = pd.DataFrame(M)
df.columns = ['obs', 'exp', 'phyloP', 'ATAC_ovary', 'ATAC_testis', 
	'DNase_ovary', 'DNase_testis', 'H3K27ac_testis', 'H3K27me3_testis', 
	'H3K36me3_ovary', 'H3K4me1_testis', 'H3K4me3_ovary', 'H3K4me3_testis', 
	'H3K9me3_ovary', 'H3K9me3_testis', 'WGBS_ovary', 'WGBS_testis']

df['scars'] = df['obs']-df['exp']

codata_columns = ['ATAC_ovary', 'ATAC_testis', 
	'DNase_ovary', 'DNase_testis', 'H3K27ac_testis', 'H3K27me3_testis', 
	'H3K36me3_ovary', 'H3K4me1_testis', 'H3K4me3_ovary', 'H3K4me3_testis', 
	'H3K9me3_ovary', 'H3K9me3_testis', 'WGBS_ovary', 'WGBS_testis']


for name in codata_columns:
	cors = get_correlations(df, name, 500, 50)
	print(name, 'Correlation O-E:', round(cors[0], 2), 'Correlation phyloP:', round(cors[1], 2))




