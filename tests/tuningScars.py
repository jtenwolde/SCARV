# script to compare scars, orion and cdts by means of scatter and ecdfs 
def load_snvs(fn):
    import pyranges as pr

    snvs = pr.read_bed(fn) 
    snvs.columns = ["Chromosome", "Start", "End", "ref", "alt", "ac", "an"]
    snvs.ac = snvs.ac.astype(int)

    return snvs


def toPercentile(scores, percentiles):
    import pandas as pd
    import numpy as np

    # extend boundaries in case min/max extends beyond sampled min/max
    percentiles[0] = -np.inf
    percentiles[-1] = np.inf 

    assert len(np.unique(percentiles)) == len(percentiles), "Percentile values are non-unique"

    percentiles = pd.cut(scores, bins=percentiles, labels=np.arange(1,101,1), include_lowest=True)
    return percentiles


import pyranges as pr
from keras.models import load_model
import keras_genomics
import joblib
import numpy as np
import pandas as pd
from scars import *
from scipy.stats import median_abs_deviation

import sys


flank = 6
sampleSize = int(1e5)
smoothingWindowFlank = int(sys.argv[1])
MAF_range = sys.argv[2]    # choice of ['singleton', 'rare', 'all', 'common', 'very_common']
Model_type = sys.argv[3]    # choice of ['singleton', 'rare']
outFn = sys.argv[4]
pop_size = 32399
codata_args = []


reference_fasta = "/rds/project/who1000-1/rds-who1000-cbrc/ref/UCSC/hg38/hg38.fa"
comSnvsFn = "/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/lil_project/common_snvs_pass.bed"
rareSnvsFn = "/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/scars_pipeline_gnomad_hg38/rare_snvs_corrected.bed"


if MAF_range in ['common', 'very_common']:
    snvs = load_snvs(comSnvsFn)
    if MAF_range == 'very_common':
        snvs = snvs[snvs.ac/snvs.an > 1/100]
elif MAF_range in ['singleton', 'rare']:
    snvs = load_snvs(rareSnvsFn)
    if MAF_range == 'singleton':
        snvs = snvs[snvs.ac == 1]
elif MAF_range == 'all':
    snvsRare = load_snvs(rareSnvsFn)
    snvsCom = load_snvs(comSnvsFn)
    snvsAllDf = snvsRare.as_df().append(snvsCom.as_df())
    snvs = pr.PyRanges(snvsAllDf).sort()


if MAF_range in ['singleton', 'rare']:
    reliable_sites = pr.read_bed("/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/scars_pipeline_gnomad_hg38/reliable_sites.bed")
else:
    reliable_sites = pr.read_bed("/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/lil_project/reliable_sites_lil_project.bed")


orion_cdts_regions = pr.read_bed("/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/competing_tracks/regions_CDTS_Orion.bed")

reliableOverlapWithCompetingScores = reliable_sites.intersect(orion_cdts_regions)
sampledLoci = scars_queries.sample_loci_from_pr(reliableOverlapWithCompetingScores, sampleSize)


if Model_type == 'singleton':
    cnn = load_model('/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/lil_project/cnn_seq_lil_proj.h5', custom_objects={'RevCompConv1D'\
     : keras_genomics.layers.RevCompConv1D, 'DenseAfterRevcompConv1D': keras_genomics.layers.DenseAfterRevcompConv1D})
    calibration_model = joblib.load("/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/lil_project/calibration_seq_lil_proj.h5")
else:
    cnn = load_model('/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/scars_pipeline_gnomad_hg38/cnn_seq.h5', custom_objects={'RevCompConv1D'\
     : keras_genomics.layers.RevCompConv1D, 'DenseAfterRevcompConv1D': keras_genomics.layers.DenseAfterRevcompConv1D})
    calibration_model = joblib.load("/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/scars_pipeline_gnomad_hg38/calibration_seq.h5")


sampledLociExtd = sampledLoci.slack(smoothingWindowFlank)
scars = scars_assess.evaluate(sampledLociExtd, reliable_sites, reference_fasta, flank, cnn, calibration_model, pop_size, snvs, codata_args, split=True)
sufficientCoverage = scars.coverage > 0.9 * (2 * smoothingWindowFlank + 1)  # 90% coverage required
scarsCovered = scars[sufficientCoverage]

obs_norm = (scarsCovered.observed - np.mean(scarsCovered.observed))/median_abs_deviation(scarsCovered.observed)
exp_norm = (scarsCovered.expected - np.mean(scarsCovered.expected))/median_abs_deviation(scarsCovered.expected)


percentiles = np.quantile(obs_norm - exp_norm, np.arange(0, 1.01, 0.01))


patho_snvs = pr.read_bed("/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/non_coding_pathogenic_variants/noncoding_pathogenic_annot_by_Orion_and_CDTS_including_scores.bed")
patho_snvs.columns = ["Chromosome", "Start", "End", "CDTS", "Orion"]


patho_snvs_extd = patho_snvs.slack(smoothingWindowFlank)
scarsPatho = scars_assess.evaluate(patho_snvs_extd, reliable_sites, reference_fasta, flank, cnn, calibration_model, pop_size, snvs, codata_args, split=True)


obs_patho_norm = (scarsPatho.observed - np.mean(scarsCovered.observed))/median_abs_deviation(scarsCovered.observed)
exp_patho_norm = (scarsPatho.expected - np.mean(scarsCovered.expected))/median_abs_deviation(scarsCovered.expected)


pathoSnvsAnnot = pr.PyRanges(patho_snvs.as_df().loc[scarsPatho.index])
percentilesPatho = toPercentile(obs_patho_norm - exp_patho_norm, percentiles)


np.savetxt(outFn, percentilesPatho.astype(int).to_numpy(), fmt='%1.1d')


