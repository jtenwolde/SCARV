# script to check the possibility of using ewens sampling fmla

def linearSmooth (X):
    return sum(X)


def match_scores_to_ranges(gr, gr_reliable_mgd, exp_reliable_flat, smoothing=linearSmooth):
    import pyranges as pr
    import pandas as pd
    import numpy as np
    import collections

    gr_reliable_spl = gr_reliable_mgd.tile(1)

    e = pd.Series(data = exp_reliable_flat, name="expected")

    gr_reliable_spl = gr_reliable_spl.insert(e)

    index = pd.Series(range(len(gr)), name="id")
    gr = gr.insert(index) 
    gr_spl = gr.tile(1)

    hits = gr_spl.join(gr_reliable_spl, suffix="_spl")
    exp = hits.as_df().groupby('id')['expected'].agg(smoothing)

    countsByIndex = collections.Counter(hits.id)
    coverage = pd.Series([countsByIndex[ix] if ix in countsByIndex.keys() else 0 for ix in exp.index], name="coverage")
    
    out = pd.concat([exp, coverage], axis=1)

    return out



def evaluate(gr, reliable_sites, reference_fasta, flank, model_bal, model_cal, pop_size, snvs, smoothing=linearSmooth):
    import numpy as np
    import pandas as pd

    reliable_snvs = snvs.intersect(reliable_sites)

    index = pd.Series(range(len(gr)), name="id")
    gr_id = gr.insert(index)

    hits = gr_id.join(reliable_snvs) 
    K = hits.as_df().groupby('id')['ac'].count()
    K_unique = np.unique(hits.as_df().groupby('id')['ac'].count())
 
    thetas_unique = {k: get_theta(k) for k in K_unique}  
 
    thetas = pd.Series([thetas_unique[k] for k in K])
    thetas.name = 'theta_hat'

    entropy = hits.as_df().groupby('id')['ac'].agg(getSkewEntropy, pop_size)
    entropy.name = 'entropy'

    gr_mgd = gr.merge()
    gr_reliable_mgd = gr_mgd.intersect(reliable_sites)

    seqs = scars_queries.query_sequence(gr_reliable_mgd, reference_fasta, flank)
    scars_queries.correct_refs(gr_reliable_mgd, snvs, seqs)    
    scars_fit.anchor_on_AC(seqs)
    
    exp = scars_assess.get_expected(seqs, np.empty((len(seqs), 0)), model_bal, model_cal, pop_size)
    exp_matched = match_scores_to_ranges(gr, gr_reliable_mgd, exp, smoothing)

    return pd.concat([thetas, entropy, exp_matched], axis=1)



def get_theta (k, pop_size=32399, theta_0 = 0.1):
    from scipy.optimize import minimize
    import collections

    theta = minimize(neg_log_likelihood, theta_0, args=(k+1, pop_size), bounds=[(0,None)]).x[0]

    return theta


def neg_log_likelihood (theta, k, pop_size):
    import math 
    import numpy as np

    loglike = k * np.log(theta) - np.sum(np.log(theta + np.arange(2*pop_size)))

    return -loglike


def getSkewEntropy (ACs, pop_size):
    pk = np.array(ACs/(2*pop_size))
    term1 = -np.nansum(pk*np.log2(pk)) 
    term2 = -np.nansum((1-pk)*np.log2(1-pk))
    return (term1+term2)/len(ACs)


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

    percentiles = pd.cut(scores, bins=percentiles, labels=np.arange(0.1,100.1,0.1), include_lowest=True)
    return percentiles


def weighted_mean(x, y, w1, w2):
    return w1*x + w2*y

def weighted_harmonic_mean(x, y, w1, w2):
    return 1/((w1/x)+(w2/y))

def fishers_method(x, y):
    import numpy as np
    return -2*(np.log(x) + np.log(y))


import pyranges as pr
from keras.models import load_model
import keras_genomics
import joblib
import numpy as np
import pandas as pd
from scars import *
from scipy.stats import median_abs_deviation
from statsmodels.distributions.empirical_distribution import ECDF


flank = 6
sampleSize = int(1e5)
smoothingWindowFlank = 250
MAF_range = 'all'    
Model_type = 'singleton'   
pop_size = 32399
codata_args = []


reference_fasta = "/rds/project/who1000-1/rds-who1000-cbrc/ref/UCSC/hg38/hg38.fa"
comSnvsFn = "/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/scars_pipeline_gnomad_hg38/nfe/variants/common_snvs_pass.bed"
rareSnvsFn = "/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/scars_pipeline_gnomad_hg38/nfe/variants/rare_snvs_corrected.bed"

snvsRare = load_snvs(rareSnvsFn)
snvsCom = load_snvs(comSnvsFn)
snvsAllDf = snvsRare.as_df().append(snvsCom.as_df())
snvs = pr.PyRanges(snvsAllDf).sort()
snvs = snvs[snvs.ac>0]

reliable_sites = pr.read_bed("/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/scars_pipeline_gnomad_hg38/nfe/quality_filtering/reliable_sites_incl_common.bed")
orion_cdts_regions = pr.read_bed("/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/competing_tracks/regions_CDTS_Orion.bed")

reliableOverlapWithCompetingScores = reliable_sites.intersect(orion_cdts_regions)
sampledLoci = scars_queries.sample_loci_from_pr(reliableOverlapWithCompetingScores, sampleSize)

cnn = load_model('/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/scars_pipeline_gnomad_hg38/nfe/cnn/cnn_singleton_input_N1e6_seqonly.h5', custom_objects={'RevCompConv1D'\
 : keras_genomics.layers.RevCompConv1D, 'DenseAfterRevcompConv1D': keras_genomics.layers.DenseAfterRevcompConv1D})
calibration_model = joblib.load("/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/scars_pipeline_gnomad_hg38/nfe/calibration/calibration_singleton_input_N1e6_seqonly.h5")

sampledLociExtd = sampledLoci.slack(smoothingWindowFlank)
scars = evaluate(sampledLociExtd, reliable_sites, reference_fasta, flank, cnn, calibration_model, pop_size, snvs)
sufficientCoverage = (scars.coverage > 0.9 * (2 * smoothingWindowFlank + 1)) & (scars.theta_hat != 0.1)  # 90% coverage required, theta==0.1 indicates failure
scarsCovered = scars[sufficientCoverage]

theta_norm = (scarsCovered.theta_hat - np.mean(scarsCovered.theta_hat))/median_abs_deviation(scarsCovered.theta_hat)
exp_norm = (scarsCovered.expected - np.mean(scarsCovered.expected))/median_abs_deviation(scarsCovered.expected)

ecdf_lift = ECDF(theta_norm - exp_norm)
ecdf_tilt = ECDF(scarsCovered.entropy)

quantiles_lift = ecdf_lift(theta_norm - exp_norm)
quantiles_tilt = ecdf_tilt(scarsCovered.entropy)

patho_snvs = pr.read_bed("/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/non_coding_pathogenic_variants/noncoding_pathogenic_annot_by_Orion_and_CDTS_including_scores.bed")
patho_snvs.columns = ["Chromosome", "Start", "End", "CDTS", "Orion"]

patho_snvs_extd = patho_snvs.slack(smoothingWindowFlank)
scarsPatho = evaluate(patho_snvs_extd, reliable_sites, reference_fasta, flank, cnn, calibration_model, pop_size, snvs)

theta_patho_norm = (scarsPatho.theta_hat - np.mean(scarsCovered.theta_hat))/median_abs_deviation(scarsCovered.theta_hat)
exp_patho_norm = (scarsPatho.expected - np.mean(scarsCovered.expected))/median_abs_deviation(scarsCovered.expected)

quantiles_lift_patho = ecdf_lift(theta_patho_norm - exp_patho_norm)
quantiles_tilt_patho = ecdf_tilt(scarsPatho.entropy)


# Evaluating various ways of combining tilt and lift (theta)
w1 = 0.9
w2 = 0.1

weighted_combo = [weighted_mean(x, y, w1, w2) for x, y in zip(quantiles_lift, quantiles_tilt)]
weighted_percentiles = np.quantile(weighted_combo, np.arange(0, 1.001, 0.001))

weighted_combo_patho = [weighted_mean(x, y, w1, w2) for x, y in zip(quantiles_lift_patho, quantiles_tilt_patho)]
weighted_percentiles_patho = toPercentile(weighted_combo_patho, weighted_percentiles)


fisher_combo = [fishers_method(x, y) for x, y in zip(quantiles_lift, quantiles_tilt)]
fisher_percentiles = np.quantile(fisher_combo, np.arange(0, 1.001, 0.001))

fisher_combo_patho = [fishers_method(x, y) for x, y in zip(quantiles_lift_patho, quantiles_tilt_patho)]
fisher_percentiles_patho = np.array([round(100 - x,1) for x in toPercentile(fisher_combo_patho, fisher_percentiles)])


v1 = 0.9
v2 = 0.1

harmonic_combo = [weighted_harmonic_mean(x, y, v1, v2) for x, y in zip(quantiles_lift, quantiles_tilt)]
harmonic_percentiles = np.quantile(harmonic_combo, np.arange(0, 1.001, 0.001))

harmonic_combo = [weighted_harmonic_mean(x, y, v1, v2) for x, y in zip(quantiles_lift_patho, quantiles_tilt_patho)]
harmonic_percentiles_patho = toPercentile(harmonic_combo, harmonic_percentiles)







