# Retry of exponential smooth 

def evaluate(gr, reliable_sites, reference_fasta, flank, model_bal, model_cal, pop_size, snvs, smoothing):
    import numpy as np
    import pandas as pd

    reliable_snvs = snvs.intersect(reliable_sites)

    gr_mgd = gr.merge()
    gr_reliable_mgd = gr_mgd.intersect(reliable_sites)

    seqs = scars_queries.query_sequence(gr_reliable_mgd, reference_fasta, flank)
    scars_queries.correct_refs(gr_reliable_mgd, snvs, seqs)    
    scars_fit.anchor_on_AC(seqs)
    
    exp = scars_assess.get_expected(seqs, np.empty((len(seqs), 0)), model_bal, model_cal, pop_size)
    K = get_k(gr_reliable_mgd, reliable_snvs)

    scores_matched = match_scores_to_ranges(gr, gr_reliable_mgd, K, exp, smoothing)

    return scores_matched


def get_k(gr, snvs_pr):
    import pandas as pd
    import numpy as np
    import pyranges as pr

    gr_spl = gr.tile(1)
    index = pd.Series(range(len(gr_spl)), name="id")
    gr_spl = gr_spl.insert(index)

    snvHits = gr_spl.join(snvs_pr, suffix="_snv").as_df()
    anyHits = (snvHits.shape[0] != 0)

    observed = np.zeros(len(gr_spl))
    if anyHits:
        ac_by_id = snvHits.groupby(['id']).agg({'ac': "count"})
        observed[ac_by_id.index] = ac_by_id['ac']

    return observed


def match_scores_to_ranges(gr, gr_reliable_mgd, k_reliable_flat, exp_reliable_flat, smoothing):
    import pyranges as pr
    import pandas as pd
    import numpy as np
    import collections

    gr_reliable_spl = gr_reliable_mgd.tile(1)

    k_i = pd.Series(data = k_reliable_flat, name="K")
    e = pd.Series(data = exp_reliable_flat, name="expected")

    gr_reliable_spl = gr_reliable_spl.insert(k_i)
    gr_reliable_spl = gr_reliable_spl.insert(e)

    index = pd.Series(range(len(gr)), name="id")
    gr = gr.insert(index) 
    gr_spl = gr.tile(1)

    # unmapped hits will be assigned K=1 & expected=1, exponentialSmooth will set this to 0
    hits = gr_spl.join(gr_reliable_spl, suffix="_spl", how='left')
    out = hits.as_df().groupby('id')[['K', 'expected']].agg(smoothing)

    countsByIndex = collections.Counter(hits[hits.K!=-1].id)
    out['coverage'] = [countsByIndex[ix] if ix in countsByIndex.keys() else 0 for ix in out.index]    

    return out


def exponentialSmooth (X, Omega=0.993):
    import numpy as np

    X = np.array(X)
    X[X==-1] = 0

    windowSize = len(X)//2

    S = 0
    for i in range(windowSize+1):
        S += (Omega**i) * (X[windowSize + i] + X[windowSize - i])/2

    S *= (1 - Omega) / (1 - (Omega**(windowSize+1)))

    return S



def betaSmooth (X, Omega=0.993):
    import numpy as np
    from scipy.stats import beta

    X = np.array(X)
    X[X==-1] = 0

    weightings = beta.pdf(np.linspace(0, 1, len(X)+2), Omega, Omega)[1:-1]
    weightings /= np.sum(weightings)

    S = np.sum(weightings * X)

    return S



def get_theta (k, pop_size=32399, theta_0 = 0.1):
    from scipy.optimize import minimize
    import collections

    theta = minimize(neg_log_likelihood, theta_0, args=(k, pop_size), bounds=[(0,None)]).x[0]

    return theta


def neg_log_likelihood (theta, k, pop_size):
    import math 
    import numpy as np

    loglike = (k+1) * np.log(theta) - np.sum(np.log(theta + np.arange(2*pop_size)))

    return -loglike


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


def percScoreToCumulativePercCountPlot(filelist, plt=None, x_upper=None):
    percScore = concatPercentiles(filelist)
    cumulativeCounts = getCumulativeCounts(percScore)
    plt = plotCumulativePercCountPlot(cumulativeCounts, plt, x_upper)
    return plt


def concatPercentiles (filelist):
    import numpy as np

    l = []
    for fn in filelist:
        arr = np.loadtxt(fn)
        l.extend(arr)

    return l


def getCumulativeCounts(percentile_values):
    import collections
    import numpy as np

    percentileTally = collections.Counter([int(10*x) for x in percentile_values])
    addMissingPercentiles(percentileTally)

    counts = [count for perc, count in sorted(percentileTally.items())]
    countsCumulative = np.cumsum(counts)

    return countsCumulative


def plotCumulativePercCountPlot(countsCumulative, plt=None, x_upper=None):
    import numpy as np

    percentiles = np.arange(0.1, 100.1, 0.1)
    plt.plot(percentiles, countsCumulative)
    
    plt.set_xlabel("Percentile")
    plt.set_ylabel("Cumulative count")

    plt.set_ylim(0, 1.1 * countsCumulative[-1])

    if x_upper is not None:
        plt.set_xlim(0, x_upper)
        plt.set_ylim(0, 1.1 * countsCumulative[x_upper*10])

    return plt


def addMissingPercentiles(percentileCountDictionary):
    import numpy as np

    percentileSetComplete = set(np.arange(1, 1001, 1))
    percentilesPresent = set(percentileCountDictionary.keys())
    
    toAdd = percentileSetComplete - percentilesPresent
    for i in list(toAdd):
        percentileCountDictionary[i] = 0

    return None


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
smoothingWindowFlank = 250  
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


for omega in np.exp(np.log([0.1, 0.35, 0.6, 0.85])/250):
    def customSmooth (x):
        return exponentialSmooth(x, Omega=omega)

    sampledLociExtd = sampledLoci.slack(smoothingWindowFlank)
    scars = evaluate(sampledLociExtd, reliable_sites, reference_fasta, flank, cnn, calibration_model, pop_size, snvs, customSmooth)
    sufficientCoverage = scars.coverage > 0.9 * (2 * smoothingWindowFlank + 1)
    scarsCovered = scars[sufficientCoverage]

    K_norm = (scarsCovered.K - np.mean(scarsCovered.K))/median_abs_deviation(scarsCovered.K)
    exp_norm = (scarsCovered.expected - np.mean(scarsCovered.expected))/median_abs_deviation(scarsCovered.expected)

    percentiles = np.quantile(K_norm - exp_norm, np.arange(0, 1.001, 0.001))

    patho_snvs = pr.read_bed("/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/non_coding_pathogenic_variants/noncoding_pathogenic_annot_by_Orion_and_CDTS_including_scores.bed")
    patho_snvs.columns = ["Chromosome", "Start", "End", "CDTS", "Orion"]

    patho_snvs_extd = patho_snvs.slack(smoothingWindowFlank)
    scarsPatho = evaluate(patho_snvs_extd, reliable_sites, reference_fasta, flank, cnn, calibration_model, pop_size, snvs, customSmooth)

    K_patho_norm = (scarsPatho.K - np.mean(scarsCovered.K))/median_abs_deviation(scarsCovered.K)
    exp_patho_norm = (scarsPatho.expected - np.mean(scarsCovered.expected))/median_abs_deviation(scarsCovered.expected)

    percentilesPatho = toPercentile(K_patho_norm - exp_patho_norm, percentiles)

    outFn = "/home/jwt44/pathoPercentiles_250_all_singleton_ExpoSmooth_" + str(round(omega, 5)) + ".txt"
    np.savetxt(outFn, percentilesPatho.astype(float).to_numpy(), fmt='%1.1f')



filenames = ["/home/jwt44/pathoPercentiles_250_all_singleton_ExpoSmooth_0.99083.txt",
            "/home/jwt44/pathoPercentiles_250_all_singleton_ExpoSmooth_0.99581.txt",
            "/home/jwt44/pathoPercentiles_250_all_singleton_ExpoSmooth_0.99796.txt",
            "/home/jwt44/pathoPercentiles_250_all_singleton_ExpoSmooth_0.99935.txt"]


import matplotlib.pyplot as plt
fig, axs = plt.subplots(1, 2, figsize=(10,5))

percScoreToCumulativePercCountPlot(filenames[:1], axs[0])
percScoreToCumulativePercCountPlot(filenames[1:2], axs[0])
percScoreToCumulativePercCountPlot(filenames[2:3], axs[0])
percScoreToCumulativePercCountPlot(filenames[3:4], axs[0])

percScoreToCumulativePercCountPlot(filenames[:1], axs[1], 5)
percScoreToCumulativePercCountPlot(filenames[1:2], axs[1], 5)
percScoreToCumulativePercCountPlot(filenames[2:3], axs[1], 5)
percScoreToCumulativePercCountPlot(filenames[3:4], axs[1], 5)

axs[0].legend(['0.99083', '0.99581', '0.99796', '0.99935'])
axs[1].legend(['0.99083', '0.99581', '0.99796', '0.99935'])

axs[0].set_xlabel("Percentile scores", fontsize=14)
axs[1].set_xlabel("Percentile scores", fontsize=14)

fig.savefig("/home/jwt44/smoothingComparison.pdf") 




for omega in np.array([0.1, 0.5, 0.9, 1, 1.1, 1.5, 2]):
    def customSmooth (x):
        return betaSmooth(x, Omega=omega)

    sampledLociExtd = sampledLoci.slack(smoothingWindowFlank)
    scars = evaluate(sampledLociExtd, reliable_sites, reference_fasta, flank, cnn, calibration_model, pop_size, snvs, customSmooth)
    sufficientCoverage = scars.coverage > 0.9 * (2 * smoothingWindowFlank + 1)
    scarsCovered = scars[sufficientCoverage]

    K_norm = (scarsCovered.K - np.mean(scarsCovered.K))/median_abs_deviation(scarsCovered.K)
    exp_norm = (scarsCovered.expected - np.mean(scarsCovered.expected))/median_abs_deviation(scarsCovered.expected)

    percentiles = np.quantile(K_norm - exp_norm, np.arange(0, 1.001, 0.001))

    patho_snvs = pr.read_bed("/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/non_coding_pathogenic_variants/noncoding_pathogenic_annot_by_Orion_and_CDTS_including_scores.bed")
    patho_snvs.columns = ["Chromosome", "Start", "End", "CDTS", "Orion"]

    patho_snvs_extd = patho_snvs.slack(smoothingWindowFlank)
    scarsPatho = evaluate(patho_snvs_extd, reliable_sites, reference_fasta, flank, cnn, calibration_model, pop_size, snvs, customSmooth)

    K_patho_norm = (scarsPatho.K - np.mean(scarsCovered.K))/median_abs_deviation(scarsCovered.K)
    exp_patho_norm = (scarsPatho.expected - np.mean(scarsCovered.expected))/median_abs_deviation(scarsCovered.expected)

    percentilesPatho = toPercentile(K_patho_norm - exp_patho_norm, percentiles)

    outFn = "/home/jwt44/pathoPercentiles_250_all_singleton_BetaSmooth_" + str(round(omega, 5)) + ".txt"
    np.savetxt(outFn, percentilesPatho.astype(float).to_numpy(), fmt='%1.1f')



filenames = ["pathoPercentiles_250_all_singleton_BetaSmooth_0.1.txt",
"pathoPercentiles_250_all_singleton_BetaSmooth_0.5.txt",
"pathoPercentiles_250_all_singleton_BetaSmooth_0.9.txt",
"pathoPercentiles_250_all_singleton_BetaSmooth_1.0.txt",
"pathoPercentiles_250_all_singleton_BetaSmooth_1.1.txt",
"pathoPercentiles_250_all_singleton_BetaSmooth_1.5.txt",
"pathoPercentiles_250_all_singleton_BetaSmooth_2.0.txt"]


import matplotlib.pyplot as plt
fig, axs = plt.subplots(1, 2, figsize=(10,5))

percScoreToCumulativePercCountPlot(filenames[:1], axs[0])
percScoreToCumulativePercCountPlot(filenames[1:2], axs[0])
percScoreToCumulativePercCountPlot(filenames[2:3], axs[0])
percScoreToCumulativePercCountPlot(filenames[3:4], axs[0])
percScoreToCumulativePercCountPlot(filenames[4:5], axs[0])
percScoreToCumulativePercCountPlot(filenames[5:6], axs[0])
percScoreToCumulativePercCountPlot(filenames[6:7], axs[0])

percScoreToCumulativePercCountPlot(filenames[:1], axs[1], 5)
percScoreToCumulativePercCountPlot(filenames[1:2], axs[1], 5)
percScoreToCumulativePercCountPlot(filenames[2:3], axs[1], 5)
percScoreToCumulativePercCountPlot(filenames[3:4], axs[1], 5)
percScoreToCumulativePercCountPlot(filenames[4:5], axs[1], 5)
percScoreToCumulativePercCountPlot(filenames[5:6], axs[1], 5)
percScoreToCumulativePercCountPlot(filenames[6:7], axs[1], 5)

axs[0].legend(['0.1', '0.5', '0.9', '1.0', '1.1', '1.5', '2.0'])
axs[1].legend(['0.1', '0.5', '0.9', '1.0', '1.1', '1.5', '2.0'])

axs[0].set_xlabel("Percentile scores", fontsize=14)
axs[1].set_xlabel("Percentile scores", fontsize=14)

fig.savefig("/home/jwt44/smoothingComparisonBeta.pdf") 








