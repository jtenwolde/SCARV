# script to compare scars, orion and cdts by means of scatter and ecdfs 

def plotCumulativePercCountPlot(countsCumulative, plt=None, x_upper=None):
    import numpy as np

    percentiles = np.arange(100 + 1)

    plt.plot(percentiles, countsCumulative)
    
    plt.set_xlabel("Percentile")
    plt.set_ylabel("Cumulative count")

    ylim_lower = -50
    ylim_upper = countsCumulative[-1] + 50

    plt.set_ylim(ylim_lower, ylim_upper)
    
    if x_upper is not None:
        plt.set_xlim(0, x_upper)
        plt.set_ylim(0, countsCumulative[x_upper] + 50)

    return plt


def getCumulativeCounts(percentile_values):
    import collections
    import numpy as np

    percentileTally = collections.Counter(percentile_values)
    addMissingPercentiles(percentileTally)

    counts = [count for perc, count in sorted(percentileTally.items())]
    countsCumulative = np.cumsum(counts)

    return countsCumulative


def addMissingPercentiles(percentileCountDictionary):
    percentileSetComplete = set(np.arange(100 + 1))
    percentilesPresent = set(percentileCountDictionary.keys())
    
    toAdd = percentileSetComplete - percentilesPresent
    for i in list(toAdd):
        percentileCountDictionary[i] = 0

    return None


def percScoreToCumulativePercCountPlot(percScore, plt=None, x_upper=None):
    cumulativeCounts = getCumulativeCounts(percScore)
    plt = plotCumulativePercCountPlot(cumulativeCounts, plt, x_upper)
    return plt


def toPercentile(scores, percentiles):
    import pandas as pd
    import numpy as np

    # extend boundaries in case min/max extends beyond sampled min/max
    percentiles[0] = -np.inf
    percentiles[-1] = np.inf 

    assert len(np.unique(percentiles)) == len(percentiles), "Percentile values are non-unique"

    percentiles = pd.cut(scores, bins=percentiles, labels=np.arange(1,101,1), include_lowest=True)
    return percentiles


def plotHist(x):
    import matplotlib.pyplot as plt

    plt.hist(x, bins = 100)

    plt.xlabel("Percentile")
    plt.ylabel("Count")

    return plt


def savePlotToPdf(plot, outFile):
    import matplotlib.pyplot as plt

    plot.savefig(outFile) 
    plot.close() 


def scoreToPercentileHistogram(score, percentiles, plot_fn):
    score_percentile = toPercentile(score, percentiles)
    histogram = plotHist(score_percentile)
    savePlotToPdf(histogram, plot_fn)
    return 


def load_snvs(fn):
    import pyranges as pr

    snvs = pr.read_bed(fn) 
    snvs.columns = ["Chromosome", "Start", "End", "ref", "alt", "ac", "an"]
    snvs.ac = snvs.ac.astype(int)

    return snvs


def loadExonsAsPr(fn):
    import pandas as pd
    import pyranges as pr

    exons_df = pd.read_csv(fn, sep='\t', header=None)
    relevantColumns = [0,1,2,10]
    exons_annot_df = exons_df.iloc[:, relevantColumns]
    exons_annot_df.columns = ["Chromosome", "Start", "End", "constr"]
    exons_annot = pr.PyRanges(exons_annot_df)

    return exons_annot



def get_codata(gr, human_ovary_bw_files=None, human_ovary_impute_vals=None, CpG_methylation_files=None, exons_annot=None):
    import numpy as np

    if human_ovary_bw_files is not None:
        human_ovary_tracks = query_codata_human_ovary(gr, human_ovary_bw_files, human_ovary_impute_vals)
    else:
        human_ovary_tracks = np.empty(shape=(gr.length,0))

    if CpG_methylation_files is not None:
        CpG_methylation = query_codata_CpG_meth(gr, CpG_methylation_files)
    else:
        CpG_methylation = np.empty(shape=(gr.length,0))
    
    if exons_annot is not None:
        exon_spline_basis = query_nearest_exon(gr, exons_annot)
    else:
        exon_spline_basis = np.empty(shape=(gr.length,0))
   

    codata = np.c_[human_ovary_tracks, CpG_methylation, exon_spline_basis]
    return codata



def query_codata_CpG_meth(gr, filename_list):
    import pyBigWig
    import numpy as np

    N = gr.length
    out = np.empty(shape=(N, 4 * len(filename_list)), dtype=np.int8)
    for ix, file in enumerate(filename_list):
        bw = pyBigWig.open(file)
        scores_nested = [bw.values(row.Chromosome, row.Start, row.End) for row in gr.as_df().itertuples()]
        scores = np.array([item for sublist in scores_nested for item in sublist])
        scores[np.isnan(scores)] = -1
        arr = one_hot_encode_methylation_signal(scores)
        out[:,(4*ix):(4*(ix+1))] = arr
        bw.close()

    return out


def one_hot_encode_methylation_signal(scores):
        import numpy as np

        N = len(scores)
        arr = np.zeros(shape=(N, 4), dtype=np.int8)

        arr[scores>60, 3] = 1
        arr[(scores>20) & (scores<=60), 2] = 1 
        arr[(scores>=0) & (scores<=20), 1] = 1
        arr[scores==-1, 0] = 1

        return arr 


# assumes bigWig files with scores that don't need binning/ohe
# furthermore it is useful to sort the filename_list to make covariates identifiable
def query_codata_human_ovary(gr, filename_list, impute_values):
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



# functions required for generating the spline bases later
# kv: knot vector, u: samples, k: index of control point, d: degree
def coxDeBoor(k, d, u, kv):
    # Test for end conditions
    if (d == 0):
        return ((u - kv[k] >= 0) & (u - kv[k + 1] <= 0))*1

    denom1 = kv[k + d] - kv[k]
    term1 = 0
    if denom1 > 0:
        term1 = ((u - kv[k]) / denom1) * coxDeBoor(k, d - 1,u,kv)

    denom2 = kv[k + d + 1] - kv[k + 1]
    term2 = 0
    if denom2 > 0:
        term2 = ((-(u - kv[k + d + 1]) / denom2) * coxDeBoor(k + 1, d - 1, u,kv))

    return term1 + term2



# order is degree + 1
def gen_spline_basis(x, order, knots_raw, intercept, i = 1):
    import numpy as np 

    knots = np.concatenate((np.repeat(knots_raw[0], order), knots_raw,
        np.repeat(knots_raw[-1], order)))

    MM = np.empty(shape=(len(x), order + len(knots_raw) - i))

    for j in range(order + len(knots_raw) - i):
        MM[:,j] = coxDeBoor(j, order, x, knots)

    if intercept:
        return MM
    else:
        return MM[:,1:]



def query_nearest_exon(gr, exons_annot):
    import pybedtools
    import pandas as pd
    import numpy as np

    # split loci into individual sites and find closest exons, all ties are included, so that strongest constraint can taken across ties
    gr_spl = gr.tile(1)
    gr_spl_with_nearest_exon = gr_spl.k_nearest(exons_annot)

    # find most constrained for all ties
    chr_pos_dist_constr = gr_spl_with_nearest_exon.as_df()[['Chromosome', 'End', 'constr', 'Distance']]

    chr_pos_dist_constr['Distance'] = np.abs(chr_pos_dist_constr['Distance'])
    chr_pos_dist_constr['Chromosome'] = chr_pos_dist_constr['Chromosome'].astype(str) # change from categorical to string type to avoid outer product blow up by groupby

    chr_pos_dist_min_constr = chr_pos_dist_constr.groupby(['Chromosome','End','Distance'], sort=False)['constr'].min().reset_index()

    dist_and_constr = chr_pos_dist_min_constr[['Distance','constr']].values
 
    # need to feed both constraint and distance to bspline function
    dist_basis = gen_spline_basis(dist_and_constr[:,0], order=3, knots_raw = [0, 1, int(2e3), int(1e5), 20083767], intercept=False)
    constr_basis = gen_spline_basis(dist_and_constr[:,1], order=3, knots_raw = [0, 0.01, 0.3, 0.99, 27.50007], intercept=False)

    r = np.repeat(range(dist_basis.shape[1]), constr_basis.shape[1])   
    c = np.tile(range(constr_basis.shape[1]), dist_basis.shape[1])

    spl_basis = dist_basis[:,c] * constr_basis[:,r]

    return spl_basis



# takes in numpy array
# returns n x 2 array with means and standard deviations
def get_normalisation (X):
    import numpy as np

    means = np.mean(X, axis=0)
    std_devs = np.std(X, axis=0)

    out = np.c_[means, std_devs]

    return out



# takes in n x 2 numpy array with means in 1st column, sd in 2nd
# returns normalised numpy array
def normalise_data (X, scaling):
    import numpy as np

    if X.dtype != 'float':
        X = X.astype('float')
    
    for i in range(X.shape[1]):
        X[:,i] = (X[:,i] - scaling[i,0])/ scaling[i,1]

    const_cols = (scaling[:,1] == 0)

    return X[:, ~const_cols]



def evaluate(gr, reliable_sites, reference_fasta, flank, model_bal, model_cal, pop_size, snvs, codata_args, codataNormalisation=None, split=False):
    import numpy as np

    gr_mgd = gr.merge()

    gr_reliable_mgd = gr_mgd.intersect(reliable_sites)

    obs = scars_assess.get_observed(gr_reliable_mgd, snvs)

    if len(codata_args) > 0: 
        assert codataNormalisation is not None, "normalisation required when providing arguments to read as codata"

        codata = get_codata(gr_reliable_mgd, *codata_args)
        codata = scars_codata.normalise_data(codata, codataNormalisation)
    else:
        codata = np.empty((len(obs),0))

    seqs = scars_queries.query_sequence(gr_reliable_mgd, reference_fasta, flank)
    scars_queries.correct_refs(gr_reliable_mgd, snvs, seqs)    
    scars_fit.anchor_on_AC(seqs)
    
    exp = scars_assess.get_expected(seqs, codata, model_bal, model_cal, pop_size)

    if not split:
        return [sum(obs), sum(exp), sum(obs)/sum(exp)]
    else:
        return scars_assess.match_scores_to_ranges(gr, gr_reliable_mgd, obs, exp)





import pyranges as pr
import pybedtools
from keras.models import load_model
import keras_genomics
import joblib
import os
import numpy as np
import pandas as pd
from scars import *
from scipy.stats import median_abs_deviation

import sys


flank = 6
sampleSize = int(1e5)
smoothingWindowFlank = 250
pop_size = 32399

chr_list = sorted(["chr" + str(i) for i in np.arange(1, 23, 1)])
chr_list.append('chrX')
chr_list.append('chrY')

reference_fasta = "/rds/project/who1000-1/rds-who1000-cbrc/ref/UCSC/hg38/hg38.fa"
genome = pybedtools.genome_registry.hg38

# missing values are imputed using the average for fold-change bigWig and using 0 for count data (eg DNase)
epidata_folder = "/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/epigenomic_tracks/"
epidata_files = os.listdir(epidata_folder)
non_methylation_bw_files = sorted([epidata_folder + f for f in epidata_files if "WGBS" not in f])
human_ovary_impute_vals_raw = list(map(lambda x: scars_queries.average_bw (x, chr_list, genome), non_methylation_bw_files))
human_ovary_impute_vals = [0 if "DNase" in fn else x for x,fn in zip(human_ovary_impute_vals_raw, non_methylation_bw_files)]
CpG_methylation_files = [epidata_folder + f for f in epidata_files if "WGBS" in f]


comSnvsFn = "/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/scars_pipeline_gnomad_hg38/nfe/variants/common_snvs_pass.bed"
rareSnvsFn = "/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/scars_pipeline_gnomad_hg38/nfe/variants/rare_snvs_corrected.bed"

rare_snvs = load_snvs(rareSnvsFn)
rare_snvs = rare_snvs[rare_snvs.ac / rare_snvs.an < 1.5/(2*pop_size)] # select probable singletons when AN is not maximum
common_snvs = load_snvs(comSnvsFn)


training_loci = pr.read_bed("/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/scars_pipeline_gnomad_hg38/nfe/training_sites/training_loci_N1e6_phyloP_min0p1_to_0p1.bed")
exons_annot = loadExonsAsPr("/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/scars_pipeline_gnomad_hg38/nfe/codata/exons_ensembl_annotated.bed")


# train model
sequence = scars_queries.query_sequence(training_loci, reference_fasta, flank)
scars_queries.correct_refs(training_loci, rare_snvs, sequence)
codata = scars_codata.get_codata(training_loci, non_methylation_bw_files, human_ovary_impute_vals, CpG_methylation_files, exons_annot)
indices, output, multiplier_split = scars_queries.split_data(training_loci, rare_snvs, sequence, pop_size)

balanced_multiplier = multiplier_split[:,0]
X_bal, codata_bal, y_bal = scars_fit.prep_data(indices, output, balanced_multiplier, sequence, codata, expand=True)
scars_fit.anchor_on_AC(X_bal, y_bal)
normalisation_values = scars_codata.get_normalisation(codata_bal)
codata_bal = scars_codata.normalise_data(codata_bal, normalisation_values)

epiColumns = np.arange(0, 12, 1)
methColumns = np.arange(12, 20, 1)
exonColumns = np.arange(20, codata.shape[1], 1)

model_bal_epi = scars_fit.fit_cnn(X_bal, y_bal, codata_bal[:, epiColumns])
model_bal_meth = scars_fit.fit_cnn(X_bal, y_bal, codata_bal[:, methColumns])
model_bal_exon = scars_fit.fit_cnn(X_bal, y_bal, codata_bal[:, exonColumns])

calibration_multiplier = multiplier_split[:,1]
X_cal, codata_cal, y_cal, weight_cal = scars_fit.prep_data(indices, output, calibration_multiplier, sequence, codata)
scars_fit.anchor_on_AC(X_cal, y_cal)
codata_cal = scars_codata.normalise_data(codata_cal, normalisation_values)

model_cal_epi = scars_fit.fit_calibration(model_bal_epi, X_cal, y_cal, weight_cal, flank, codata_cal[:, epiColumns])
model_cal_meth = scars_fit.fit_calibration(model_bal_meth, X_cal, y_cal, weight_cal, flank, codata_cal[:, methColumns])
model_cal_exon = scars_fit.fit_calibration(model_bal_exon, X_cal, y_cal, weight_cal, flank, codata_cal[:, exonColumns])


# get percentiles

orion_cdts_regions = pr.read_bed("/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/competing_tracks/regions_CDTS_Orion.bed")

reliable_sites = pr.read_bed("/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/scars_pipeline_gnomad_hg38/nfe/quality_filtering/reliable_sites_incl_common.bed")
reliableOverlapWithCompetingScores = reliable_sites.intersect(orion_cdts_regions)
sampledLoci = scars_queries.sample_loci_from_pr(reliableOverlapWithCompetingScores, sampleSize)


sampledLociExtd = sampledLoci.slack(smoothingWindowFlank)
scars_epi = evaluate(sampledLociExtd, reliable_sites, reference_fasta, flank, model_bal_epi, model_cal_epi, pop_size, common_snvs, [non_methylation_bw_files, human_ovary_impute_vals, None, None], normalisation_values[epiColumns],  split=True)
scars_meth = evaluate(sampledLociExtd, reliable_sites, reference_fasta, flank, model_bal_meth, model_cal_meth, pop_size, common_snvs, [None, None, CpG_methylation_files, None], normalisation_values[methColumns],  split=True)
scars_exon = evaluate(sampledLociExtd, reliable_sites, reference_fasta, flank, model_bal_exon, model_cal_exon, pop_size, common_snvs, [None, None, None, exons_annot], normalisation_values[exonColumns],  split=True)

cnn = load_model('/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/scars_pipeline_gnomad_hg38/nfe/cnn/cnn_singleton_input_N1e6_seqonly.h5', custom_objects={'RevCompConv1D'\
 : keras_genomics.layers.RevCompConv1D, 'DenseAfterRevcompConv1D': keras_genomics.layers.DenseAfterRevcompConv1D})
calibration_model = joblib.load("/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/scars_pipeline_gnomad_hg38/nfe/calibration/calibration_singleton_input_N1e6_seqonly.h5")

scars = scars_assess.evaluate(sampledLociExtd, reliable_sites, reference_fasta, flank, cnn, calibration_model, pop_size, common_snvs, [],  split=True)


obs_norm_epi = (scars_epi.observed - np.mean(scars_epi.observed))/median_abs_deviation(scars_epi.observed)
exp_norm_epi = (scars_epi.expected - np.mean(scars_epi.expected))/median_abs_deviation(scars_epi.expected)
percentiles_epi = np.quantile(obs_norm_epi - exp_norm_epi, np.arange(0, 1.01, 0.01))

obs_norm_meth = (scars_meth.observed - np.mean(scars_meth.observed))/median_abs_deviation(scars_meth.observed)
exp_norm_meth = (scars_meth.expected - np.mean(scars_meth.expected))/median_abs_deviation(scars_meth.expected)
percentiles_meth = np.quantile(obs_norm_meth - exp_norm_meth, np.arange(0, 1.01, 0.01))

obs_norm_exons = (scars_exon.observed - np.mean(scars_exon.observed))/median_abs_deviation(scars_exon.observed)
exp_norm_exons = (scars_exon.expected - np.mean(scars_exon.expected))/median_abs_deviation(scars_exon.expected)
percentiles_exons = np.quantile(obs_norm_exons - exp_norm_exons, np.arange(0, 1.01, 0.01))

obs_norm = (scars.observed - np.mean(scars.observed))/median_abs_deviation(scars.observed)
exp_norm = (scars.expected - np.mean(scars.expected))/median_abs_deviation(scars.expected)
percentiles = np.quantile(obs_norm - exp_norm, np.arange(0, 1.01, 0.01))



# evaluate pathos
patho_snvs = pr.read_bed("/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/scars_pipeline_gnomad_hg38/non_coding_pathogenic_variants/annotated/noncoding_pathogenic_HGMD_ClinVar_hg38_annot_with_Orion_CDTS_prototypeNfeSCARS.bed")
patho_snvs.columns = ["Chromosome", "Start", "End", "Orion", "CDTS", "SCARS"]

patho_snvs_extd = patho_snvs.slack(smoothingWindowFlank)
scarsPatho_epi = evaluate(patho_snvs_extd, reliable_sites, reference_fasta, flank, model_bal_epi, model_cal_epi, pop_size, common_snvs, [non_methylation_bw_files, human_ovary_impute_vals, None, None], normalisation_values[epiColumns], split=True)
scarsPatho_meth = evaluate(patho_snvs_extd, reliable_sites, reference_fasta, flank, model_bal_meth, model_cal_meth, pop_size, common_snvs, [None, None, CpG_methylation_files, None], normalisation_values[methColumns], split=True)
scarsPatho_exon = evaluate(patho_snvs_extd, reliable_sites, reference_fasta, flank, model_bal_exon, model_cal_exon, pop_size, common_snvs, [None, None, None, exons_annot], normalisation_values[exonColumns], split=True)
scarsPatho = scars_assess.evaluate(patho_snvs_extd, reliable_sites, reference_fasta, flank, cnn, calibration_model, pop_size, common_snvs, [], split=True)


obs_patho_norm_epi = (scarsPatho_epi.observed - np.mean(scars_epi.observed))/median_abs_deviation(scars_epi.observed)
exp_patho_norm_epi = (scarsPatho_epi.expected - np.mean(scars_epi.expected))/median_abs_deviation(scars_epi.expected)

obs_patho_norm_meth = (scarsPatho_meth.observed - np.mean(scars_meth.observed))/median_abs_deviation(scars_meth.observed)
exp_patho_norm_meth = (scarsPatho_meth.expected - np.mean(scars_meth.expected))/median_abs_deviation(scars_meth.expected)

obs_patho_norm_exon = (scarsPatho_exon.observed - np.mean(scars_exon.observed))/median_abs_deviation(scars_exon.observed)
exp_patho_norm_exon = (scarsPatho_exon.expected - np.mean(scars_exon.expected))/median_abs_deviation(scars_exon.expected)

obs_patho_norm = (scarsPatho.observed - np.mean(scars.observed))/median_abs_deviation(scars.observed)
exp_patho_norm = (scarsPatho.expected - np.mean(scars.expected))/median_abs_deviation(scars.expected)



pathoSnvsAnnot = pr.PyRanges(patho_snvs.as_df().loc[scarsPatho.index])
percentilesPatho_epi = toPercentile(obs_patho_norm_epi - exp_patho_norm_epi, percentiles_epi)
percentilesPatho_meth = toPercentile(obs_patho_norm_meth - exp_patho_norm_meth, percentiles_meth)
percentilesPatho_exon = toPercentile(obs_patho_norm_exon - exp_patho_norm_exon, percentiles_exons)
percentilesPatho = toPercentile(obs_patho_norm - exp_patho_norm, percentiles)


fig, axs = plt.subplots(1, 2, figsize=(10,5))

percScoreToCumulativePercCountPlot(percentilesPatho_epi, axs[0])
percScoreToCumulativePercCountPlot(percentilesPatho_meth, axs[0])
percScoreToCumulativePercCountPlot(percentilesPatho_exon, axs[0])
percScoreToCumulativePercCountPlot(percentilesPatho, axs[0])

percScoreToCumulativePercCountPlot(percentilesPatho_epi, axs[1], 5)
percScoreToCumulativePercCountPlot(percentilesPatho_meth, axs[1], 5)
percScoreToCumulativePercCountPlot(percentilesPatho_exon, axs[1], 5)
percScoreToCumulativePercCountPlot(percentilesPatho, axs[1], 5)

axs[0].legend(['Non-methylation epigenetic data', 'Methylation status', 'Exon correction', 'None'])
axs[1].legend(['Non-methylation epigenetic data', 'Methylation status', 'Exon correction', 'None'])

comparePlotFile = "/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/scars_pipeline_gnomad_hg38/plots/compareCodataByCategory.pdf"
fig.savefig(comparePlotFile)




