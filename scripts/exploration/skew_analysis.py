# script to check the possibility of using ewens sampling fmla

def evaluate(gr, reliable_sites, snvs):
    import numpy as np
    import pandas as pd

    reliable_snvs = snvs.intersect(reliable_sites)

    index = pd.Series(range(len(gr)), name="id")
    gr_id = gr.insert(index)

    hits = gr_id.join(reliable_snvs)

    K = hits.as_df().groupby('id')['ac'].count()
    K.name = 'K'

    entropy = hits.as_df().groupby('id')['ac'].agg(getSkewEntropy, pop_size)
    entropy.name = 'entropy'

    homozygosity = hits.as_df().groupby('id')['ac'].agg(getSkewHomozygosity)
    homozygosity.name = 'homozygosity'

    return pd.concat([K, entropy, homozygosity], axis=1)


def getSkewEntropy (ACs, pop_size):
    pk = np.array(ACs/(2*pop_size))
    term1 = -np.nansum(pk*np.log2(pk)) 
    term2 = -np.nansum((1-pk)*np.log2(1-pk))
    return (term1+term2)/len(ACs)


def getSkewHomozygosity (ACs):
    return sum(ACs**2)/(sum(ACs)**2)


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

    percentiles = pd.cut(scores, bins=percentiles, labels=np.arange(1,1001,1), include_lowest=True)
    return percentiles


def getRelFreqHist(hits):
    import collections
    import numpy as np

    labels = ['1', '2', '3', '4-5', '6-10', '11-50', '51-32399']
    binned = pd.cut(hits.ac, [0,1,2,3,5,10,50,pop_size], labels=labels)
    counts = collections.Counter(binned)
    return np.array([counts[lab]/len(binned) for lab in labels])


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

sampledLociExtd = sampledLoci.slack(smoothingWindowFlank)
skewScores = evaluate(sampledLociExtd, reliable_sites, snvs)

percentileEntropy = np.quantile(skewScores.entropy, np.arange(0, 1.001, 0.001))
percentileHomozygosity = np.quantile(skewScores.homozygosity, np.arange(0, 1.001, 0.001))

patho_snvs = pr.read_bed("/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/non_coding_pathogenic_variants/noncoding_pathogenic_annot_by_Orion_and_CDTS_including_scores.bed")
patho_snvs.columns = ["Chromosome", "Start", "End", "CDTS", "Orion"]

patho_snvs_extd = patho_snvs.slack(smoothingWindowFlank)
skewScoresPatho = evaluate(patho_snvs_extd, reliable_sites, snvs)

percentilesPathoEntropy = toPercentile(skewScoresPatho.entropy, percentileEntropy)
percentilesPathoHomozygosity = toPercentile(skewScoresPatho.homozygosity, percentileHomozygosity)

np.savetxt(outFnEntropy, percentilesPathoEntropy.astype(int).to_numpy(), fmt='%1.1d')
np.savetxt(outFnHomozygosity, percentilesPathoHomozygosity.astype(int).to_numpy(), fmt='%1.1d')




import matplotlib.pyplot as plt

hitsPatho = patho_snvs_extd.join(snvs.intersect(reliable_sites))
hits = sampledLociExtd.join(snvs.intersect(reliable_sites))

arr1 = getRelFreqHist(hitsPatho)
arr2 = getRelFreqHist(hits)

nReps=1000
nPathos = 1166
M = np.empty((nReps, len(labels)))
for i in range(nReps):
    gr_narrow = pr.PyRanges(patho_snvs.as_df().iloc[np.random.choice(range(nPathos),nPathos)])
    gr = gr_narrow.slack(smoothingWindowFlank)
    index = pd.Series(range(len(gr)), name="id")
    gr_id = gr.insert(index)
    hits = gr_id.join(reliable_snvs_near_snvs)
    arr1 = getRelFreqHist(hits)
    M[i,:] = arr1/arr2


fig, axs = plt.subplots(2, 2, figsize=(15, 15))

axs[0,0].bar(labels, arr1)
axs[0,1].bar(labels, arr2)

axs[0,0].set_xlabel("AC bin", fontsize=12)
axs[0,1].set_xlabel("AC bin", fontsize=12)

axs[0,0].set_ylabel("Relative frequency pathogenic non-coding", fontsize=12)
axs[0,1].set_ylabel("Relative frequency random sites", fontsize=12)

axs[1,0].bar(labels, arr1/arr2)
axs[1,0].set_ylim(0.9 * np.min(np.median(M, axis=0)), 1.1 * np.max(np.median(M, axis=0)))
axs[1,0].set_xlabel("AC bin", fontsize=12)
axs[1,0].set_ylabel("Median ratio of relative frequency", fontsize=12)
axs[1,0].axhline(1, color='k')

axs[1,1].boxplot([skewScores.entropy,skewScoresPatho.entropy],positions=[1,2], labels=['Random sites', 'Non-coding pathogenic'])
axs[1,1].set_ylim(0, 0.1)
axs[1,1].set_ylabel("Shannon entropy", fontsize=12)

fig.savefig("/home/jwt44/skewAssessment.pdf")











