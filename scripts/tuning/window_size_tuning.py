
def get_observed_and_expected_singletons(gr, snvs, cnn, calibration_model, reference_fasta):
    import numpy as np
    import pandas as pd
    import pyranges as pr
    from scarv import scarv_assess

    max_window_size = gr.lengths()[0]
    flank = cnn.input_shape[1]//2
    observed_singletons = np.zeros(shape=(len(gr), max_window_size), dtype=np.uint8)
    counter = 0 

    chr_list = gr.Chromosome.cat.categories
    for chrom in chr_list:
        min_ix = gr[chrom].Start.iloc[0]
        max_ix = gr[chrom].End.iloc[-1]

        snv_indicator = np.zeros(max_ix-min_ix, dtype=bool)

        snvs_in_chrom = snvs[chrom][min_ix:max_ix]
        if snvs_in_chrom.length == 0:
            continue

        snv_present = snvs[chrom][min_ix:max_ix].End
        snv_indicator[snv_present - min_ix - 1] = True

        for loci in gr[chrom].as_df().itertuples():
            observed_singletons[counter] = snv_indicator[(loci.Start - min_ix):(loci.End - min_ix)]
            counter += 1

    print('done with observed')

    expected_singletons = np.zeros(shape=(len(gr), max_window_size))
    gr_extd = gr.slack(flank)
    counter2 = 0

    for chrom in chr_list:
        min_ix = gr_extd[chrom].Start.iloc[0]
        max_ix = gr_extd[chrom].End.iloc[-1]

        chrom_gr = pr.from_dict({'Chromosome': [chrom], 'Start': [min_ix], 'End': [max_ix]})
        chrom_seq = pr.get_fasta(chrom_gr, reference_fasta)[0].upper()

        nucs_ohe = np.array(pd.get_dummies(pd.Series(list(chrom_seq)))[['A','C','G','T']])
        k_mers = np.lib.index_tricks.as_strided(nucs_ohe, shape=(len(nucs_ohe)-2*flank, 2*flank+1, 4), strides=(4,4,1))

        print('loaded sequence')

        relevant_indices = []
        for row in gr[chrom].merge().slack(flank).as_df().itertuples():
            relevant_indices.extend([*range(row.Start, row.End-2*flank)])
        relevant_indices = np.array(relevant_indices)

        relevant_k_mers = k_mers[relevant_indices - min_ix]
        relevant_preds_balanced = cnn.predict(relevant_k_mers)
        relevant_preds_calibrated = scarv_assess.calibrate(np.log(relevant_preds_balanced), calibration_model, relevant_k_mers)
        relevant_p_alts = np.sum(relevant_preds_calibrated * (1 - relevant_k_mers[:,flank,:]), axis=1)

        print('done with predicting')

        old_start = min_ix
        ix = 0
        for loci in gr_extd[chrom].as_df().itertuples():
            if loci.Start - old_start < max_window_size:
                ix += loci.Start - old_start
            else:
                ix += max_window_size

            expected_singletons[counter2] = 2 * pop_size * relevant_p_alts[ix:(ix+max_window_size)]
            counter2 += 1
            old_start = loci.Start

    return observed_singletons, expected_singletons


import pyranges as pr
from scarv import *
import numpy as np
import pandas as pd
from keras.models import load_model
import keras_genomics
import joblib
from scipy import stats
import pybedtools

ancestry = "nfe"
pop_size = 32399
max_window_size = 1001
half_ws = max_window_size//2
n_loci = int(1e5)

reference_fasta = "/rds/project/who1000-1/rds-who1000-cbrc/ref/UCSC/hg38/hg38.fa"
cnn_file = "/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/scarv_pipeline_gnomad_hg38/" + ancestry + "/cnn/cnn_25.h5"
calibration_file = "/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/scarv_pipeline_gnomad_hg38/" + ancestry + "/calibration/calibration_25.h5"

cnn = load_model(cnn_file, custom_objects={'RevCompConv1D': keras_genomics.layers.RevCompConv1D})
calibration_model = joblib.load(calibration_file)

reliable_sites = pr.read_bed("/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/scarv_pipeline_gnomad_hg38/" + ancestry + "/quality_filtering/reliable_sites.bed")
training_loci = pr.read_bed("/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/scarv_pipeline_gnomad_hg38/" + ancestry + "/training_sites/training_loci.bed")
singletons = scarv_queries.load_variants("/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/scarv_pipeline_gnomad_hg38/" + ancestry + "/variants/pass_singleton_snvs.bed")

training_loci_subsample = scarv_queries.sample_loci_from_pr(training_loci, n_loci)
training_loci_extd = training_loci_subsample.slack(half_ws)

obs_by_site, exp_by_site = get_observed_and_expected_singletons(training_loci_extd, singletons, cnn, calibration_model, reference_fasta)

medians = []
discrep_dists = []
for w in np.arange(half_ws+1):
    o = np.sum(obs_by_site[:, (half_ws-w):(half_ws+w+1)], axis=1)/(2*w+1)
    e = np.sum(exp_by_site[:, (half_ws-w):(half_ws+w+1)], axis=1)/(2*w+1)
    
    o = o[~np.isnan(e)]
    e = e[~np.isnan(e)]

    medians.append(np.median(o-e))
    discrep_dists.append(list(o-e))

