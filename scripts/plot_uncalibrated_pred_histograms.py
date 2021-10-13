from scarv import *
import numpy as np
import pandas as pd
import itertools 
import joblib
from keras.models import load_model
import keras_genomics
import sys

k = 13
ancestry = "nfe"

project_folder = "/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/scars_pipeline_gnomad_hg38/" + ancestry + "/"

cnn_file = project_folder + "cnn/cnn_" + str(k) + ".h5"
cnn = load_model(cnn_file, custom_objects={'RevCompConv1D': keras_genomics.layers.RevCompConv1D})

k_mers = ['X'.join(tup) for tup in itertools.product(['A','C','G','T'], repeat = k)]
nucs = [list(seq) for seq in k_mers]
nucs_df = pd.DataFrame(nucs)

nucs_cat = nucs_df.apply(lambda x: pd.Categorical(x, categories = ['A', 'C', 'X', 'G', 'T']))
nucs_bin = pd.get_dummies(nucs_cat)
k_mers_ohe = np.array(nucs_bin, dtype=np.int8).reshape(nucs_bin.shape[0], 2*k-1, 5) 

preds_uncalibrated_snv_del = cnn.predict(k_mers_ohe)
preds_by_kmer_snv_del = pd.DataFrame(preds_uncalibrated_snv_del, index=k_mers, columns=['p_A', 'p_C', 'p_X', 'p_G', 'p_T'])


k_mers = ['X'.join(tup) for tup in itertools.product(['A','C','G','T'], repeat = k-1)]
k_mers = ["X" + kmer + "X" for kmer in k_mers]
nucs = [list(seq) for seq in k_mers]
nucs_df = pd.DataFrame(nucs)

nucs_cat = nucs_df.apply(lambda x: pd.Categorical(x, categories = ['A', 'C', 'X', 'G', 'T']))
nucs_bin = pd.get_dummies(nucs_cat)
k_mers_ohe = np.array(nucs_bin, dtype=np.int8).reshape(nucs_bin.shape[0], 2*k-1, 5) 

preds_uncalibrated_ins = cnn.predict(k_mers_ohe)
preds_by_kmer_ins = pd.DataFrame(preds_uncalibrated_ins, index=k_mers, columns=['p_A', 'p_C', 'p_X', 'p_G', 'p_T'])

# save in case of hpc problems
outFile = "/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/scars_pipeline_gnomad_hg38/nfe/preds_by_kmer_uncalibrated.pkl"
preds_by_kmer = pd.concat([preds_by_kmer_snv_del, preds_by_kmer_ins])
preds_by_kmer.to_pickle(outFile)









import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt

preds_by_kmer_file = "/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/scars_pipeline_gnomad_hg38/nfe/preds_by_kmer.pkl"
preds_by_kmer = pd.read_pickle(preds_by_kmer_file)

preds_by_kmer_snv_del = preds_by_kmer.loc[[x[0] != "X" for x in preds_by_kmer.index], :]
preds_by_kmer_ins = preds_by_kmer.loc[[x[0] == "X" for x in preds_by_kmer.index], :]


p_ins_bin = preds_by_kmer_ins[['p_A', 'p_C', 'p_G', 'p_T']].values.flatten()
p_del_bin = preds_by_kmer_snv_del.p_X.values

p_snv_A = preds_by_kmer_snv_del.loc[[x[12]!='A' for x in preds_by_kmer_snv_del.index], 'p_A'].values
p_snv_C = preds_by_kmer_snv_del.loc[[x[12]!='C' for x in preds_by_kmer_snv_del.index], 'p_C'].values
p_snv_G = preds_by_kmer_snv_del.loc[[x[12]!='G' for x in preds_by_kmer_snv_del.index], 'p_G'].values
p_snv_T = preds_by_kmer_snv_del.loc[[x[12]!='T' for x in preds_by_kmer_snv_del.index], 'p_T'].values

p_snv = np.concatenate([p_snv_A, p_snv_C, p_snv_G, p_snv_T])

##computing the bin properties (same for both distributions)
num_bin = 300000

bin_lims_snv = np.linspace(1e-8, np.max(p_snv), num_bin + 1)
bin_centers_snv = 0.5*(bin_lims_snv[:-1]+bin_lims_snv[1:])

bin_lims_ins = np.linspace(1e-14, 1e-8, num_bin + 1)
bin_centers_ins = 0.5*(bin_lims_ins[:-1]+bin_lims_ins[1:])

bin_lims_del = np.linspace(1e-10, 1e-6, num_bin + 1)
bin_centers_del = 0.5*(bin_lims_del[:-1]+bin_lims_del[1:])

##computing the histograms
hist_snv, _ = np.histogram(p_snv, bins=bin_lims_snv)
hist_ins, _ = np.histogram(p_ins_bin, bins=bin_lims_ins)
hist_del, _ = np.histogram(p_del_bin, bins=bin_lims_del)

##normalizing
hist_snv = hist_snv/np.max(hist_snv)
hist_ins = hist_ins/np.max(hist_ins)
hist_del = hist_del/np.max(hist_del)


##plotting
fig, axes = plt.subplots(1, 2, figsize=(10, 5))

axes[0].step(bin_centers_snv, hist_snv, color='red')
axes[0].set_xscale('log')

axes[1].step(bin_centers_ins, hist_ins, color='blue')
axes[1].step(bin_centers_del, hist_del, color='green')
axes[1].legend(['Insertion', 'Deletion'])
axes[1].set_xscale('log')

fig.savefig('/home/jwt44/calibrated_prob_hists.pdf')
