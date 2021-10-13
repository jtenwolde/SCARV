from scarv import *
import numpy as np
import pandas as pd
import itertools 
import joblib
from keras.models import load_model
import keras_genomics
import sys

k = 13
ancestry = sys.argv[1]
outFile = sys.argv[2]

project_folder = "/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/scars_pipeline_gnomad_hg38/" + ancestry + "/"
cnn_file = project_folder + "cnn/cnn_" + str(k) + ".h5"
calibration_file_snv_del = project_folder + "calibration/calibration_snv_del_" + str(k) + ".h5"
calibration_file_ins = project_folder + "calibration/calibration_ins_" + str(k) + ".h5"

cnn = load_model(cnn_file, custom_objects={'RevCompConv1D': keras_genomics.layers.RevCompConv1D})
calibration_model_snv_del = joblib.load(calibration_file_snv_del)
calibration_model_ins = joblib.load(calibration_file_ins)

k_mers = ['X'.join(tup) for tup in itertools.product(['A','C','G','T'], repeat = k)]
nucs = [list(seq) for seq in k_mers]
nucs_df = pd.DataFrame(nucs)

nucs_cat = nucs_df.apply(lambda x: pd.Categorical(x, categories = ['A', 'C', 'X', 'G', 'T']))
nucs_bin = pd.get_dummies(nucs_cat)
k_mers_ohe = np.array(nucs_bin, dtype=np.int8).reshape(nucs_bin.shape[0], 2*k-1, 5) 

preds_uncalibrated_snv_del = cnn.predict(k_mers_ohe)
preds_calibrated_snv_del = scarv_assess.calibrate(preds_uncalibrated_snv_del, calibration_model_snv_del, k_mers_ohe)
preds_by_kmer_snv_del = pd.DataFrame(preds_calibrated_snv_del, index=k_mers, columns=['p_A', 'p_C', 'p_X', 'p_G', 'p_T'])


k_mers = ['X'.join(tup) for tup in itertools.product(['A','C','G','T'], repeat = k-1)]
k_mers = ["X" + kmer + "X" for kmer in k_mers]
nucs = [list(seq) for seq in k_mers]
nucs_df = pd.DataFrame(nucs)

nucs_cat = nucs_df.apply(lambda x: pd.Categorical(x, categories = ['A', 'C', 'X', 'G', 'T']))
nucs_bin = pd.get_dummies(nucs_cat)
k_mers_ohe = np.array(nucs_bin, dtype=np.int8).reshape(nucs_bin.shape[0], 2*k-1, 5) 

preds_uncalibrated_ins = cnn.predict(k_mers_ohe)
preds_calibrated_ins = scarv_assess.calibrate(preds_uncalibrated_ins, calibration_model_ins, k_mers_ohe)
preds_by_kmer_ins = pd.DataFrame(preds_calibrated_ins, index=k_mers, columns=['p_A', 'p_C', 'p_X', 'p_G', 'p_T'])

preds_by_kmer = pd.concat([preds_by_kmer_snv_del, preds_by_kmer_ins])
preds_by_kmer.to_pickle(outFile)



