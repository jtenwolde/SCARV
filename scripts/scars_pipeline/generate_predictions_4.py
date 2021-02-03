from scars import scars_assess
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
calibration_file = project_folder + "calibration/calibration_" + str(k) + ".h5"

cnn = load_model(cnn_file, custom_objects={'RevCompConv1D': keras_genomics.layers.RevCompConv1D})
calibration_model = joblib.load(calibration_file)

k_mers = [''.join(tup) for tup in itertools.product(['A','C','G','T'], repeat = k)]
nucs = [list(seq) for seq in k_mers]
nucs_df = pd.DataFrame(nucs)

nucs_cat = nucs_df.apply(lambda x: pd.Categorical(x, categories = ['A', 'C', 'G', 'T']))
nucs_bin = pd.get_dummies(nucs_cat)
k_mers_ohe = np.array(nucs_bin, dtype=np.int8).reshape(nucs_bin.shape[0], k, 4) 

preds_uncalibrated = cnn.predict(k_mers_ohe)
preds_calibrated = scars_assess.calibrate(preds_uncalibrated, calibration_model, k_mers_ohe)

preds_by_kmer = pd.DataFrame(preds_calibrated, index=k_mers, columns=['p_A', 'p_C', 'p_G', 'p_T'])
preds_by_kmer.to_pickle(outFile)
