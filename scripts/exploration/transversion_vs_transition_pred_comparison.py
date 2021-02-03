# script to verify relative mutation probabilities transversion vs. transition
import numpy as np
from keras.models import load_model
import keras_genomics
import matplotlib.pyplot as plt
import joblib

# generate random sequences
n_samples = int(1e5)
seqs = np.random.multinomial(1, [1/4.]*4, size=33 * n_samples).reshape(n_samples, 33, 4)

# predict
cnn_file = "/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/scars_pipeline_gnomad_hg38/nfe/cnn/cnn_33.h5"
cnn = load_model(cnn_file, custom_objects={'RevCompConv1D': keras_genomics.layers.RevCompConv1D})

preds = cnn.predict(seqs)

# get p_ts and p_tv
flank = seqs.shape[1]//2
p_ref = np.sum(preds * seqs[:, flank, :], axis=1)
p_ts = np.sum(preds * np.roll(seqs[:, flank, :], 2), axis=1)

p_tv_expanded = (preds * (1 - (np.roll(seqs[:, flank, :], 2) + seqs[:, flank, :]))).reshape(4*n_samples)
p_tv = p_tv_expanded[p_tv_expanded > 0]

# boxplot
data = [p_tv, p_ts]

fig, ax = plt.subplots()
ax.boxplot(data, showmeans=True)

ax.set_xticklabels(['Transversion', 'Transition'])

ax.set_xlabel('Type of SNV')
ax.set_ylabel('Predicted probability')

fig.savefig("/home/jwt44/boxplot_preds_transition_vs_transversion.pdf")
