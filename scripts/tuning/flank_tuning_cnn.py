def create_cnn(flank):
    import keras
    import keras_genomics
    import math
   
    model = keras.models.Sequential()

    model.add(keras_genomics.layers.RevCompConv1D(filters=50, kernel_size=5, input_shape=(2*flank+1, 4), 
        padding="same", activation="relu", kernel_initializer="he_normal"))

    model.add(keras_genomics.layers.RevCompConv1D(filters=100, kernel_size=model.layers[-1].output_shape[1],
        activation="relu", kernel_initializer="he_normal"))

    model.add(keras_genomics.layers.RevCompConv1D(filters=2, kernel_size=model.layers[-1].output_shape[1],
        activation="softmax", kernel_initializer="he_normal"))
    
    model.add(keras.layers.Flatten())

    return model


from scars import *

import numpy as np
import pandas as pd
import pyranges as pr
import pybedtools
import keras
from sklearn.metrics import roc_auc_score
import sys

outfile = sys.argv[1]
f = open(outfile, 'w')

pop_size = 32399

genome = pybedtools.genome_registry.hg38
reference_fasta = "/rds/project/who1000-1/rds-who1000-cbrc/ref/UCSC/hg38/hg38.fa"

pass_singletons = scars_queries.load_variants("/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/scars_pipeline_gnomad_hg38/nfe/variants/pass_singleton_snvs.bed")
pass_deletions = scars_queries.load_variants("/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/scars_pipeline_gnomad_hg38/nfe/variants/pass_deletions.bed")
training_loci = pr.read_bed("/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/scars_pipeline_gnomad_hg38/nfe/training_sites/training_loci.bed")

training_loci_subsample_size = int(1e8)
training_loci_subsample = scars_queries.sample_loci_from_pr(training_loci, training_loci_subsample_size)

max_flank = 20
sequence_max_flank = scars_queries.query_sequence(training_loci_subsample, max_flank, genome, reference_fasta)
scars_queries.correct_refs(training_loci_subsample, pass_singletons, sequence_max_flank)

for flank in range(4, 20 + 1, 4):
    sequence = sequence_max_flank[:, (max_flank-flank):(max_flank+flank+1), :]
    indices, output, multiplier_split = scars_queries.split_data(training_loci_subsample, pass_singletons, pass_deletions, sequence, pop_size)

    balanced_multiplier = multiplier_split[:,0]
    X, y = scars_fit.prep_data(indices, output, balanced_multiplier, sequence, expand=True)
    
    train_index = np.random.choice(a = np.arange(len(X)), size = len(X)//2, replace = False)
    X_train_val, y_train_val = X[train_index], y[train_index]
    X_test, y_test = np.delete(X, train_index, 0), np.delete(y, train_index, 0)

    val_index = np.random.choice(a = np.arange(X_train_val.shape[0]), size = round(0.1*X_train_val.shape[0]), replace=False)
    X_train, y_train = np.delete(X_train_val, val_index, axis=0), np.delete(y_train_val, val_index, axis=0)
    X_val, y_val = X_train_val[val_index], y_train_val[val_index]

    model_bal = create_cnn(flank)
    optim = keras.optimizers.Adam(lr = 0.0001)
    model_bal.compile(optimizer=optim,loss='categorical_crossentropy', metrics=['accuracy'])
    ival = scars_diagnostics.IntervalEvaluation(validation_data=(X_val, y_val), flank=flank, patience=5, interval=1)
    model_bal.fit(x = X_train, y = y_train, epochs=100, batch_size = 64, callbacks = [ival])

    preds = model_bal.predict(X_test)
    preds_alt = 1-np.sum(preds * X_test[:, flank, :], axis=1)
    y_alt = 1-np.sum(y_test * X_test[:, flank, :], axis=1)
    y_alt_sampled = np.random.binomial(n=1, p=preds_alt)

    auc_ht = scars_diagnostics.hand_till_auc(preds, y_test)
    max_auc_ht = scars_diagnostics.max_auc(preds)

    auc_binary = roc_auc_score(y_alt, preds_alt)
    max_auc_binary = roc_auc_score(y_alt_sampled, preds_alt)

    print('%d\t%.3f\t%.3f\t%.3f\t%.3f' % (2*flank+1, auc_ht, max_auc_ht, auc_binary, max_auc_binary), file=f)

f.close()






