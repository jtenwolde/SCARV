def simple_model (k, n):
    import keras
    import keras_genomics

    model = keras.models.Sequential()

    if n==0:
        model.add(keras_genomics.layers.RevCompConv1D(filters=2, kernel_size=k, input_shape=(k,4), activation="softmax", kernel_initializer="he_normal"))
        model.add(keras.layers.Flatten())
        return model
    else:
        model.add(keras_genomics.layers.RevCompConv1D(filters=100, kernel_size=k, input_shape=(k,4), activation="relu", kernel_initializer="he_normal")) 
        for i in range(n-1):
            model.add(keras_genomics.layers.RevCompConv1D(filters=100, kernel_size=model.layers[-1].output_shape[1], activation="relu", kernel_initializer="he_normal")) 
        model.add(keras_genomics.layers.RevCompConv1D(filters=2, kernel_size=model.layers[-1].output_shape[1], activation="softmax", kernel_initializer="he_normal"))
        model.add(keras.layers.Flatten())
        return model


from scars import *

import numpy as np
import pandas as pd
import pyranges as pr
import pybedtools
import keras
from sklearn.metrics import roc_auc_score
from sklearn.metrics import average_precision_score

outfile = sys.argv[1]
f = open(outfile, 'w')

pop_size = 32399

genome = pybedtools.genome_registry.hg38
reference_fasta = "/rds/project/who1000-1/rds-who1000-cbrc/ref/UCSC/hg38/hg38.fa"

pass_singletons = scars_queries.load_variants("/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/scars_pipeline_gnomad_hg38/nfe/variants/pass_singleton_snvs.bed")
pass_deletions = scars_queries.load_variants("/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/scars_pipeline_gnomad_hg38/nfe/variants/pass_deletions.bed")
training_loci = pr.read_bed("/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/scars_pipeline_gnomad_hg38/nfe/training_sites/training_loci.bed")

################## RUN SEARCH OVER FLANK PARAMETER #######################
training_loci_subsample_size = int(1e8)
training_loci_subsample = scars_queries.sample_loci_from_pr(training_loci, training_loci_subsample_size)

max_flank = 20
sequence_max_flank = scars_queries.query_sequence(training_loci_subsample, max_flank, genome, reference_fasta)
scars_queries.correct_refs(training_loci_subsample, pass_singletons, sequence_max_flank)

for flank in range(2, 20 + 1, 2):
    sequence = sequence_max_flank[:, (max_flank-flank):(max_flank+flank+1), :]
    indices, output, multiplier_split = scars_queries.split_data(training_loci_subsample, pass_singletons, pass_deletions, sequence, pop_size)

    balanced_multiplier = multiplier_split[:,0]
    X, y = scars_fit.prep_data(indices, output, balanced_multiplier, sequence, expand=True)
    
    index = np.random.choice(a = np.arange(len(X)), size = len(X)//2, replace = False)

    X_train, y_train = X[index], y[index]
    X_test, y_test = np.delete(X, index, 0), np.delete(y, index, 0)

    ix = np.random.choice(range(X.shape[0]), round(0.9*X.shape[0]), replace=False)

    X_val, y_val = np.delete(X, ix, axis=0), np.delete(y, ix, axis=0)
    X_train, y_train = X[ix], y[ix]

    for n in [0,3]:
        model_bal = simple_model(2*flank+1, n)
        optim = keras.optimizers.Adam(lr = 0.0001)
        model_bal.compile(optimizer=optim,loss='categorical_crossentropy', metrics=['accuracy'])
        ival = scars_diagnostics.IntervalEvaluation(validation_data=(X_val, y_val), flank=flank, patience=5, interval=1)
        model_bal.fit(x = X_train, y = y_train, epochs=100, batch_size = 64, callbacks = [ival])

        preds = model_bal.predict(X_test)
        auc_mult = scars_diagnostics.hand_till_auc(preds, y_test)

        p_mut = np.sum((1-X_test[:,flank,:]) * preds, axis=1)
        y_mut = np.sum(y_test * (1-X_test[:,flank,:]), axis=1)

        auc = roc_auc_score(y_mut, p_mut)
        ap = average_precision_score(y_mut, p_mut)
        brier = np.mean((y_mut-p_mut)**2)

        print(flank, n, auc_mult, auc, ap, brier, file=f)
        
f.close()


