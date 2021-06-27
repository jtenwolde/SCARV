def create_cnn(flank):
    import keras
    import keras_genomics

    model = keras.models.Sequential()

    model.add(keras_genomics.layers.RevCompConv1D(filters=100, kernel_size=5, input_shape=(2*flank+1, 5),
        padding="same", activation="relu", kernel_initializer="he_normal"))

    model.add(keras_genomics.layers.RevCompConv1D(filters=100, kernel_size=5, padding="same",
        activation="relu", kernel_initializer="he_normal"))

    model.add(keras_genomics.layers.RevCompConv1D(filters=100, kernel_size=model.layers[-1].output_shape[1],
        activation="relu", kernel_initializer="he_normal"))

    model.add(keras_genomics.layers.RevCompConv1D(filters=3, kernel_size=model.layers[-1].output_shape[1],
                activation="softmax", kernel_initializer="he_normal"))

    model.add(keras.layers.Flatten())

    model.add(keras.layers.Lambda(sum_middle_elements))

    return model



def sum_middle_elements(x):
    from keras import backend as K

    A_and_C = x[:, 0:2,]
    X = K.sum(x[:, 2:4], axis=1, keepdims=True)
    G_and_T = x[:, 4:]
    return K.concatenate([A_and_C, X, G_and_T], axis=1)




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

ancestry = "nfe"
pop_size_male = 13596
pop_size_female = 18703

pop_split = {'XY': pop_size_male, 'XX': pop_size_female}
pop_size = pop_split['XY'] + pop_split['XX']


genome = pybedtools.genome_registry.hg38
chrXnonPAR = pr.from_dict({"Chromosome": ['chrX', 'chrX'], "Start": [0, 2781479], "End": [10001, 155701383]})
reference_fasta = "/rds/project/who1000-1/rds-who1000-cbrc/ref/UCSC/hg38/hg38.fa"

project_folder = "/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/scars_pipeline_gnomad_hg38/" + ancestry + "/"
singleton_snvs = scars_queries.load_variants(project_folder + "variants/pass_singleton_snvs.bed")
singleton_deletions = scars_queries.load_variants(project_folder + "variants/pass_singleton_1bp_deletions.bed")
singleton_snvs_and_deletions = pr.PyRanges(pd.concat([singleton_snvs.as_df(), singleton_deletions.as_df()], ignore_index=True))
singleton_insertions = scars_queries.load_variants(project_folder + "variants/pass_singleton_1bp_insertions.bed")

training_loci = pr.read_bed(project_folder + "training_sites/training_loci.bed")

training_loci_subsample_size = int(5e7)
training_loci_subsample_snvs_and_dels = scars_queries.sample_loci_from_pr(training_loci, training_loci_subsample_size)

max_flank = 20
sequence_snvs_and_dels_max_flank = scars_queries.query_sequence(training_loci_subsample_snvs_and_dels, max_flank, genome, reference_fasta)
scars_queries.correct_refs(training_loci_subsample_snvs_and_dels, singleton_snvs, sequence_snvs_and_dels_max_flank)


training_loci_insertions = training_loci[training_loci.lengths()>1]
training_loci_insertions.End -= 1
training_loci_subsample_insertions = scars_queries.sample_loci_from_pr(training_loci_insertions, training_loci_subsample_size)

training_loci_subsample_insertions.End += 1
sequence_ins_max_flank = scars_queries.query_sequence(training_loci_subsample_insertions, max_flank, genome, reference_fasta)
training_loci_subsample_insertions.End -= 1

sequence_ins_max_flank = sequence_ins_max_flank[sequence_ins_max_flank[:, 2*max_flank,2]==1]


for flank in range(2, 13, 1):
    sequence_snvs_and_dels = sequence_snvs_and_dels_max_flank[:, (2*max_flank-2*flank):(2*max_flank+2*flank+1), :]
    indices, output, multiplier_split = scars_queries.split_data(training_loci_subsample_snvs_and_dels, chrXnonPAR, singleton_snvs_and_deletions, sequence_snvs_and_dels, pop_split)
    balanced_multiplier = multiplier_split[:,0]
    X_bal_snv_del, y_bal_snv_del = scars_fit.prep_data(indices, output, balanced_multiplier, sequence_snvs_and_dels, expand=True)

    sequence_ins = sequence_ins_max_flank[:, (2*max_flank-2*flank):(2*max_flank+2*flank+1), :]
    indices_ins, output_ins, multiplier_split_ins = scars_queries.split_data(training_loci_subsample_insertions, chrXnonPAR, singleton_insertions, sequence_ins, pop_split)
    balanced_multiplier_ins = multiplier_split_ins[:,0]
    X_bal_ins, y_bal_ins = scars_fit.prep_data(indices_ins, output_ins, balanced_multiplier_ins, sequence_ins, expand=True)

    X = np.concatenate([X_bal_snv_del, X_bal_ins])
    y = np.concatenate([y_bal_snv_del, y_bal_ins])

    train_index = np.random.choice(a = np.arange(len(X)), size = len(X)//2, replace = False)
    X_train_val, y_train_val = X[train_index], y[train_index]
    X_test, y_test = np.delete(X, train_index, 0), np.delete(y, train_index, 0)

    val_index = np.random.choice(a = np.arange(X_train_val.shape[0]), size = round(0.1*X_train_val.shape[0]), replace=False)
    X_train, y_train = np.delete(X_train_val, val_index, axis=0), np.delete(y_train_val, val_index, axis=0)
    X_val, y_val = X_train_val[val_index], y_train_val[val_index]

    model_bal = create_cnn(2*flank)
    optim = keras.optimizers.Adam(lr = 0.0001)
    model_bal.compile(optimizer=optim,loss='categorical_crossentropy', metrics=['accuracy'])
    ival = scars_diagnostics.IntervalEvaluation(validation_data=(X_val, y_val), flank=flank, patience=5, interval=1)
    model_bal.fit(x = X_train, y = y_train, epochs=100, batch_size = 64, callbacks = [ival])

    preds = model_bal.predict(X_test)
    preds_alt = 1-np.sum(preds * X_test[:, 2*flank, :], axis=1)
    y_alt = 1-np.sum(y_test * X_test[:, 2*flank, :], axis=1)

    auc_ht = scars_diagnostics.hand_till_auc(preds, y_test)
    auc_binary = roc_auc_score(y_alt, preds_alt)

    print(2*flank+1, auc_ht, auc_binary, file=f)

f.close()


