from scarv import *

import numpy as np
import pandas as pd
import pyranges as pr
import pybedtools
import joblib
import sys

flank = int(sys.argv[1])
ancestry = sys.argv[2]
pop_size_male = int(sys.argv[3])
pop_size_female = int(sys.argv[4])

pop_split = {'XY': pop_size_male, 'XX': pop_size_female}

genome = pybedtools.genome_registry.hg38
chrXnonPAR = pr.from_dict({"Chromosome": ['chrX', 'chrX'], "Start": [0, 2781479], "End": [10001, 155701383]})
reference_fasta = "/rds/project/who1000-1/rds-who1000-cbrc/ref/UCSC/hg38/hg38.fa"

project_folder = "/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/scars_pipeline_gnomad_hg38/" + ancestry + "/"
singleton_snvs = scarv_queries.load_variants(project_folder + "variants/pass_singleton_snvs.bed")
singleton_deletions = scarv_queries.load_variants(project_folder + "variants/pass_singleton_1bp_deletions.bed")
singleton_snvs_and_deletions = pr.PyRanges(pd.concat([singleton_snvs.as_df(), singleton_deletions.as_df()], ignore_index=True))
singleton_insertions = scarv_queries.load_variants(project_folder + "variants/pass_singleton_1bp_insertions.bed")

training_loci = pr.read_bed(project_folder + "training_sites/training_loci.bed")

training_loci_subsample_size = int(5e7)
training_loci_subsample_snvs_and_dels = scarv_queries.sample_loci_from_pr(training_loci, training_loci_subsample_size)


sequence_snvs_and_dels = scarv_queries.query_sequence(training_loci_subsample_snvs_and_dels, flank, genome, reference_fasta)
scarv_queries.correct_refs(training_loci_subsample_snvs_and_dels, singleton_snvs, sequence_snvs_and_dels)
indices, output, multiplier_split = scarv_queries.split_data(training_loci_subsample_snvs_and_dels, chrXnonPAR, singleton_snvs_and_deletions, sequence_snvs_and_dels, pop_split)

balanced_multiplier = multiplier_split[:,0]
X_bal_snv_del, y_bal_snv_del = scarv_fit.prep_data(indices, output, balanced_multiplier, sequence_snvs_and_dels, expand=True)

calibration_multiplier = multiplier_split[:,1]
X_cal_snv_del, y_cal_snv_del, weight_cal_snv_del = scarv_fit.prep_data(indices, output, calibration_multiplier, sequence_snvs_and_dels)


training_loci_insertions = training_loci[training_loci.lengths()>1]
training_loci_insertions.End -= 1
training_loci_subsample_insertions = scarv_queries.sample_loci_from_pr(training_loci_insertions, training_loci_subsample_size)
    
training_loci_subsample_insertions.End += 1
sequence_ins = scarv_queries.query_sequence(training_loci_subsample_insertions, flank, genome, reference_fasta)
training_loci_subsample_insertions.End -= 1

sequence_ins = sequence_ins[sequence_ins[:, 2*flank,2]==1]
indices_ins, output_ins, multiplier_split_ins = scarv_queries.split_data(training_loci_subsample_insertions, chrXnonPAR, singleton_insertions, sequence_ins, pop_split)

balanced_multiplier_ins = multiplier_split_ins[:,0]
X_bal_ins, y_bal_ins = scarv_fit.prep_data(indices_ins, output_ins, balanced_multiplier_ins, sequence_ins, expand=True)

calibration_multiplier_ins = multiplier_split_ins[:,1]
X_cal_ins, y_cal_ins, weight_cal_ins = scarv_fit.prep_data(indices_ins, output_ins, calibration_multiplier_ins, sequence_ins)


X_bal = np.concatenate([X_bal_snv_del, X_bal_ins])
y_bal = np.concatenate([y_bal_snv_del, y_bal_ins])

cnn = scarv_fit.fit_cnn(X_bal, y_bal)
cnn.save(project_folder + "cnn/cnn_" + str(2*flank+1) + ".h5")

calibration_model_snvs_and_dels = scarv_fit.fit_calibration(cnn, X_cal_snv_del, y_cal_snv_del, weight_cal_snv_del)
joblib.dump(calibration_model_snvs_and_dels, project_folder + "calibration/calibration_snv_del_" + str(2*flank+1) + ".h5")

calibration_model_ins = scarv_fit.fit_calibration(cnn, X_cal_ins, y_cal_ins, weight_cal_ins)
joblib.dump(calibration_model_ins, project_folder + "calibration/calibration_ins_" + str(2*flank+1) + ".h5")

