from scars import *

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
pass_singletons = scars_queries.load_variants(project_folder + "variants/pass_singleton_snvs.bed")
pass_deletions = scars_queries.load_variants(project_folder + "variants/pass_deletions.bed")
training_loci = pr.read_bed(project_folder + "training_sites/training_loci.bed")

training_loci_subsample_size = int(5e7)
training_loci_subsample = scars_queries.sample_loci_from_pr(training_loci, training_loci_subsample_size)

sequence = scars_queries.query_sequence(training_loci_subsample, flank, genome, reference_fasta)
scars_queries.correct_refs(training_loci_subsample, pass_singletons, sequence)
indices, output, multiplier_split = scars_queries.split_data(training_loci_subsample, chrXnonPAR, pass_singletons, pass_deletions, sequence, pop_split)

balanced_multiplier = multiplier_split[:,0]
X_bal, y_bal = scars_fit.prep_data(indices, output, balanced_multiplier, sequence, expand=True)
cnn = scars_fit.fit_cnn(X_bal, y_bal)
cnn.save(project_folder + "cnn/cnn_" + str(2*flank+1) + ".h5")

calibration_multiplier = multiplier_split[:,1]
X_cal, y_cal, weight_cal = scars_fit.prep_data(indices, output, calibration_multiplier, sequence)
calibration_model = scars_fit.fit_calibration(cnn, X_cal, y_cal, weight_cal)
joblib.dump(calibration_model, project_folder + "calibration/calibration_" + str(2*flank+1) + ".h5")

