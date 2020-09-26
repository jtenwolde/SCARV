from scars import *

import numpy as np
import pandas as pd
import pyranges as pr
import pybedtools
import os
import joblib


################## SET PARAMETERS ##################

# set chromosomes under consideration, ethnicity, population size and number of flanking nucleotides for model covariates
chr_list = sorted(["chr" + str(i) for i in np.arange(1, 23, 1)])
chr_list.append('chrX')
chr_list.append('chrY')

ethnicity = "nfe"   #afr
pop_size = 32399    #21042
flank = 6
training_loci_subsample_size = int(1e6)


################## SET FILENAMES FOR VARIANTS, QUALITY SCORES AND REFERENCE GENOME ##################

intermediate_file_folder = "/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/scars_pipeline_gnomad_hg38/"

# vcf and coverage tsv downloaded from gnomad 
gnomad_vcf = "/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/gnomad_v3_hg38/gnomad.genomes.r3.0.sites.vcf.bgz"
gnomad_coverage_tsv = "/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/gnomad_v3_hg38/gnomad.genomes.r3.0.coverage.summary.tsv.bgz"

# phyloP from reference on, minOPR calculated myself and lifted over hg38, ensembl ftp to access annotations
phyloP_bw = "/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/hg38.phyloP100way.bw"
ensembl_ftp = "ftp://ftp.ensembl.org/pub/release-100/gtf/homo_sapiens/Homo_sapiens.GRCh38.100.gtf.gz"

# set genome for knowledge on chromosome size, set fasta filename for sequence queries
genome = pybedtools.genome_registry.hg38
reference_fasta = "/rds/project/who1000-1/rds-who1000-cbrc/ref/UCSC/hg38/hg38.fa"


################## SET FILENAMES AND IMPUTATION VALUES FOR EPIGENETIC CODATA ##################

# filenames to query codata for use as model covariates
# missing values are imputed using the average for fold-change bigWig and using 0 for count data (eg DNase)
epidata_folder = "/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/epigenomic_tracks/"
epidata_files = os.listdir(epidata_folder)
non_methylation_bw_files = sorted([epidata_folder + f for f in epidata_files if "WGBS" not in f])
human_ovary_impute_vals_raw = list(map(lambda x: scars_queries.average_bw (x, chr_list, genome), non_methylation_bw_files))
human_ovary_impute_vals = [0 if "DNase" in fn else x for x,fn in zip(human_ovary_impute_vals_raw, non_methylation_bw_files)]
CpG_methylation_files = [epidata_folder + f for f in epidata_files if "WGBS" in f]


################## QUERY VARIANTS FROM VCF FILE ##################

# query required variants, multiple queries in one function call for efficiency
rare_snvs_query = {'variant_type': "snv", 'ethnicity': ethnicity, 'maf_range': [0, 0.001], 'PASS': True}
fail_snvs_query = {'variant_type': "snv", 'ethnicity': ethnicity, 'maf_range': [0, 1], 'PASS': False}
common_vars_query = {'variant_type': None, 'ethnicity': ethnicity, 'maf_range': [0.001, 0.999], 'PASS': None}

queries = [rare_snvs_query, fail_snvs_query, common_vars_query]
data_dict = scars_queries.query_vcf(gnomad_vcf, queries)

rare_snvs = coordlist_to_pyranges(data_dict[0], entryNames=["Chromosome", "Start", "End", "ref", "alt", "ac", "an"])
fail_snvs = coordlist_to_pyranges(data_dict[1], entryNames=["Chromosome", "Start", "End", "ref", "alt", "ac", "an"])
common_vars = coordlist_to_pyranges(data_dict[2], entryNames=["Chromosome", "Start", "End", "ref", "alt", "ac", "an"])


################## QUERY RELIABLE SITES AND TRAINING LOCI ##################

reliable_sites = scars_queries.get_reliable_sites(gnomad_coverage_tsv, genome, fail_snvs, common_vars)
training_loci = scars_queries.get_training_loci(reliable_sites, phyloP_bw, ensembl_ftp)


################## SPLIT TRAINING LOCI IN TWO, FOR TRAINING SEQ ONLY MODEL AND FULL MODEL ##################

# sequence only model will be used to assess exon constraint, which will provide an additional covariate in the full model
seqonly_training_loci = scars_queries.sample_loci_from_pr(training_loci, training_loci_subsample_size)
remaining_loci = training_loci.subtract(seqonly_training_loci)
full_training_loci = scars_queries.sample_loci_from_pr(remaining_loci, training_loci_subsample_size)

################## FIT AND ASSESS MODEL WITH NO CODATA ##################

# fit model without codata to assess exon constraint
sequence = scars_queries.query_sequence(seqonly_training_loci, reference_fasta, flank)
scars_queries.correct_refs(seqonly_training_loci, rare_snvs, sequence)
indices, output, multiplier_split = scars_queries.split_data(seqonly_training_loci, rare_snvs, sequence, pop_size)

balanced_multiplier = multiplier_split[:,0]
X_bal, y_bal = scars_fit.prep_data(indices, output, balanced_multiplier, sequence, expand=True)
scars_fit.anchor_on_AC(X_bal, y_bal)
model_bal = scars_fit.fit_cnn(X_bal, y_bal)

calibration_multiplier = multiplier_split[:,1]
X_cal, y_cal, weight_cal = scars_fit.prep_data(indices, output, calibration_multiplier, sequence)
scars_fit.anchor_on_AC(X_cal, y_cal)
model_cal = scars_fit.fit_calibration(model_bal, X_cal, y_cal, weight_cal, flank)

test_multiplier = multiplier_split[:,2]
X_test, y_test, weight_test = scars_fit.prep_data(indices, output, test_multiplier, sequence)
scars_fit.anchor_on_AC(X_test, y_test)

preds_bal = model_bal.predict([X_test, np.empty((X_test.shape[0], 0))])
preds_cal = scars_assess.predict_and_calibrate_platt_iso(model_cal, np.log(preds_bal), X_test) 

n_bins_ece = 20
auc = scars_diagnostics.hand_till_auc(preds_cal, y_test, weight_test)
ece = scars_diagnostics.expected_calibration_error(n_bins_ece, preds_cal, y_test, weight_test, plot=True)


################## QUERY AND EVALUATE CONSTRAINT FOR EXONS ##################

# impose expected > 5 to avoid extreme o/e values (often 0 or >>1)
exons = scars_queries.query_ensembl(ensembl_ftp, "exon")
codata_args_for_exon_evaluation = []
exon_constr = scars_assess.evaluate(exons, reliable_sites, reference_fasta, flank, model_bal, model_cal, pop_size, rare_snvs, codata_args_for_exon_evaluation, split=True)
exon_constr = exon_constr.loc[exon_constr['expected'] > 5]

# set up dictionary {exon_ensembl_id: o/e}
# exons with no score (eg not covered, too small) are removed from the analysis 
exon_ids = [r.attrs['exon_id'] for r in exons.at(exon_constr.index)]
exon_constr_dict = {exon_id: constr for exon_id, constr in zip(exon_ids, exon_constr['obs_exp'])}
exons_annot = exons.each(lambda x: scars_codata.annot_with_constr(x, exon_constr_dict, '-1'))\
                   .filter(lambda x: x.attrs['constr']!='-1')\
                   .saveas(intermediate_file_folder + 'exons_annotated.bed')
 

################## FIT AND ASSESS FULL MODEL ##################

sequence = scars_queries.query_sequence(full_training_loci, reference_fasta, flank, genome)
scars_queries.correct_refs(full_training_loci, rare_snvs, sequence)
codata = scars_codata.get_codata(full_training_loci, non_methylation_bw_files, human_ovary_impute_vals, CpG_methylation_files, exons_annot, chr_list)
indices, output, multiplier_split = scars_queries.split_data(full_training_loci, rare_snvs, sequence, pop_size)

balanced_multiplier = multiplier_split[:,0]
X_bal, codata_bal, y_bal = scars_fit.prep_data(indices, output, balanced_multiplier, sequence, codata, expand=True)
scars_fit.anchor_on_AC(X_bal, y_bal)
normalisation_values = scars_codata.get_normalisation(codata_bal)
codata_bal = scars_codata.normalise_data(codata_bal, normalisation_values)
model_bal_full = scars_fit.fit_cnn(X_bal, y_bal, codata_bal)

calibration_multiplier = multiplier_split[:,1]
X_cal, codata_cal, y_cal, weight_cal = scars_fit.prep_data(indices, output, calibration_multiplier, sequence, codata)
scars_fit.anchor_on_AC(X_cal, y_cal)
codata_cal = scars_codata.normalise_data(codata_cal, normalisation_values)
model_cal_full = scars_fit.fit_calibration(model_bal_full, X_cal, y_cal, weight_cal, flank, codata_cal)

test_multiplier = multiplier_split[:,2]
X_test, codata_test, y_test, weight_test = scars_fit.prep_data(indices, output, test_multiplier, sequence, codata)
scars_fit.anchor_on_AC(X_test, y_test)
codata_test = scars_codata.normalise_data(codata_test, normalisation_values)

preds_bal = model_bal_full.predict([X_test, codata_test])
preds_cal = scars_assess.predict_and_calibrate_platt_iso(model_cal_full, np.log(preds_bal), X_test) 

n_bins_ece = 20
auc_full = scars_diagnostics.hand_till_auc(preds_cal, y_test, weight_test)
ece_full = scars_diagnostics.expected_calibration_error(n_bins_ece, preds_cal, y_test, weight_test, plot=True)

# save full models
model_bal_full.save(intermediate_file_folder + "cnn_full.h5")
joblib.dump(model_cal_full, intermediate_file_folder + "calibration_full.h5")

################## CREATE GENOME-WIDE VECTOR WITH OBS_EXP VALUES ##################

#model_bal_full = load_model(intermediate_file_folder + 'cnn_full.h5', custom_objects={'RevCompConv1D'\
# : keras_genomics.layers.RevCompConv1D, 'DenseAfterRevcompConv1D': keras_genomics.layers.DenseAfterRevcompConv1D})
#model_cal_full = joblib.load(intermediate_file_folder + "calibration_full.h5")

resolution=10
window_size=550
out_bw_fn = intermediate_file_folder + "obs_exp_sliding_window_res" + str(resolution) + "_width" + str(window_size) + "_hg38.bw"
scars_assess.genome_wide_score_to_bigWig(genome, out_bw_fn, rare_snvs, reliable_sites_bw_fn, reference_fasta,\
    flank, model_bal_full, model_cal_full, pop_size, window_size, resolution, chr_list,\
    non_methylation_bw_files, human_ovary_impute_vals, CpG_methylation_files, exons_annot, chr_list)


################## QUERY AND ASSESS NON-CODING PATHOGENIC VARIANTS ##################

# parameters to query non-coding pathogenic variants from hgmd and clinvar using an exon_flank to exclude splice sites
exon_flank = 10
hgmd_vcf = "/home/km369/.random_stuff/2019-4/hgmd_pro_2019.4_hg38.vcf"
clinvar_vcf = "/rds/project/who1000-1/rds-who1000-cbrc/ref/ClinVar/20200615/GRCh38/clinvar_20200615.vcf.gz"

non_coding_patho_snvs = scars_noncoding.get_non_coding_pathogenic(exon_flank, chr_list, genome, ensembl_ftp, hgmd_vcf, clinvar_vcf)
non_coding_patho_scores = scars_assess.query_score_from_bigWig(non_coding_patho_snvs, out_bw_fn)




