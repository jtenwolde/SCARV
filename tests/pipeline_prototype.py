from scars import *

import pybedtools
import numpy as np
import pandas as pd
import os


################## SET PARAMETERS ##################

# set chromosomes under consideration, ethnicity, population size and number of flanking nucleotides for model covariates
chr_list = sorted(["chr" + str(i) for i in np.arange(1,23,1)]) 
ethnicity = "nfe"   #afr
pop_size = 32399    #21042
flank = 6


################## SET FILENAMES FOR VARIANTS, QUALITY SCORES AND REFERENCE GENOME ##################

# set directory, contains required files, but limited to include only chr1:0-1,000,000
prototyping_dir = "/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/prototyping_scars_files_first_Mb_chr1/"

# vcf and coverage tsv downloaded from gnomad 
gnomad_vcf = prototyping_dir + "first_Mb.vcf.bgz"
gnomad_coverage_tsv = prototyping_dir + "first_Mb_coverage.tsv.bgz"

# phyloP from reference on, minOPR calculated myself and lifted over hg38, ensembl ftp to access annotations
phyloP_bw = prototyping_dir + "first_MB_hg38.phyloP100way.bw"
minOPR_file = "/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/opr/minOPR_gt0.99_hg38_liftedOver_sorted.bed.gz"
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
common_vars_query = {'variant_type': None, 'ethnicity': ethnicity, 'maf_range': [0.001, 1], 'PASS': None}

queries = [rare_snvs_query, fail_snvs_query, common_vars_query]
data_dict = scars_queries.query_vcf(gnomad_vcf, queries)

rare_snvs = scars_queries.convert_to_bed(data_dict[0])
fail_snvs = scars_queries.convert_to_bed(data_dict[1])
common_vars = scars_queries.convert_to_bed(data_dict[2])


################## QUERY RELIABLE SITES AND TRAINING LOCI ##################

# save the reliable sites (covered, high minOPR, no overlap with common/FAIL variants) as bigWig for quick access
# subsequently query training loci (reliable, neutral phyloP, non-exonic)
reliable_sites_bw_fn = "/home/jwt44/reliable_sites.bw"
reliable_sites = scars_queries.get_reliable_sites(minOPR_file, gnomad_coverage_tsv, chr_list, genome, fail_snvs, common_vars, bw_file=reliable_sites_bw_fn)
training_loci = scars_queries.get_training_loci(reliable_sites_bw_fn, phyloP_bw, ensembl_ftp, chr_list, genome)


################## SPLIT TRAINING LOCI IN TWO, FOR TRAINING SEQ ONLY MODEL AND FULL MODEL ##################

# sequence only model will be used to assess exon constraint, which will provide an additional covariate in the full model
N = len(training_loci)
n = N // 2

split_ix = sorted(np.random.choice(range(N), n, replace=False))

seqonly_training_loci = training_loci.at(split_ix)
full_training_loci = training_loci.at(sorted(set(range(N)) - set(split_ix)))


################## FIT AND ASSESS MODEL WITH NO CODATA ##################

# fit model without codata to assess exon constraint
sequence = scars_queries.query_sequence(seqonly_training_loci, reference_fasta, flank, genome)
indices, output, multiplier_split = scars_queries.split_data(seqonly_training_loci, rare_snvs, genome, reference_fasta, pop_size)

balanced_multiplier = multiplier_split[:,0]
X_bal, y_bal = scars_fit.prep_data(indices, output, balanced_multiplier, sequence, expand=True)
model_bal = scars_fit.fit_cnn_seq_only(X_bal, y_bal, flank)

calibration_multiplier = multiplier_split[:,1]
X_cal, y_cal, weight_cal = scars_fit.prep_data(indices, output, calibration_multiplier, sequence)
model_cal = scars_fit.fit_calibration(model_bal, X_cal, y_cal, weight_cal, flank)

test_multiplier = multiplier_split[:,2]
X_test, y_test, weight_test = scars_fit.prep_data(indices, output, test_multiplier, sequence)

preds_bal = model_bal.predict(X_test)
preds_cal = scars_assess.predict_and_calibrate_platt_iso(model_cal, np.log(preds_bal), X_test) 

n_bins_ece = 20
auc = scars_diagnostics.hand_till_auc(preds_cal, y_test, weight_test)
ece = scars_diagnostics.expected_calibration_error(n_bins_ece, preds_cal, y_test, weight_test, plot=True)


################## QUERY AND EVALUATE CONSTRAINT FOR EXONS ##################

# impose expected > 5 to avoid extreme o/e values (often 0 or >>1)
exons = scars_queries.query_ensembl(ensembl_ftp, "exon", chr_list)
exon_constr = scars_assess.evaluate(exons, reliable_sites_bw_fn, reference_fasta, flank, genome, model_bal, model_cal, pop_size, rare_snvs, split=True)
exon_constr = exon_constr.loc[exon_constr['expected'] > 5]

# set up dictionary {exon_ensembl_id: o/e}
# exons with no score (eg not covered, too small) are removed from the analysis 
exon_ids = [r.attrs['exon_id'] for r in exons.at(exon_constr.index)]
exon_constr_dict = {exon_id: constr for exon_id, constr in zip(exon_ids, exon_constr['obs_exp'])}
exons_annot = exons.each(lambda x: scars_codata.annot_with_constr(x, exon_constr_dict, '-1'))\
                   .filter(lambda x: x.attrs['constr']!='-1')\
                   .saveas()


################## FIT AND ASSESS FULL MODEL ##################

sequence = scars_queries.query_sequence(full_training_loci, reference_fasta, flank, genome)
codata = scars_codata.get_codata(full_training_loci, non_methylation_bw_files, human_ovary_impute_vals, CpG_methylation_files, exons_annot, chr_list)
indices, output, multiplier_split = scars_queries.split_data(full_training_loci, rare_snvs, genome, reference_fasta, pop_size)

balanced_multiplier = multiplier_split[:,0]
X_bal, codata_bal, y_bal = scars_fit.prep_data(indices, output, balanced_multiplier, sequence, codata, expand=True)
model_bal_full = scars_fit.fit_cnn(X_bal, y_bal, flank, codata_bal)

calibration_multiplier = multiplier_split[:,1]
X_cal, codata_cal, y_cal, weight_cal = scars_fit.prep_data(indices, output, calibration_multiplier, sequence, codata)
model_cal_full = scars_fit.fit_calibration(model_bal_full, X_cal, y_cal, weight_cal, flank, codata)

test_multiplier = multiplier_split[:,2]
X_test, codata_test, y_test, weight_test = scars_fit.prep_data(indices, output, test_multiplier, sequence, codata)

preds_bal = model_bal_full.predict([X_test, codata_test])
preds_cal = scars_assess.predict_and_calibrate_platt_iso(model_cal_full, np.log(preds_bal), X_test) 

n_bins_ece = 20
auc_full = scars_diagnostics.hand_till_auc(preds_cal, y_test, weight_test)
ece_full = scars_diagnostics.expected_calibration_error(n_bins_ece, preds_cal, y_test, weight_test, plot=True)


################## CREATE GENOME-WIDE VECTOR WITH OBS_EXP VALUES ##################

resolution=10
window_size=550
out_bw_fn = prototyping_dir + "testy_test.bw"
scars_assess.genome_wide_score_to_bigWig(genome, out_bw_fn, rare_snvs, reliable_sites_bw_fn, reference_fasta,\
    flank, model_bal_full, model_cal_full, pop_size, window_size, resolution, chr_list, *codata_args)


################## QUERY AND ASSESS NON-CODING PATHOGENIC VARIANTS ##################

# parameters to query non-coding pathogenic variants from hgmd and clinvar using an exon_flank to exclude splice sites
exon_flank = 10
hgmd_vcf = "/rds/project/who1000-1/rds-who1000-cbrc/ref/HGMD/20191111/hgmd_pro_2019.3_hg38.vcf.gz"
clinvar_vcf = "/rds/project/who1000-1/rds-who1000-cbrc/ref/ClinVar/20200615/GRCh38/clinvar_20200615.vcf.gz"

non_coding_patho_snvs = scars_noncoding.get_non_coding_pathogenic(exon_flank, chr_list, genome, ensembl_ftp, hgmd_vcf, clinvar_vcf)
non_coding_patho_scores = scars_assess.query_score_from_bigWig(non_coding_patho_snvs, out_bw_fn)





