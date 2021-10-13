# script to fit classification model

from scarv import scarv_classifier, scarv_queries, scarv_assess
from sklearn.metrics import roc_auc_score, roc_curve
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import pyranges as pr
from xgboost import XGBClassifier
import os

os.chdir("/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/classifier")

CTCF = pr.read_bed('regulatory_build/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20190329_CTCF_simplified.bed')
Prom = pr.read_bed('regulatory_build/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20190329_Promoter_simplified.bed')
PromFlanking = pr.read_bed('regulatory_build/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20190329_Promoter_Flanking_simplified.bed')
Enhancer = pr.read_bed('regulatory_build/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20190329_Enhancer_simplified.bed')
OpenChrom = pr.read_bed('regulatory_build/homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20190329_Open_Chromatin_Region_simplified.bed')

IntronDist = pr.read_bed('gencode_v27/gencode_v27_introns_distal_simplified.bed')
IntronCis = pr.read_bed('gencode_v27/gencode_v27_introns_cis_simplified.bed')
CDS = pr.read_bed('gencode_v27/gencode_v27_CDS_simplified.bed')
UTR = pr.read_bed('gencode_v27/gencode_v27_UTR_simplified.bed')
ncRNA = pr.read_bed('gencode_v27/gencode_v27_ncRNAs_simplified.bed')

SpliceSites = pr.read_bed('gencode_v27/gencode_splice_sites_hg38_v27.bed')

# declare in increasing order of priority
genome_segmentation = {'CTCF': CTCF, 'Prom': Prom, 'PromFlanking': PromFlanking,
                        'Enhancer': Enhancer, 'OpenChrom': OpenChrom, 'IntronDist': IntronDist,
                        'UTR': UTR, 'ncRNA': ncRNA, 'IntronCis': IntronCis, 'CDS': CDS}

scars_fn = "/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/scars_pipeline_gnomad_hg38/nfe/scars_tracks/scars_hg38_incl_1bp_indels_raw.bed.gz"
linsight_fn = "/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/classifier/covariates/LINSIGHT/LINSIGHT_hg38.bed.gz"
remm_fn = "/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/classifier/covariates/ReMM/ReMM.v0.3.1_hg38.bed.gz"
funseq_fn = "/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/classifier/covariates/funseq/hg38_NCscore_funseq216.bed.gz"
expression_annotation_LogStdExp_fn = "/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/classifier/covariates/GTEx/expression_values_by_gene_LogStdExp.bed"

score_annotation_fns = {'LINSIGHT': linsight_fn, 'SCARS': scars_fn, 'ReMM': remm_fn, 'funseq': funseq_fn}
gene_annotation_fns = {'LogStdExp': expression_annotation_LogStdExp_fn}


patho_SNVs_fn = "/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/scars_pipeline_gnomad_hg38/non_coding_pathogenic_variants/noncoding_pathogenic_HGMD_Regulatory_DM_DM?_8bpFromSplice_hg38_sorted.bed"
patho_SNVs = pr.read_bed(patho_SNVs_fn)
patho_SNVs.columns = ['Chromosome', 'Start', 'End', 'Alt']

patho_SNVs_annotated = scarv_classifier.load_data(patho_SNVs, patho_SNVs_fn, score_annotation_fns, gene_annotation_fns)

benign_SNVs_fn = "/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/classifier/variants/gnomAD_hg38_covered_SNVs_AFgt0p01_AFRandNFE.bed"
benign_SNVs = pr.read_bed(benign_SNVs_fn)
benign_SNVs.columns = ['Chromosome', 'Start', 'End', 'Ref', 'Alt']


ratio = 15

matched_benign_SNVs = scarv_classifier.query_matching_variants(patho_SNVs, benign_SNVs, ratio, genome_segmentation, SpliceSites)
matched_benign_SNVs_fn = "/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/classifier/variants/matched_benign_SNVs_ratio" + str(ratio) + ".bed"
scarv_queries.writeToBed(matched_benign_SNVs, matched_benign_SNVs_fn)
matched_benign_SNVs_annotated = scarv_classifier.load_data(matched_benign_SNVs, matched_benign_SNVs_fn, score_annotation_fns, gene_annotation_fns)

patho_SNV_data = patho_SNVs_annotated.as_df().drop(['Chromosome', 'Start', 'End', 'Alt'], 1)
patho_SNV_data['pathogenic'] = 1

benign_SNV_data = matched_benign_SNVs_annotated.as_df().drop(['Chromosome', 'Start', 'End'], 1)
benign_SNV_data['pathogenic'] = 0

data = pd.concat([patho_SNV_data, benign_SNV_data], axis=0).reset_index(drop=True)
data_shuffled = data.sample(frac=1)


Y = data_shuffled['pathogenic'].to_numpy()
X = data_shuffled[['LINSIGHT', 'SCARS', 'ReMM', 'funseq', 'LogStdExp']].to_numpy()


clf = XGBClassifier(eval_metric='logloss', use_label_encoder=False)
clf.fit(X, Y)
clf_fn = "/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/classifier/models/scars_clf.model"
clf.save_model(clf_fn)

