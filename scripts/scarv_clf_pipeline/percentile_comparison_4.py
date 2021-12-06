from scarv import scarv_classifier, scarv_queries, scarv_assess
from sklearn.metrics import roc_auc_score, roc_curve
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import pyranges as pr
from xgboost import XGBClassifier
import os
import math

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

scarv_fn = "/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/scarv_pipeline_gnomad_hg38/nfe/scarv_tracks/scarv_hg38_incl_1bp_indels_raw.bed.gz"
linsight_fn = "/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/classifier/covariates/LINSIGHT/LINSIGHT_hg38.bed.gz"
remm_fn = "/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/classifier/covariates/ReMM/ReMM.v0.3.1_hg38.bed.gz"
funseq_fn = "/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/classifier/covariates/funseq/hg38_NCscore_funseq216.bed.gz"

score_annotation_fns = {'LINSIGHT': linsight_fn, 'SCARV': scarv_fn,
                        'ReMM': remm_fn, 'funseq': funseq_fn}

expression_annotation_LogStdExp_fn = "/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/classifier/covariates/GTEx/expression_values_by_gene_LogStdExp.bed"
gene_annotation_fns = {'LogStdExp': expression_annotation_LogStdExp_fn}

ncER_perc_fn = "/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/competing_tracks/ncER/ncER_1bp/ncER_perc_1bp_sorted_hg38.bed.gz"
competing_score_annotation_fns = {'ncER': ncER_perc_fn}


patho_SNVs_fn = "/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/scarv_pipeline_gnomad_hg38/non_coding_pathogenic_variants/noncoding_pathogenic_HGMD_Regulatory_DM_DM?_8bpFromSplice_hg38_sorted.bed"
patho_SNVs = pr.read_bed(patho_SNVs_fn)
patho_SNVs.columns = ['Chromosome', 'Start', 'End', 'Alt']

patho_SNVs_annotated = scarv_classifier.load_data(patho_SNVs, patho_SNVs_fn, score_annotation_fns, gene_annotation_fns)
patho_SNVs_competing_scores = scarv_classifier.load_data(patho_SNVs, patho_SNVs_fn, competing_score_annotation_fns)

benign_SNVs_fn = "/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/classifier/variants/gnomAD_hg38_covered_SNVs_AFgt0p01_AFRandNFE.bed"
benign_SNVs = pr.read_bed(benign_SNVs_fn)
benign_SNVs.columns = ['Chromosome', 'Start', 'End', 'Ref', 'Alt']


ncER_training_SNVs = pr.read_bed("/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/classifier/variants/ncER_SNVs/ncER_training_SNVs_hg38.bed")[[]]

ratio = 15

matched_benign_SNVs = scarv_classifier.query_matching_variants(patho_SNVs, benign_SNVs, ratio, genome_segmentation, SpliceSites)
matched_benign_SNVs_fn = "/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/classifier/variants/matched_benign_SNVs_ratio" + str(ratio) + ".bed"
scarv_queries.writeToBed(matched_benign_SNVs, matched_benign_SNVs_fn)
matched_benign_SNVs_annotated = scarv_classifier.load_data(matched_benign_SNVs, matched_benign_SNVs_fn, score_annotation_fns, gene_annotation_fns)

ncER_train_overlap = patho_SNVs_annotated.join(ncER_training_SNVs, how='left')
ncER_train_indicator = (ncER_train_overlap.Start_b != -1).astype(np.int8).to_numpy()


patho_SNV_data = patho_SNVs_annotated.as_df().drop(['Start', 'Alt'], 1)
patho_SNV_data['pathogenic'] = 1

benign_SNV_data = matched_benign_SNVs_annotated.as_df().drop(['Start'], 1)
benign_SNV_data['pathogenic'] = 0

data = pd.concat([patho_SNV_data, benign_SNV_data], axis=0).reset_index(drop=True)

data['ncER_training'] = 0 
data.loc[np.where(ncER_train_indicator==1)[0], 'ncER_training'] = 1

data_shuffled = data.sample(frac=1)


Y = data_shuffled['pathogenic'].to_numpy()
X = data_shuffled.drop(['Chromosome', 'End', 'pathogenic'], axis=1).to_numpy()
Pos = data_shuffled[['Chromosome', 'End']].reset_index(drop=True)

random_loci_fn = "/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/competing_tracks/ncER/ncER_1bp/ncER_random_loci_1_in_30k.bed"
random_loci = pr.read_bed(random_loci_fn)
random_loci = random_loci[random_loci.Chromosome!="chrY"]
random_loci_data = scarv_classifier.load_data(random_loci, random_loci_fn, score_annotation_fns, gene_annotation_fns)
X_random = random_loci_data.as_df().drop(['Chromosome', 'Start', 'End'], 1)


k = 10
folds = np.random.choice(np.arange(k), size=X.shape[0], replace=True)
df = pd.DataFrame()

for i in range(k):
    X_train, X_test = X[folds!=i], X[folds==i]
    Y_train, Y_test = Y[folds!=i], Y[folds==i]
    Pos_test = Pos.loc[folds==i]

    clf = XGBClassifier(eval_metric='logloss', use_label_encoder=False)
    clf.fit(X_train[:,:-1], Y_train)

    preds_random = clf.predict_proba(X_random)[:,1] + np.random.normal(0, 1e-9, X_random.shape[0])
    
    X_non_ncER_training = X_test[X_test[:,-1]!=1, :-1]
    Y_non_ncER_training = Y_test[X_test[:,-1]!=1] 
    Pos_non_ncER_training = Pos_test.loc[X_test[:,-1]!=1]

    preds_test = clf.predict_proba(X_non_ncER_training[Y_non_ncER_training==1])[:,1]
    percentiles = np.quantile(preds_random, np.arange(0, 1.001, 0.001))
    percentiles_test = scarv_assess.toPercentile(preds_test, percentiles)

    df_fold = pd.concat([Pos_non_ncER_training.loc[Y_non_ncER_training==1].reset_index(drop=True), pd.DataFrame(percentiles_test)],1)
    df_fold.columns = ['Chromosome', 'End', 'SCARV_clf']
    df = df.append(df_fold)


df['Start'] = df['End'] - 1


############################################################
CADD_track = "/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/CADD_v1.6/whole_genome_SNVs.tsv.gz"
query_SNV_CADD_cmd = "tabix " + CADD_track + " -R <(awk '{gsub(\"chr\", \"\"); print $1,$3,$4}' " + patho_SNVs_fn + \
    ") | awk '{print \"chr\"$1,$2-1,$2,$3,$4,$5}' OFS='\t' - | bedmap --echo --echo-map --skip-unmapped - " + patho_SNVs_fn + \
    " | sed 's/|/\t/g' - | sed 's/;/\t/g' - | awk '$5==$10 || $5==$14 || $5==$18 {print $1,$2,$3,$6}' OFS='\t' - | sort | uniq | sort-bed - "
CADD_scores_raw = subprocess.Popen(query_SNV_CADD_cmd, shell=True, executable='/bin/bash', stdout=subprocess.PIPE).communicate()

chromosomes = [x.decode("utf-8") for x in CADD_scores_raw[0].split()[::4]]
starts = [int(x) for x in CADD_scores_raw[0].split()[1::4]]
ends = [int(x) for x in CADD_scores_raw[0].split()[2::4]]
scores = [float(x) for x in CADD_scores_raw[0].split()[3::4]]

CADD_samples_df = pd.read_csv("/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/CADD_v1.6/CADD_samples_1_in_2500.txt", header=None)
CADD_samples =  CADD_samples_df.iloc[:,0].tolist()
CADD_percentiles = np.quantile(CADD_samples, np.arange(0, 1.001, 0.001))

CADD_percentile_scores = scarv_assess.toPercentile(scores, CADD_percentiles)
CADD_percentile_scores_inv = [100-float(x) for x in CADD_percentile_scores]

CADD_df = pd.DataFrame({'Chromosome': chromosomes,
                   'Start': starts,
                   'End': ends,
                   'CADD': CADD_percentile_scores_inv})
min_CADD_by_pos = CADD_df.groupby(['Chromosome', 'Start', 'End'])['CADD'].min().reset_index()
CADD_gr = pr.PyRanges(min_CADD_by_pos)

patho_SNVs_competing_scores_full = patho_SNVs_competing_scores.join(CADD_gr, how='left').drop(['Start_b', 'End_b'])
############################################################


gr = pr.PyRanges(df)
gr = gr.join(patho_SNVs_competing_scores_full)

df_full = gr[['ncER', 'SCARV_clf', 'CADD']].as_df().drop(['Chromosome', 'Start', 'End'], 1)
df_full_noNaN = df_full.dropna()

df_full_noNaN.loc[:, 'CADD'] = (100 - df_full_noNaN['CADD'].astype(float))

df_full_noNaN.loc[:, 'ncER'] = [math.ceil(10*x)/10 for x in df_full_noNaN['ncER']]
df_full_noNaN.loc[:, 'SCARV_clf'] = [math.ceil(10*x)/10 for x in df_full_noNaN['SCARV_clf']]
df_full_noNaN.loc[:, 'CADD'] = [math.ceil(10*x)/10 for x in df_full_noNaN['CADD']]


import matplotlib.pyplot as plt

fig, axs = plt.subplots(1, 2, figsize=(10,5))

scarv_assess.percScoreToCumulativePercCountPlot(df_full_noNaN['SCARV_clf'], axs[0])
scarv_assess.percScoreToCumulativePercCountPlot(df_full_noNaN['ncER'], axs[0])
scarv_assess.percScoreToCumulativePercCountPlot(df_full_noNaN['CADD'], axs[0])

scarv_assess.percScoreToCumulativePercCountPlot(df_full_noNaN['SCARV_clf'], axs[1], 95)
scarv_assess.percScoreToCumulativePercCountPlot(df_full_noNaN['ncER'], axs[1], 95)
scarv_assess.percScoreToCumulativePercCountPlot(df_full_noNaN['CADD'], axs[1], 95)

fig.tight_layout()
axs[0].legend(['SCARV-clf', 'ncER', 'CADD'], loc="lower center")

fig.savefig(outFile) 















