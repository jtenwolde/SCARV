import os
import subprocess
import sys
import math
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scarv import *

outFile = sys.argv[1]

# query pathogenic SNVs
exon_flank = 8
ensembl_ftp = "ftp://ftp.ensembl.org/pub/release-100/gtf/homo_sapiens/Homo_sapiens.GRCh38.100.gtf.gz"
hgmd_vcf = "/home/et341/rds-et341/shared/2019-4/hgmd_pro_2019.4_hg38.vcf"
hgmd_annotation = "/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/scarv_pipeline_gnomad_hg38/non_coding_pathogenic_variants/hgmd_hg38.txt.gz"

hgmd_pathogenic_snv_file = "/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/scarv_pipeline_gnomad_hg38/non_coding_pathogenic_variants/noncoding_pathogenic_HGMD_Regulatory_DM_DM?_8bpFromSplice_hg38_sorted.bed"
snvs_annot_by_all_file = "/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/scarv_pipeline_gnomad_hg38/non_coding_pathogenic_variants/noncoding_pathogenic_positions_HGMD_Regulatory_DM_DM?_8bpFromSplice_hg38_sorted_annot_by_all_scores.bed"

non_coding_patho_snvs = scarv_noncoding.getNonCodingPathogenic(exon_flank, ensembl_ftp, hgmd_vcf, hgmd_annotation)
scarv_queries.writeToBed(non_coding_patho_snvs, hgmd_pathogenic_snv_file)

SCARV_track = "/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/scarv_pipeline_gnomad_hg38/nfe/scarv_tracks/scarv_hg38_incl_1bp_indels_autosomes.bed.gz"
CDTS_track = "/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/reimplementing_CDTS_gnomad_hg38/CDTS/CDTS_diff_perc_coordsorted_autosomes_gnomAD_N32299_hg38.bed.gz"
Orion_track = "/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/reimplementing_Orion_gnomad_hg38/Orion/Orion_autosomes_hg38.bed.gz"
gwRVIS_track = "/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/gwRVIS/gwrvis_single_nt.hg38.bed.gz"

SCARV_annotated_cmd = "tabix " + SCARV_track + " -R " + hgmd_pathogenic_snv_file + " | cut -f1-3" 
CDTS_annotated_cmd = "<(tabix " + CDTS_track + " -R " + hgmd_pathogenic_snv_file + ")"
Orion_annotated_cmd = "<(tabix " + Orion_track + " -R " + hgmd_pathogenic_snv_file + ")"
gwRVIS_annotated_cmd = "<(tabix " + gwRVIS_track + " -R " + hgmd_pathogenic_snv_file + ")"

annotated_SNVs_cmd = SCARV_annotated_cmd + " | bedmap --echo --skip-unmapped - " + CDTS_annotated_cmd + " | bedmap --echo --skip-unmapped - " + \
        Orion_annotated_cmd + " | bedmap --echo --skip-unmapped - " + gwRVIS_annotated_cmd + " > " + snvs_annot_by_all_file
subprocess.Popen(annotated_SNVs_cmd, shell=True, executable='/bin/bash')

SCARV_scores = list(map(float, os.popen("tabix " + SCARV_track + " -R " + snvs_annot_by_all_file).read().split()[3::4]))
CDTS_scores = list(map(float, os.popen("tabix " + CDTS_track + " -R " + snvs_annot_by_all_file).read().split()[6::7]))
Orion_scores = list(map(float, os.popen("tabix " + Orion_track + " -R " + snvs_annot_by_all_file).read().split()[3::4]))

gwRVIS_samples = np.loadtxt("/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/gwRVIS/approx_1e5k_samples.txt")
gwRVIS_raw_scores = list(map(float, os.popen("tabix " + gwRVIS_track + " -R " + snvs_annot_by_all_file).read().split()[3::4]))
gwRVIS_perc_boundaries = np.quantile(gwRVIS_samples, np.arange(0,1.001,0.001))
gwRVIS_scores = scarv_assess.toPercentile(gwRVIS_raw_scores, gwRVIS_perc_boundaries)

SCARV_scores_corrected = [100-x for x in SCARV_scores]
CDTS_scores_corrected = [100-float(x) for x in pd.cut(CDTS_scores, np.arange(0, 100.1, 0.1), labels=np.arange(0.1, 100.1, 0.1))]
Orion_scores_corrected = [100-x for x in Orion_scores]
gwRVIS_scores_corrected = [100-float(x) for x in gwRVIS_scores]

print(np.median(SCARV_scores_corrected), np.median(CDTS_scores_corrected), np.median(Orion_scores_corrected), np.median(gwRVIS_scores_corrected))


fig, axs = plt.subplots(1, 2, figsize=(10,5))

axs[0].text(-19, 463, "a", size=30)
axs[1].text(-0.7, 260, "b", size=30)

axs[0].tick_params(axis='both', which='major', labelsize=15)
axs[1].tick_params(axis='both', which='major', labelsize=15)

axs[0].grid()
axs[1].grid()

scarv_assess.percScoreToCumulativePercCountPlot([100-x for x in SCARV_scores], axs[0])
scarv_assess.percScoreToCumulativePercCountPlot([100-x for x in CDTS_scores], axs[0])
scarv_assess.percScoreToCumulativePercCountPlot([100-x for x in Orion_scores], axs[0])
scarv_assess.percScoreToCumulativePercCountPlot([100-x for x in gwRVIS_scores], axs[0])

scarv_assess.percScoreToCumulativePercCountPlot([100-x for x in SCARV_scores], axs[1], 95)
scarv_assess.percScoreToCumulativePercCountPlot([100-x for x in CDTS_scores], axs[1], 95)
scarv_assess.percScoreToCumulativePercCountPlot([100-x for x in Orion_scores], axs[1], 95)
scarv_assess.percScoreToCumulativePercCountPlot([100-x for x in gwRVIS_scores], axs[1], 95)

fig.tight_layout()
fig.subplots_adjust(bottom=0.28)
axs[0].legend(['SCARV', 'CDTS', 'Orion', 'gwRVIS'], loc="lower center", bbox_to_anchor=(1.1, -0.45), ncol=5, prop={'size': 20}, frameon=False)

fig.savefig(outFile) 
