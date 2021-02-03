def percScoreToCumulativePercCountPlot(percScore, plt=None, x_upper=None, dashed=False):
    cumulativeCounts = getCumulativeCounts(percScore)
    plt = plotCumulativePercCountPlot(cumulativeCounts, plt, x_upper, dashed)

    return plt


def getCumulativeCounts(percentile_values):
    import collections
    import numpy as np

    percentileTally = collections.Counter([int(10*x) for x in percentile_values])
    addMissingPercentiles(percentileTally)

    counts = [count for perc, count in sorted(percentileTally.items())]
    countsCumulative = np.cumsum(counts)

    return countsCumulative


def plotCumulativePercCountPlot(countsCumulative, plt=None, x_upper=None, dashed=False):
    import numpy as np

    percentiles = np.arange(0.1, 100.1, 0.1)

    if dashed:
        plt.plot(percentiles, countsCumulative, '--', linewidth=3)
    else:
        plt.plot(percentiles, countsCumulative, linewidth=3)

    plt.set_xlabel("Percentile", fontsize=16)
    plt.set_ylabel("Cumulative count", fontsize=16)

    plt.set_ylim(0, 1.1 * countsCumulative[-1])

    if x_upper is not None:
        plt.set_xlim(0, x_upper)
        plt.set_ylim(0, 1.1 * max([line.get_ydata()[x_upper * 10 - 1] for line in axs[1].lines]))

    return plt


def addMissingPercentiles(percentileCountDictionary):
    import numpy as np

    percentileSetComplete = set(np.arange(1, 1001, 1))
    percentilesPresent = set(percentileCountDictionary.keys())
    
    toAdd = percentileSetComplete - percentilesPresent
    for i in list(toAdd):
        percentileCountDictionary[i] = 0

    return None


import os
import subprocess
import sys
import math
import matplotlib.pyplot as plt
from scars import *

outFile = sys.argv[1]

# query pathogenic SNVs
exon_flank = 8
ensembl_ftp = "ftp://ftp.ensembl.org/pub/release-100/gtf/homo_sapiens/Homo_sapiens.GRCh38.100.gtf.gz"
hgmd_vcf = "/home/km369/.random_stuff/2019-4/hgmd_pro_2019.4_hg38.vcf"
hgmd_annotation = "/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/scars_pipeline_gnomad_hg38/non_coding_pathogenic_variants/hgmd_hg38.txt.gz"

hgmd_pathogenic_snv_filename = "/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/scars_pipeline_gnomad_hg38/non_coding_pathogenic_variants/noncoding_pathogenic_HGMD_Regulatory_DM_DM?_8bpFromSplice_hg38_sorted.bed"
snvs_annot_by_all_filename = "/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/scars_pipeline_gnomad_hg38/non_coding_pathogenic_variants/noncoding_pathogenic_positions_HGMD_Regulatory_DM_DM?_8bpFromSplice_hg38_sorted_annot_by_all_scores.bed"

non_coding_patho_snvs = scars_noncoding.getNonCodingPathogenic(exon_flank, ensembl_ftp, hgmd_vcf, hgmd_annotation)
scars_queries.writeToBed(non_coding_patho_snvs, hgmd_pathogenic_snv_filename)

SCARS_track = "/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/scars_pipeline_gnomad_hg38/nfe/scars_tracks/scars_autosomes_hg38.bed.gz"
CDTS_track = "/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/reimplementing_CDTS_gnomad_hg38/CDTS/CDTS_diff_perc_coordsorted_autosomes_gnomAD_N32299_hg38.bed.gz"
Orion_track = "/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/reimplementing_Orion_gnomad_hg38/Orion/Orion_autosomes_hg38.bed.gz"
CADD_track = "/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/CADD_v1.6/whole_genome_SNVs.tsv.gz"

SCARS_annotated_cmd = "tabix " + SCARS_track + " -R " + hgmd_pathogenic_snv_file + " | cut -f1-3" 
CDTS_annotated_cmd = "<(tabix " + CDTS_track + " -R " + hgmd_pathogenic_snv_file + ")"
Orion_annotated_cmd = "<(tabix " + Orion_track + " -R " + hgmd_pathogenic_snv_file + ")"
CADD_annotated_cmd = "<(tabix " + CADD_track + " -R <(awk '{gsub(\"chr\",\"\"); print $1,$3}' " + hgmd_pathogenic_snv_file + ") | awk '{print \"chr\"$1,$2-1,$2}' OFS='\t' - | bedops -m - )"

annotated_SNVs_cmd = SCARS_annotated_cmd + " | bedmap --echo --skip-unmapped - " + CDTS_annotated_cmd + " | bedmap --echo --skip-unmapped - " + \
        Orion_annotated_cmd + " | bedmap --echo --skip-unmapped - " + CADD_annotated_cmd + " > " + snvs_annot_by_all_file
subprocess.Popen(annotated_SNVs_cmd, shell=True, executable='/bin/bash')

SCARS_scores = list(map(float, os.popen("tabix " + SCARS_track + " -R " + snvs_annot_by_all_file).read().split()[3::4]))
CDTS_scores = list(map(float, os.popen("tabix " + CDTS_track + " -R " + snvs_annot_by_all_file).read().split()[6::7]))
Orion_scores = list(map(float, os.popen("tabix " + Orion_track + " -R " + snvs_annot_by_all_file).read().split()[3::4]))


query_SNV_CADD_cmd = "tabix " + CADD_track + " -R <(awk '{gsub(\"chr\", \"\"); print $1,$3,$4}' " + snvs_annot_by_all_filename + \
    ") | awk '{print \"chr\"$1,$2-1,$2,$3,$4,$5}' OFS='\t' - | bedmap --echo --echo-map --skip-unmapped - " + hgmd_pathogenic_snv_file + \
    " | sed 's/|/\t/g' - | sed 's/;/\t/g' - | awk '$5==$10 || $5==$14 || $5==$18 {print $1,$2,$3,$6}' OFS='\t' - | sort | uniq | sort-bed - | cut -f4"
CADD_scores_raw = subprocess.Popen(query_SNV_CADD_cmd, shell=True, executable='/bin/bash', stdout=subprocess.PIPE).communicate()
CADD_scores = list(map(float,x[0].split()))

CADD_samples_raw = os.popen("zcat " + CADD_track + " | grep -v '#' - | awk -v seed=$RANDOM 'BEGIN{srand(seed)} {x=rand(); if (x<2500) {print $5}}' - ").read().split()
CADD_samples = list(map(float, CADD_samples_raw))
CADD_percentiles = np.quantile(CADD_samples, np.arange(0, 1.001, 0.001))

CADD_percentile_scores = scars_assess.toPercentile(CADD_scores, CADD_percentiles)
CADD_percentile_scores_inv = [100-float(x) for x in CADD_percentile_scores]


# Ensure one digit after comma, ranging from 0.1 to 100.0
SCARS_corrected = [round(x, 1) for x in SCARS_scores]
CDTS_corrected = [math.ceil(10*x)/10 for x in CDTS_scores]
Orion_corrected = [round(x, 1) for x in Orion_scores]
CADD_corrected = [round(x, 1) for x in CADD_percentile_scores_inv]


fig, axs = plt.subplots(1, 2, figsize=(10,5))

axs[0].text(-19, 463, "a", size=18)
axs[1].text(-0.7, 260, "b", size=18)

percScoreToCumulativePercCountPlot(SCARS_corrected, axs[0])
percScoreToCumulativePercCountPlot(CDTS_corrected, axs[0])
percScoreToCumulativePercCountPlot(Orion_corrected, axs[0])
percScoreToCumulativePercCountPlot(CADD_corrected, axs[0], dashed=True)

percScoreToCumulativePercCountPlot(SCARS_corrected, axs[1], 5)
percScoreToCumulativePercCountPlot(CDTS_corrected, axs[1], 5)
percScoreToCumulativePercCountPlot(Orion_corrected, axs[1], 5)
percScoreToCumulativePercCountPlot(CADD_corrected, axs[1], 5, dashed=True)

fig.tight_layout()
fig.subplots_adjust(bottom=0.2)
axs[0].legend(['SCARS', 'CDTS', 'Orion', 'CADD'], loc="lower center", bbox_to_anchor=(1.1, -0.28), ncol=4, prop={'size': 12})

fig.savefig(outFile) 







