# MK regions exploration
from scars import *

import matplotlib.pyplot as plt
import pyranges as pr
import pandas as pd
import numpy as np
import os

MK_regions_fn = "/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/classifier/humpy_analysis/MK-humpy400_-1_47-expansion_hg38.bed"
MK_regions = pr.read_bed(MK_regions_fn)
MK_regions.columns = ["Chromosome", "Start", "End", "nearest_min_cov", "max_atac"]

MK_regions.index = np.arange(len(MK_regions))
MK_regions_split = MK_regions.tile(1)

SCARS_clf_track_fn = "/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/classifier/tracks/scars_clf.bed.gz"

scores_raw = os.popen("tabix " + SCARS_clf_track_fn + " -R " + MK_regions_fn).read().split()
chr_indices = np.where(["chr" in x for x in scores_raw])[0]
scores_raw_split = [scores_raw[chr_indices[i]:chr_indices[i+1]] for i in np.arange(len(chr_indices)-1)]
scores = [sublist + [np.nan] if len(sublist)==3 else sublist for sublist in scores_raw_split]

scores_gr = scars_queries.coordlist_to_pyranges(scores, ["Chromosome", "Start", "End", "SCARS_clf"])
scores_gr.SCARS_clf = scores_gr.SCARS_clf.astype(float)

MK_regions_split_annot = MK_regions_split.join(scores_gr)
df = MK_regions_split_annot.as_df()[["index", "nearest_min_cov", "max_atac", "SCARS_clf"]]

summary = df.groupby("index").agg({"nearest_min_cov": "first", "max_atac": "first", "SCARS_clf": np.median})
summary['max_atac_bin'] = pd.qcut(summary.max_atac, 8)
print(summary.corr())




BPD_genes = pd.read_csv("/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/classifier/humpy_analysis/BPD_diagnostic_grade_genes.tsv", sep="\t", header=None)
BPD_genes = BPD_genes.rename(columns={1: "gene_name"})

ensembl_ftp = "ftp://ftp.ensembl.org/pub/release-100/gtf/homo_sapiens/Homo_sapiens.GRCh38.100.gtf.gz"
ensembl_genes = scars_queries.query_ensembl(ensembl_ftp, "gene")

ensembl_genes.gene_name = [str_spl[5] for str_spl in ensembl_genes.attributes.str.split('"')]
ensembl_genes = ensembl_genes[['gene_name']]

ensembl_genes.BPD_indicator = [name in BPD_genes.gene_name.values for name in ensembl_genes.gene_name]


BPD_indicator_by_index_raw = MK_regions.k_nearest(ensembl_genes).as_df()[["index", "BPD_indicator"]]
BPD_indicator_by_index = BPD_indicator_by_index_raw.groupby("index").any()

summary = pd.concat([summary, BPD_indicator_by_index], axis=1)




fig, axs = plt.subplots(1, 2, figsize=(10,5))

summary.boxplot(column='SCARS_clf', by='max_atac_bin', rot=45, fontsize=12, ax=axs[0])
axs[0].set_ylabel("Median SCARS-clf", fontsize=13)
axs[0].set_xlabel("ATAC-seq peak", fontsize=13)
axs[0].set_title("")


summary.boxplot(column='SCARS_clf', by='BPD_indicator', fontsize=12, ax=axs[1])
axs[1].set_ylabel("Median SCARS-clf", fontsize=13)
axs[1].set_xlabel("Diagnostic-grade BPD gene", fontsize=13)
axs[1].set_title("")

fig.suptitle("")
fig.tight_layout()
fig.savefig("/home/jwt44/mk_humpy_analysis.pdf")













