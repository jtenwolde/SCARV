# script to get genome-wide percentiles for classifier

def query_annotation_chromosome (annot_fn, name, chromosome):
    import os
    import pyranges as pr
    import numpy as np

    query = os.popen("tabix " + annot_fn + " " + chromosome).read().split()

    x = np.where(["chr" in x for x in query])[0]
    ncols = x[1] - x[0]

    gr_annot = pr.from_dict({'Chromosome': query[::ncols], 'Start': list(map(int, query[1::ncols])), 'End': list(map(int, query[2::ncols])), name: list(map(float, query[3::ncols]))})
    
    df_annot = gr_annot.tile(1).as_df()
    df_annot['Chromosome'] = df_annot['Chromosome'].astype("object")

    gr_annot_max = pr.PyRanges(df_annot.groupby(['Chromosome','Start','End'])[name].max().dropna().reset_index())

    return gr_annot_max



def query_nearest_chromosome (annot_fn, chrom):
	import numpy as np
	import pyranges as pr

	annot = pr.read_bed(annot_fn)[chrom].as_df()
	annot.columns = ['Chromosome', 'Start', 'End', "LogStdExp"]

	score = np.empty(chr_lengths[chrom])
	row_ix = 0

	row_val = annot.loc[row_ix, 'LogStdExp']
	row_val_next = annot.loc[row_ix + 1, 'LogStdExp']

	row_start = annot.loc[row_ix, 'Start']
	row_start_next = annot.loc[row_ix + 1, 'Start']

	row_end = annot.loc[row_ix, 'End']
	row_end_next = annot.loc[row_ix + 1, 'End']

	old_boundary = 0
	while row_ix < annot.shape[0]:
		new_boundary = max(row_end, (row_end + row_start_next)//2)
		score[old_boundary:new_boundary] = row_val
		old_boundary = new_boundary

		if row_ix == annot.shape[0] - 1:
			score[old_boundary:] = row_val_next
			break

		if row_end_next > row_end:
			row_val = row_val_next 
			row_start = row_start_next
			row_end = row_end_next

		row_val_next = annot.loc[row_ix + 1, 'LogStdExp']
		row_start_next = annot.loc[row_ix + 1, 'Start']	
		row_end_next = annot.loc[row_ix + 1, 'End']

		row_ix += 1

	return score


import os
import numpy as np
import pandas as pd
import random
from scarv import scarv_assess
import sys
import pyranges as pr
import pybedtools
from xgboost import XGBClassifier
import scipy.stats as ss

ancestry = "nfe" 

chr_list = ["chr" + str(i) for i in range(1, 23)]
chr_list.extend(["chrX"])

chr_lengths_raw = [248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 
145138636, 138394717, 133797422, 135086622, 133275309, 114364328, 107043718, 101991189, 
90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 50818468, 156040895]
chr_lengths = dict(zip(chr_list, chr_lengths_raw))

gencode_exons = pr.read_bed("/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/classifier/gencode_v27/gencode_exons_hg38_v27.bed").merge()

scars_fn = "/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/scars_pipeline_gnomad_hg38/nfe/scars_tracks/scars_hg38_incl_1bp_indels_raw.bed.gz"
linsight_fn = "/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/classifier/covariates/LINSIGHT/LINSIGHT_hg38.bed.gz"
expression_annotation_LogStdExp_fn = "/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/classifier/covariates/GTEx/expression_values_by_gene_LogStdExp.bed"
remm_fn = "/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/classifier/covariates/ReMM/ReMM.v0.3.1_hg38.bed.gz"
funseq_fn = "/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/classifier/covariates/funseq/hg38_NCscore_funseq216.bed.gz"

preds = np.empty(sum(chr_lengths_raw))
counter = 0

clf_fn = "/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/classifier/models/scars_clf.model"
clf = XGBClassifier() 
clf.load_model(clf_fn)

for chrom in chr_list:
	print("doing", chrom)

	chrom_gr = pr.from_dict({'Chromosome': [chrom], 'Start': [0], 'End': [chr_lengths[chrom]]})

	X = np.empty(shape=(chr_lengths[chrom], 5))
	X[:] = np.nan

	LINSIGHT = query_annotation_chromosome (linsight_fn, "LINSIGHT", chrom)
	X[LINSIGHT.Start, 0] = LINSIGHT.LINSIGHT

	SCARS = query_annotation_chromosome (scars_fn, "SCARS", chrom)
	X[SCARS.Start, 1] = SCARS.SCARS

	funseq = query_annotation_chromosome (funseq_fn, "funseq", chrom)
	X[funseq.Start, 2] = funseq.funseq

	ReMM = query_annotation_chromosome (remm_fn, "ReMM", chrom)
	X[ReMM.Start, 3] = ReMM.ReMM

	expression_annotation_LogStdExp = query_nearest_chromosome (expression_annotation_LogStdExp_fn, chrom)
	X[:, 4] = expression_annotation_LogStdExp

	chrom_preds = clf.predict_proba(X)[:,1]
	chrom_preds[gencode_exons[chrom].tile(1).Start] = np.nan

	preds[counter:(counter+chr_lengths[chrom])] = chrom_preds
	counter += chr_lengths[chrom]

percentiles = (100 * pd.Series(preds).rank()/np.sum(~np.isnan(preds))).to_numpy()


# write to file
outFolder = "/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/classifier/tracks/"
write_counter = 0
for chrom in chr_list:
	chrom_length = chr_lengths[chrom]
	score = percentiles[write_counter:(write_counter + chrom_length)]
	
	df = pd.DataFrame({'Chromosome': chrom, 'Start': np.arange(chrom_length), 'End': np.arange(chrom_length) + 1, 'Score': score})
	
	outFile = outFolder + "scars_clf_" + chrom + ".bed.gz"
	df.to_csv(outFile, sep='\t', columns=["Chromosome", "Start", "End", "Score"], index=False, header=False, compression='gzip')
	
	write_counter += chr_lengths[chrom]






