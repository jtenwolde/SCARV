import os
import numpy as np
import pandas as pd
import random
from scars import scars_assess
import sys

ancestry = sys.argv[1] 
window_size = 575

chr_list = ["chr" + str(i) for i in range(1, 23)]
chr_list.extend(["chrXnonPAR", "chrXPAR"])

chr_lengths_raw = [248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 
145138636, 138394717, 133797422, 135086622, 133275309, 114364328, 107043718, 101991189, 
90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 50818468, 156040895, 156040895]
chr_lengths = dict(zip(chr_list, chr_lengths_raw))

observed_entropy_sum = 0
expected_entropy_sum = 0
n_sites = 0

track_folder = "/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/scars_pipeline_gnomad_hg38/" + ancestry + "/scars_tracks/"

for chrom in chr_list:
	print("doing", chrom)
	inFile = track_folder + "entropies_" + chrom + ".bed.gz"
	df = pd.read_csv(inFile, sep='\t', header=None, usecols=[3,4,5], names=["H_Obs", "H_Exp", "Coverage"])
	df_covered = df.loc[df.Coverage > 0.9 * window_size]

	obs_sum_chrom, exp_sum_chrom = df_covered[['H_Obs', 'H_Exp']].sum()

	observed_entropy_sum += obs_sum_chrom
	expected_entropy_sum += exp_sum_chrom
	n_sites += df_covered.shape[0]

average_observed_entropy = observed_entropy_sum/n_sites
average_expected_entropy = expected_entropy_sum/n_sites



# sample proportion of sites to be able to approximate median absolute deviation
observed_entropy_deviations = []
expected_entropy_deviations = []

for chrom in chr_list:
	print("doing", chrom)
	inFile = track_folder + "entropies_" + chrom + ".bed.gz"
	skip = random.sample(range(chr_lengths[chrom]),
		k = chr_lengths[chrom] - chr_lengths[chrom]//300)

	df = pd.read_csv(inFile, sep='\t', header=None, usecols=[3,4,5], 
		names=["H_Obs", "H_Exp", "Coverage"], skiprows=skip)
	df_covered = df.loc[df.Coverage > 0.9 * window_size]
	
	observed_entropy_extension = abs(df_covered.H_Obs - average_observed_entropy).tolist()
	expected_entropy_extension = abs(df_covered.H_Exp - average_expected_entropy).tolist()

	observed_entropy_deviations.extend(observed_entropy_extension)
	expected_entropy_deviations.extend(expected_entropy_extension)

	
median_absolute_deviation_H_obs = np.median(observed_entropy_deviations)
median_absolute_deviation_H_exp = np.median(expected_entropy_deviations)



for chrom in chr_list:
	print("doing", chrom)

	inFile = track_folder + "entropies_" + chrom + ".bed.gz"
	outFile = track_folder + "scars_" + chrom + ".bed.gz"

	df = pd.read_csv(inFile, sep='\t', header=None, names=["Chromosome", "Start", "End", "H_Obs", "H_Exp", "Coverage"])
	
	df['Scars'] = (df.H_Obs - average_observed_entropy)/median_absolute_deviation_H_obs - (df.H_Exp - average_expected_entropy)/median_absolute_deviation_H_exp
	df_output = df[df.Coverage> 0.9 * window_size]

	df_output.to_csv(outFile, sep='\t', columns=["Chromosome", "Start", "End", "Scars", "Coverage"], index=False, header=False, compression='gzip')



SCARS_samples = os.popen("zcat " + track_folder + "scars_chr*.gz | awk -v seed=$RANDOM 'BEGIN{srand(seed)} {x=rand(); if (x<1/2500) {print $4}}' - ").read().split()
SCARS_percentiles =  np.quantile(list(map(float, SCARS_samples)), np.arange(0, 1.001, 0.001))

for chrom in chr_list:
	print("doing", chrom)
	inFile = track_folder + "scars_" + chrom + ".bed.gz"
	outFile = track_folder + "scars_percentiles_" + chrom + ".bed.gz"

	df = pd.read_csv(inFile, sep='\t', header=None, names=["Chromosome", "Start", "End", "Scars", "Coverage"])	
	df['Percentile'] = scars_assess.toPercentile(df['Scars'], SCARS_percentiles)

	df.to_csv(outFile, sep='\t', columns=["Chromosome", "Start", "End", "Percentile"], index=False, header=False, compression='gzip')


