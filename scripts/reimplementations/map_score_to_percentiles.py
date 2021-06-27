from scars import scars_assess
import pandas as pd
import numpy as np
import os

chr_list = ["chr" + str(i) for i in range(1, 23)]

Orion_samples = os.popen("cat Orion/Orion_chr*.bed | awk -v seed=$RANDOM 'BEGIN{srand(seed)} {x=rand(); if (x<1/2500) {print $4}}' - ").read().split()
Orion_percentiles =  np.quantile(list(map(float, Orion_samples)), np.arange(0, 1.001, 0.001))

for chrom in chr_list:
	print("doing", chrom)

	inFile = "Orion/Orion_" + chrom + ".bed"
	outFile = "Orion/Orion_percentiles_" + chrom + ".bed.gz"

	df = pd.read_csv(inFile, sep='\t', header=None, names=["Chromosome", "Start", "End", "Orion"])	
	df['Percentile'] = scars_assess.toPercentile(df['Orion'], Orion_percentiles)

	df.to_csv(outFile, sep='\t', columns=["Chromosome", "Start", "End", "Percentile"], index=False, header=False, compression='gzip')

