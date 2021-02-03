from scars import *

import numpy as np
import pandas as pd
import pyranges as pr
import sys

ancestry = sys.argv[1]
pop_size = int(sys.argv[2])				  # N_XX + N_XY
pop_size_chrXnonPAR = int(sys.argv[3])    # 2 * N_XX + N_XY

inFileSNVs = sys.argv[4]
inFileFAILVariants = sys.argv[5]

outFileReliableSites = sys.argv[6]
outFileTrainingLoci = sys.argv[7]

phyloP_bw = "/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/hg38.phyloP100way.bw"
gnomad_coverage_tsv = "/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/gnomad_v3_hg38/gnomad.genomes.r3.0.1.coverage.summary.tsv.bgz"
ensembl_ftp = "ftp://ftp.ensembl.org/pub/release-100/gtf/homo_sapiens/Homo_sapiens.GRCh38.100.gtf.gz"

pass_snvs = scars_queries.load_variants(inFileSNVs)

chrXnonPAR = pr.from_dict({"Chromosome": ['chrX', 'chrX'], 
	"Start": [0, 2781479], "End": [10001, 155701383]})
pass_snvs_autosomes_and_chrXPAR = pass_snvs.subtract(chrXnonPAR)
pass_snvs_chrXnonPAR = pass_snvs.join(chrXnonPAR, how='right').drop(['Start_b', 'End_b'])

pass_nonsingletons_autosomes_and_chrXPAR = pass_snvs_autosomes_and_chrXPAR[pass_snvs_autosomes_and_chrXPAR.ac/pass_snvs_autosomes_and_chrXPAR.an>1.5/(2*pop_size)]
pass_nonsingletons_chrXnonPAR = pass_snvs_chrXnonPAR[pass_snvs_chrXnonPAR.ac/pass_snvs_chrXnonPAR.an>1.5/pop_size_chrXnonPAR]
pass_nonsingleton_snvs = pr.PyRanges(pd.concat([pass_nonsingletons_autosomes_and_chrXPAR.as_df(), pass_nonsingletons_chrXnonPAR.as_df()], ignore_index=True))

fail_variants = scars_queries.load_variants(inFileFAILVariants)

reliable_sites = scars_queries.get_reliable_sites(gnomad_coverage_tsv, fail_variants, pop_split)
reliable_sites = reliable_sites[reliable_sites.Chromosome != "chrY"]
scars_queries.write_to_bed(reliable_sites, outFileReliableSites)

training_loci = scars_queries.get_training_loci(reliable_sites, phyloP_bw, ensembl_ftp, pass_nonsingleton_snvs)
scars_queries.write_to_bed(training_loci, outFileTrainingLoci)
