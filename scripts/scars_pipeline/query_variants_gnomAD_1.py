from scars import *

import numpy as np
import pandas as pd
import pyranges as pr
import pybedtools
import sys

ancestry = sys.argv[1]
outFileSNVs = sys.argv[2]
outFileDeletions = sys.argv[3]
outFileFAILVariants = sys.argv[4]

gnomad_vcf = "/rds/project/who1000-1/rds-who1000-cbrc/user/jwt44/gnomad_v3_hg38/gnomad.genomes.r3.0.sites.vcf.bgz"

# multiple simultaneous queries so that vcf file is read only once
pass_snvs_query = {'variant_type': "snv", 'ancestry': ancestry, 'PASS': True}
pass_indels_query = {'variant_type': "indel", 'ancestry': ancestry, 'PASS': True}
fail_vars_query = {'variant_type': None, 'ancestry': ancestry, 'PASS': False}

queries = [pass_snvs_query, pass_indels_query, fail_vars_query]
data_dict = scars_queries.query_vcf(gnomad_vcf, queries)

pass_snvs = scars_queries.coordlist_to_pyranges(data_dict[0], entryNames=["Chromosome", "Start", "End", "ref", "alt", "ac", "an"])
pass_indels = scars_queries.coordlist_to_pyranges(data_dict[1], entryNames=["Chromosome", "Start", "End", "ref", "alt", "ac", "an"])
fail_vars = scars_queries.coordlist_to_pyranges(data_dict[2], entryNames=["Chromosome", "Start", "End", "ref", "alt", "ac", "an"])

# the vcf query flips ref & alt for minor alleles with MAF > 0.5, which requires adjustments in the case of multiallelicism
corrected_inverted_multiallelics = scars_snv_correction.correct_inverted_multiallelics(pass_snvs)
corrected_pass_snvs = scars_snv_correction.insert_corrected_snvs(corrected_inverted_multiallelics, pass_snvs)

corrected_pass_snvs.ac = corrected_pass_snvs.ac.astype(int)
corrected_pass_snvs.an = corrected_pass_snvs.an.astype(int)
scars_queries.write_to_bed(corrected_pass_snvs, outFileSNVs)

pass_deletions = pass_indels[(pass_indels.ref.str.len() > 1) & (pass_indels.alt.str.len() == 1)]
pass_deletions.ac = pass_deletions.ac.astype(int)
pass_deletions.an = pass_deletions.an.astype(int)
scars_queries.write_to_bed(pass_deletions, outFileDeletions)

fail_snvs_and_deletions = fail_vars[fail_vars.alt.str.len() == 1]
scars_queries.write_to_bed(fail_snvs_and_deletions, outFileFAILVariants)


