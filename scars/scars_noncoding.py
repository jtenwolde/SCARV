from scars import scars_queries

# function that converts HGMD vcf to pybedtool object
# filters for: 1) snvs, 2) CLASS=="DM" or "DM?"
def VCFtoBedTool_HGMD(VCF_file, chr_list):
    import vcf
    import pybedtools

    HGMD_reader = vcf.Reader(filename=VCF_file, encoding='utf-8')
    HGMD_out_l = []
    for record in HGMD_reader:
        is_snv = (len(record.REF)==1) and (len(record.ALT[0])==1)       # filter for snvs
        is_dm = (record.INFO['CLASS'] in ['DM', 'DM?'])                 # filter for class
        if (is_snv and is_dm):
            HGMD_out_l += [(record.CHROM, record.POS-1, record.POS)]    
    HGMD_bed = pybedtools.BedTool(HGMD_out_l)\
                         .each(prepend_chr)

    out = HGMD_bed.filter(lambda x: x.chrom in chr_list)\
                  .sort()

    return out



# function that converts ClinVar vcf to pybedtool object
# filters for: 1) snvs, 2) CLNSIG=="Pathogenic" or "Likely_pathogenic"
def VCFtoBedTool_ClinVar(VCF_file, chr_list, clnsig_l):
    import vcf
    import pybedtools
    
    ClinVar_reader = vcf.Reader(filename=VCF_file, encoding='utf-8')
    ClinVar_out_l = []
    for record in ClinVar_reader:
        is_snv = (record.INFO['CLNVC'] == "single_nucleotide_variant")                      # filter for snvs
        is_p_or_lp = False if 'CLNSIG' not in record.INFO.keys()\
                else (record.INFO['CLNSIG'][0] in clnsig_l)      # filter for clinical significance
        if (is_snv and is_p_or_lp):
            ClinVar_out_l += [(record.CHROM, record.POS-1, record.POS)]
    ClinVar_bed = pybedtools.BedTool(ClinVar_out_l)\
                            .each(prepend_chr)

    out = ClinVar_bed.filter(lambda x: x.chrom in chr_list)\
                     .sort()

    return out



def get_non_coding_pathogenic(exon_flank, chr_list, genome, ensembl_ftp, hgmd_vcf, clinvar_vcf):
    import pybedtools

    exons = scars_queries.query_ensembl(ensembl_ftp, "exon", chr_list)
    transcripts = scars_queries.query_ensembl(ensembl_ftp, "transcript", chr_list)
    introns = transcripts.subtract(exons)

    exons_extensions = exons.slop(b=exon_flank, genome=genome)\
                            .intersect(introns)

    hgmd_snvs = VCFtoBedTool_HGMD(hgmd_vcf, chr_list)
    clinvar_snvs_p_or_lp = VCFtoBedTool_ClinVar(clinvar_vcf, chr_list, ['Pathogenic', 'Likely_pathogenic'])
    clinvar_snvs_b_or_lb = VCFtoBedTool_ClinVar(clinvar_vcf, chr_list, ['Benign', 'Likely_benign'])

    exclusion = exons.cat(*[exons_extensions, clinvar_snvs_b_or_lb])

    out_mgd = hgmd_snvs.cat(clinvar_snvs_p_or_lp)\
                       .subtract(exclusion)
    out = pybedtools.BedTool().window_maker(b=out_mgd, w=1)

    return out
