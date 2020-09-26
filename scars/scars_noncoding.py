from scars import scars_queries

# function that converts HGMD vcf to pybedtool object
# filters for: 1) snvs, 2) CLASS=="DM" or "DM?"
def VCFtoBedTool_HGMD(vcf_file, chr_list):
    from pysam import VariantFile
    import pybedtools

    HGMD_out_l = []
    bcf_in = VariantFile(vcf_file)

    for record in bcf_in.fetch():
        is_snv = (len(record.ref)==1 and len(record.alts[0])==1)
        is_dm = record.info['CLASS'] in ['DM', 'DM?']
        if (is_snv and is_dm):
            HGMD_out_l += [(record.contig, record.start, record.stop)]

    HGMD_bed = pybedtools.BedTool(HGMD_out_l)\
                         .each(scars_queries.prepend_chr)

    out = HGMD_bed.filter(lambda x: x.chrom in chr_list)\
                  .sort()
     
    return out



# function that converts ClinVar vcf to pybedtool object
# filters for: 1) snvs, 2) CLNSIG=="Pathogenic" or "Likely_pathogenic"
def VCFtoBedTool_ClinVar(VCF_file, chr_list, clnsig_l):
    from pysam import VariantFile
    import pybedtools
    
    bcf_in = VariantFile(VCF_file)
    ClinVar_out_l = []
    for record in bcf_in.fetch():
        is_snv = (record.info['CLNVC'] == "single_nucleotide_variant")                      # filter for snvs
        is_p_or_lp = False if 'CLNSIG' not in record.info.keys()\
                else (record.info['CLNSIG'][0] in clnsig_l)      # filter for clinical significance
        if (is_snv and is_p_or_lp):
            ClinVar_out_l += [(record.contig, record.start, record.stop)]
    ClinVar_bed = pybedtools.BedTool(ClinVar_out_l)\
                            .each(scars_queries.prepend_chr)

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
