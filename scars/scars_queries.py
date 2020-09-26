from scars import scars_assess


def coordlist_to_pyranges (pos_l, entryNames = ["Chromosome", "Start", "End"]):
    assert len(pos_l[0]) == len(entryNames), "not enough entry names provided"

    import pyranges as pr 

    values_by_entry = zip(*pos_l)
    pos_dict = {entry: values for entry, values in zip(entryNames, values_by_entry)}
    pos_pr = pr.from_dict(pos_dict)

    return pos_pr 



def sample_loci_from_pr (gr, n_loci):
    import pyranges as pr
    import numpy as np

    N = gr.length
    gr_spl = gr.tile(1)

    ix = np.random.choice(range(N), n_loci, replace=False)
    sample_pr = pr.PyRanges(gr_spl.as_df().iloc[ix,:])

    return sample_pr



def process_record(record, query):
    AC_str = "AC_" + query['ethnicity']
    AF_str = "AF_" + query['ethnicity']
    AN_str = "AN_" + query['ethnicity']

    if query['PASS'] is not None:
            if ("PASS" in list(record.filter)) != query['PASS']:
                return ()

    if query['variant_type'] is not None:
        if record.info['variant_type'] != query['variant_type']:
            return ()
    if record.info[AC_str][0] == 0:
        return ()

    if ((record.info[AF_str][0] < query['maf_range'][0]) or (record.info[AF_str][0] > query['maf_range'][1])) and\
        ((1-record.info[AF_str][0] < query['maf_range'][0]) or (1-record.info[AF_str][0] > query['maf_range'][1])):
        return ()

    if (record.info[AF_str][0] <= 0.5):
        out = (record.contig, record.start, record.stop, record.ref, record.alts[0],\
            str(record.info[AC_str][0]), str(record.info[AN_str]))
    else:
        out = (record.contig, record.start, record.stop, record.alts[0], record.ref,\
            str(record.info[AN_str] - record.info[AC_str][0]), str(record.info[AN_str]))

    return out



def query_vcf(vcf_file, queries):
    from pysam import VariantFile

    data_dict = dict()
    for ix in range(len(queries)):
        data_dict[ix] = list()

    bcf_in = VariantFile(vcf_file)

    for record in bcf_in.fetch():
        for ix, query in enumerate(queries):
            to_add = process_record(record, query)
            if len(to_add) > 0:
                data_dict[ix].append(to_add)

    return(data_dict)



# assumes format chrX:pos, mean, median, over_1, over_5, over_10, over_15, over_20, over_25, over_30, over_50, over_100
def filter_coverage(coverage_file):
    import gzip

    f = gzip.open(coverage_file, 'r') 
    next(f) #remove header

    out_l = []
    for line in f:
        line_spl = line.decode().split('\t')
        if (float(line_spl[6]) >= 0.9) and (float(line_spl[11]) <= 0.1):
            pos = line_spl[0].split(':')
            out_l.append((pos[0], int(pos[1])-1, int(pos[1])))

    return out_l



def filter_phyloP(phyloP_bw, phyloP_range=[-3,0]):
    import pyBigWig

    bw = pyBigWig.open(phyloP_bw)

    out_l = []
    for chrom in bw.chroms().keys():
        vals = bw.intervals(chrom)
        l = [(chrom, record[0], record[1]) for record in vals\
            if ((record[2] >= phyloP_range[0]) and (record[2] <= phyloP_range[1]))]
        out_l.extend(l)

    return out_l



def get_reliable_sites(coverage_file, genome, snvs_fail_pr, common_vars_pr):
    import pyranges as pr

    covered_loci_l = filter_coverage(coverage_file)
    covered_loci_pr = coordlist_to_pyranges(covered_loci_l)
    covered_loci_mgd_pr = covered_loci_pr.merge()

    ranges = covered_loci_mgd_pr.subtract(snvs_fail_pr)\
                                .subtract(common_vars_pr)\

    return ranges




def get_training_loci(reliable_sites_pr, phyloP_bw, ensembl_ftp):
    import pyranges as pr
    
    phylop_filtered_l = filter_phyloP(phyloP_bw)
    phylop_filtered_pr = coordlist_to_pyranges(phylop_filtered_l)
    phylop_filtered_mgd_pr = phylop_filtered_pr.merge()

    exons = query_ensembl(ensembl_ftp, "exon")

    training_loci_pr = reliable_sites_pr.intersect(phylop_filtered_mgd_pr)\
                                        .subtract(exons)
    
    return training_loci_pr



# gtf format uses inclusive range 
def query_ensembl(ensembl_ftp, annot):
    import pandas as pd
    import pyranges as pr

    data = pd.read_csv(ensembl_ftp, sep='\t', header=None,\
        names=['Chromosome', 'source', 'type', 'Start', 'End', 'score', 'strand',\
        'phase', 'attributes'], skiprows=[i for i in range(5)], compression="gzip")
    data['Start'] = data['Start'] - 1 

    data_correct_type = data.loc[data['type']==annot]
    data_correct_type["Chromosome"] = ["chr" + str(chrom) for chrom in data_correct_type["Chromosome"]]

    data_correct_type_pr = pr.PyRanges(data_correct_type)

    return data_correct_type_pr



def query_sequence(gr, reference_fasta, flank):
    import pandas as pd
    import numpy as np
    import pyranges as pr

    gr_extd = gr.slack(flank)
    seqs = pr.get_fasta(gr_extd, reference_fasta)

    seqs_nested = [[stri[i:j] for i, j in zip(range(len(stri)-2*flank), range(2*flank+1, len(stri)+1))] for stri in seqs]
    seqs_flat = [item for sublist in seqs_nested for item in sublist]

    seqs_flat_upper = [seq.upper() for seq in seqs_flat]

    nucs = [list(seq) for seq in seqs_flat_upper]
    nucs_df = pd.DataFrame(nucs)

    # one hot encode the sequences
    nucs_cat = nucs_df.apply(lambda x: pd.Categorical(x, categories = ['A', 'C', 'G', 'T']))
    nucs_bin = pd.get_dummies(nucs_cat)
    nucs_rshpd = np.array(nucs_bin, dtype=float).reshape(nucs_bin.shape[0], 2*flank+1, 4) 

    # address unidentified nucleotides by replacing them with [1/4,1/4,1/4,1/4]
    indices_N = np.where(np.sum(nucs_rshpd, axis=2) == 0)
    nucs_rshpd[indices_N[0], indices_N[1]] = np.repeat([[0.25, 0.25, 0.25, 0.25]], len(indices_N[0]), axis=0)

    return nucs_rshpd


def correct_refs(gr, snvs_pr, seq):
    import numpy as np
    import pandas as pd
    import pyranges as pr

    flank = seq.shape[1]//2

    gr_spl = gr.tile(1)
    row_ix = pd.Series(range(len(gr_spl)), name="id")
    gr_spl = gr_spl.insert(row_ix)

    gr_spl_on_snvs_hits = gr_spl.join(snvs_pr)
    
    anyHits = (gr_spl_on_snvs_hits.length != 0)
    
    if anyHits:
        gr_spl_hit_ids = gr_spl_on_snvs_hits.id

        nucs_df = pd.DataFrame(gr_spl_on_snvs_hits.ref)
        nucs_cat = nucs_df.apply(lambda x: pd.Categorical(x, categories = ['A', 'C', 'G', 'T']))
        
        refs_from_vcf_ohe = np.array(pd.get_dummies(nucs_cat))
        seq[gr_spl_hit_ids, flank] = refs_from_vcf_ohe
        


def balance_data(AN, AC):
    import collections
    import numpy as np

    samples = np.random.choice(len(AN), sum(AC), p=AN/sum(AN))
    counts = collections.Counter(samples)
    
    keys = [int(key) for key in counts.keys()]
    vals = [int(val) for val in counts.values()]

    AN_out = np.repeat(0, len(AN))
    AN_out[keys] = vals

    return AN_out



# does not require gr to be width 1
def split_data(gr, snvs_pr, seq, pop_size):
    import numpy as np
    import pandas as pd
    import pyranges as pr

    n_sites = sum(gr.lengths())
    flank = seq.shape[1]//2

    gr_spl = gr.tile(1)
    row_ix = pd.Series(range(n_sites), name="id")
    gr_spl = gr_spl.insert(row_ix)

    neutral_snvs = snvs_pr.join(gr_spl)

    AC = np.array(neutral_snvs.ac, dtype=int)

    AC_skew = np.random.binomial(n=AC, p=0.6)
    AC_bal = AC - AC_skew

    AN_skew = np.random.binomial(n=2*pop_size, p=0.6, size=n_sites)
    AN_rem = np.repeat(2*pop_size, n_sites) - AN_skew
    AN_bal = balance_data(AN_rem, AC_bal)

    N_bal = np.concatenate([AN_bal, AC_bal])
    N_skew = np.concatenate([AN_skew, AC_skew])

    N_cal = np.random.binomial(N_skew, p=2/3)
    N_test = N_skew - N_cal

    alts_ohe = np.array(pd.get_dummies(neutral_snvs.alt))
    refs_ohe = seq[:,flank,:]
    nuc = np.concatenate([refs_ohe, alts_ohe])

    indices = np.concatenate((np.arange(n_sites), neutral_snvs.id), axis=0)

    return indices, nuc, np.c_[N_bal, N_cal, N_test]



def average_bw (bw_fn, chr_list, genome):
    import pyBigWig
    import numpy as np

    bw = pyBigWig.open(bw_fn)

    mean_by_chr = np.empty(len(chr_list))
    chr_lengths = [genome[chrom][1] for chrom in chr_list]
    for ix, chrom in enumerate(chr_list):
        mean_by_chr[ix] = bw.stats(chrom)[0]

    bw.close()

    return np.average(mean_by_chr, weights=chr_lengths)



