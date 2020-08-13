#script requires: pyVCF, pybedtools, shutil, urllib, contextlib, numpy, tempfile
from scars import scars_assess


# function to prepend "chr" in front of chromosome names to make all ranges constistent
def prepend_chr(x):
    x.chrom = 'chr' + x.chrom
    return x



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

    if (record.info[AF_str][0] < query['maf_range'][0]) or (record.info[AF_str][0] > query['maf_range'][1]):
        return ()

    out = (record.contig, record.start, record.stop, record.ref, record.alts[0],\
        str(record.info[AC_str][0]), str(record.info[AN_str]))

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



def convert_to_bed(data_list):
    import pybedtools

    bed = pybedtools.BedTool(data_list).sort()

    return(bed)



# assumes format chrX:pos, mean, median, over_1, over_5, over_10, over_15, over_20, over_25, over_30, over_50, over_100
def filter_coverage(coverage_file):
    import pybedtools
    import gzip

    f = gzip.open(coverage_file, 'r') 
    next(f)

    out_l = []
    for line in f:
        line_spl = line.decode().split('\t')
        if (float(line_spl[6]) >= 0.9) and (float(line_spl[11]) <= 0.1):
            pos = line_spl[0]
            out_l += [(pos[:4], int(pos[5:])-1, int(pos[5:]))]

    bed = pybedtools.BedTool(out_l)\
                    .sort()\
                    .merge()
    return bed



def filter_phyloP(phyloP_bw, chr_list, range=[-3,0]):
    import pyBigWig
    import pybedtools
    import numpy as np

    bw = pyBigWig.open(phyloP_bw)

    bed_l = []
    chr_present = list(bw.chroms().keys())
    for chrom in (chr_list and chr_present):
        vals = bw.intervals(chrom)
        l = [(chrom, r[0], r[1]) for r in vals if ((r[2] >= range[0]) and (r[2] <= range[1]))]
        bed_l += [pybedtools.BedTool(l).merge()]

    if len(bed_l) > 1:
        out = bed_l[0].cat(*bed_l[1:], postmerge=False)
    else:
        out = bed_l[0]

    return out


def filter_sites(bed, reliable_sites_bw):
    import pyBigWig
    import pybedtools

    out = []
    bw = pyBigWig.open(reliable_sites_bw)
    for r in bed:
        vals = bw.intervals(r.chrom, r.start, r.end)
        if vals is not None:
            for pos in vals:
                out+=[(r.chrom, max(pos[0], r.start), min(pos[1], r.end))]
    bed_reliable = pybedtools.BedTool(out)

    return bed_reliable



def get_reliable_sites(minOPR_file, coverage_file, chr_list, genome, snvs_fail, common_vars, bw_file=None):
    import pybedtools
    import pyBigWig

    minOPR_pass = pybedtools.BedTool(minOPR_file)
    coverage = filter_coverage(coverage_file)

    ranges = minOPR_pass.intersect(coverage)\
                        .subtract(snvs_fail)\
                        .subtract(common_vars)\
                        .sort()\
                        .merge()\
                        .saveas()

    # save reliable sites in form of bw for time efficient reuse
    chroms, starts, ends = [r.chrom for r in ranges], [r.start for r in ranges], [r.end for r in ranges]
    bw = pyBigWig.open(bw_file, "w")
    header = [(chrom, genome[chrom][1]) for chrom in chr_list]
    bw.addHeader(header)
    bw.addEntries(chroms, starts, ends=ends, values=[1.0] * len(chroms))
    bw.close()



def get_training_loci(reliable_sites_bw, phyloP_bw, ensembl_ftp, chr_list, genome):
    import pybedtools

    genome_bed = pybedtools.BedTool([(chrom, 0, genome[chrom][1]) for chrom in chr_list])

    reliable_sites = filter_sites(genome_bed, reliable_sites_bw)
    phylop = filter_phyloP(phyloP_bw, chr_list)
    exons = query_ensembl(ensembl_ftp, "exon", chr_list)

    out = reliable_sites.intersect(phylop)\
                        .subtract(exons)\
                        .sort()\
                        .merge()
    return out



# function that downloads the ensembl annotation
# filters out specific annotation using annot argument
# assumes chr_list contains "chr"
def query_ensembl(ensembl_ftp, annot, chr_list):
    import pybedtools
    import shutil
    import urllib.request as request
    from contextlib import closing
    import numpy as np
    import tempfile

    tmpdir = tempfile.TemporaryDirectory()

    with closing(request.urlopen(ensembl_ftp)) as r:
        with open(tmpdir.name + 'ensembl_data', 'wb') as f:
           shutil.copyfileobj(r, f)

    data = pybedtools.BedTool(tmpdir.name + 'ensembl_data')
    if "chr" not in data[0].chrom:
        data = data.each(prepend_chr)

    out = data.filter(lambda x: x.chrom in chr_list)\
              .filter(lambda x: x.fields[2]==annot)\
              .sort()  #no merge to preserve attributes

    return out



def query_sequence(bed, reference_fasta, flank, genome):
    import pybedtools
    import pandas as pd
    import numpy as np

    bed_annot = bed.slop(b=flank, genome=genome)\
                   .sequence(fi=reference_fasta)

    seqs = []
    for line in open(bed_annot.seqfn):
        if '>' in line: continue
        seqs += [line.strip()]

    seqs_nested = [[stri[i:j].upper() for i, j in zip(range(len(stri)-2*flank), range(2*flank+1, len(stri)+1))] for stri in seqs]
    seqs_flat = [item for sublist in seqs_nested for item in sublist]

    nucs = [list(seq) for seq in seqs_flat]
    nucs_df = pd.DataFrame(nucs)

    # one hot encode the sequences
    nucs_cat = nucs_df.apply(lambda x: pd.Categorical(x, categories = ['A', 'C', 'G', 'T']))
    nucs_bin = pd.get_dummies(nucs_cat)
    nucs_rshpd = np.array(nucs_bin, dtype=float).reshape(nucs_bin.shape[0], 2*flank+1, 4) 

    # address unidentified nucleotides by replacing them with [1/4,1/4,1/4,1/4]
    indices_N = np.where(np.sum(nucs_rshpd, axis=2) == 0)
    nucs_rshpd[indices_N[0], indices_N[1]] = np.repeat([[0.25, 0.25, 0.25, 0.25]], len(indices_N[0]), axis=0)

    return nucs_rshpd



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



def split_and_annot_data(loci, snvs_rare, genome, reference_fasta, pop_size):
    import pybedtools
    import numpy as np
    import pandas as pd

    n_sites = sum([len(r) for r in loci])

    neutral_snvs = snvs_rare.intersect(loci)

    AC = np.array([int(snv.fields[5]) for snv in neutral_snvs])

    AC_skew = np.random.binomial(n=AC, p=0.6)
    AC_bal = AC - AC_skew

    AN_skew = np.random.binomial(n=2*pop_size, p=0.6, size=n_sites)
    AN_rem = np.repeat(2*pop_size, n_sites) - AN_skew
    AN_bal = balance_data(AN_rem, AC_bal)

    N_bal = np.concatenate([AN_bal, AC_bal])
    N_skew = np.concatenate([AN_skew, AC_skew])

    N_cal = np.random.binomial(N_skew, p=2/3)
    N_test = N_skew - N_cal

    alts_ohe = np.array(pd.get_dummies([snv.fields[4] for snv in neutral_snvs]))
    refs_ohe = query_sequence(loci, reference_fasta, 0, genome)[:,0,:]
    nuc = np.concatenate([refs_ohe, alts_ohe])

    loci_spl = pybedtools.BedTool().window_maker(b=loci,w=1)
    hits = findOverlaps(neutral_snvs, loci_spl)
    indices = np.concatenate((np.arange(len(loci_spl)), hits[:,1]), axis=0)

    return indices, nuc, np.c_[N_bal, N_cal, N_test]



def findOverlaps(bed1, bed2):
    import pyranges as pr
    import pybedtools
    import pandas as pd
    import numpy as np

    gr1 = bed_to_pr(bed1)
    gr2 = bed_to_pr(bed2)

    gr1.id = np.arange(len(gr1))
    gr2.id = np.arange(len(gr2))

    out = gr1.join(gr2, suffix="_2").as_df()[['id', 'id_2']].to_numpy()

    return out



def bed_to_pr (bed):
    import pyranges as pr
    import pybedtools

    chroms = [int(r.chrom.replace("chr", "")) for r in bed]
    starts = [r.start for r in bed]
    ends = [r.end for r in bed]

    gr = pr.from_dict({"Chromosome": chroms, "Start": starts, "End": ends})

    return gr



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



