from scars import scars_queries

def query_codata_CpG_meth(loci, filename_list):
    import pyBigWig
    import numpy as np

    N = sum([len(r) for r in loci])
    out = np.empty(shape=(N, 4 * len(filename_list)), dtype=np.int8)
    for ix, file in enumerate(filename_list):
        bw = pyBigWig.open(file)
        scores_nested = [bw.values(r.chrom, r.start, r.end) for r in loci]
        scores = np.array([item for sublist in scores_nested for item in sublist])
        scores[np.isnan(scores)] = -1
        arr = np.zeros(shape=(N, 4), dtype=np.int8)
        arr[scores>0.6, 3] = 1
        arr[(scores>0.2) & (scores<=0.6), 2] = 1 
        arr[(scores>=0) & (scores<=0.2), 1] = 1
        arr[scores==-1, 0] = 1
        out[:,(4*ix):(4*(ix+1))] = arr

    return out



# assumes bigWig files with scores that don't need binning/ohe
# furthermore it is useful to sort the filename_list to make covariates identifiable
def query_codata_human_ovary(loci, filename_list, impute_values):
    import pyBigWig
    import numpy as np

    N = sum([len(r) for r in loci])
    out = np.empty(shape=(N, len(filename_list)))
    for ix, file in enumerate(filename_list):
        bw = pyBigWig.open(file)
        scores_nested = [bw.values(r.chrom, r.start, r.end) for r in loci]
        scores = np.array([item for sublist in scores_nested for item in sublist])
        scores[scores==np.nan] = impute_values[ix]
        out[:,ix] = scores

    return out



# functions required for generating the spline bases later
# kv: knot vector, u: samples, k: index of control point, d: degree
def coxDeBoor(k, d, u, kv):
    # Test for end conditions
    if (d == 0):
        return ((u - kv[k] >= 0) & (u - kv[k + 1] <= 0))*1

    denom1 = kv[k + d] - kv[k]
    term1 = 0
    if denom1 > 0:
        term1 = ((u - kv[k]) / denom1) * coxDeBoor(k, d - 1,u,kv)

    denom2 = kv[k + d + 1] - kv[k + 1]
    term2 = 0
    if denom2 > 0:
        term2 = ((-(u - kv[k + d + 1]) / denom2) * coxDeBoor(k + 1, d - 1, u,kv))

    return term1 + term2



# order is degree + 1
def gen_spline_basis(x, order, knots_raw, intercept, i = 1):
    import numpy as np 

    knots = np.concatenate((np.repeat(knots_raw[0], order), knots_raw,
        np.repeat(knots_raw[-1], order)))

    MM = np.empty(shape=(len(x), order + len(knots_raw) - i))

    for j in range(order + len(knots_raw) - i):
        MM[:,j] = coxDeBoor(j, order, x, knots)

    if intercept:
        return MM
    else:
        return MM[:,1:]



def annot_with_constr(r, constr_by_ID_dict, impute):
    if r.attrs['exon_id'] in constr_by_ID_dict.keys():
        r.attrs['constr'] = str(constr_by_ID_dict[r.attrs['exon_id']])
    else:
        r.attrs['constr'] = impute
    return r



def query_nearest_exon(loci, exons_annot, chr_list):
    import pybedtools
    import pandas as pd
    import numpy as np

    # split loci into individual sites
    loci_spl = pybedtools.BedTool().window_maker(b=loci, w=1)

    # find closest exons, all ties are included, so that maximum pLI can taken across ties
    loci_with_nearest_exon = loci_spl.closest(exons_annot, d=True)

    # find maximum pLI for all ties
    chr_pos_dist_constr_l = [[r[0], r[2], int(r[-1]), float(r[11].split(';')[-2][7:])] for r in loci_with_nearest_exon]
    chr_pos_dist_constr_df = pd.DataFrame(chr_pos_dist_constr_l, columns=['chr','pos','dist','constr'])
    chr_pos_dist_min_constr = chr_pos_dist_constr_df.groupby(['chr','pos','dist'], sort=False)['constr'].min().reset_index()

    dist_and_constr = chr_pos_dist_min_constr[['dist','constr']].values
    constr_unknown = (dist_and_constr[:,1]==-1)
 
    # need to feed both constraint and distance to bspline function
    dist_basis = gen_spline_basis(dist_and_constr[:,0], order=3, knots_raw = [0, 1, int(2e3), int(1e5), 20083767], intercept=False)
    constr_basis = gen_spline_basis(dist_and_constr[:,1], order=3, knots_raw = [0, 0.01, 0.3, 0.99, 27.50007], intercept=False)

    r = np.repeat(range(dist_basis.shape[1]), constr_basis.shape[1])   
    c = np.tile(range(constr_basis.shape[1]), dist_basis.shape[1])

    spl_basis = dist_basis[:,c] * constr_basis[:,r]
    spl_basis[constr_unknown] = np.nan

    return spl_basis




def get_codata(loci, human_ovary_bw_files, human_ovary_impute_vals, CpG_methylation_files, exons_annot, chr_list):
    import numpy as np

    human_ovary_tracks = query_codata_human_ovary(loci, human_ovary_bw_files, human_ovary_impute_vals)
    CpG_methylation = query_codata_CpG_meth(loci, CpG_methylation_files)
    exon_spline_basis = query_nearest_exon(loci, exons_annot, chr_list)

    codata = np.c_[human_ovary_tracks, CpG_methylation, exon_spline_basis]
    return codata




