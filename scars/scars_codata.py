from scars import scars_queries

def get_codata(gr, human_ovary_bw_files, human_ovary_impute_vals, CpG_methylation_files, exons_annot):
    import numpy as np

    human_ovary_tracks = query_codata_human_ovary(gr, human_ovary_bw_files, human_ovary_impute_vals)
    CpG_methylation = query_codata_CpG_meth(gr, CpG_methylation_files)
    exon_spline_basis = query_nearest_exon(gr, exons_annot)

    codata = np.c_[human_ovary_tracks, CpG_methylation, exon_spline_basis]
    return codata


def query_codata_CpG_meth(gr, filename_list):
    import pyBigWig
    import numpy as np

    N = gr.length
    out = np.empty(shape=(N, 4 * len(filename_list)), dtype=np.int8)
    for ix, file in enumerate(filename_list):
        bw = pyBigWig.open(file)
        scores_nested = [bw.values(row.Chromosome, row.Start, row.End) for row in gr.as_df().itertuples()]
        scores = np.array([item for sublist in scores_nested for item in sublist])
        scores[np.isnan(scores)] = -1
        arr = one_hot_encode_methylation_signal(scores)
        out[:,(4*ix):(4*(ix+1))] = arr
        bw.close()

    return out


def one_hot_encode_methylation_signal(scores):
        import numpy as np

        N = len(scores)
        arr = np.zeros(shape=(N, 4), dtype=np.int8)

        arr[scores>60, 3] = 1
        arr[(scores>20) & (scores<=60), 2] = 1 
        arr[(scores>=0) & (scores<=20), 1] = 1
        arr[scores==-1, 0] = 1

        return arr 


# assumes bigWig files with scores that don't need binning/ohe
# furthermore it is useful to sort the filename_list to make covariates identifiable
def query_codata_human_ovary(gr, filename_list, impute_values):
    import pyranges as pr
    import pyBigWig
    import numpy as np

    N = gr.length
    out = np.empty(shape=(N, len(filename_list)))
    for ix, file in enumerate(filename_list):
        bw = pyBigWig.open(file)
        scores_nested = [bw.values(row.Chromosome, row.Start, row.End) for row in gr.as_df().itertuples()]
        scores = np.array([item for sublist in scores_nested for item in sublist])
        scores[scores==np.nan] = impute_values[ix]
        out[:,ix] = scores
        bw.close()

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



def query_nearest_exon(gr, exons_annot):
    import pybedtools
    import pandas as pd
    import numpy as np

    # split loci into individual sites and find closest exons, all ties are included, so that strongest constraint can taken across ties
    gr_spl = gr.tile(1)
    gr_spl_with_nearest_exon = gr_spl.k_nearest(exons_annot)

    # find most constrained for all ties
    chr_pos_dist_constr = gr_spl_with_nearest_exon.as_df()[['Chromosome', 'End', 'constr', 'Distance']]

    chr_pos_dist_constr['Distance'] = np.abs(chr_pos_dist_constr['Distance'])
    chr_pos_dist_constr['Chromosome'] = chr_pos_dist_constr['Chromosome'].astype(str) # change from categorical to string type to avoid outer product blow up by groupby

    chr_pos_dist_min_constr = chr_pos_dist_constr.groupby(['Chromosome','End','Distance'], sort=False)['constr'].min().reset_index()

    dist_and_constr = chr_pos_dist_min_constr[['Distance','constr']].values
 
    # need to feed both constraint and distance to bspline function
    dist_basis = gen_spline_basis(dist_and_constr[:,0], order=3, knots_raw = [0, 1, int(2e3), int(1e5), 20083767], intercept=False)
    constr_basis = gen_spline_basis(dist_and_constr[:,1], order=3, knots_raw = [0, 0.01, 0.3, 0.99, 27.50007], intercept=False)

    r = np.repeat(range(dist_basis.shape[1]), constr_basis.shape[1])   
    c = np.tile(range(constr_basis.shape[1]), dist_basis.shape[1])

    spl_basis = dist_basis[:,c] * constr_basis[:,r]

    return spl_basis



# takes in numpy array
# returns n x 2 array with means and standard deviations
def get_normalisation (X):
    import numpy as np

    means = np.mean(X, axis=0)
    std_devs = np.std(X, axis=0)

    out = np.c_[means, std_devs]

    return out



# takes in n x 2 numpy array with means in 1st column, sd in 2nd
# returns normalised numpy array
def normalise_data (X, scaling):
    import numpy as np

    if X.dtype != 'float':
        X = X.astype('float')
    
    for i in range(X.shape[1]):
        X[:,i] = (X[:,i] - scaling[i,0])/ scaling[i,1]

    const_cols = (scaling[:,1] == 0)

    return X[:, ~const_cols]




