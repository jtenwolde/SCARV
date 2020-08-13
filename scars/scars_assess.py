from scars import scars_queries, scars_codata

def predict_and_calibrate_platt_iso(model_calibration, Z, sequence_features):
    import numpy as np

    flank = int((sequence_features.shape[1]-1)/2)

    Z_ref = np.nansum(Z * sequence_features[:, flank, :], axis=1)
    neg_inf_class_scores = ~np.isfinite(Z_ref) #capture which class scores are neg infinity

    pred_ref_clbt = np.empty(shape=(Z_ref.shape[0]))
    pred_ref_clbt[~neg_inf_class_scores] = model_calibration.transform(Z_ref[~neg_inf_class_scores])
    pred_ref_clbt[neg_inf_class_scores] = model_calibration.transform([-1000]) #replace with clipped lowest value

    #temporary matrix storing inf for the references alleles, so that softmax can be applied to alts only
    tmp = np.inf * sequence_features[:,flank,:]
    tmp[np.isnan(tmp)] = 0

    class_scores_alts = Z * (1-sequence_features[:, flank,:]) - tmp
    class_preds_alts = np.apply_along_axis(softmax, 1, class_scores_alts)

    class_preds_alts_scaled = np.apply_along_axis(lambda x: x*(1-pred_ref_clbt), 0, class_preds_alts)
    class_preds_ref = np.apply_along_axis(lambda x: x*pred_ref_clbt, 0, sequence_features[:, flank,:])

    class_preds = class_preds_alts_scaled + class_preds_ref
    return class_preds


def softmax(v):
    import numpy as np
    val = np.exp(v)/np.sum(np.exp(v))
    return val



def get_expected(bed, reference_fasta, flank, genome, codata, model_bal, model_cal, pop_size):
    import numpy as np
    import pybedtools

    bed_spl = pybedtools.BedTool().window_maker(b=bed, w=1)

    seqs = scars_queries.query_sequence(bed_spl, reference_fasta, flank, genome)
    
    if codata is not None:
        preds_bal = model_bal.predict([seqs, codata])
    else:
        preds_bal = model_bal.predict(seqs)

    preds = predict_and_calibrate_platt_iso(model_cal, np.log(preds_bal), seqs)
    expected = (1-np.sum(preds * seqs[:, flank, :], axis=1)) * 2 * pop_size

    return expected


def get_observed(bed, snvs_rare):
    import pandas as pd
    import pybedtools
    import numpy as np

    bed_spl = pybedtools.BedTool().window_maker(b=bed, w=1)
    hits = bed_spl.intersect(snvs_rare, wao=True)

    counts = []
    for record in hits:
        counts += [[record.fields[0], record.fields[2], record.fields[8]]]

    df = pd.DataFrame(counts, columns=['chr', 'pos', 'ac'])
    df_filled = df.replace('.', '0').astype({'chr': 'str','pos': 'int32', 'ac': 'int32'})
    observed = df_filled.groupby(['chr', 'pos']).sum()['ac'].to_numpy()

    return observed


def match_scores_to_ranges(bed, bed_reliable_mgd, obs_reliable_flat, exp_reliable_flat):
    import pyranges as pr
    import pybedtools
    import pandas as pd
    import numpy as np

    gr_reliable = scars_queries.bed_to_pr(bed_reliable_mgd)
    gr_reliable_spl = gr_reliable.tile(1)

    gr = scars_queries.bed_to_pr(bed)

    o = pd.Series(data = obs_reliable_flat, name="observed")
    e = pd.Series(data = exp_reliable_flat, name="expected")

    gr_reliable_spl = gr_reliable_spl.insert(o)
    gr_reliable_spl = gr_reliable_spl.insert(e)

    gr.id = np.arange(len(gr))
    out = gr.join(gr_reliable_spl, suffix="_spl").as_df().groupby(['id']).agg({'observed': "sum", 'expected': "sum"})
    out['obs_exp'] = out.apply(lambda row: row.observed / row.expected, axis=1)

    return out



def evaluate(bed, reliable_sites_bw, reference_fasta, flank, genome, model_bal, model_cal, pop_size, snvs, *codata_args, split=False):
    import numpy as np

    bed_mgd = bed.sort()\
                 .merge()

    bed_reliable_mgd = scars_queries.filter_sites(bed_mgd, reliable_sites_bw)

    obs = get_observed(bed_reliable_mgd, snvs)

    if len(codata_args) > 0: 
        codata = scars_codata.get_codata(bed_reliable_mgd, *codata_args)
    else:
        codata = None
    exp = get_expected(bed_reliable_mgd, reference_fasta, flank, genome, codata, model_bal, model_cal, pop_size)

    if not split:
        return [sum(obs), sum(exp), sum(obs)/sum(exp)]
    else:
        return match_scores_to_ranges(bed, bed_reliable_mgd, obs, exp)



def genome_wide_score_to_bigWig(genome, out_bw, snvs, reliable_sites_bw, reference_fasta, flank, model_bal, model_cal, pop_size, window_size, resolution, chr_list, *codata_args):
    import pyBigWig
    import pybedtools
    import numpy as np
    import pandas as pd

    bw = pyBigWig.open(out_bw, "w")
    header = [(chrom, genome[chrom][1]) for chrom in chr_list]
    bw.addHeader(header)
    for chrom in chr_list:
        l = [tuple([chrom] + list(genome[chrom]))]
        chrom_bed = pybedtools.BedTool(l)
        chrom_bed_reliable = scars_queries.filter_sites(chrom_bed, reliable_sites_bw)
        obs = get_observed(chrom_bed_reliable, snvs)
        
        codata = scars_codata.get_codata(chrom_bed_reliable, *codata_args)
        exp = get_expected(chrom_bed_reliable, reference_fasta, flank, genome, codata, model_bal, model_cal, pop_size)

        l2 = [(chrom, str(start), str(start + window_size)) for start in np.arange(0, genome[chrom][1] - window_size, resolution)]
        chrom_bed_windows = pybedtools.BedTool(l2)
        matched = match_scores_to_ranges(chrom_bed_windows, chrom_bed_reliable, obs, exp)

        vals = np.repeat(np.nan, len(l2))
        vals[matched.index] = matched['obs_exp']
        bw.addEntries(chrom, 0, values=vals, span=resolution, step=resolution)
    bw.close()

    return None


def query_score_from_bigWig (bed, bw_fn):
    import pyBigWig
    import pybedtools

    bw = pyBigWig.open(bw_fn)
    scores = [bw.values(r.chrom, r.end - 1, r.end)[0] for r in bed]
    bw.close()

    return scores







