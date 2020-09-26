from scars import scars_queries, scars_codata, scars_fit

def predict_and_calibrate_platt_iso(model_calibration, Z, sequence_features):
    import numpy as np

    flank = int((sequence_features.shape[1]-1)/2)

    Z_ref = np.nansum(Z * sequence_features[:, flank, :], axis=1)
    neg_inf_class_scores = ~np.isfinite(Z_ref) #capture which class scores are neg infinity
    allWeightOnReference = np.where(Z_ref == 0)

    pred_ref_clbt = np.empty(shape=(Z_ref.shape[0]))
    pred_ref_clbt[~neg_inf_class_scores] = model_calibration.transform(Z_ref[~neg_inf_class_scores])
    pred_ref_clbt[neg_inf_class_scores] = model_calibration.transform([-1000]) #replace with clipped lowest value

    #temporary matrix storing inf for the references alleles, so that softmax can be applied to alts only
    tmp = np.inf * sequence_features[:,flank,:]
    tmp[np.isnan(tmp)] = 0

    class_scores_alts = Z * (1-sequence_features[:, flank,:]) - tmp
    class_scores_alts[np.where(np.isnan(class_scores_alts))] = -np.inf # situation occurs when alt is assigned p=1 by cnn
    class_preds_alts = np.apply_along_axis(softmax, 1, class_scores_alts)

    class_preds_alts_scaled = np.apply_along_axis(lambda x: x*(1-pred_ref_clbt), 0, class_preds_alts)
    class_preds_ref = np.apply_along_axis(lambda x: x*pred_ref_clbt, 0, sequence_features[:, flank,:])

    class_preds = class_preds_alts_scaled + class_preds_ref
    class_preds[allWeightOnReference] = sequence_features[allWeightOnReference, flank]

    return class_preds



def softmax(v):
    import numpy as np
    val = np.exp(v)/np.sum(np.exp(v))
    return val



def get_expected(seqs, codata, model_bal, model_cal, pop_size):
    import numpy as np
    import pybedtools

    flank = seqs.shape[1]//2

    preds_bal = model_bal.predict([seqs, codata])
    preds = predict_and_calibrate_platt_iso(model_cal, np.log(preds_bal), seqs)
    
    expected = (1-np.sum(preds * seqs[:, flank, :], axis=1)) * 2 * pop_size

    return expected


#snvs_pr must be annotated with ac column (int) if provided
def get_observed(gr, snvs_pr):
    import pandas as pd
    import numpy as np
    import pyranges as pr

    gr_spl = gr.tile(1)
    index = pd.Series(range(len(gr_spl)), name="id")
    gr_spl = gr_spl.insert(index)

    snvHits = gr_spl.join(snvs_pr, suffix="_snv").as_df()
    anyHits = (snvHits.shape[0] != 0)

    observed = np.zeros(len(gr_spl))
    if anyHits:
        ac_by_id = snvHits.groupby(['id']).agg({'ac': "sum"})
        observed[ac_by_id.index] = ac_by_id['ac']

    return observed


def match_scores_to_ranges(gr, gr_reliable_mgd, obs_reliable_flat, exp_reliable_flat):
    import pyranges as pr
    import pandas as pd
    import numpy as np
    import collections

    gr_reliable_spl = gr_reliable_mgd.tile(1)

    o = pd.Series(data = obs_reliable_flat, name="observed")
    e = pd.Series(data = exp_reliable_flat, name="expected")

    gr_reliable_spl = gr_reliable_spl.insert(o)
    gr_reliable_spl = gr_reliable_spl.insert(e)

    index = pd.Series(range(len(gr)), name="id")
    gr = gr.insert(index) 

    hits = gr.join(gr_reliable_spl, suffix="_spl")
    out = hits.as_df().groupby(['id']).agg({'observed': "sum", 'expected': "sum"})

    countsByIndex = collections.Counter(hits.id)
    out['coverage'] = [countsByIndex[ix] if ix in countsByIndex.keys() else 0 for ix in out.index]    

    return out



def evaluate(gr, reliable_sites, reference_fasta, flank, model_bal, model_cal, pop_size, snvs, codata_args, codataNormalisation=None, split=False):
    import numpy as np

    gr_mgd = gr.merge()

    gr_reliable_mgd = gr_mgd.intersect(reliable_sites)

    obs = get_observed(gr_reliable_mgd, snvs)

    if len(codata_args) > 0: 
        assert codataNormalisation is not None, "normalisation required when providing arguments to read as codata"

        codata = scars_codata.get_codata(gr_reliable_mgd, *codata_args)
        codata = scars_codata.normalise_data(codata, codataNormalisation)
    else:
        codata = np.empty((len(obs),0))

    seqs = scars_queries.query_sequence(gr_reliable_mgd, reference_fasta, flank)
    scars_queries.correct_refs(gr_reliable_mgd, snvs, seqs)    
    scars_fit.anchor_on_AC(seqs)
    
    exp = get_expected(seqs, codata, model_bal, model_cal, pop_size)

    if not split:
        return [sum(obs), sum(exp), sum(obs)/sum(exp)]
    else:
        return match_scores_to_ranges(gr, gr_reliable_mgd, obs, exp)





