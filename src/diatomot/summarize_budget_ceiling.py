from pathlib import Path
from io import BytesIO
from zipfile import ZipFile
import argparse
import os
if (__name__ == '__main__'):
    for variable in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS']:
        os.environ[variable] = '1'
import numpy as np
import pandas as pd
from verification_inputs import Resource
from verification_inputs import output_directory
KEYS = ['analysis_id', 'mask_id', 'reveal_replicate']
SEED = 20261111
NBOOT = 5000

def summary_draws(frame, numerator, denominator, cluster_ids, bootstrap):
    weighted = pd.DataFrame({
        'cluster_id' : frame.catchment_id,
        'top'        : frame.state_weight * numerator,
        'bottom'     : frame.state_weight * denominator,
    })
    cluster = weighted.groupby('cluster_id')[['top', 'bottom']].sum().reindex(cluster_ids, fill_value = 0.0)
    cluster_arrays = (cluster.top.to_numpy(), cluster.bottom.to_numpy())
    top = cluster_arrays[0]
    bottom = cluster_arrays[1]
    boot_bottom = bootstrap @ bottom
    draws = np.divide(bootstrap @ top, boot_bottom, out = np.full(len(boot_bottom), np.nan), where = boot_bottom > 0)
    estimate = float(top.sum() / bottom.sum())
    interval = np.nanquantile(draws, [0.025, 0.975])
    lower = interval[0]
    upper = interval[1]
    return ({
        'estimate'           : estimate,
        'ci_lower'           : float(lower),
        'ci_upper'           : float(upper),
        'denominator_weight' : float(bottom.sum()),
        'valid_bootstraps'   : int(np.isfinite(draws).sum()),
    })

def read_ceiling(path):
    if (path.is_dir()):
        path = path / 'oracle_states.parquet'
    if (path.suffix.lower() == '.zip')       :
        with ZipFile(path) as archive:
            frame = pd.read_csv(BytesIO(archive.read('oracle_state_requirements_and_costs.csv')), float_precision = 'round_trip')
    elif (path.suffix.lower() == '.parquet') :
        frame = pd.read_parquet(path)
    else                                     :
        frame = pd.read_csv(path, float_precision = 'round_trip')
    return (frame.rename(columns = {'minimum_reveals_single_valve' : 'minimum_reveals'}))

def assemble_states(resource, oracle, partition):
    columns = KEYS + ['catchment_id', 'site_id', 'selection_policy', 'verification_fraction', 'verification_budget', 'initial_class_resolved', 'class_resolved']
    files = []
    for name in resource.names(partition + '/terminal/'):
        if (name.endswith('.parquet')):
            files.append(name)
    frames = []
    for name in files:
        frames.append(resource.read_parquet(name, columns = columns))
    frame = pd.concat(frames, ignore_index = True)
    if (partition == 'development') :
        expected = 1952
    else                            :
        expected = 3374
    if (frame.analysis_id.nunique() != expected or len(frame) != expected * 450):
        raise ValueError('Population summary requires the complete archived partition')
    if (frame.duplicated(KEYS + ['selection_policy', 'verification_fraction']).any()):
        raise ValueError('Duplicate verification states')
    states = frame.loc[frame.selection_policy == 'BAV'].drop(columns = ['selection_policy', 'class_resolved']).copy()
    for policy in ['BAV', 'B3', 'BAV_uniform']:
        values = frame.loc[frame.selection_policy == policy, KEYS + ['verification_fraction', 'class_resolved']]
        states = states.merge(values.rename(columns = {'class_resolved' : policy}), on = KEYS + ['verification_fraction'], validate = 'one_to_one')
    subset = oracle.loc[oracle.partition == partition, KEYS + ['initial_resolved', 'minimum_reveals', 'minimum_reveals_same_batch']]
    if (len(subset) != expected * 25 or subset.duplicated(KEYS).any()):
        raise ValueError('Incomplete or duplicate ceiling population')
    states = states.merge(subset, on = KEYS, validate = 'many_to_one')
    if (not (states.initial_resolved.astype(bool) == states.initial_class_resolved.astype(bool)).all()):
        raise ValueError('Initial class states differ')
    meta = states[['catchment_id', 'site_id', 'analysis_id']].drop_duplicates()
    ncatchments = meta.catchment_id.nunique()
    nsites = meta.groupby('catchment_id').site_id.nunique()
    nrecords = meta.groupby(['catchment_id', 'site_id']).analysis_id.transform('count')
    meta['state_weight'] = 1.0 / (ncatchments * meta.catchment_id.map(nsites) * nrecords * 25)
    states = states.merge(meta, on = ['catchment_id', 'site_id', 'analysis_id'], validate = 'many_to_one')
    states['oracle_same_batch'] = (states.minimum_reveals_same_batch >= 0) & (states.minimum_reveals_same_batch <= states.verification_budget)
    states['oracle_unit'] = (states.minimum_reveals >= 0) & (states.minimum_reveals <= states.verification_budget)
    for policy in ['BAV', 'B3', 'BAV_uniform']:
        if ((states[policy].astype(bool) & ~states.oracle_same_batch).any()):
            raise ValueError(policy + ' exceeds the same-batch ceiling')
    states['ambiguous'] = (~states.initial_resolved.astype(bool)).astype(float)
    return (states)

def main():
    parser = argparse.ArgumentParser(description = 'Summarize verification budget ceilings.')
    parser.add_argument('--resource', type = Path, required = True, help = 'Online Resource 16')
    parser.add_argument('--ceiling', type = Path, required = True, help = 'Online Resource 17 ZIP or recomputed oracle_states.parquet')
    parser.add_argument('--out', type = Path, required = True)
    args = parser.parse_args()
    args.out = output_directory(args.out)
    oracle = read_ceiling(args.ceiling)
    source = Resource(args.resource)
    summaries = []
    try:
        source.check_inputs()
        for partition in ['development', 'confirmation']:
            states = assemble_states(source, oracle, partition)
            ids = sorted(states.catchment_id.unique())
            size = len(ids)
            weights = np.random.default_rng(SEED).multinomial(size, np.full(size, 1 / size), size = NBOOT).astype(float) / size
            for fraction, group in states.groupby('verification_fraction'):
                if (not np.isclose(group.state_weight.sum(), 1.0, atol = 1e-10)):
                    raise ValueError('State weights do not sum to one')
                ambiguity = group.ambiguous.to_numpy()
                metrics = {}
                for policy in ['BAV', 'B3', 'BAV_uniform', 'oracle_same_batch', 'oracle_unit']:
                    metrics[policy + '_recovery'] = group[policy].astype(float).to_numpy()
                metrics['same_batch_minus_BAV'] = metrics['oracle_same_batch_recovery'] - metrics['BAV_recovery']
                metrics['unit_minus_same_batch'] = metrics['oracle_unit_recovery'] - metrics['oracle_same_batch_recovery']
                for metric, outcome in metrics.items():
                    row = {
                        'partition'             : partition,
                        'verification_fraction' : fraction,
                        'metric'                : metric,
                    }
                    row.update(summary_draws(group, outcome * ambiguity, ambiguity, ids, weights))
                    summaries.append(row)
    finally:
        source.close()
    args.out.mkdir(parents = True, exist_ok = True)
    pd.DataFrame(summaries).to_csv(args.out / 'budget_ceiling_summary.csv', index = False)

if (__name__ == '__main__'):
    main()
