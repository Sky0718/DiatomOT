import json
import numpy as np
import pandas as pd
SOURCE = None
TERMINALS = None
OUT = None
POLICIES = ['B0', 'B1', 'B2', 'B3', 'BAV', 'BAV_uniform']
FRACTIONS = [0.1, 0.25, 0.5]
GROUP = ['selection_policy', 'verification_fraction']
IDS = ['catchment_id', 'site_id', 'analysis_id']
NBOOT = 5000
SEED = 20261111
METRICS = {
    'class_recovery'                   : ('recovered', 'ambiguous'),
    'good_moderate_recovery'           : ('gm_recovered', 'gm_ambiguous'),
    'class_resolution'                 : ('class_resolved', None),
    'good_moderate_resolution'         : ('good_moderate_resolved', None),
    'mean_width'                       : ('tdi_width', None),
    'mean_width_reduction'             : ('tdi_width_reduction', None),
    'mean_actual_cost'                 : ('verified_valves', None),
    'mean_budget'                      : ('verification_budget', None),
    'mean_cost_initially_ambiguous'    : ('ambiguous_cost', 'ambiguous'),
    'mean_width_initially_ambiguous'   : ('ambiguous_width', 'ambiguous'),
    'mean_cost_initially_gm_ambiguous' : ('gm_ambiguous_cost', 'gm_ambiguous'),
    'recoveries_per_valve'             : ('recovered', 'verified_valves'),
    'initial_class_resolution'         : ('initial_class_resolved', None),
    'initial_good_moderate_resolution' : ('initial_good_moderate_resolved', None),
}

def evaluate(array, denominator, weights):
    numerator = weights @ array
    if (denominator is None):
        return (numerator)
    divisor = weights @ denominator
    return (np.divide(numerator, divisor, out = np.full_like(numerator, np.nan), where = divisor > 0))

def analyse(partition):
    columns = IDS + GROUP + ['mask_id', 'reveal_replicate', 'initial_class_resolved', 'class_resolved', 'initial_good_moderate_resolved', 'good_moderate_resolved', 'tdi_width', 'tdi_width_reduction', 'verified_valves', 'verification_budget', 'initial_unresolved_valves']
    if (SOURCE is not None) :
        files = []
        for name in SOURCE.names(partition + '/terminal/'):
            if (name.endswith('.parquet')):
                files.append(name)
        frames = []
        for name in files:
            frames.append(SOURCE.read_parquet(name, columns = columns))
        frame = pd.concat(frames, ignore_index = True)
    else                    :
        files = sorted((TERMINALS / partition / 'terminal').glob('part_*.parquet'))
        frames = []
        for path in files:
            frames.append(pd.read_parquet(path, columns = columns))
        frame = pd.concat(frames, ignore_index = True)
    if (partition == 'development') :
        expected = 1952
    else                            :
        expected = 3374
    assert (frame['analysis_id'].nunique() == expected)
    assert (len(frame) == expected * 5 * 5 * 6 * 3)
    assert (not frame.duplicated(['analysis_id', 'mask_id', 'reveal_replicate'] + GROUP).any())
    assert (set(frame['selection_policy']) == set(POLICIES))
    assert (set(frame['verification_fraction']) == set(FRACTIONS))
    assert (frame.groupby(['analysis_id'] + GROUP).size().eq(25).all())
    for column in ['initial_class_resolved', 'initial_good_moderate_resolved']:
        assert (frame.groupby(['analysis_id', 'mask_id'])[column].nunique().eq(1).all())
        frame[column] = frame[column].astype(float)
    frame['class_resolved'] = frame['class_resolved'].astype(float)
    frame['good_moderate_resolved'] = frame['good_moderate_resolved'].astype(float)
    assert ((frame['class_resolved'] >= frame['initial_class_resolved']).all())
    assert ((frame['good_moderate_resolved'] >= frame['initial_good_moderate_resolved']).all())
    frame['ambiguous'] = 1 - frame['initial_class_resolved']
    frame['gm_ambiguous'] = 1 - frame['initial_good_moderate_resolved']
    frame['recovered'] = frame['ambiguous'] * frame['class_resolved']
    frame['gm_recovered'] = frame['gm_ambiguous'] * frame['good_moderate_resolved']
    frame['ambiguous_cost'] = frame['ambiguous'] * frame['verified_valves']
    frame['ambiguous_width'] = frame['ambiguous'] * frame['tdi_width']
    frame['gm_ambiguous_cost'] = frame['gm_ambiguous'] * frame['verified_valves']
    numerator_columns = []
    for pair in METRICS.values():
        numerator_columns.append(pair[0])
    denominator_columns = []
    for pair in METRICS.values():
        if (pair[1] is not None):
            denominator_columns.append(pair[1])
    value_columns = sorted(set(numerator_columns + denominator_columns))
    records = frame.groupby(GROUP + IDS, observed = True)[value_columns].mean().reset_index()
    sites = records.groupby(GROUP + ['catchment_id', 'site_id'], observed = True)[value_columns].mean().reset_index()
    catchments = sites.groupby(GROUP + ['catchment_id'], observed = True)[value_columns].mean().reset_index()
    destination = OUT / partition
    destination.mkdir(parents = True, exist_ok = True)
    records.to_parquet(destination / 'record_estimates.parquet', index = False)
    sites.to_parquet(destination / 'site_estimates.parquet', index = False)
    catchments.to_csv(destination / 'catchment_basis.csv', index = False)
    ids = sorted(catchments['catchment_id'].unique())
    size = len(ids)
    if (partition == 'development') :
        expected_catchments = 144
    else                            :
        expected_catchments = 226
    assert (size == expected_catchments)
    rng = np.random.default_rng(SEED)
    weights = rng.multinomial(size, np.full(size, 1 / size), size = NBOOT).astype(float) / size
    np.savez_compressed(destination / 'paired_bootstrap_weights.npz', weights = weights, catchment_ids = np.asarray(ids))
    point_weights = np.full((1, size), 1 / size)
    summaries = []
    distributions = {}
    points = {}
    for group_entry in catchments.groupby(GROUP, observed = True):
        group_key = group_entry[0]
        policy = group_key[0]
        fraction = group_key[1]
        data = group_entry[1]
        data = data.set_index('catchment_id').loc[ids]
        for metric_entry in METRICS.items():
            metric = metric_entry[0]
            metric_columns = metric_entry[1]
            numerator = metric_columns[0]
            denominator = metric_columns[1]
            array = data[numerator].to_numpy(dtype = float)
            if (denominator is None) :
                divisor = None
            else                     :
                divisor = data[denominator].to_numpy(dtype = float)
            samples = evaluate(array, divisor, weights)
            estimate = float(evaluate(array, divisor, point_weights)[0])
            assert (np.isfinite(samples).all())
            interval = np.quantile(samples, [0.025, 0.975])
            lower = interval[0]
            upper = interval[1]
            key = (policy, fraction, metric)
            distributions[key] = samples
            points[key] = estimate
            if (divisor is not None) :
                supported_catchments = int((divisor > 0).sum())
            else                     :
                supported_catchments = size
            summaries.append({
                'partition'             : partition,
                'selection_policy'      : policy,
                'verification_fraction' : fraction,
                'metric'                : metric,
                'estimate'              : estimate,
                'ci_lower'              : lower,
                'ci_upper'              : upper,
                'catchments'            : size,
                'supported_catchments'  : supported_catchments,
            })
    summary = pd.DataFrame(summaries)
    summary.to_csv(destination / 'summary.csv', index = False)
    errors = []
    for policy in POLICIES:
        for fraction in FRACTIONS:
            r0 = points[policy, fraction, 'initial_class_resolution']
            recovery = points[policy, fraction, 'class_recovery']
            rall = points[policy, fraction, 'class_resolution']
            errors.append(abs(rall - (r0 + (1 - r0) * recovery)))
    assert (max(errors) < 1e-12)
    baseline_path = OUT / 'development_baseline_lock.json'
    if (partition == 'development') :
        assert (not baseline_path.exists()), 'Do not replace an existing baseline lock'
        tie_order = {
            'B3' : 0,
            'B2' : 1,
            'B1' : 2,
            'B0' : 3,
        }
        candidates = []
        for policy in ['B0', 'B1', 'B2', 'B3']:
            candidates.append({
                'policy'      : policy,
                'recovery'    : points[policy, 0.25, 'class_recovery'],
                'actual_cost' : points[policy, 0.25, 'mean_actual_cost'],
            })

        def baseline_priority(row):
            return ((-row['recovery'], row['actual_cost'], tie_order[row['policy']]))
        candidates.sort(key = baseline_priority)
        selected = candidates[0]['policy']
        lock = {
            'selected_baseline'       : selected,
            'primary_budget_fraction' : 0.25,
            'selection_population'    : 'development',
            'selection_metric'        : 'class_recovery',
            'candidates'              : candidates,
            'conditional_weighting'   : 'ratio under the original equal-catchment/equal-site/equal-record/equal-mask/equal-reveal population measure',
        }
        baseline_path.write_text(json.dumps(lock, indent = 2), encoding = 'utf-8')
    else                            :
        lock = json.loads(baseline_path.read_text())
        selected = lock['selected_baseline']
    comparisons = []
    pairs = []
    for policy in POLICIES:
        if (policy != selected):
            pairs.append((policy, selected))
    for pair in [('BAV', 'B3'), ('BAV', 'BAV_uniform')]:
        if (pair not in pairs):
            pairs.append(pair)
    for policy, baseline in pairs:
        for fraction in FRACTIONS:
            for metric in METRICS:
                differences = distributions[policy, fraction, metric] - distributions[baseline, fraction, metric]
                interval = np.quantile(differences, [0.025, 0.975])
                lower = interval[0]
                upper = interval[1]
                comparisons.append({
                    'partition'             : partition,
                    'selection_policy'      : policy,
                    'comparator'            : baseline,
                    'verification_fraction' : fraction,
                    'metric'                : metric,
                    'difference'            : points[policy, fraction, metric] - points[baseline, fraction, metric],
                    'ci_lower'              : lower,
                    'ci_upper'              : upper,
                    'paired_catchments'     : size,
                })
    pd.DataFrame(comparisons).to_csv(destination / 'paired_comparisons.csv', index = False)
    bootstrap_arrays = {}
    for distribution_entry in distributions.items():
        distribution_key = distribution_entry[0]
        p = distribution_key[0]
        f = distribution_key[1]
        m = distribution_key[2]
        value = distribution_entry[1]
        bootstrap_arrays[f'{p}_{f:.2f}_{m}'] = value
    np.savez_compressed(destination / 'bootstrap_distributions.npz', ** bootstrap_arrays)
