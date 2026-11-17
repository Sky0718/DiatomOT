import numpy as np
from verification_engine import POLICIES
from verification_engine import seed_for
from verification_engine import simulate
from verification_engine import index_value
from verification_engine import class_index
from verification_engine import quality_ratio

def make_seeds(record, mask, policy):
    decisions = int(np.ceil(mask['counts'] / 5).sum())
    seeds = np.zeros((max(decisions, 1), len(record['genera'])), dtype = np.uint64)
    for decision in range(decisions):
        for genus_index, genus in enumerate(record['genera']):
            seeds[decision, genus_index] = seed_for(f"{record['analysis_id']}|{mask['mask_id']}|{policy}|{decision}|{genus}")
    return (seeds)

def make_queue(record, mask, replicate):
    scores = []
    taxa = []
    offsets = []
    rows = []
    for genus_index, genus in enumerate(record['genera']):
        offsets.append(len(scores))
        members = np.where(record['taxon_genus'] == genus_index)[0]
        repeated = np.repeat(members, mask['hidden'][members])
        seed = seed_for(f"{record['analysis_id']}|{mask['mask_id']}|{replicate}|{genus}|DiatomOT-reveal-v1")
        queue = np.random.default_rng(seed).permutation(repeated)
        genus_scores = record['taxon_scores'][queue]
        genus_taxa = record['taxon_codes'][queue]
        scores.extend(genus_scores.tolist())
        taxa.extend(genus_taxa.tolist())
        if (len(queue) > 0):
            rows.append({
                'analysis_id'      : record['analysis_id'],
                'mask_id'          : mask['mask_id'],
                'reveal_replicate' : replicate,
                'genus'            : genus,
                'queue_seed'       : str(seed),
                'taxon_codes'      : genus_taxa.tolist(),
                'scores_tenths'    : genus_scores.tolist(),
            })
    return ((np.asarray(scores, dtype = np.int64), np.asarray(taxa, dtype = np.int64), np.asarray(offsets, dtype = np.int64), rows))

def process_record(record, policies, fractions):
    terminals = []
    actions = []
    queues = []
    N = record['N']
    divisor = record['divisor']
    truth = index_value(record['truth_num'], N)
    for mask in record['masks']:
        initial_lower = index_value(mask['low_num'], N)
        initial_upper = index_value(mask['high_num'], N)
        initial_count = class_index(initial_lower, divisor) - class_index(initial_upper, divisor) + 1
        initial_gm = (quality_ratio(initial_lower, divisor) >= 0.6) == (quality_ratio(initial_upper, divisor) >= 0.6)
        hidden = int(mask['hidden'].sum())
        max_budget = min(hidden, int(np.floor(max(fractions) * hidden + 0.5)))
        seed_arrays = {}
        for policy in policies:
            if (policy in ('BAV', 'BAV_uniform') and initial_count > 1) :
                seed_arrays[policy] = make_seeds(record, mask, policy)
            else                                                        :
                seed_arrays[policy] = np.zeros((max(max_budget, 1), len(record['genera'])), dtype = np.uint64)
        for replicate in range(1, 6):
            queue = make_queue(record, mask, replicate)
            queue_scores = queue[0]
            queue_taxa = queue[1]
            queue_offsets = queue[2]
            queue_rows = queue[3]
            queues.extend(queue_rows)
            for policy in policies:
                if (policy == 'BAV_uniform') :
                    prior = record['uniform_cdf']
                else                         :
                    prior = record['prior_cdf']
                random_seed = np.uint64(seed_for(f"{record['analysis_id']}|{mask['mask_id']}|{replicate}|B0|random-verification"))
                for fraction in fractions:
                    budget = min(hidden, int(np.floor(fraction * hidden + 0.5)))
                    result = simulate(POLICIES.index(policy), mask['counts'], record['low_support'], record['high_support'], prior, record['support_scores'], record['score_lengths'], mask['low_num'], mask['high_num'], N, divisor, record['thresholds'], budget, record['fixed_priority'], seed_arrays[policy], random_seed, queue_scores, queue_taxa, queue_offsets, record['truth_num'])
                    low_num = result[0]
                    high_num = result[1]
                    actual = result[2]
                    decisions = result[3]
                    batches = result[4]
                    batch_seeds = result[5]
                    lower = index_value(low_num, N)
                    upper = index_value(high_num, N)
                    eqr_lower = quality_ratio(upper, divisor)
                    eqr_upper = quality_ratio(lower, divisor)
                    count = class_index(lower, divisor) - class_index(upper, divisor) + 1
                    common = {
                        'analysis_id'           : record['analysis_id'],
                        'mask_id'               : mask['mask_id'],
                        'reveal_replicate'      : replicate,
                        'selection_policy'      : policy,
                        'verification_fraction' : fraction,
                    }
                    sequence = f"{record['analysis_id']}:{mask['mask_id']}:{replicate}:{policy}:{fraction:.2f}"
                    terminal = {}
                    terminal.update(common)
                    terminal['site_id'] = record['site_id']
                    terminal['catchment_id'] = record['catchment_id']
                    terminal['verification_budget'] = budget
                    terminal['N'] = N
                    terminal['initial_unresolved_valves'] = hidden
                    terminal['verified_valves'] = actual
                    terminal['remaining_unresolved_valves'] = hidden - actual
                    terminal['initial_tdi_lower'] = initial_lower
                    terminal['initial_tdi_upper'] = initial_upper
                    terminal['initial_tdi_width'] = initial_upper - initial_lower
                    terminal['tdi_lower'] = lower
                    terminal['tdi_upper'] = upper
                    terminal['tdi_width'] = upper - lower
                    terminal['tdi_width_reduction'] = initial_upper - initial_lower - upper + lower
                    terminal['initial_possible_class_count'] = initial_count
                    terminal['possible_class_count'] = count
                    terminal['initial_class_resolved'] = initial_count == 1
                    terminal['class_resolved'] = count == 1
                    terminal['initial_good_moderate_resolved'] = initial_gm
                    terminal['good_moderate_resolved'] = (eqr_lower >= 0.6) == (eqr_upper >= 0.6)
                    terminal['eqr_lower'] = eqr_lower
                    terminal['eqr_upper'] = eqr_upper
                    terminal['true_tdi'] = truth
                    terminal['true_class'] = class_index(truth, divisor)
                    if (count == 1) :
                        assessment = 'class_resolved'
                    else            :
                        assessment = 'class_ambiguous'
                    terminal['qa_status'] = assessment
                    terminal['verification_sequence'] = sequence
                    terminal['actions_count'] = decisions
                    terminals.append(terminal)
                    for decision in range(decisions):
                        row = batches[decision]
                        genus_index = int(row[0])
                        queue_start = int(row[1])
                        taken = int(row[2])
                        global_start = queue_offsets[genus_index] + queue_start
                        actual_taxa = queue_taxa[global_start:global_start + taken]
                        action = {}
                        action.update(common)
                        action['verification_budget'] = budget
                        action['decision_index'] = decision
                        action['genus'] = record['genera'][genus_index]
                        action['queue_start'] = queue_start
                        action['queue_end'] = queue_start + taken
                        action['requested_batch'] = int(row[8])
                        action['verified_in_batch'] = taken
                        action['verified_cumulative'] = int(row[9])
                        action['score_sum_tenths'] = int(row[3])
                        action['tdi_lower'] = index_value(row[4], N)
                        action['tdi_upper'] = index_value(row[5], N)
                        action['expected_resolution_gain_per_valve'] = row[6]
                        action['expected_depth_reduction_per_valve'] = row[7]
                        action['candidate_seed'] = str(int(batch_seeds[decision]))
                        action['verified_taxon_codes'] = actual_taxa.tolist()
                        action['verification_sequence'] = sequence
                        actions.append(action)
    return ((terminals, actions, queues))
