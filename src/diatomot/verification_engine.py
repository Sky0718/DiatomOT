import hashlib
import numpy as np
from numba import njit
POLICIES = ('B0', 'B1', 'B2', 'B3', 'BAV', 'BAV_uniform')
BUDGET_FRACTIONS = (0.1, 0.25, 0.5)

def seed_for(value):
    return (int.from_bytes(hashlib.sha256(value.encode('utf-8')).digest()[:8], 'big'))

@njit(cache = False)
def next_uniform(state):
    state = state + np.uint64(11400714819323198485)
    value = state
    value = (value ^ value >> np.uint64(30)) * np.uint64(13787848793156543929)
    value = (value ^ value >> np.uint64(27)) * np.uint64(10723151780598845931)
    value = value ^ value >> np.uint64(31)
    uniform = float(value >> np.uint64(11)) / 9007199254740992.0
    return ((state, uniform))

@njit(cache = False)
def index_value(score_numerator, total):
    return (25.0 * (score_numerator / 10.0) / total - 25.0)

@njit(cache = False)
def quality_ratio(tdi, divisor):
    return (min(1.0, 0.8 * (100.0 - tdi) / divisor))

@njit(cache = False)
def class_index(tdi, divisor):
    q = quality_ratio(tdi, divisor)
    value = 0
    for threshold in [0.2, 0.4, 0.6, 0.8]:
        if (q >= threshold):
            value += 1
    return (value)

@njit(cache = False)
def boundary_depth(lower, upper, thresholds):
    total = 0.0
    for threshold in thresholds:
        if (lower < threshold < upper):
            total += min(threshold - lower, upper - threshold)
    return (total)

@njit(cache = False)
def choose_genus(policy, remaining, low_support, high_support, prior_cdf, support_scores, score_lengths, low_num, high_num, total, divisor, thresholds, budget, fixed_priority, candidate_seeds, random_state):
    selected = -1
    if (policy == 0):
        random_draw = next_uniform(random_state)
        random_state = random_draw[0]
        uniform = random_draw[1]
        target = int(uniform * remaining.sum())
        cumulative = 0
        for genus in range(len(remaining)):
            cumulative += remaining[genus]
            if (target < cumulative):
                return ((genus, 1, random_state, 0.0, 0.0, np.uint64(0)))
    if (policy == 2):
        for genus in fixed_priority:
            if (genus >= 0 and remaining[genus] > 0):
                return ((genus, min(5, remaining[genus], budget), random_state, 0.0, 0.0, np.uint64(0)))
    if (policy == 1):
        maximum = -1
        for genus in range(len(remaining)):
            if (remaining[genus] > maximum):
                maximum = remaining[genus]
                selected = genus
        return ((selected, min(5, remaining[selected], budget), random_state, 0.0, 0.0, np.uint64(0)))
    if (policy == 2 or policy == 3):
        best_span = -1
        best_count = -1
        for genus in range(len(remaining)):
            if (remaining[genus] == 0):
                continue
            span = high_support[genus] - low_support[genus]
            if (span > best_span or (span == best_span and remaining[genus] > best_count)):
                selected = genus
                best_span = span
                best_count = remaining[genus]
        return ((selected, min(5, remaining[selected], budget), random_state, 0.0, 0.0, np.uint64(0)))
    current_lower = index_value(low_num, total)
    current_upper = index_value(high_num, total)
    current_depth = boundary_depth(current_lower, current_upper, thresholds)
    best_gain = -1.0
    best_depth = -1e+100
    best_span = -1
    best_h = 0
    best_seed = np.uint64(0)
    for genus in range(len(remaining)):
        if (remaining[genus] == 0):
            continue
        h = min(5, remaining[genus], budget)
        state = candidate_seeds[genus]
        resolutions = 0
        depth_sum = 0.0
        for draw in range(64):
            hypothetical_sum = 0
            for valve in range(h):
                candidate_draw = next_uniform(state)
                state = candidate_draw[0]
                uniform = candidate_draw[1]
                score_index = 0
                while (score_index < score_lengths[genus] - 1 and uniform >= prior_cdf[genus, score_index]):
                    score_index += 1
                hypothetical_sum += support_scores[genus, score_index]
            next_low = low_num + hypothetical_sum - h * low_support[genus]
            next_high = high_num - h * high_support[genus] + hypothetical_sum
            lower = index_value(next_low, total)
            upper = index_value(next_high, total)
            if (class_index(lower, divisor) == class_index(upper, divisor)):
                resolutions += 1
            depth_sum += boundary_depth(lower, upper, thresholds)
        gain = resolutions / (64.0 * h)
        depth_gain = (current_depth - depth_sum / 64.0) / h
        span = high_support[genus] - low_support[genus]
        replace = gain > best_gain
        if (gain == best_gain):
            replace = depth_gain > best_depth
            if (depth_gain == best_depth):
                replace = span > best_span
        if (replace):
            selected = genus
            best_gain = gain
            best_depth = depth_gain
            best_span = span
            best_h = h
            best_seed = candidate_seeds[genus]
    return ((selected, best_h, random_state, best_gain, best_depth, best_seed))

@njit(cache = False)
def simulate(policy, initial_counts, low_support, high_support, prior_cdf, support_scores, score_lengths, low_num, high_num, total, divisor, thresholds, budget, fixed_priority, seeds, random_seed, queue_scores, queue_taxa, queue_offsets, truth_num):
    remaining = initial_counts.copy()
    consumed = np.zeros(len(remaining), dtype = np.int64)
    actions = np.zeros((max(budget, 1), 10), dtype = np.float64)
    action_seeds = np.zeros(max(budget, 1), dtype = np.uint64)
    actual = 0
    decisions = 0
    random_state = np.uint64(random_seed)
    initial_lower = index_value(low_num, total)
    initial_upper = index_value(high_num, total)
    truth = index_value(truth_num, total)
    initial_class_count = class_index(initial_lower, divisor) - class_index(initial_upper, divisor) + 1
    while (actual < budget):
        lower = index_value(low_num, total)
        upper = index_value(high_num, total)
        if (class_index(lower, divisor) == class_index(upper, divisor)):
            break
        if (remaining.sum() == 0):
            break
        selection = choose_genus(policy, remaining, low_support, high_support, prior_cdf, support_scores, score_lengths, low_num, high_num, total, divisor, thresholds, budget - actual, fixed_priority, seeds[decisions], random_state)
        genus = selection[0]
        h = selection[1]
        random_state = selection[2]
        gain = selection[3]
        depth_gain = selection[4]
        action_seed = selection[5]
        assert (genus >= 0 and h > 0)
        start_offset = consumed[genus]
        sum_score = 0
        for step in range(h):
            position = queue_offsets[genus] + consumed[genus]
            score = queue_scores[position]
            assert (remaining[genus] > 0)
            previous_low = low_num
            previous_high = high_num
            low_num += score - low_support[genus]
            high_num -= high_support[genus] - score
            consumed[genus] += 1
            remaining[genus] -= 1
            actual += 1
            sum_score += score
            assert (low_num >= previous_low)
            assert (high_num <= previous_high)
            assert (low_num <= truth_num <= high_num)
            lower = index_value(low_num, total)
            upper = index_value(high_num, total)
            if (class_index(lower, divisor) == class_index(upper, divisor)):
                assert (class_index(lower, divisor) == class_index(truth, divisor))
                break
        actions[decisions, 0] = genus
        actions[decisions, 1] = start_offset
        actions[decisions, 2] = consumed[genus] - start_offset
        actions[decisions, 3] = sum_score
        actions[decisions, 4] = low_num
        actions[decisions, 5] = high_num
        actions[decisions, 6] = gain
        actions[decisions, 7] = depth_gain
        actions[decisions, 8] = h
        actions[decisions, 9] = actual
        action_seeds[decisions] = action_seed
        decisions += 1
    assert (0 <= actual <= budget)
    assert (remaining.sum() + actual == initial_counts.sum())
    lower = index_value(low_num, total)
    upper = index_value(high_num, total)
    low_class = class_index(lower, divisor)
    high_class = class_index(upper, divisor)
    assert (initial_lower - 1e-10 <= lower <= truth + 1e-10)
    assert (truth - 1e-10 <= upper <= initial_upper + 1e-10)
    assert (low_class - high_class + 1 <= initial_class_count)
    return ((low_num, high_num, actual, decisions, actions[:decisions], action_seeds[:decisions]))
