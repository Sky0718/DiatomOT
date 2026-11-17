import hashlib
import math
import numpy as np
from numba import njit

def class_from_numerator(numerator, total, divisor):
    value = 25.0 * (numerator / 10.0) / total - 25.0
    ratio = min(1.0, 0.8 * (100.0 - value) / divisor)
    result = 0
    for threshold in [0.2, 0.4, 0.6, 0.8]:
        if (ratio >= threshold):
            result += 1
    return (result)

def required_gains(low_num, high_num, truth_num, total, divisor):
    if (not low_num <= truth_num <= high_num):
        raise ValueError('Reference numerator is outside initial bounds')
    label = class_from_numerator(truth_num, total, divisor)
    left = int(low_num)
    right = int(truth_num)
    while (left < right):
        middle = (left + right) // 2
        if (class_from_numerator(middle, total, divisor) <= label) :
            right = middle
        else                                                       :
            left = middle + 1
    minimum_numerator = left
    left = int(truth_num)
    right = int(high_num)
    while (left < right):
        middle = (left + right + 1) // 2
        if (class_from_numerator(middle, total, divisor) >= label) :
            left = middle
        else                                                       :
            right = middle - 1
    maximum_numerator = left
    return ((minimum_numerator - low_num, high_num - maximum_numerator, label))

def frozen_queues(record, mask, replicate):
    result = []
    for group, genus in enumerate(record['genera']):
        members = np.flatnonzero(record['taxon_genus'] == group)
        repeated = np.repeat(members, mask['hidden'][members])
        payload = f"{record['analysis_id']}|{mask['mask_id']}|{replicate}|{genus}|DiatomOT-reveal-v1"
        seed = int.from_bytes(hashlib.sha256(payload.encode('utf-8')).digest()[:8], 'big')
        order = np.random.default_rng(seed).permutation(repeated)
        result.append(np.asarray(record['taxon_scores'][order], dtype = np.int64))
    return (result)

@njit(cache = False)
def prefix_cost_kernel(prefix_a, prefix_b, lengths, required_a, required_b, budget):
    if (required_a == 0 and required_b == 0):
        return (0)
    previous = np.full((budget + 1, required_a + 1), -1, dtype = np.int64)
    previous[0, 0] = 0
    used_max = 0
    best = budget + 1
    for group in range(len(lengths)):
        current = np.full((budget + 1, required_a + 1), -1, dtype = np.int64)
        length = lengths[group]
        for used in range(min(used_max, best - 1) + 1):
            for gain_a in range(required_a + 1):
                gain_b = previous[used, gain_a]
                if (gain_b < 0):
                    continue
                limit = min(length, budget - used, best - 1 - used)
                for take in range(limit + 1):
                    next_used = used + take
                    next_a = min(required_a, gain_a + prefix_a[group, take])
                    next_b = min(required_b, gain_b + prefix_b[group, take])
                    if (next_a == required_a and next_b == required_b) :
                        best = min(best, next_used)
                    elif (next_b > current[next_used, next_a])         :
                        current[next_used, next_a] = next_b
        previous = current
        used_max = min(budget, used_max + length)
    if (best <= budget):
        return (best)
    return (-1)

@njit(cache = False)
def batch_prefix_cost_kernel(prefix_a, prefix_b, lengths, full_lengths, required_a, required_b, budget):
    if (required_a == 0 and required_b == 0):
        return (0)
    previous = np.full((budget + 1, required_a + 1, 2), -1, dtype = np.int64)
    previous[0, 0, 0] = 0
    used_max = 0
    best = budget + 1
    for group in range(len(lengths)):
        current = np.full((budget + 1, required_a + 1, 2), -1, dtype = np.int64)
        length = lengths[group]
        for used in range(min(used_max, best - 1) + 1):
            for gain_a in range(required_a + 1):
                for partial in range(2):
                    gain_b = previous[used, gain_a, partial]
                    if (gain_b < 0):
                        continue
                    limit = min(length, budget - used, best - 1 - used)
                    for take in range(limit + 1):
                        next_partial = partial
                        if (take % 5 != 0 and take != full_lengths[group]):
                            next_partial += 1
                        if (next_partial > 1):
                            continue
                        next_used = used + take
                        next_a = min(required_a, gain_a + prefix_a[group, take])
                        next_b = min(required_b, gain_b + prefix_b[group, take])
                        if (next_a == required_a and next_b == required_b)       :
                            best = min(best, next_used)
                        elif (next_b > current[next_used, next_a, next_partial]) :
                            current[next_used, next_a, next_partial] = next_b
        previous = current
        used_max = min(budget, used_max + length)
    if (best <= budget):
        return (best)
    return (-1)

def minimum_prefix_cost(queues, lows, highs, required_a, required_b, budget, same_batch = False):
    if (required_a < 0 or required_b < 0 or budget < 0):
        raise ValueError('Requirements and budget must be nonnegative')
    budget = int(budget)
    step = 0
    for queue, low, high in zip(queues, lows, highs):
        if (np.any((queue < low) | (queue > high))):
            raise ValueError('A queue score is outside its support')
        for score in queue[:budget]:
            step = math.gcd(step, int(score - low))
            step = math.gcd(step, int(high - score))
    if (step == 0):
        if (required_a == 0 and required_b == 0):
            return (0)
        return (-1)
    required_a = (int(required_a) + step - 1) // step
    required_b = (int(required_b) + step - 1) // step
    prefix_a = np.zeros((len(queues), budget + 1), dtype = np.int64)
    prefix_b = np.zeros_like(prefix_a)
    lengths = np.zeros(len(queues), dtype = np.int64)
    full_lengths = np.zeros(len(queues), dtype = np.int64)
    for group, queue in enumerate(queues):
        full_lengths[group] = len(queue)
        length = min(len(queue), budget)
        if (highs[group] == lows[group]):
            length = 0
        lengths[group] = length
        for position in range(length):
            prefix_a[group, position + 1] = prefix_a[group, position] + (int(queue[position]) - int(lows[group])) // step
            prefix_b[group, position + 1] = prefix_b[group, position] + (int(highs[group]) - int(queue[position])) // step
    if (required_a > required_b):
        swapped_prefixes = (prefix_b, prefix_a)
        prefix_a = swapped_prefixes[0]
        prefix_b = swapped_prefixes[1]
        swapped_requirements = (required_b, required_a)
        required_a = swapped_requirements[0]
        required_b = swapped_requirements[1]
    if (same_batch):
        return (int(batch_prefix_cost_kernel(prefix_a, prefix_b, lengths, full_lengths, required_a, required_b, budget)))
    return (int(prefix_cost_kernel(prefix_a, prefix_b, lengths, required_a, required_b, budget)))

def record_minimum(record, mask, replicate, maximum_budget, same_batch = False):
    requirements = required_gains(int(mask['low_num']), int(mask['high_num']), int(record['truth_num']), int(record['N']), float(record['divisor']))
    required_a = requirements[0]
    required_b = requirements[1]
    label = requirements[2]
    queues = frozen_queues(record, mask, replicate)
    minimum = minimum_prefix_cost(queues, record['low_support'], record['high_support'], required_a, required_b, maximum_budget, same_batch = same_batch)
    return ({
        'analysis_id'       : int(record['analysis_id']),
        'mask_id'           : int(mask['mask_id']),
        'reveal_replicate'  : int(replicate),
        'initial_resolved'  : required_a == 0 and required_b == 0,
        'true_class'        : int(label),
        'lower_requirement' : int(required_a),
        'upper_requirement' : int(required_b),
        'minimum_reveals'   : minimum,
        'maximum_budget'    : int(maximum_budget),
        'same_batch'        : bool(same_batch),
    })

def prefix_witness(queues, lows, highs, required_a, required_b, budget):
    previous = {(0, 0, 0) : ()}
    for group, queue in enumerate(queues):
        current = {}
        gain_a = 0
        gain_b = 0
        options = [(0, 0, 0)]
        for take in range(1, min(len(queue), budget) + 1):
            gain_a += int(queue[take - 1]) - int(lows[group])
            gain_b += int(highs[group]) - int(queue[take - 1])
            options.append((take, gain_a, gain_b))
        for previous_entry in previous.items():
            previous_key = previous_entry[0]
            used = previous_key[0]
            old_a = previous_key[1]
            old_b = previous_key[2]
            allocation = previous_entry[1]
            for take, gain_a, gain_b in options:
                if (used + take > budget):
                    break
                key = (used + take, min(required_a, old_a + gain_a), min(required_b, old_b + gain_b))
                if (key not in current):
                    current[key] = allocation + (take,)
        previous = current
    solutions = []
    for final_entry in previous.items():
        final_key = final_entry[0]
        used = final_key[0]
        gain_a = final_key[1]
        gain_b = final_key[2]
        allocation = final_entry[1]
        if (gain_a >= required_a and gain_b >= required_b):
            solutions.append((used, allocation))
    if (solutions):
        return (min(solutions))
    return (None)
