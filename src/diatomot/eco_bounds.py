import hashlib

import numpy as np

def stable_int(value : str) -> int:
    encoded = value.encode('utf-8')
    digest = hashlib.sha256(encoded).hexdigest()
    return (int(digest[:16], 16))

def split_catchment(value : str) -> str:
    if ((not value) or (value == 'MISSING_CATCHMENT')):
        return ('missing_catchment')
    encoded = ('DiatomOT-v2-catchment-20261105|' + value).encode('utf-8')
    digest = hashlib.sha256(encoded).hexdigest()
    bucket = (int(digest[:8], 16) % 100)
    if (bucket < 40):
        return ('development')
    return ('confirmation')

def make_mask(counts, indicators, keys, fraction, mechanism, seed = 1701, analysis_id = ''):
    supplied_counts = np.asarray(counts, dtype = float)
    scores = np.asarray(indicators, dtype = float)
    taxon_keys = np.asarray(keys, dtype = str)
    if ((supplied_counts.ndim != 1) or (supplied_counts.size == 0)):
        raise ValueError('Counts must be a nonempty one-dimensional array.')
    if ((scores.shape != supplied_counts.shape) or (taxon_keys.shape != supplied_counts.shape)):
        raise ValueError('Counts, scores and taxon keys must have the same shape.')
    if ((not np.isfinite(supplied_counts).all()) or (not np.isfinite(scores).all())):
        raise ValueError('Counts and scores must be finite.')
    if (np.any(supplied_counts != np.floor(supplied_counts))):
        raise ValueError('Counts must be integers.')
    if (np.any(supplied_counts < 0)):
        raise ValueError('Counts must be nonnegative.')
    count_vector = np.asarray(counts, dtype = np.int64)
    total = int(count_vector.sum())
    if ((total <= 0) or (not (0 <= fraction < 1))):
        raise ValueError('The count total must be positive and the fraction must be in [0, 1).')
    missing_count = int(np.floor((float(fraction) * total) + 0.5))
    if (mechanism == 'mcar') :
        seed_text = f'{analysis_id}|{fraction:.8f}|{seed}|DiatomOT-v2-mask'
        actual_seed = stable_int(seed_text)
        generator = np.random.default_rng(actual_seed)
        masked = generator.multivariate_hypergeometric(count_vector, missing_count)
    else                     :
        actual_seed = 0
        if (mechanism == 'rarest_first')           :
            priority = count_vector
        elif (mechanism == 'low_indicator_first')  :
            priority = scores
        elif (mechanism == 'high_indicator_first') :
            priority = (-scores)
        else                                       :
            raise ValueError(f'Unknown masking mechanism: {mechanism}')
        order = np.lexsort((taxon_keys, priority))
        masked = np.zeros_like(count_vector)
        remaining = missing_count
        for index in order:
            taken = min(int(count_vector[index]), remaining)
            masked[index] = taken
            remaining = (remaining - taken)
            if (remaining == 0):
                break
    if ((int(masked.sum()) != missing_count) or np.any(masked > count_vector)):
        raise ArithmeticError('The mask does not preserve the requested count.')
    return (masked, actual_seed)

def score_bounds(counts, masked, scores, lower_scores, upper_scores):
    count_vector = np.asarray(counts, dtype = float)
    masked_vector = np.asarray(masked, dtype = float)
    score_vector = np.asarray(scores, dtype = float)
    lower_vector = np.asarray(lower_scores, dtype = float)
    upper_vector = np.asarray(upper_scores, dtype = float)
    arrays = (count_vector, masked_vector, score_vector, lower_vector, upper_vector)
    if ((count_vector.ndim != 1) or (count_vector.size == 0)):
        raise ValueError('Counts must be a nonempty one-dimensional array.')
    for array in arrays:
        if (array.shape != count_vector.shape):
            raise ValueError('All inputs must have the same shape.')
        if (not np.isfinite(array).all()):
            raise ValueError('All inputs must be finite.')
    if ((np.any(count_vector < 0)) or (np.any(masked_vector < 0))):
        raise ValueError('Counts must be nonnegative.')
    if ((np.any(masked_vector > count_vector)) or (count_vector.sum() <= masked_vector.sum())):
        raise ValueError('Masked counts must leave a positive identified count.')
    if ((np.any(lower_vector > score_vector)) or (np.any(upper_vector < score_vector))):
        raise ValueError('A reference score is outside its genus support.')
    if ((np.any(lower_vector < 1)) or (np.any(upper_vector > 5))):
        raise ValueError('Score supports must lie within [1, 5].')
    total = count_vector.sum()
    masked_total = masked_vector.sum()
    known_score = np.dot((count_vector - masked_vector), score_vector)
    truth = ((25 * np.dot(count_vector, score_vector) / total) - 25)
    point = ((25 * known_score / (total - masked_total)) - 25)
    lower = ((25 * (known_score + np.dot(masked_vector, lower_vector)) / total) - 25)
    upper = ((25 * (known_score + np.dot(masked_vector, upper_vector)) / total) - 25)
    blind_lower = ((25 * (known_score + masked_total) / total) - 25)
    blind_upper = ((25 * (known_score + (5 * masked_total)) / total) - 25)
    width = (upper - lower)
    if (not ((blind_lower - 1e-9) <= lower <= (truth + 1e-9))):
        raise ArithmeticError('The lower endpoint is inconsistent with the reference score.')
    if (not ((truth - 1e-9) <= upper <= (blind_upper + 1e-9))):
        raise ArithmeticError('The upper endpoint is inconsistent with the reference score.')
    expected_width = (25 * np.dot(masked_vector, (upper_vector - lower_vector)) / total)
    if (abs(width - expected_width) >= 1e-8):
        raise ArithmeticError('The interval width is inconsistent with its contributions.')
    result = {
        'truth'       : float(truth),
        'point'       : float(point),
        'lower'       : float(lower),
        'upper'       : float(upper),
        'width'       : float(width),
        'blind_lower' : float(blind_lower),
        'blind_upper' : float(blind_upper),
        'blind_width' : float(blind_upper - blind_lower),
    }
    return (result)
