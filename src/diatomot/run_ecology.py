import argparse
import gzip
import hashlib
import json
import platform
import re
from pathlib import Path

import numpy as np
import pandas as pd

from eco_bounds import make_mask
from eco_bounds import score_bounds
from eco_bounds import split_catchment

SEEDS = [1701, 1702, 1703, 1704, 1705]
FRACTIONS = [0.05, 0.10, 0.20, 0.30, 0.50]
DEFINITIONS = [('orig', 'TDI5LM.orig'), ('current', 'TDI5LM')]
INPUT_FILES = [
    'candidate_sample_metadata.csv',
    'candidate_sample_taxon_long.csv.gz',
    'darleq3_rda_coefficients.csv',
    'OPEN_DATA_TAXON_INFO.csv',
]

def sha(path):
    digest = hashlib.sha256()
    with (Path(path).open('rb')) as stream:
        while (True):
            block = stream.read(1048576)
            if (not block):
                break
            digest.update(block)
    return (digest.hexdigest())

def genus_token(name, valid_genera):
    match = re.match(r'^([A-Z][a-z]+)(?:\s|$)', str(name).strip())
    if (match is not None):
        token = match.group(1)
        if (token in valid_genera):
            return (token)
    return ('UNRESOLVED_GENUS')

def load_dictionary(input_dir):
    dictionary = pd.read_csv(input_dir / 'darleq3_rda_coefficients.csv', low_memory = False)
    taxon_info = pd.read_csv(input_dir / 'OPEN_DATA_TAXON_INFO.csv', dtype = str)
    genus_rows = taxon_info.loc[taxon_info['TAXON_RANK'].eq('Genus'), 'TAXON_NAME']
    valid_genera = set(genus_rows.str.strip())
    accepted_names = dictionary['ValidName'].fillna(dictionary['TaxonName'])
    assigned_genera = []
    for name in accepted_names:
        assigned_genera.append(genus_token(name, valid_genera))
    dictionary['genus'] = assigned_genera
    support_rows = []
    groups = dictionary.groupby('genus', sort = True)
    for group in groups:
        genus = group[0]
        entries = group[1]
        support = {
            'genus'           : genus,
            'dictionary_rows' : len(entries),
        }
        for definition in DEFINITIONS:
            prefix = definition[0]
            column = definition[1]
            values = pd.to_numeric(entries[column], errors = 'coerce').dropna()
            if ((genus == 'UNRESOLVED_GENUS') or (len(values) == 0)) :
                support[prefix + '_low'] = 1.0
                support[prefix + '_high'] = 5.0
            else                                                     :
                if ((not np.isfinite(values).all()) or (values.min() < 1) or (values.max() > 5)):
                    raise ValueError(f'Invalid coefficient support for genus {genus}.')
                support[prefix + '_low'] = float(values.min())
                support[prefix + '_high'] = float(values.max())
        support_rows.append(support)
    supports = pd.DataFrame(support_rows)
    keyed = dictionary.dropna(subset = ['NBSCode']).copy()
    if (keyed['NBSCode'].duplicated().any()):
        raise ValueError('The coefficient dictionary contains duplicate NBS keys.')
    return (keyed, supports)

def load_inputs(input_dir):
    input_dir = Path(input_dir).resolve()
    for name in INPUT_FILES:
        if (not (input_dir / name).is_file()):
            raise FileNotFoundError(input_dir / name)
    metadata = pd.read_csv(input_dir / 'candidate_sample_metadata.csv', low_memory = False)
    taxon_counts = pd.read_csv(input_dir / 'candidate_sample_taxon_long.csv.gz', low_memory = False)
    dictionary_result = load_dictionary(input_dir)
    dictionary = dictionary_result[0]
    supports = dictionary_result[1]
    taxon_counts = taxon_counts.merge(
        dictionary[['NBSCode', 'genus']],
        left_on = 'TAXON_LIST_ITEM_KEY',
        right_on = 'NBSCode',
        how = 'left',
        validate = 'many_to_one',
    )
    taxon_counts['genus'] = taxon_counts['genus'].fillna('UNRESOLVED_GENUS')
    taxon_counts = taxon_counts.merge(supports, on = 'genus', how = 'left', validate = 'many_to_one')
    for definition in DEFINITIONS:
        prefix = definition[0]
        taxon_counts[prefix + '_low'] = taxon_counts[prefix + '_low'].fillna(1.0)
        taxon_counts[prefix + '_high'] = taxon_counts[prefix + '_high'].fillna(5.0)
    original_scored = taxon_counts['tdi5lm_orig_indicator'].notna()
    current_scored = taxon_counts['tdi5lm_indicator'].notna()
    if (not original_scored.equals(current_scored)):
        raise ValueError('The coefficient definitions have different scored taxon sets.')
    scoring = taxon_counts.loc[original_scored].copy()
    counts = scoring['count'].to_numpy(dtype = float)
    if ((not np.isfinite(counts).all()) or np.any(counts < 0) or np.any(counts != np.floor(counts))):
        raise ValueError('Scored valve counts must be finite nonnegative integers.')
    if (scoring.duplicated(['ANALYSIS_ID', 'TAXON_LIST_ITEM_KEY']).any()):
        raise ValueError('An analysis contains duplicate taxon keys.')
    score_columns = [('orig', 'tdi5lm_orig_indicator'), ('current', 'tdi5lm_indicator')]
    for definition in score_columns:
        prefix = definition[0]
        column = definition[1]
        scores = scoring[column]
        inside_support = ((scoring[prefix + '_low'] <= scores) & (scores <= scoring[prefix + '_high']))
        if ((not np.isfinite(scores).all()) or (not inside_support.all())):
            raise ValueError('A scored taxon lies outside its dictionary genus support.')
    if ((metadata['ANALYSIS_ID'].duplicated().any()) or (metadata['SAMPLE_ID'].duplicated().any())):
        raise ValueError('Analysis and sample identifiers must be unique.')
    if (metadata.empty):
        raise ValueError('The prepared metadata contain no eligible samples.')
    metadata['catchment'] = metadata['CATCHMENT'].fillna('MISSING_CATCHMENT')
    metadata['catchment'] = metadata['catchment'].astype(str).str.strip()
    metadata['catchment'] = metadata['catchment'].replace('', 'MISSING_CATCHMENT')
    metadata['partition'] = metadata['catchment'].map(split_catchment)
    metadata['fully_scored'] = np.isclose(
        metadata['raw_count_sum'], metadata['scoring_count_sum'], atol = 0, rtol = 0,
    )
    metadata['N_gt_1000'] = (metadata['raw_count_sum'] > 1000)
    if (set(scoring['ANALYSIS_ID']) != set(metadata['ANALYSIS_ID'])):
        raise ValueError('Scored counts and metadata must contain the same analysis identifiers.')
    return (metadata, scoring, supports)

def output_directory(path):
    destination = Path(path).resolve()
    repository = Path(__file__).resolve().parents[2]
    if (destination.is_relative_to(repository)):
        raise ValueError('Write analysis outputs outside the code repository.')
    return (destination)

def require_empty_directory(path):
    destination = output_directory(path)
    if (destination.exists()):
        if (not destination.is_dir()):
            raise FileExistsError(destination)
        for existing in destination.iterdir():
            raise FileExistsError(f'Output directory is not empty: {destination}')
    return (destination)

def masking_regimes():
    regimes = []
    for seed in SEEDS:
        regimes.append(('mcar', seed))
    directional = ['rarest_first', 'low_indicator_first', 'high_indicator_first']
    for mechanism in directional:
        regimes.append((mechanism, 0))
    return (regimes)

def sample_identity(analysis_id, metadata, total):
    identity = {
        'ANALYSIS_ID'   : int(analysis_id),
        'SITE_ID'       : int(metadata['SITE_ID']),
        'catchment'     : metadata['catchment'],
        'fully_scored'  : bool(metadata['fully_scored']),
        'N'             : total,
        'raw_N'         : int(metadata['raw_count_sum']),
        'non_scoring_N' : (int(metadata['raw_count_sum']) - total),
        'N_gt_1000'     : bool(metadata['N_gt_1000']),
    }
    return (identity)

def contribution_records(identity, taxon_rows, masked, original_low, original_high, seed):
    contributions = pd.DataFrame({
        'genus'              : taxon_rows['genus'].to_numpy(str),
        'width_contribution' : (25 * masked * (original_high - original_low) / identity['N']),
        'masked_valves'      : masked,
    })
    records = []
    groups = contributions.groupby('genus', sort = False)
    for group in groups:
        genus = group[0]
        values = group[1]
        record = identity.copy()
        record['replicate_seed'] = seed
        record['genus'] = genus
        record['width_contribution'] = float(values['width_contribution'].sum())
        record['masked_valves'] = int(values['masked_valves'].sum())
        records.append(record)
    return (records)

def evaluate_sample(analysis_id, metadata, taxon_rows, regimes):
    taxon_rows = taxon_rows.sort_values('TAXON_LIST_ITEM_KEY')
    counts = taxon_rows['count'].to_numpy(np.int64)
    taxon_keys = taxon_rows['TAXON_LIST_ITEM_KEY'].to_numpy(str)
    original_scores = taxon_rows['tdi5lm_orig_indicator'].to_numpy(float)
    current_scores = taxon_rows['tdi5lm_indicator'].to_numpy(float)
    original_low = taxon_rows['orig_low'].to_numpy(float)
    original_high = taxon_rows['orig_high'].to_numpy(float)
    current_low = taxon_rows['current_low'].to_numpy(float)
    current_high = taxon_rows['current_high'].to_numpy(float)
    total = int(counts.sum())
    if (total != int(metadata['scoring_count_sum'])):
        raise ValueError(f'Scored count total differs from metadata for analysis {analysis_id}.')
    original_value = float((25 * np.dot(counts, original_scores) / total) - 25)
    current_value = float((25 * np.dot(counts, current_scores) / total) - 25)
    original_difference = abs(original_value - metadata['recomputed_TDI5LM_orig'])
    current_difference = abs(current_value - metadata['recomputed_TDI5LM'])
    maximum_error = max(original_difference, current_difference)
    if (maximum_error >= 1e-8):
        raise ValueError(f'Prepared reference index differs for analysis {analysis_id}.')
    identity = sample_identity(analysis_id, metadata, total)
    version = identity.copy()
    version['SAMPLE_DATE'] = metadata['SAMPLE_DATE']
    version['site_alkalinity'] = metadata['site_alkalinity_numeric']
    version['archive_TDI5LM'] = metadata['TDI5LM']
    version['orig'] = original_value
    version['current'] = current_value
    version['delta'] = (current_value - original_value)
    version['orig_archive_error'] = (original_value - metadata['TDI5LM'])
    version['current_archive_error'] = (current_value - metadata['TDI5LM'])
    mask_rows = []
    genus_rows = []
    for fraction in FRACTIONS:
        for regime in regimes:
            mechanism = regime[0]
            seed = regime[1]
            mask_result = make_mask(
                counts, original_scores, taxon_keys, fraction, mechanism, seed, str(analysis_id),
            )
            masked = mask_result[0]
            actual_seed = mask_result[1]
            original_bounds = score_bounds(counts, masked, original_scores, original_low, original_high)
            current_bounds = score_bounds(counts, masked, current_scores, current_low, current_high)
            record = identity.copy()
            record['fraction_requested'] = fraction
            record['mechanism'] = mechanism
            record['replicate_seed'] = seed
            record['actual_seed'] = str(actual_seed)
            record['m'] = int(masked.sum())
            record['fraction_actual'] = float(masked.sum() / total)
            for key in original_bounds:
                record['orig_' + key] = original_bounds[key]
            for key in current_bounds:
                record['current_' + key] = current_bounds[key]
            record['definition_union_lower'] = min(original_bounds['lower'], current_bounds['lower'])
            record['definition_union_upper'] = max(original_bounds['upper'], current_bounds['upper'])
            record['definition_union_hull_width'] = (record['definition_union_upper'] - record['definition_union_lower'])
            inner_lower = max(original_bounds['lower'], current_bounds['lower'])
            inner_upper = min(original_bounds['upper'], current_bounds['upper'])
            record['definition_union_gap'] = max(0.0, (inner_lower - inner_upper))
            mask_rows.append(record)
            if ((fraction == 0.2) and (mechanism == 'mcar')):
                contributions = contribution_records(
                    identity, taxon_rows, masked, original_low, original_high, seed,
                )
                genus_rows.extend(contributions)
    return (mask_rows, version, genus_rows, maximum_error)

def run_partition(split, metadata, scoring, supports, destination, input_dir, input_hashes):
    selected = metadata.loc[metadata['partition'].eq(split)].set_index('ANALYSIS_ID')
    if (selected.empty):
        raise ValueError(f'No eligible analyses in partition {split}.')
    destination = require_empty_directory(destination)
    destination.mkdir(parents = True, exist_ok = True)
    metadata.to_csv(destination / 'partition_manifest.csv', index = False)
    supports.to_csv(destination / 'genus_score_support.csv', index = False)
    selected_counts = scoring.loc[scoring['ANALYSIS_ID'].isin(selected.index)]
    groups = selected_counts.groupby('ANALYSIS_ID', sort = True)
    regimes = masking_regimes()
    pending = []
    versions = []
    contributions = []
    written = 0
    write_header = True
    maximum_error = 0.0
    result_path = (destination / 'mask_results.csv.gz')
    with (gzip.open(result_path, 'wt', encoding = 'utf-8', newline = '')) as stream:
        for group in groups:
            analysis_id = group[0]
            taxon_rows = group[1]
            sample_result = evaluate_sample(analysis_id, selected.loc[analysis_id], taxon_rows, regimes)
            pending.extend(sample_result[0])
            versions.append(sample_result[1])
            contributions.extend(sample_result[2])
            maximum_error = max(maximum_error, sample_result[3])
            if (len(pending) >= 10000):
                pd.DataFrame(pending).to_csv(stream, index = False, header = write_header)
                written = (written + len(pending))
                write_header = False
                pending.clear()
        if (pending):
            pd.DataFrame(pending).to_csv(stream, index = False, header = write_header)
            written = (written + len(pending))
    expected = (len(selected) * len(FRACTIONS) * len(regimes))
    if ((written != expected) or (len(versions) != len(selected))):
        raise ArithmeticError('The experiment produced an unexpected number of records.')
    pd.DataFrame(versions).to_csv(destination / 'coefficient_definition_results.csv', index = False)
    pd.DataFrame(contributions).to_csv(destination / 'primary_genus_contributions.csv.gz', index = False)
    for name in INPUT_FILES:
        if (sha(input_dir / name) != input_hashes[name]):
            raise RuntimeError(f'Input file changed during the run: {name}')
    output_hashes = {}
    for path in sorted(destination.iterdir()):
        if (path.is_file()):
            output_hashes[path.name] = sha(path)
    receipt = {
        'status'                                     : 'COMPLETE',
        'split'                                      : split,
        'analyses'                                   : len(versions),
        'case_rows'                                  : written,
        'expected_case_rows'                         : expected,
        'fully_scored_analyses'                      : int(selected['fully_scored'].sum()),
        'sample_reference_max_absolute_disagreement' : maximum_error,
        'bound_invariant_failures'                   : 0,
        'python'                                     : platform.python_version(),
        'numpy'                                      : np.__version__,
        'pandas'                                     : pd.__version__,
        'fractions'                                  : FRACTIONS,
        'replicate_seeds'                            : SEEDS,
        'inputs'                                     : input_hashes,
        'outputs'                                    : output_hashes,
    }
    receipt_path = (destination / 'run_receipt.json')
    receipt_path.write_text((json.dumps(receipt, indent = 2) + '\n'), encoding = 'utf-8')
    return (receipt)

def main():
    parser = argparse.ArgumentParser(description = 'Calculate genus-constrained trophic-index intervals.')
    parser.add_argument('--inputs', type = Path, required = True)
    parser.add_argument('--output', type = Path, required = True)
    parser.add_argument('--split', choices = ['development', 'confirmation', 'both'], default = 'both')
    arguments = parser.parse_args()
    input_dir = arguments.inputs.resolve()
    output = output_directory(arguments.output)
    splits = ['development', 'confirmation']
    if (arguments.split != 'both'):
        splits = [arguments.split]
    for split in splits:
        require_empty_directory(output / split)
    input_hashes = {}
    for name in INPUT_FILES:
        input_hashes[name] = sha(input_dir / name)
    loaded = load_inputs(input_dir)
    metadata = loaded[0]
    scoring = loaded[1]
    supports = loaded[2]
    for split in splits:
        run_partition(split, metadata, scoring, supports, (output / split), input_dir, input_hashes)

if (__name__ == '__main__'):
    main()
