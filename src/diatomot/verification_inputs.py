from io import BytesIO
from pathlib import Path
from zipfile import ZipFile
import hashlib
import json
import numpy as np
import pandas as pd
from verification_engine import class_index
from verification_engine import index_value
INPUT_NAMES = ('development_prior.csv', 'hidden_taxon_counts.parquet', 'initial_states.csv.gz', 'records.csv', 'taxon_codes.csv')

class Resource:

    def __init__(self, location):
        self.path = Path(location).resolve()
        self.archive = None
        if (self.path.is_file())      :
            self.archive = ZipFile(self.path)
        elif (not self.path.is_dir()) :
            raise FileNotFoundError(self.path)
        self.metadata = json.loads(self.read_bytes('metadata.json'))

    def read_bytes(self, name):
        if (self.archive is not None):
            return (self.archive.read(name))
        target = (self.path / name).resolve()
        if (not target.is_relative_to(self.path)):
            raise ValueError(name)
        return (target.read_bytes())

    def read_csv(self, name, ** options):
        return (pd.read_csv(BytesIO(self.read_bytes(name)), float_precision = 'round_trip', ** options))

    def read_parquet(self, name, ** options):
        return (pd.read_parquet(BytesIO(self.read_bytes(name)), ** options))

    def names(self, prefix):
        if (self.archive is not None):
            members = []
            for name in self.archive.namelist():
                if (name.startswith(prefix)):
                    members.append(name)
            return (sorted(members))
        folder = (self.path / prefix).resolve()
        members = []
        for path in folder.rglob('*'):
            if (path.is_file()):
                members.append(path.relative_to(self.path).as_posix())
        return (sorted(members))

    def check_inputs(self):
        expected = self.metadata['files_sha256']
        for name in INPUT_NAMES:
            member = 'inputs/' + name
            digest = hashlib.sha256(self.read_bytes(member)).hexdigest()
            if (digest != expected[member]):
                raise ValueError('Input checksum mismatch: ' + member)

    def close(self):
        if (self.archive is not None):
            self.archive.close()

def reconstruct(resource, record_ids = None, partition = None):
    resource.check_inputs()
    metadata = resource.read_csv('inputs/records.csv')
    prior_table = resource.read_csv('inputs/development_prior.csv')
    taxon_table = resource.read_csv('inputs/taxon_codes.csv', dtype = {'TAXON_LIST_ITEM_KEY' : str})
    initial = resource.read_csv('inputs/initial_states.csv.gz', compression = 'gzip', dtype = {'mask_seed' : str})
    hidden = resource.read_parquet('inputs/hidden_taxon_counts.parquet')
    if (partition is not None):
        metadata = metadata.loc[metadata['partition'] == partition]
    if (record_ids is not None):
        chosen = set(record_ids)
        metadata = metadata.loc[metadata['analysis_id'].isin(chosen)]
        if (set(metadata['analysis_id']) != chosen):
            raise ValueError('Requested records are absent from the chosen population')
    metadata = metadata.sort_values('analysis_id').set_index('analysis_id')
    hidden = hidden.loc[hidden['analysis_id'].isin(metadata.index)]
    initial = initial.set_index(['analysis_id', 'mask_id'])
    taxon_table = taxon_table.set_index('taxon_code')
    support = {}
    priors = {}
    for genus, group in prior_table.groupby('genus', sort = True):
        group = group.sort_values('score')
        scores = np.rint(group['score'].to_numpy() * 10).astype(np.int64)
        counts = group['development_valves'].to_numpy(dtype = np.float64)
        support[genus] = scores
        priors[genus] = (counts + 1.0) / (counts.sum() + len(scores))
    records = {}
    for aid, all_masks in hidden.groupby('analysis_id', sort = True):
        meta = metadata.loc[aid]
        first = all_masks.loc[all_masks['mask_id'] == 1701].sort_values('taxon_code')
        codes = first['taxon_code'].to_numpy(dtype = np.int64)
        n = first['total_valves'].to_numpy(dtype = np.int64)
        s10 = np.rint(first['score'].to_numpy() * 10).astype(np.int64)
        genera = sorted(first['genus'].unique())
        genus_index = dict(zip(genera, range(len(genera))))
        genus_codes = []
        for name in first['genus']:
            genus_codes.append(genus_index[name])
        g = np.asarray(genus_codes, dtype = np.int64)
        lookup = taxon_table.loc[codes]
        if (not np.array_equal(s10, lookup['score_tenths'].to_numpy())):
            raise ValueError('Taxon score mismatch')
        if (not np.array_equal(first['genus'].to_numpy(), lookup['genus'].to_numpy())):
            raise ValueError('Taxon genus mismatch')
        N = int(n.sum())
        if (N != int(meta.N)):
            raise ValueError('Record count mismatch')
        divisor = float(100 - 9.933 * np.exp(0.81 * np.log10(meta.alkalinity)))
        thresholds = 100 - np.asarray([0.2, 0.4, 0.6, 0.8]) * divisor / 0.8
        lower_support = []
        for genus in genera:
            lower_support.append(support[genus].min())
        lo = np.asarray(lower_support, dtype = np.int64)
        upper_support = []
        for genus in genera:
            upper_support.append(support[genus].max())
        hi = np.asarray(upper_support, dtype = np.int64)
        support_lengths = []
        for genus in genera:
            support_lengths.append(len(support[genus]))
        support_size = max(support_lengths)
        scores = np.zeros((len(genera), support_size), dtype = np.int64)
        cdf = np.ones((len(genera), support_size), dtype = np.float64)
        uniform_cdf = cdf.copy()
        lengths = np.zeros(len(genera), dtype = np.int64)
        for index, genus in enumerate(genera):
            k = len(support[genus])
            lengths[index] = k
            scores[index, :k] = support[genus]
            cdf[index, :k] = np.cumsum(priors[genus])
            uniform_cdf[index, :k] = np.arange(1, k + 1) / k
            cdf[index, k - 1] = 1.0
        truth_num = int(np.dot(n, s10))
        masks = []
        for mask_id in [1701, 1702, 1703, 1704, 1705]:
            part = all_masks.loc[all_masks['mask_id'] == mask_id].sort_values('taxon_code')
            if (not np.array_equal(codes, part['taxon_code'].to_numpy())):
                raise ValueError('Inconsistent taxon ordering across masks')
            masked = part['hidden_valves'].to_numpy(dtype = np.int64)
            counts = np.bincount(g, weights = masked, minlength = len(genera)).astype(np.int64)
            known_num = int(np.dot(n - masked, s10))
            low_num = known_num + int(np.dot(counts, lo))
            high_num = known_num + int(np.dot(counts, hi))
            stored = initial.loc[aid, mask_id]
            lower = index_value(low_num, N)
            upper = index_value(high_num, N)
            truth = index_value(truth_num, N)
            if (not np.allclose([lower, upper, truth], [stored.tdi_lower, stored.tdi_upper, stored.true_tdi], rtol = 0, atol = 1e-12)):
                raise ValueError('Initial assessment mismatch')
            if (int(masked.sum()) != int(stored.initial_unresolved_valves)):
                raise ValueError('Initial hidden count mismatch')
            resolved = class_index(lower, divisor) == class_index(upper, divisor)
            if (resolved != bool(stored.initial_class_resolved)):
                raise ValueError('Initial class mismatch')
            masks.append({
                'mask_id'   : mask_id,
                'mask_seed' : int(stored.mask_seed),
                'hidden'    : masked,
                'counts'    : counts,
                'low_num'   : low_num,
                'high_num'  : high_num,
            })
        priority = []
        for name in ['Navicula', 'Gomphonema', 'Nitzschia']:
            priority.append(genus_index.get(name, -1))
        records[int(aid)] = {
            'analysis_id'    : int(aid),
            'site_id'        : int(meta.site_id),
            'catchment_id'   : str(meta.catchment_id),
            'partition'      : meta.partition,
            'alkalinity'     : float(meta.alkalinity),
            'N'              : N,
            'genera'         : genera,
            'taxon_genus'    : g,
            'taxon_codes'    : codes,
            'taxon_scores'   : s10,
            'low_support'    : lo,
            'high_support'   : hi,
            'prior_cdf'      : cdf,
            'uniform_cdf'    : uniform_cdf,
            'support_scores' : scores,
            'score_lengths'  : lengths,
            'divisor'        : divisor,
            'thresholds'     : thresholds,
            'truth_num'      : truth_num,
            'fixed_priority' : np.asarray(priority, dtype = np.int64),
            'masks'          : masks,
        }
    return (records)

def output_directory(location):
    destination = Path(location).resolve()
    repository = Path(__file__).resolve().parents[2]
    if (destination.is_relative_to(repository)):
        raise ValueError('Output directories must be outside the repository.')
    if (destination.exists()):
        if (not destination.is_dir() or any(destination.iterdir())):
            raise FileExistsError('Output directory must be absent or empty: ' + str(destination))
    return (destination)
