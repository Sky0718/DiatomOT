from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import argparse
import math
import os
if (__name__ == '__main__'):
    for variable in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS']:
        os.environ[variable] = '1'
import pandas as pd
from budget_ceiling import record_minimum
from verification_inputs import Resource
from verification_inputs import reconstruct
from verification_inputs import output_directory

def budget_count(hidden_count, fraction):
    if (hidden_count < 0 or not 0 <= fraction <= 1):
        raise ValueError('Invalid count or budget fraction')
    return (min(hidden_count, int(math.floor(fraction * hidden_count + 0.5))))

def one_chunk(records):
    rows = []
    for record in records:
        for mask in record['masks']:
            hidden = int(sum(mask['counts']))
            cap = budget_count(hidden, 0.5)
            for replicate in range(1, 6):
                row = record_minimum(record, mask, replicate, cap, same_batch = False)
                batch = record_minimum(record, mask, replicate, cap, same_batch = True)
                row['minimum_reveals_same_batch'] = batch['minimum_reveals']
                for name in ['partition', 'site_id', 'catchment_id']:
                    row[name] = record[name]
                row['hidden_valves'] = hidden
                for label, fraction in [('10', 0.1), ('25', 0.25), ('50', 0.5)]:
                    row['budget_' + label] = budget_count(hidden, fraction)
                if (row['minimum_reveals_same_batch'] >= 0):
                    if (not 0 <= row['minimum_reveals'] <= row['minimum_reveals_same_batch']):
                        raise ValueError('Single-valve ceiling cost exceeds same-batch ceiling cost')
                rows.append(row)
    return (rows)

def main():
    parser = argparse.ArgumentParser(description = 'Calculate exact verification budget ceilings.')
    parser.add_argument('--resource', type = Path, required = True, help = 'Online Resource 16 ZIP or extracted folder')
    parser.add_argument('--partition', choices = ('development', 'confirmation', 'both'), default = 'both')
    parser.add_argument('--record-ids', nargs = '+', type = int)
    parser.add_argument('--out', type = Path, required = True)
    parser.add_argument('--workers', type = int, default = 4)
    args = parser.parse_args()
    if (args.workers < 1):
        raise ValueError('Workers must be positive')
    destination = output_directory(args.out.resolve())
    source = Resource(args.resource)
    try:
        if (args.partition == 'both') :
            partition = None
        else                          :
            partition = args.partition
        records = reconstruct(source, args.record_ids, partition)
    finally:
        source.close()

    def record_order(r):
        return ((r['partition'], r['analysis_id']))
    values = sorted(records.values(), key = record_order)
    if (not values):
        raise ValueError('No selected records')
    chunks = []
    for start in range(0, len(values), 12):
        chunks.append(values[start:start + 12])
    destination.mkdir(parents = True, exist_ok = True)
    if (args.workers > 1) :
        pool = ProcessPoolExecutor(max_workers = args.workers)
    else                  :
        pool = None
    if (pool) :
        results = pool.map(one_chunk, chunks)
    else      :
        results = map(one_chunk, chunks)
    frames = []
    try:
        for number, rows in enumerate(results):
            frame = pd.DataFrame(rows)
            frame.to_parquet(destination / f'part_{number:04d}.parquet', index = False)
            frames.append(frame)
    finally:
        if (pool):
            pool.shutdown()
    frame = pd.concat(frames, ignore_index = True)
    keys = ['partition', 'analysis_id', 'mask_id', 'reveal_replicate']
    if (frame.duplicated(keys).any() or len(frame) != len(values) * 25):
        raise ValueError('Duplicate or incomplete ceiling states')
    frame.to_parquet(destination / 'oracle_states.parquet', index = False)

if (__name__ == '__main__'):
    main()
