from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import argparse
import pandas as pd
from verification_engine import POLICIES
from verification_engine import BUDGET_FRACTIONS
from verification_replay import process_record
from verification_inputs import Resource
from verification_inputs import reconstruct
from verification_inputs import output_directory

def run_chunk(records):
    terminals = []
    actions = []
    queues = []
    for record in records:
        record_results = process_record(record, list(POLICIES), list(BUDGET_FRACTIONS))
        t = record_results[0]
        a = record_results[1]
        q = record_results[2]
        terminals.extend(t)
        actions.extend(a)
        queues.extend(q)
    return ((terminals, actions, queues))

def save_chunk(destination, index, result):
    counts = {}
    for label, rows in zip(('terminal', 'actions', 'reveal_queues'), result):
        counts[label] = len(rows)
        if (rows):
            path = destination / label / f'part_{index:04d}.parquet'
            pd.DataFrame(rows).to_parquet(path, index = False, compression = 'zstd')
    return (counts)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--resource', type = Path, required = True)
    parser.add_argument('--partition', choices = ['development', 'confirmation'], required = True)
    parser.add_argument('--out', type = Path, required = True)
    parser.add_argument('--record-ids', nargs = '+', type = int)
    parser.add_argument('--workers', type = int, default = 4)
    args = parser.parse_args()
    if (args.workers < 1):
        raise ValueError('Workers must be positive')
    destination = output_directory(args.out.resolve() / args.partition)
    source = Resource(args.resource)
    try:
        records = reconstruct(source, args.record_ids, args.partition)
    finally:
        source.close()
    if (not records):
        raise ValueError('No selected records')
    values = list(records.values())
    chunks = []
    for index in range(0, len(values), 16):
        chunks.append(values[index:index + 16])
    for label in ['terminal', 'actions', 'reveal_queues']:
        (destination / label).mkdir(parents = True, exist_ok = False)
    totals = {
        'terminal'      : 0,
        'actions'       : 0,
        'reveal_queues' : 0,
    }
    if (args.workers == 1) :
        for index, chunk in enumerate(chunks):
            counts = save_chunk(destination, index, run_chunk(chunk))
            for key in totals:
                totals[key] += counts[key]
    else                   :
        with ProcessPoolExecutor(max_workers = args.workers) as pool:
            for index, result in enumerate(pool.map(run_chunk, chunks)):
                counts = save_chunk(destination, index, result)
                for key in totals:
                    totals[key] += counts[key]
    if (totals['terminal'] != len(records) * 5 * 5 * 6 * 3):
        raise ValueError('Terminal count mismatch')

if (__name__ == '__main__'):
    main()
