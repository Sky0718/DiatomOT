from pathlib import Path
import argparse
import verification_summary
from verification_inputs import Resource
from verification_inputs import output_directory

def main():
    parser = argparse.ArgumentParser()
    source = parser.add_mutually_exclusive_group(required = True)
    source.add_argument('--resource', type = Path)
    source.add_argument('--terminals', type = Path)
    parser.add_argument('--partition', choices = ['development', 'confirmation'], required = True)
    parser.add_argument('--out', type = Path, required = True)
    args = parser.parse_args()
    verification_summary.OUT = args.out.resolve()
    output_directory(verification_summary.OUT / args.partition)
    if (args.partition == 'development'):
        if ((verification_summary.OUT / 'development_baseline_lock.json').exists()):
            raise FileExistsError('Choose a new output directory for the development comparison.')
    if (args.partition == 'confirmation'):
        if (not (verification_summary.OUT / 'development_baseline_lock.json').is_file()):
            raise FileNotFoundError('Aggregate the development population first')
    if (args.resource is not None) :
        verification_summary.SOURCE = Resource(args.resource)
    else                           :
        verification_summary.TERMINALS = args.terminals.resolve()
    try:
        if (verification_summary.SOURCE is not None):
            verification_summary.SOURCE.check_inputs()
        verification_summary.analyse(args.partition)
    finally:
        if (verification_summary.SOURCE is not None):
            verification_summary.SOURCE.close()

if (__name__ == '__main__'):
    main()
