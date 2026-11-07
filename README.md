# DiatomOT

Analysis code for *Taxonomic concentration of diatom index uncertainty*.

## Research question

How much trophic information remains when a diatom assemblage retains valve counts and genus identities but some finer taxon labels are unresolved? These scripts calculate the Trophic Diatom Index values compatible with that information and quantify how individual genera contribute to the remaining range.

The analysis starts from Environment Agency river diatom counts and the DARLEQ3 coefficient dictionary. Recorded taxon labels are masked at specified fractions while their counts and correct genera are retained. Each unresolved count can take the lowest or highest coefficient assigned to its genus in the complete dictionary. Combining those limits with the identified counts gives lower and upper index endpoints. The scripts also calculate an unrestricted range, in which unresolved counts can take any coefficient from 1 to 5, and a point estimate based on identified counts.

This repository contains the numerical analysis code and instructions. Input datasets are obtained from their original providers. Numerical results, microscopy images and figure-generation code are excluded from this repository. All input, output and environment directories in the commands below are outside the repository.

## Repository structure

```text
GitHub Repo/
  README.md
  .gitignore
  docs/
    data.md
    usage.md
  environment/
    requirements.txt
  src/
    diatomot/
      prepare_inputs.py
      eco_bounds.py
      run_ecology.py
      summarize_ecology.py
      precision.py
      concentration.py
      coverage.py
```

## Requirements and data

Use Python 3.11 with NumPy 2.4.2, pandas 2.3.3 and pyreadr 0.5.9, as specified in [requirements.txt](environment/requirements.txt). pyreadr reads the coefficient dictionary directly; an R installation is not required.

Obtain `DIAT_OPEN_DATA_METRICS.csv`, `DIAT_OPEN_DATA_SITE.csv` and `DIAT_OPEN_DATA_TAXA.csv` from the Environment Agency [diatom archive](https://environment.data.gov.uk/ecology/explorer/downloads/DIAT_OPEN_DATA.zip), and `OPEN_DATA_TAXON_INFO.csv` from its [taxon archive](https://environment.data.gov.uk/ecology/explorer/downloads/OPEN_DATA_TAXON_INFO.zip). The fifth input is [darleq3_data.rda](https://raw.githubusercontent.com/nsj3/darleq3/e750ae0511b96b815aacb836d839f979bce980eb/data/darleq3_data.rda) from the specified DARLEQ3 commit.

Place these five files in `../diatomot-work/raw/`. The [data-source notes](docs/data.md) identify the study snapshot, original archive checksums and upstream licence terms. Taxon matching uses the original BIOSYS identifiers and DARLEQ3 NBS codes. Preserve the source filenames and columns.

## Installation

Run all commands from the repository root. Create an environment outside the repository:

```text
python -m venv ../diatomot-env
```

Activate it with `../diatomot-env/Scripts/Activate.ps1` in PowerShell, or `source ../diatomot-env/bin/activate` in a POSIX shell. Install the dependencies:

```text
python -m pip install -r environment/requirements.txt
```

## Analysis commands

Prepare the inputs, calculate the masking experiments, and derive the statistical summaries in this order:

```text
python -B src/diatomot/prepare_inputs.py --data-dir ../diatomot-work/raw --coefficients ../diatomot-work/raw/darleq3_data.rda --output ../diatomot-work/inputs
python -B src/diatomot/run_ecology.py --inputs ../diatomot-work/inputs --output ../diatomot-work/experiments --split both
python -B src/diatomot/summarize_ecology.py --experiments ../diatomot-work/experiments --output ../diatomot-work/summaries --split both
python -B src/diatomot/precision.py --experiments ../diatomot-work/experiments --output ../diatomot-work/precision
python -B src/diatomot/concentration.py --experiments ../diatomot-work/experiments --summaries ../diatomot-work/summaries --output ../diatomot-work/concentration
python -B src/diatomot/coverage.py --inputs ../diatomot-work/inputs --output ../diatomot-work/coverage
```

Each output directory must be absent or empty. For experiments and summaries, this requirement applies to each partition subdirectory. The `-B` option keeps Python bytecode caches out of the repository. Both partitions are needed for the precision and concentration analyses. Coverage requires only the prepared inputs. See the [usage reference](docs/usage.md) for arguments and output files.

## Analysis components

| Module | Role |
| --- | --- |
| `prepare_inputs.py` | Match taxon identifiers, establish sample eligibility and prepare count and coefficient tables |
| `eco_bounds.py` | Assign catchment partitions, construct masks and calculate feasible index endpoints |
| `run_ecology.py` | Apply random and selective masking under both coefficient definitions |
| `summarize_ecology.py` | Aggregate errors, interval widths, threshold attainment and derived sample classes |
| `precision.py` | Calculate weighted width quantiles and threshold probabilities |
| `concentration.py` | Summarize genus contributions, abundance shares and catchment-omission estimates |
| `coverage.py` | Describe records by scoring completeness, year, site, catchment and partition |

## Calculation conventions

The primary coefficient definition is `TDI5LM.orig`; `TDI5LM` supplies the alternative. The primary population contains records in which every counted valve has a coefficient. Additional summaries use all eligible scored subcommunities and the fully scored subset with at most 1,000 recorded valves. Catchment allocation uses a fixed hash, keeping observations from the same catchment together.

Masking fractions are 5%, 10%, 20%, 30% and 50%. Random loss samples valves without replacement, using five replicates with seeds 1701 through 1705. Selective loss removes the lowest-scoring, highest-scoring or least-abundant taxa first. Each selective condition uses one mask. Requested fractions are converted to integer valve counts before scoring.

Statistical summaries average masks within records, records within sites and sites within catchments, then give catchments equal weight. Confidence intervals use 5,000 catchment bootstrap resamples with seed 20261106. Paired quantities retain the same catchment draws. Width quantiles use the inverse weighted empirical distribution; attainment probabilities use thresholds of 5, 10, 15 and 20 index points with a numerical tolerance of 1e-9.

Concentration summaries include the combined contributions of *Navicula*, *Gomphonema* and *Nitzschia*. Overall shares divide summed catchment contributions by summed catchment totals. Masked-abundance shares use masked scoring counts as the denominator. Coverage descriptions give records equal weight. Derived sample classes use the specified TDI calibration and site alkalinity from 5 to 250 mg/L as CaCO3. The manuscript and Supplementary Information provide the scientific interpretation and full definitions.
