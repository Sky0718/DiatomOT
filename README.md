# DiatomOT

Analysis code for *Taxonomic verification priorities and budget limits in river diatom monitoring*.

## Research question

Incomplete diatom counts raise two monitoring questions: whether the retained taxonomic information already supports an assessment, and which additional identifications would be most useful. DiatomOT calculates compatible Trophic Diatom Index (TDI) ranges, locates their genus contributions, and evaluates targeted identification against finite-budget recovery ceilings.

The first analysis starts from Environment Agency river diatom counts and the DARLEQ3 coefficient dictionary. Finer taxon labels are masked while counts and correct genera are retained. Each unresolved count can take the lowest or highest coefficient assigned to its genus in the complete dictionary. Combining those limits with the identified counts gives lower and upper index endpoints. The scripts also calculate an unrestricted range, in which unresolved counts can take any coefficient from 1 to 5, and a point estimate based on identified counts.

The verification analysis replays additional identification on the stored masked records. Boundary-aware verification (BAV) uses development-derived ecological score distributions to direct identification towards unresolved class boundaries. Four simple policies and a uniform-prior BAV implementation use the same stored masks and reveal queues. A retrospective calculation with advance knowledge of hidden scores determines the largest class recovery permitted by the specified queues, action rules and budgets.

This repository contains numerical analysis code and instructions. Raw inputs come from their original providers; verification inputs and archived numerical results accompany the article as Online Resources 16 and 17. Input data, numerical results, microscopy images and figure-generation code remain outside this repository. All data, output and environment paths in the commands below are external to it.

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
      verification_inputs.py
      verification_engine.py
      verification_replay.py
      verification_summary.py
      run_verification.py
      summarize_verification.py
      budget_ceiling.py
      run_budget_ceiling.py
      summarize_budget_ceiling.py
```

## Requirements and data

Use Python 3.11 and the versions in [requirements.txt](environment/requirements.txt): NumPy 2.4.2, pandas 2.3.3, pyreadr 0.5.9, Numba 0.64.0, llvmlite 0.46.0 and PyArrow 23.0.1. pyreadr reads the DARLEQ3 coefficient dictionary directly; an R installation is not required.

For the archive analyses, obtain `DIAT_OPEN_DATA_METRICS.csv`, `DIAT_OPEN_DATA_SITE.csv` and `DIAT_OPEN_DATA_TAXA.csv` from the Environment Agency [diatom archive](https://environment.data.gov.uk/ecology/explorer/downloads/DIAT_OPEN_DATA.zip), and `OPEN_DATA_TAXON_INFO.csv` from its [taxon archive](https://environment.data.gov.uk/ecology/explorer/downloads/OPEN_DATA_TAXON_INFO.zip). The fifth input is [darleq3_data.rda](https://raw.githubusercontent.com/nsj3/darleq3/e750ae0511b96b815aacb836d839f979bce980eb/data/darleq3_data.rda) from the specified DARLEQ3 commit. Place these five files in `../diatomot-work/raw/`.

For verification, place **Online Resource 16 — Complete targeted verification results** and **Online Resource 17 — Verification budget ceilings** in `../diatomot-work/resources/`. The examples call these archives `OR16.zip` and `OR17.zip`; substitute their actual filenames if different. OR16 includes the population tables, original masks, development score counts, terminal outcomes and reveal queues. OR17 contains the archived ceiling states and their summaries. The verification calculations can use these companion resources independently of the raw-archive preparation commands.

The [data-source notes](docs/data.md) describe the source snapshot, attribution, required members and how the two input routes relate. Preserve source filenames, identifiers and columns.

## Installation

Run commands from the repository root. Create an environment outside the repository:

```text
python -m venv ../diatomot-env
```

Activate it with `../diatomot-env/Scripts/Activate.ps1` in PowerShell, or `source ../diatomot-env/bin/activate` in a POSIX shell, then install the dependencies:

```text
python -m pip install -r environment/requirements.txt
```

## Archive analyses

Prepare the inputs, calculate the masking experiments and derive the statistical summaries in this order:

```text
python -B src/diatomot/prepare_inputs.py --data-dir ../diatomot-work/raw --coefficients ../diatomot-work/raw/darleq3_data.rda --output ../diatomot-work/inputs
python -B src/diatomot/run_ecology.py --inputs ../diatomot-work/inputs --output ../diatomot-work/experiments --split both
python -B src/diatomot/summarize_ecology.py --experiments ../diatomot-work/experiments --output ../diatomot-work/summaries --split both
python -B src/diatomot/precision.py --experiments ../diatomot-work/experiments --output ../diatomot-work/precision
python -B src/diatomot/concentration.py --experiments ../diatomot-work/experiments --summaries ../diatomot-work/summaries --output ../diatomot-work/concentration
python -B src/diatomot/coverage.py --inputs ../diatomot-work/inputs --output ../diatomot-work/coverage
```

Both partitions are needed for precision and concentration. Coverage requires only the prepared inputs. Each destination must be absent or empty; for partitioned analyses this applies to the requested partition subdirectory. The `-B` option keeps Python bytecode caches out of the repository.

## Verification and budget ceilings

The archived terminal outcomes can be summarized directly. Aggregate development first, then confirmation using the comparator selected in development:

```text
python -B src/diatomot/summarize_verification.py --resource ../diatomot-work/resources/OR16.zip --partition development --out ../diatomot-work/verification-statistics
python -B src/diatomot/summarize_verification.py --resource ../diatomot-work/resources/OR16.zip --partition confirmation --out ../diatomot-work/verification-statistics
python -B src/diatomot/summarize_budget_ceiling.py --resource ../diatomot-work/resources/OR16.zip --ceiling ../diatomot-work/resources/OR17.zip --out ../diatomot-work/ceiling-statistics
```

To calculate the verification trajectories and exact ceilings themselves, use:

```text
python -B src/diatomot/run_verification.py --resource ../diatomot-work/resources/OR16.zip --partition development --out ../diatomot-work/verification --workers 4
python -B src/diatomot/run_verification.py --resource ../diatomot-work/resources/OR16.zip --partition confirmation --out ../diatomot-work/verification --workers 4
python -B src/diatomot/run_budget_ceiling.py --resource ../diatomot-work/resources/OR16.zip --out ../diatomot-work/ceilings --workers 4
```

The full ceiling calculation searches 133,150 record-mask-queue states and can be computationally intensive. `--record-ids` restricts either calculation to selected monitoring records. Full-population summaries require all records, masks and reveal replicates. The [usage reference](docs/usage.md) gives the subset commands and commands for summarizing newly calculated outputs.

## Analysis components

| Modules | Role |
| --- | --- |
| `prepare_inputs.py`, `eco_bounds.py` | Match taxon identifiers, define eligible records, assign catchment partitions, construct masks and calculate index endpoints |
| `run_ecology.py`, `summarize_ecology.py` | Calculate random and selective masking results, population estimates and bootstrap intervals |
| `precision.py`, `concentration.py`, `coverage.py` | Summarize precision thresholds, genus contributions and monitoring-record coverage |
| `verification_inputs.py` | Read the companion input tables and reconstruct original masked records |
| `verification_engine.py`, `verification_replay.py` | Apply identification policies using shared reveal queues and deterministic candidate draws |
| `run_verification.py`, `verification_summary.py`, `summarize_verification.py` | Calculate trajectories, select the development comparator and summarize paired recovery and effort |
| `budget_ceiling.py`, `run_budget_ceiling.py`, `summarize_budget_ceiling.py` | Find exact queue-prefix costs and summarize conditional budget ceilings |

## Calculation conventions

The primary coefficient definition is `TDI5LM.orig`; `TDI5LM` supplies the alternative. The primary archive population contains records in which every counted valve has a coefficient. Additional summaries use eligible scored subcommunities and the fully scored subset with at most 1,000 recorded valves. Catchment allocation uses a fixed hash, keeping observations from the same catchment together.

Archive masking fractions are 5%, 10%, 20%, 30% and 50%. Random loss samples valves without replacement, using five replicates with seeds 1701 through 1705. Selective loss removes the lowest-scoring, highest-scoring or least-abundant taxa first. Requested fractions become integer valve counts through nearest-integer rounding with half-up ties. Archive summaries average masks within records, records within sites and sites within catchments, giving catchments equal weight. Their intervals use 5,000 catchment bootstrap resamples with seed 20261106. Width quantiles use the inverse weighted empirical distribution; attainment thresholds are 5, 10, 15 and 20 TDI points with tolerance 1e-9.

Verification uses the original 20% random masks, five matched reveal replicates and budgets of 10%, 25% and 50% of hidden valves. The six policies are B0, B1, B2, B3, BAV and BAV_uniform. BAV evaluates 64 hypothetical batches per candidate genus. The development prior adds one count to each permitted score. The within-genus reveal queue is common to all policies and budgets. Score sums use integer tenths; EQR is capped at one, and equality at a class boundary belongs to the higher class. Identification stops as soon as the full compatible range supports one class, including a check after every valve within a batch.

Verification and ceiling summaries use 5,000 catchment bootstrap resamples with seed 20261111 and preserve the paired draws. Conditional recovery divides weighted recovered-and-initially-ambiguous mass by the original weighted initial-ambiguity mass. The simple comparator is selected among B0–B3 by development recovery at the 25% budget; ties favour lower actual effort, then B3, B2, B1 and B0.

The exact retrospective ceilings condition on known hidden scores and fixed within-genus queues. One permits the same five-valve selections as BAV, with shorter terminal selections; the other permits switching genus after each valve. Both search up to the rounded 50% budget. A stored minimum cost of `-1` identifies states without a resolving allocation within that cap. The manuscript and Supplementary Information give the full scientific definitions and interpretation.
