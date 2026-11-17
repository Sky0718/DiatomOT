# Usage reference

Run all commands from the repository root. Paths may be absolute or relative to the current working directory. Keep raw inputs, companion resources, calculated outputs and Python environments outside the repository.

## Environment and inputs

Use Python 3.11 and the versions in [environment/requirements.txt](../environment/requirements.txt). The [README](../README.md) contains environment setup and the main command sequence. The [data-source notes](data.md) describe the raw-archive inputs and the companion resources.

The original archive analyses read the five provider files. Verification and ceiling analyses read Online Resource 16 directly and can run independently of the archive-analysis commands. Online Resource 17 supplies archived ceiling states for direct aggregation. In the examples, these files are saved as `../diatomot-work/resources/OR16.zip` and `../diatomot-work/resources/OR17.zip`.

## Archive-analysis arguments

| Script | Arguments |
| --- | --- |
| [prepare_inputs.py](../src/diatomot/prepare_inputs.py) | `--data-dir` for the four BIOSYS CSV files, `--coefficients` for `darleq3_data.rda`, and `--output` for prepared tables |
| [run_ecology.py](../src/diatomot/run_ecology.py) | `--inputs` for prepared tables, `--output` for experiment tables, and optional `--split` |
| [summarize_ecology.py](../src/diatomot/summarize_ecology.py) | `--experiments` for experiment tables, `--output` for summaries, and optional `--split` |
| [precision.py](../src/diatomot/precision.py) | `--experiments` for both experiment partitions and `--output` for precision summaries |
| [concentration.py](../src/diatomot/concentration.py) | `--experiments` and `--summaries` for both partitions, plus `--output` for concentration summaries |
| [coverage.py](../src/diatomot/coverage.py) | `--inputs` for prepared tables and `--output` for coverage summaries |

The two partition-aware scripts accept `--split development`, `--split confirmation` or `--split both`; the default is `both`. Their output parents contain matching `development/` and `confirmation/` subdirectories. Use the same parent directory for separate partition commands. Precision and concentration require both partition subdirectories.

## Verification arguments

| Script | Arguments |
| --- | --- |
| [run_verification.py](../src/diatomot/run_verification.py) | Required `--resource`, `--partition` and `--out`; optional `--record-ids` and `--workers` |
| [summarize_verification.py](../src/diatomot/summarize_verification.py) | Either `--resource` for archived outcomes or `--terminals` for a newly calculated output parent; required `--partition` and `--out` |
| [run_budget_ceiling.py](../src/diatomot/run_budget_ceiling.py) | Required `--resource` and `--out`; optional `--partition`, `--record-ids` and `--workers` |
| [summarize_budget_ceiling.py](../src/diatomot/summarize_budget_ceiling.py) | Required `--resource`, `--ceiling` and `--out` |

`--resource` accepts the OR16 ZIP or the directory containing its extracted members. Verification replay and aggregation require `--partition development` or `--partition confirmation`. The ceiling calculation uses both partitions unless `--partition` selects one. Replay and ceiling calculation default to four workers; use a positive integer for `--workers`.

`--ceiling` accepts the OR17 ZIP, its `oracle_state_requirements_and_costs.csv`, a newly calculated `oracle_states.parquet`, or a directory containing that Parquet file. Both complete populations and all three budgets are required by the ceiling summary. The reader accepts `minimum_reveals_single_valve` from OR17 or `minimum_reveals` from the calculation output for the single-valve cost.

All output paths must be outside the repository. Destinations must be absent or empty. For verification replay and aggregation, this restriction applies to the requested partition subdirectory; the two partitions share one output parent. A summary parent may therefore contain its development outputs when confirmation is added. Choose a new external destination to retain earlier results.

## Subset calculations

To calculate the four monitoring examples using one worker:

```text
python -B src/diatomot/run_verification.py --resource ../diatomot-work/resources/OR16.zip --partition confirmation --record-ids 677036 686137 654561 662592 --out ../diatomot-work/verification-examples --workers 1
```

To calculate exact ceiling costs for one confirmation record:

```text
python -B src/diatomot/run_budget_ceiling.py --resource ../diatomot-work/resources/OR16.zip --partition confirmation --record-ids 677036 --out ../diatomot-work/ceiling-example --workers 1
```

These subsets describe selected records. Full-population aggregation requires the complete development or confirmation population. Each full replay record contributes 450 terminal rows: six policies, three budgets, five masks and five reveal replicates.

## Summarizing newly calculated outputs

After calculating both verification partitions with the commands in the README, aggregate development first and then confirmation:

```text
python -B src/diatomot/summarize_verification.py --terminals ../diatomot-work/verification --partition development --out ../diatomot-work/replayed-statistics
python -B src/diatomot/summarize_verification.py --terminals ../diatomot-work/verification --partition confirmation --out ../diatomot-work/replayed-statistics
```

Development aggregation writes `development_baseline_lock.json` in the output parent. Confirmation reads this selection from the same parent. The comparator is selected among B0–B3 at the 25% budget, first by class recovery, then lower actual effort, then B3, B2, B1 and B0. The archived selection is B3.

After calculating both ceiling populations, combine their costs with the archived verification outcomes:

```text
python -B src/diatomot/summarize_budget_ceiling.py --resource ../diatomot-work/resources/OR16.zip --ceiling ../diatomot-work/ceilings/oracle_states.parquet --out ../diatomot-work/recomputed-ceiling-statistics
```

The ceiling summary rebuilds the population weights from OR16, uses the same 5,000 catchment draws for paired comparisons, and reports recovery for BAV, B3, BAV_uniform and both ceilings. It also reports same-batch-minus-BAV and single-valve-minus-same-batch differences. The `budget_ceiling_summary.csv` fields and metric identifiers correspond to the summary in OR17.

## Output tables

| Destination | Main files and contents |
| --- | --- |
| Prepared inputs | `candidate_sample_metadata.csv`, `candidate_sample_taxon_long.csv.gz`, `darleq3_rda_coefficients.csv` and `OPEN_DATA_TAXON_INFO.csv` form the experiment inputs; inventory and eligibility tables describe sample selection |
| Experiment partition | `mask_results.csv.gz` contains mask-level endpoints; `coefficient_definition_results.csv` contains full-count indices; `primary_genus_contributions.csv.gz` contains genus contributions at the primary condition |
| Archive-summary partition | `all_condition_summary.csv` contains condition estimates and intervals; `primary_endpoint.json` contains primary summaries; additional tables retain record, site and catchment quantities, coefficient comparisons and genus contributions |
| Precision | `threshold_readouts.csv`, `weighted_quantiles.csv` and `catchment_threshold_cdf.csv` contain threshold probabilities and weighted width summaries |
| Concentration | `catchment_shares_all.csv` contains catchment shares; partition-specific tables contain bootstrap ratios and omission estimates; `distribution_summary.csv` and `overall_ratios_and_influence.csv` contain aggregate quantities |
| Coverage | `group_summary.csv`, `year_coverage.csv`, `site_coverage.csv`, `catchment_coverage.csv`, `coverage_intersections.csv` and `partition_coverage.csv` describe the retained archive |
| Verification partition | Parquet chunks in `terminal/`, `actions/` and `reveal_queues/` contain terminal outcomes, policy choices and shared reveal sequences |
| Verification-summary partition | `record_estimates.parquet`, `site_estimates.parquet`, `catchment_basis.csv`, `summary.csv`, `paired_comparisons.csv` and bootstrap NPZ files; the output parent also contains `development_baseline_lock.json` |
| Ceiling calculation | Chunk Parquet files and `oracle_states.parquet` contain the minimum resolving cost under both action protocols for each record-mask-queue state |
| Ceiling summary | `budget_ceiling_summary.csv` contains recovery estimates, percentile intervals and paired differences for both partitions and all budgets |

Index errors and interval widths use TDI points. Archive threshold probabilities and genus shares are stored as fractions. Verification budgets are fractions of the initially hidden valves. Recovery calculations retain initially ambiguous cases as the denominator; all-case class resolution is a separate quantity. Consult column names and the corresponding source calculation when combining tables.

Keep prepared tables, experiment partitions and summaries from the same source snapshot together. Some steps read input hashes stored by an earlier step, so renaming or editing intermediate tables changes their relationship to those inputs. Verification input integrity is checked against the companion metadata before calculation. Code and calculated outputs occupy separate directories throughout.
