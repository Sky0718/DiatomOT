# Usage reference

Run all commands from the repository root.

## Environment

Use Python 3.11 and the versions in [environment/requirements.txt](../environment/requirements.txt). The [README](../README.md) contains environment setup and the complete command sequence. The [data-source notes](data.md) describe the five required raw files.

## Command arguments

| Script | Arguments |
| --- | --- |
| [prepare_inputs.py](../src/diatomot/prepare_inputs.py) | `--data-dir` for the four BIOSYS CSV files, `--coefficients` for `darleq3_data.rda`, and `--output` for prepared tables |
| [run_ecology.py](../src/diatomot/run_ecology.py) | `--inputs` for prepared tables, `--output` for experiment tables, and optional `--split` |
| [summarize_ecology.py](../src/diatomot/summarize_ecology.py) | `--experiments` for experiment tables, `--output` for summaries, and optional `--split` |
| [precision.py](../src/diatomot/precision.py) | `--experiments` for both experiment partitions and `--output` for precision summaries |
| [concentration.py](../src/diatomot/concentration.py) | `--experiments` and `--summaries` for both partitions, plus `--output` for concentration summaries |
| [coverage.py](../src/diatomot/coverage.py) | `--inputs` for prepared tables and `--output` for coverage summaries |

The two partition-aware scripts accept `--split development`, `--split confirmation` or `--split both`; the default is `both`. Their output parents contain matching `development/` and `confirmation/` subdirectories. Use the same parent directory for the two separate partition commands. Precision and concentration require both partition subdirectories.

All paths may be absolute or relative to the current working directory. Output paths must be outside the repository. Each destination must be absent or empty; for experiment and summary parents, this applies to the requested partition subdirectories. Choose a new external destination to retain an earlier set of results.

## Output tables

| Destination | Main files and contents |
| --- | --- |
| Prepared inputs | `candidate_sample_metadata.csv`, `candidate_sample_taxon_long.csv.gz`, `darleq3_rda_coefficients.csv` and `OPEN_DATA_TAXON_INFO.csv` form the experiment inputs; inventory and eligibility tables describe sample selection |
| Experiment partition | `mask_results.csv.gz` contains mask-level endpoints; `coefficient_definition_results.csv` contains full-count indices; `primary_genus_contributions.csv.gz` contains genus contributions at the primary condition |
| Summary partition | `all_condition_summary.csv` contains condition estimates and intervals; `primary_endpoint.json` contains primary summaries; additional tables retain record, site and catchment quantities, coefficient comparisons and genus contributions |
| Precision | `threshold_readouts.csv`, `weighted_quantiles.csv` and `catchment_threshold_cdf.csv` contain threshold probabilities and weighted width summaries |
| Concentration | `catchment_shares_all.csv` contains catchment shares; partition-specific tables contain bootstrap ratios and omission estimates; `distribution_summary.csv` and `overall_ratios_and_influence.csv` contain aggregate quantities |
| Coverage | `group_summary.csv`, `year_coverage.csv`, `site_coverage.csv`, `catchment_coverage.csv`, `coverage_intersections.csv` and `partition_coverage.csv` describe the retained archive |

CSV files contain numerical results and identifiers. JSON files contain numerical summaries or input and calculation metadata. Index errors and interval widths use TDI points; threshold probabilities and genus shares are stored as fractions. Consult column names and the corresponding calculation in the source when combining tables.

Keep the prepared tables, experiment partitions and summaries from the same source snapshot together. Some steps read file hashes recorded by an earlier step, so renaming or editing intermediate tables changes their relationship to those inputs. The source code is in [src/diatomot](../src/diatomot); `eco_bounds.py` supplies shared calculations to the command scripts.
