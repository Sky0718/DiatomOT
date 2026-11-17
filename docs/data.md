# Data sources

All paths below are relative to the repository root. Download the two Environment Agency archives and extract the following four files into `../diatomot-work/raw/`:

| Source | Required files |
| --- | --- |
| [Diatom records](https://environment.data.gov.uk/ecology/explorer/downloads/DIAT_OPEN_DATA.zip) | `DIAT_OPEN_DATA_METRICS.csv`, `DIAT_OPEN_DATA_SITE.csv`, `DIAT_OPEN_DATA_TAXA.csv` |
| [Taxon information](https://environment.data.gov.uk/ecology/explorer/downloads/OPEN_DATA_TAXON_INFO.zip) | `OPEN_DATA_TAXON_INFO.csv` |

The study used the export updated on 20 October 2026 and downloaded on 6 November 2026 (UTC). These download addresses can change as the archive is updated. The SHA-256 checksums identify the ZIP files used in the study:

```text
DIAT_OPEN_DATA.zip
4e1061463fd2b0b1abb2c967b74fef34e056ff4b9b189c3bd082042e4f6d3321
OPEN_DATA_TAXON_INFO.zip
2d4ce3ec5107d07ebac165a4126347f68ce50695726454ae09bb352c9ddf24fa
```

Download [darleq3_data.rda](https://raw.githubusercontent.com/nsj3/darleq3/e750ae0511b96b815aacb836d839f979bce980eb/data/darleq3_data.rda) from the pinned [DARLEQ3 source version](https://github.com/nsj3/darleq3/tree/e750ae0511b96b815aacb836d839f979bce980eb) and place it in the same raw directory. The scripts read the `darleq3_taxa` table directly from this file; an R installation is not required.

The Environment Agency supplies the monitoring data under the Open Government Licence, as stated in its [download directory](https://environment.data.gov.uk/ecology/explorer/downloads/). DARLEQ3 was developed by Steve Juggins and Martyn Kelly; its [package description](https://github.com/nsj3/darleq3/blob/e750ae0511b96b815aacb836d839f979bce980eb/DESCRIPTION) specifies GPL-2. Obtain the source files from their providers and retain their attribution and licence terms.

## Verification inputs and archived results

The verification analyses read **Online Resource 16 — Complete targeted verification results**, supplied with the article. Keep its ZIP or extracted directory outside the repository. The examples use `../diatomot-work/resources/OR16.zip`; this short name is a local convenience and is not a different dataset.

| OR16 member | Content |
| --- | --- |
| `inputs/records.csv` | Record identifiers, monitoring partition, site, catchment, alkalinity and total valve counts |
| `inputs/taxon_codes.csv` | Stable taxon codes, genera and scores in integer tenths |
| `inputs/hidden_taxon_counts.parquet` | Full and hidden taxon counts for the five original masks |
| `inputs/development_prior.csv` | Permitted scores and development valve counts used to form the ecological prior |
| `inputs/initial_states.csv.gz` | Initial interval endpoints, mask seeds and class states |
| `metadata.json` | Dataset identifiers and member checksums |
| Partition outcome directories | Stored terminal outcomes, action histories and reveal queues |

The reader reconstructs the replay arrays from these tables. The development prior adds one count to each permitted score. The original masks are read from the stored hidden counts. These analyses do not require the raw-archive preparation output: the companion resource fixes the exact eligible population and masked states used for verification. It contains 1,952 development and 3,374 confirmation records with eligible alkalinity.

**Online Resource 17 — Verification budget ceilings** contains the archived exact-ceiling states and numerical summaries. The examples use `../diatomot-work/resources/OR17.zip`. Its `oracle_state_requirements_and_costs.csv` contains 133,150 record-mask-queue states across both populations. `minimum_reveals_single_valve` and `minimum_reveals_same_batch` give the resolving costs under the two action protocols. A value of `-1` indicates that no resolving allocation was found within the rounded 50% hidden-valve search cap. Other CSV files provide budget-specific states, population summaries, paired differences and bootstrap quantities.

The ceiling summary can read this ZIP directly or use the `oracle_states.parquet` produced by `run_budget_ceiling.py`. Newly calculated states use `minimum_reveals` for the single-valve cost; the summary reader accepts both column names. Preserve the record, mask and reveal-replicate identifiers when relating either source to OR16.

The companion resources contain data and numerical outputs. The calculation code is in this repository. Their scientific definitions, weighting and assessment conditions are described in the manuscript and Supplementary Information. The [usage reference](usage.md) gives the command arguments and output tables.
