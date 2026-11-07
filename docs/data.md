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

See the [usage reference](usage.md) for command arguments and output tables.
