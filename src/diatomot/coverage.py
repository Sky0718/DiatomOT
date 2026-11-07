import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from eco_bounds import split_catchment

GROUPS = ["complete_score", "contains_unscored_or_unmapped"]
METRICS = [
    "scoring_fraction_all_valves",
    "tdi_orig_scoring_denominator",
    "raw_count_sum",
    "scoring_count_sum",
    "unscored_count",
]

def output_directory(value):
    destination = Path(value).expanduser().resolve()
    repository = Path(__file__).resolve().parents[2]
    if (destination.is_relative_to(repository)):
        raise ValueError("The output directory must be outside this repository.")
    if (destination.exists()):
        if (not destination.is_dir()):
            raise NotADirectoryError(destination)
        if (any(destination.iterdir())):
            raise FileExistsError(destination)
    destination.mkdir(parents = True, exist_ok = True)
    return (destination)

def read_records(source):
    records = pd.read_csv(source / "candidate_sample_metadata.csv", low_memory = False)
    if (records.empty):
        raise ValueError("No eligible monitoring records.")
    if ((not records["ANALYSIS_ID"].is_unique) or (not records["SAMPLE_ID"].is_unique)):
        raise ValueError("Sample and analysis identifiers must be unique.")
    if (not records["data_integrity_ok"].all()):
        raise ValueError("The input contains records that failed integrity checks.")
    if (not records["candidate_standard_counts"].all()):
        raise ValueError("The input contains ineligible records.")
    fully_scored = records["raw_count_sum"].eq(records["scoring_count_sum"])
    if (not np.array_equal(fully_scored, records["scoring_fraction_all_valves"].eq(1))):
        raise ValueError("Scoring fractions disagree with scored counts.")
    expected_fraction = (records["scoring_count_sum"] / records["raw_count_sum"])
    np.testing.assert_allclose(records["scoring_fraction_all_valves"], expected_fraction, atol = 1e-14, rtol = 0)
    if (not np.array_equal(records["scoring_count_sum"], records["orig_scoring_count_sum"])):
        raise ValueError("The coefficient definitions have different scored counts.")
    catchments = records["CATCHMENT"].fillna("MISSING_CATCHMENT").astype(str).str.strip()
    records["catchment"] = catchments.replace("", "MISSING_CATCHMENT")
    records["fully_scored"] = fully_scored
    records["coverage_group"] = np.where(fully_scored, GROUPS[0], GROUPS[1])
    records["year"] = pd.to_datetime(records["sample_date_parsed"], errors = "raise").dt.year
    records["partition"] = records["catchment"].map(split_catchment)
    records["unscored_count"] = (records["raw_count_sum"] - records["scoring_count_sum"])
    if (not np.array_equal(records["unscored_count"], records["non_scoring_count"])):
        raise ValueError("Unscored counts disagree with total minus scored counts.")
    if (not records["dictionary_unmapped_count"].le(records["non_scoring_count"]).all()):
        raise ValueError("Unmapped counts exceed unscored counts.")
    records["tdi_orig_scoring_denominator"] = records["recomputed_TDI5LM_orig"]
    return (records)

def group_summaries(records):
    populations = [("all_candidates", records)]
    for group in GROUPS:
        populations.append((group, records.loc[records["coverage_group"].eq(group)]))
    rows = []
    for population in populations:
        group = population[0]
        subset = population[1]
        if (subset.empty):
            continue
        row = {
            "coverage_group"                  : group,
            "records"                         : len(subset),
            "sites"                           : int(subset["SITE_ID"].nunique()),
            "catchments"                      : int(subset["catchment"].nunique()),
            "record_fraction_of_candidates"   : (len(subset) / len(records)),
            "coordinate_missing_records"      : int(subset[["FULL_EASTING", "FULL_NORTHING"]].isna().any(axis = 1).sum()),
            "raw_N_gt1000_records"            : int(subset["raw_count_sum"].gt(1000).sum()),
            "any_dictionary_unmapped_records" : int(subset["dictionary_unmapped_count"].gt(0).sum()),
            "count_pooled_scoring_fraction"   : float(subset["scoring_count_sum"].sum() / subset["raw_count_sum"].sum()),
        }
        for metric in METRICS:
            values = subset[metric]
            row[metric + "_mean"] = float(values.mean())
            row[metric + "_min"] = float(values.min())
            row[metric + "_max"] = float(values.max())
            for quantile in [0.10, 0.25, 0.50, 0.75, 0.90]:
                name = metric + "_q" + str(int(quantile * 100)).zfill(2)
                row[name] = float(values.quantile(quantile, interpolation = "linear"))
        rows.append(row)
    return (pd.DataFrame(rows))

def annual_coverage(records):
    annual = records.groupby("year").agg(
        all_records = ("ANALYSIS_ID", "size"),
        complete_records = ("fully_scored", "sum"),
        all_sites = ("SITE_ID", "nunique"),
        all_catchments = ("catchment", "nunique"),
    )
    first_year = int(records["year"].min())
    last_year = int(records["year"].max())
    annual = annual.reindex(range(first_year, (last_year + 1)), fill_value = 0)
    annual["other_records"] = (annual["all_records"] - annual["complete_records"])
    annual["complete_record_inclusion_fraction"] = (annual["complete_records"] / annual["all_records"])
    for group in GROUPS:
        grouped = records.loc[records["coverage_group"].eq(group)].groupby("year")
        annual[group + "_sites"] = grouped["SITE_ID"].nunique().reindex(annual.index, fill_value = 0)
        annual[group + "_catchments"] = grouped["catchment"].nunique().reindex(annual.index, fill_value = 0)
    return (annual.reset_index())

def unit_coverage(records, unit):
    grouped = records.groupby(unit)
    summary = grouped.agg(
        all_records = ("ANALYSIS_ID", "size"),
        complete_records = ("fully_scored", "sum"),
        first_year = ("year", "min"),
        last_year = ("year", "max"),
    )
    summary["other_records"] = (summary["all_records"] - summary["complete_records"])
    summary["complete_record_inclusion_fraction"] = (summary["complete_records"] / summary["all_records"])
    summary["coverage_membership"] = "both_groups"
    summary.loc[summary["other_records"].eq(0), "coverage_membership"] = "complete_only"
    summary.loc[summary["complete_records"].eq(0), "coverage_membership"] = "other_only"
    if (unit == "catchment") :
        summary["sites"] = grouped["SITE_ID"].nunique()
    else                     :
        summary["catchment"] = grouped["catchment"].first()
    return (summary.reset_index())

def coverage_intersections(records):
    rows = []
    for unit in ["SITE_ID", "catchment"]:
        complete = set(records.loc[records["fully_scored"], unit])
        other = set(records.loc[~records["fully_scored"], unit])
        row = {
            "unit"          : unit,
            "complete"      : len(complete),
            "other"         : len(other),
            "intersection"  : len(complete & other),
            "complete_only" : len(complete - other),
            "other_only"    : len(other - complete),
            "union"         : len(complete | other),
        }
        rows.append(row)
    return (pd.DataFrame(rows))

def main():
    parser = argparse.ArgumentParser(description = "Summarize monitoring coverage and coefficient completeness.")
    parser.add_argument("--inputs", type = Path, required = True)
    parser.add_argument("--output", type = Path, required = True)
    arguments = parser.parse_args()
    source = arguments.inputs.expanduser().resolve()
    destination = output_directory(arguments.output)
    records = read_records(source)
    group_summaries(records).to_csv(destination / "group_summary.csv", index = False)
    annual_coverage(records).to_csv(destination / "year_coverage.csv", index = False)
    unit_coverage(records, "catchment").to_csv(destination / "catchment_coverage.csv", index = False)
    unit_coverage(records, "SITE_ID").to_csv(destination / "site_coverage.csv", index = False)
    coverage_intersections(records).to_csv(destination / "coverage_intersections.csv", index = False)
    partitions = records.groupby(["partition", "coverage_group"]).agg(
        records = ("ANALYSIS_ID", "size"),
        sites = ("SITE_ID", "nunique"),
        catchments = ("catchment", "nunique"),
    )
    partitions.reset_index().to_csv(destination / "partition_coverage.csv", index = False)

if (__name__ == "__main__"):
    main()
