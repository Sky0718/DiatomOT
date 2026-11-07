import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

METRICS = [
    "point_abs_error",
    "point_signed_error",
    "genus_width",
    "unrestricted_width",
    "width_saved",
    "genus_usable10",
    "unrestricted_usable10",
    "usable10_difference",
    "genus_contains",
    "unrestricted_contains",
    "definition_union_hull_width",
    "definition_union_usable10",
    "point_sample_class_disagreement",
    "point_good_boundary_disagreement",
    "genus_sample_class_resolved",
    "unrestricted_sample_class_resolved",
    "definition_union_sample_class_resolved",
    "genus_false_class_certificate",
]
DEFINITION_METRICS = [
    "delta",
    "absolute_delta",
    "delta_at_least5",
    "delta_at_least10",
    "nonzero_delta",
    "definition_sample_class_diff",
    "definition_good_boundary_diff",
]
POPULATIONS = [
    "complete_counts",
    "all_scoring_subcommunities",
    "complete_counts_raw_N_le1000",
]

def eqr_and_class(tdi, alkalinity):
    broadcast = np.broadcast_arrays(
        np.asarray(tdi, dtype = float),
        np.asarray(alkalinity, dtype = float),
    )
    index_values = broadcast[0]
    alkalinity_values = broadcast[1]
    eligible = (
        np.isfinite(alkalinity_values)
        & (alkalinity_values >= 5)
        & (alkalinity_values <= 250)
        & np.isfinite(index_values)
    )
    ratios = np.full(index_values.shape, np.nan)
    expected = (9.933 * np.exp(0.81 * np.log10(alkalinity_values[eligible])))
    raw_ratios = (0.8 * (100 - index_values[eligible]) / (100 - expected))
    ratios[eligible] = np.minimum(1.0, raw_ratios)
    classes = np.full(index_values.shape, -1, dtype = int)
    classes[eligible] = np.searchsorted(
        [0.2, 0.4, 0.6, 0.8], ratios[eligible], side = "right"
    )
    return (ratios, classes)

def metrics(frame, versions):
    result = frame.copy()
    version_lookup = versions.set_index("ANALYSIS_ID")
    result["site_alkalinity"] = result["ANALYSIS_ID"].map(
        version_lookup["site_alkalinity"]
    )
    signed_error = (result["orig_point"] - result["orig_truth"])
    result["point_abs_error"] = abs(signed_error)
    result["point_signed_error"] = signed_error
    result["genus_width"] = result["orig_width"]
    result["unrestricted_width"] = result["orig_blind_width"]
    result["width_saved"] = (result["orig_blind_width"] - result["orig_width"])
    result["genus_usable10"] = (result["orig_width"] <= (10 + 1e-9)).astype(float)
    result["unrestricted_usable10"] = (
        result["orig_blind_width"] <= (10 + 1e-9)
    ).astype(float)
    result["usable10_difference"] = (
        result["genus_usable10"] - result["unrestricted_usable10"]
    )
    genus_contains = (
        (result["orig_lower"] <= (result["orig_truth"] + 1e-9))
        & (result["orig_upper"] >= (result["orig_truth"] - 1e-9))
    )
    unrestricted_contains = (
        (result["orig_blind_lower"] <= (result["orig_truth"] + 1e-9))
        & (result["orig_blind_upper"] >= (result["orig_truth"] - 1e-9))
    )
    result["genus_contains"] = genus_contains.astype(float)
    result["unrestricted_contains"] = unrestricted_contains.astype(float)
    result["definition_union_usable10"] = (
        result["definition_union_hull_width"] <= (10 + 1e-9)
    ).astype(float)
    alkalinity = result["site_alkalinity"].to_numpy(dtype = float)
    reference = eqr_and_class(result["orig_truth"], alkalinity)
    point = eqr_and_class(result["orig_point"], alkalinity)
    lower_class = eqr_and_class(result["orig_lower"], alkalinity)[1]
    upper_class = eqr_and_class(result["orig_upper"], alkalinity)[1]
    unrestricted_lower = eqr_and_class(result["orig_blind_lower"], alkalinity)[1]
    unrestricted_upper = eqr_and_class(result["orig_blind_upper"], alkalinity)[1]
    alternative_lower = eqr_and_class(result["current_lower"], alkalinity)[1]
    alternative_upper = eqr_and_class(result["current_upper"], alkalinity)[1]
    eligible = (reference[1] >= 0)
    genus_resolved = (lower_class == upper_class)
    union_resolved = (
        genus_resolved
        & (alternative_lower == alternative_upper)
        & (lower_class == alternative_lower)
    )
    false_resolution = (eligible & genus_resolved & (lower_class != reference[1]))
    if (false_resolution.any()):
        raise ValueError("A resolved interval disagrees with its reference class.")
    result["eqr_eligible"] = eligible
    result["point_sample_class_disagreement"] = np.where(
        eligible, (point[1] != reference[1]), np.nan
    )
    result["point_good_boundary_disagreement"] = np.where(
        eligible, ((point[0] >= 0.6) != (reference[0] >= 0.6)), np.nan
    )
    result["genus_sample_class_resolved"] = np.where(
        eligible, genus_resolved, np.nan
    )
    result["unrestricted_sample_class_resolved"] = np.where(
        eligible, (unrestricted_lower == unrestricted_upper), np.nan
    )
    result["definition_union_sample_class_resolved"] = np.where(
        eligible, union_resolved, np.nan
    )
    result["genus_false_class_certificate"] = np.where(
        eligible, false_resolution, np.nan
    )
    return (result)

def hierarchical(frame, columns, extra = None):
    grouping = []
    if (extra is not None):
        grouping = list(extra)
    analysis_keys = (grouping + ["catchment", "SITE_ID", "ANALYSIS_ID"])
    site_keys = (grouping + ["catchment", "SITE_ID"])
    catchment_keys = (grouping + ["catchment"])
    analysis_groups = frame.groupby(analysis_keys, observed = True, sort = True)
    analyses = analysis_groups[columns].mean()
    sites = analyses.groupby(level = site_keys, sort = True).mean()
    catchments = sites.groupby(level = catchment_keys, sort = True).mean()
    return (analyses.reset_index(), sites.reset_index(), catchments.reset_index())

def bootstrap(values, seed = 20261106):
    finite_values = np.asarray(values, dtype = float)
    finite_values = finite_values[np.isfinite(finite_values)]
    count = len(finite_values)
    if (count == 0):
        return (np.nan, np.nan, np.nan, 0)
    random = np.random.default_rng(seed)
    indices = random.integers(0, count, size = (5000, count))
    resampled_means = finite_values[indices].mean(axis = 1)
    estimate = float(finite_values.mean())
    lower = float(np.quantile(resampled_means, 0.025))
    upper = float(np.quantile(resampled_means, 0.975))
    return (estimate, lower, upper, count)

def population_rows(frame, population):
    if (population == "complete_counts")                :
        selected = frame["fully_scored"]
    elif (population == "all_scoring_subcommunities")   :
        selected = np.ones(len(frame), dtype = bool)
    elif (population == "complete_counts_raw_N_le1000") :
        selected = (frame["fully_scored"] & ~frame["N_gt_1000"])
    else                                                :
        raise ValueError(f"Unknown population: {population}")
    return (frame.loc[selected].copy())

def metric_intervals(catchments, columns):
    summaries = {}
    for column in columns:
        estimate = bootstrap(catchments[column])
        summaries[column] = {
            "estimate"   : estimate[0],
            "ci95_lower" : estimate[1],
            "ci95_upper" : estimate[2],
            "catchments" : estimate[3],
        }
    return (summaries)

def write_json(path, value):
    content = json.dumps(value, indent = 2, allow_nan = False)
    Path(path).write_text((content + "\n"), encoding = "utf-8")

def check_output_directory(path):
    destination = Path(path).resolve()
    repository = Path(__file__).resolve().parents[2]
    if (destination.is_relative_to(repository)):
        raise ValueError("Choose an output directory outside the repository.")
    if (destination.exists()):
        if (not destination.is_dir()):
            raise FileExistsError(destination)
        if (any(destination.iterdir())):
            raise FileExistsError(f"Output directory is not empty: {destination}")
    return (destination)

def file_hash(path):
    digest = hashlib.sha256()
    with (Path(path).open("rb")) as stream:
        while (True):
            block = stream.read(1048576)
            if (not block):
                break
            digest.update(block)
    return (digest.hexdigest())

def check_run_receipt(source):
    path = (source / "run_receipt.json")
    if (not path.exists()):
        raise FileNotFoundError(path)
    receipt = json.loads(path.read_text(encoding = "utf-8"))
    if (receipt.get("status") != "COMPLETE"):
        raise ValueError(f"Experiment is incomplete: {source}")
    outputs = receipt.get("outputs", {})
    for name in outputs:
        path = (source / name).resolve()
        if (not path.is_relative_to(source.resolve())):
            raise ValueError("An experiment receipt contains an external path.")
        if (file_hash(path) != outputs[name]):
            raise ValueError(f"Experiment file changed: {name}")

def summarize_masks(frame, destination):
    records = []
    primary = None
    for population in POPULATIONS:
        subset = population_rows(frame, population)
        levels = hierarchical(subset, METRICS, ["fraction_requested", "mechanism"])
        analyses = levels[0]
        sites = levels[1]
        catchments = levels[2]
        catchments.to_csv(
            (destination / f"{population}_catchment_metrics.csv"), index = False
        )
        if (population == "complete_counts"):
            analyses.to_csv(
                (destination / "complete_counts_analysis_metrics.csv.gz"), index = False
            )
            sites.to_csv(
                (destination / "complete_counts_site_metrics.csv.gz"), index = False
            )
            primary_analyses = analyses.loc[
                (analyses["fraction_requested"] == 0.2)
                & (analyses["mechanism"] == "mcar")
            ]
            primary_catchments = catchments.loc[
                (catchments["fraction_requested"] == 0.2)
                & (catchments["mechanism"] == "mcar")
            ]
            primary_analyses.to_csv(
                (destination / "primary_analysis_metrics.csv"), index = False
            )
            primary_catchments.to_csv(
                (destination / "primary_catchment_metrics.csv"), index = False
            )
            primary = metric_intervals(primary_catchments, METRICS)
            primary["analyses"] = int(primary_analyses["ANALYSIS_ID"].nunique())
            primary["sites"] = int(primary_analyses["SITE_ID"].nunique())
            primary["mask_replicates_per_analysis"] = 5
            primary["scope"] = (
                "complete scoring counts; original score definition; correct genus "
                "supplied; controlled 20% label loss"
            )
        conditions = catchments.groupby(["fraction_requested", "mechanism"], sort = True)
        for condition, condition_catchments in conditions:
            fraction = condition[0]
            mechanism = condition[1]
            analysis_selection = (
                (analyses["fraction_requested"] == fraction)
                & (analyses["mechanism"] == mechanism)
            )
            site_selection = (
                (sites["fraction_requested"] == fraction)
                & (sites["mechanism"] == mechanism)
            )
            row = {
                "population" : population,
                "fraction"   : fraction,
                "mechanism"  : mechanism,
                "analyses"   : int(analyses.loc[analysis_selection, "ANALYSIS_ID"].nunique()),
                "sites"      : int(sites.loc[site_selection, "SITE_ID"].nunique()),
                "catchments" : len(condition_catchments),
            }
            for column in METRICS:
                interval = bootstrap(condition_catchments[column])
                row[column] = interval[0]
                row[(column + "_ci95_lower")] = interval[1]
                row[(column + "_ci95_upper")] = interval[2]
                row[(column + "_groups")] = interval[3]
            records.append(row)
    if (primary is None):
        raise ValueError("No primary population was summarized.")
    pd.DataFrame(records).to_csv((destination / "all_condition_summary.csv"), index = False)
    write_json((destination / "primary_endpoint.json"), primary)
    return (primary)

def summarize_definitions(versions, destination):
    values = versions.copy()
    values["absolute_delta"] = abs(values["delta"])
    values["delta_at_least5"] = (values["absolute_delta"] >= (5 - 1e-9)).astype(float)
    values["delta_at_least10"] = (values["absolute_delta"] >= (10 - 1e-9)).astype(float)
    values["nonzero_delta"] = (values["absolute_delta"] > 1e-9).astype(float)
    original = eqr_and_class(values["orig"], values["site_alkalinity"])
    alternative = eqr_and_class(values["current"], values["site_alkalinity"])
    values["orig_eqr"] = original[0]
    values["current_eqr"] = alternative[0]
    values["orig_sample_class"] = original[1]
    values["current_sample_class"] = alternative[1]
    eligible = (original[1] >= 0)
    values["definition_sample_class_diff"] = np.where(
        eligible, (original[1] != alternative[1]), np.nan
    )
    values["definition_good_boundary_diff"] = np.where(
        eligible, ((original[0] >= 0.6) != (alternative[0] >= 0.6)), np.nan
    )
    values.to_csv((destination / "coefficient_definition_sample_results.csv"), index = False)
    summaries = {}
    for population in POPULATIONS:
        subset = population_rows(values, population)
        levels = hierarchical(subset, DEFINITION_METRICS)
        catchments = levels[2]
        catchments.to_csv(
            (destination / f"{population}_definition_catchment_metrics.csv"), index = False
        )
        summary = metric_intervals(catchments, DEFINITION_METRICS)
        summary["analyses"] = len(subset)
        summary["sites"] = int(subset["SITE_ID"].nunique())
        summary["eqr_eligible_analyses"] = int((subset["orig_sample_class"] >= 0).sum())
        summary["maximum_absolute_delta"] = float(subset["absolute_delta"].max())
        summaries[population] = summary
    write_json((destination / "coefficient_definition_summary.json"), summaries)

def summarize_genera(source, destination):
    contributions = pd.read_csv(source / "primary_genus_contributions.csv.gz")
    contributions = contributions.loc[contributions["fully_scored"]]
    grouping = ["catchment", "SITE_ID", "ANALYSIS_ID", "genus"]
    analysis_means = contributions.groupby(grouping)["width_contribution"].mean()
    analysis_matrix = analysis_means.unstack("genus", fill_value = 0)
    site_matrix = analysis_matrix.groupby(level = ["catchment", "SITE_ID"]).mean()
    catchment_matrix = site_matrix.groupby(level = "catchment").mean()
    contributions_mean = catchment_matrix.mean().sort_values(ascending = False)
    summary = pd.DataFrame(
        {
            "genus"                              : contributions_mean.index,
            "catchment_equal_width_contribution" : contributions_mean.to_numpy(),
        }
    )
    summary["share_of_width"] = (
        contributions_mean.to_numpy() / contributions_mean.sum()
    )
    summary["cumulative_share"] = summary["share_of_width"].cumsum()
    summary.to_csv((destination / "genus_width_contributions.csv"), index = False)
    catchment_matrix.to_csv(destination / "genus_width_catchment_matrix.csv")

def summarize_partition(experiments, destination, partition):
    source = (experiments / partition)
    check_run_receipt(source)
    destination.mkdir(parents = True, exist_ok = True)
    raw = pd.read_csv(source / "mask_results.csv.gz", dtype = {"actual_seed" : str})
    versions = pd.read_csv(source / "coefficient_definition_results.csv")
    values = metrics(raw, versions)
    summarize_masks(values, destination)
    summarize_definitions(versions, destination)
    summarize_genera(source, destination)
    inputs = {}
    for name in ["mask_results.csv.gz", "coefficient_definition_results.csv", "primary_genus_contributions.csv.gz"]:
        inputs[name] = file_hash(source / name)
    outputs = {}
    for path in sorted(destination.iterdir()):
        if (path.is_file()):
            outputs[path.name] = file_hash(path)
    report = {
        "status"                            : "COMPLETE",
        "split"                             : partition,
        "raw_cases"                         : len(raw),
        "analyses"                          : len(versions),
        "genus_bound_containment_failures"  : int((values["genus_contains"] != 1).sum()),
        "invalid_sample_class_certificates" : int(
            values["genus_false_class_certificate"].fillna(0).sum()
        ),
        "bootstrap_resamples"               : 5000,
        "bootstrap_seed"                    : 20261106,
        "inputs"                            : inputs,
        "outputs"                           : outputs,
    }
    write_json((destination / "summary_receipt.json"), report)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--experiments", type = Path, required = True)
    parser.add_argument("--output", type = Path, required = True)
    parser.add_argument(
        "--split", choices = ["development", "confirmation", "both"], default = "both"
    )
    arguments = parser.parse_args()
    experiments = arguments.experiments.resolve()
    if (arguments.split == "both") :
        partitions = ["development", "confirmation"]
    else                           :
        partitions = [arguments.split]
    destinations = {}
    for partition in partitions:
        destinations[partition] = check_output_directory(arguments.output / partition)
    for partition in partitions:
        summarize_partition(experiments, destinations[partition], partition)

if (__name__ == "__main__"):
    main()
