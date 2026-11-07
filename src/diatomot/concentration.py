import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from summarize_ecology import check_output_directory
from summarize_ecology import check_run_receipt
from summarize_ecology import file_hash
from summarize_ecology import write_json

GENERA = ["Navicula", "Gomphonema", "Nitzschia"]
PARTITIONS = ["development", "confirmation"]
SEED = 20261106
RESAMPLES = 5000

def hierarchy(analysis_values):
    site_values = analysis_values.groupby(level = ["catchment", "SITE_ID"], sort = True).mean()
    catchment_values = site_values.groupby(level = "catchment", sort = True).mean()
    return (catchment_values)

def distribution(values):
    values = np.asarray(values, dtype = float)
    if (not np.isfinite(values).all()):
        raise ValueError("Catchment shares must be finite.")
    quantiles = np.quantile(values, [0.1, 0.25, 0.5, 0.75, 0.9], method = "linear")
    majority = (values > 0.5)
    summary = {
        "n"                        : len(values),
        "min"                      : float(values.min()),
        "q10"                      : float(quantiles[0]),
        "q25"                      : float(quantiles[1]),
        "q50"                      : float(quantiles[2]),
        "q75"                      : float(quantiles[3]),
        "q90"                      : float(quantiles[4]),
        "max"                      : float(values.max()),
        "mean_of_catchment_ratios" : float(values.mean()),
        "n_gt_50pct"               : int(majority.sum()),
        "fraction_gt_50pct"        : float(majority.mean()),
    }
    return (summary)

def maximum_difference(left, right):
    difference = np.abs(np.asarray(left) - np.asarray(right))
    return (float(difference.max()))

def require_close(difference, tolerance, message):
    if ((not np.isfinite(difference)) or (difference >= tolerance)):
        raise ValueError(message)

def check_contributions(contributions):
    if (contributions.empty or (not (contributions["N"] > 0).all())):
        raise ValueError("Complete scoring counts must be nonempty and positive.")
    if (not (contributions["N"] == contributions["raw_N"]).all()):
        raise ValueError("Concentration analysis requires complete scoring counts.")
    if (not (contributions["non_scoring_N"] == 0).all()):
        raise ValueError("Complete scoring counts cannot include unscored taxa.")
    expected_seeds = {1701, 1702, 1703, 1704, 1705}
    if (set(contributions["replicate_seed"].unique()) != expected_seeds):
        raise ValueError("The primary condition requires all five mask seeds.")
    duplicates = contributions.duplicated(["ANALYSIS_ID", "replicate_seed", "genus"])
    if (duplicates.any()):
        raise ValueError("Repeated analysis, seed and genus combinations were found.")
    repetitions = contributions.groupby(["ANALYSIS_ID", "genus"]).size()
    if (not (repetitions == 5).all()):
        raise ValueError("Each analysis and genus must have five masks.")
    site_counts = contributions.groupby("ANALYSIS_ID")["SITE_ID"].nunique()
    if (site_counts.max() != 1):
        raise ValueError("Each analysis must belong to one site.")

def paired_ratios(catchments, partition):
    pairs = {
        "width_share"       : ["three_genus_width", "total_width"],
        "masked_mass_share" : ["three_genus_masked_mass", "total_masked_mass"],
    }
    if (len(catchments) < 2):
        raise ValueError("At least two catchments are required for omission analysis.")
    random = np.random.default_rng(SEED)
    sample_indices = random.integers(
        0, len(catchments), size = (RESAMPLES, len(catchments))
    )
    bootstrap = pd.DataFrame(
        {
            "replicate" : np.arange(RESAMPLES),
            "partition" : partition,
        }
    )
    omissions = catchments[["catchment", "partition"]].copy()
    summaries = {}
    for name in pairs:
        numerator = catchments[pairs[name][0]].to_numpy(dtype = float)
        denominator = catchments[pairs[name][1]].to_numpy(dtype = float)
        estimate = float(numerator.mean() / denominator.mean())
        resampled_numerator = numerator[sample_indices].mean(axis = 1)
        resampled_denominator = denominator[sample_indices].mean(axis = 1)
        resampled_ratios = (resampled_numerator / resampled_denominator)
        omitted_numerator = (numerator.sum() - numerator)
        omitted_denominator = (denominator.sum() - denominator)
        omitted_ratios = (omitted_numerator / omitted_denominator)
        bootstrap[name] = resampled_ratios
        omissions[name] = omitted_ratios
        omissions[(name + "_change_pp")] = (100 * (omitted_ratios - estimate))
        influential_index = int(np.argmax(np.abs(omitted_ratios - estimate)))
        minimum_index = int(np.argmin(omitted_ratios))
        maximum_index = int(np.argmax(omitted_ratios))
        influential_change = (omitted_ratios[influential_index] - estimate)
        summaries[name] = {
            "ratio_of_means"             : estimate,
            "numerator_catchment_mean"   : float(numerator.mean()),
            "denominator_catchment_mean" : float(denominator.mean()),
            "ci95_percentile"            : np.quantile(resampled_ratios, [0.025, 0.975]).tolist(),
            "leave_one_out_min"          : float(omitted_ratios.min()),
            "leave_one_out_min_removed"  : catchments.iloc[minimum_index]["catchment"],
            "leave_one_out_max"          : float(omitted_ratios.max()),
            "leave_one_out_max_removed"  : catchments.iloc[maximum_index]["catchment"],
            "largest_absolute_change_pp" : float(100 * abs(influential_change)),
            "largest_signed_change_pp"   : float(100 * influential_change),
            "most_influential_removed"   : catchments.iloc[influential_index]["catchment"],
        }
    return (bootstrap, omissions, summaries)

def process(experiments, summaries, destination, partition):
    experiment_directory = (experiments / partition)
    summary_directory = (summaries / partition)
    check_run_receipt(experiment_directory)
    raw = pd.read_csv(experiment_directory / "primary_genus_contributions.csv.gz")
    support = pd.read_csv(experiment_directory / "genus_score_support.csv")
    support = support.set_index("genus")
    saved = pd.read_csv(summary_directory / "genus_width_catchment_matrix.csv")
    saved = saved.set_index("catchment")
    primary = pd.read_csv(summary_directory / "primary_catchment_metrics.csv")
    primary = primary.set_index("catchment")
    genus_summary = pd.read_csv(summary_directory / "genus_width_contributions.csv")
    genus_summary = genus_summary.set_index("genus")
    contributions = raw.loc[raw["fully_scored"]].copy()
    check_contributions(contributions)
    keys = ["catchment", "SITE_ID", "ANALYSIS_ID"]
    genus_keys = (keys + ["genus"])
    width_means = contributions.groupby(genus_keys, sort = True)["width_contribution"].mean()
    width_matrix = width_means.unstack("genus", fill_value = 0)
    reconstruction = hierarchy(width_matrix)
    reconstruction = reconstruction.reindex(
        index = saved.index, columns = saved.columns, fill_value = 0
    )
    reconstruction_difference = maximum_difference(reconstruction, saved)
    require_close(
        reconstruction_difference, 1e-10, "Genus contributions differ from the summary."
    )
    primary = primary.loc[saved.index]
    if (not (primary["fraction_requested"] == 0.2).all()):
        raise ValueError("Concentration requires the 20% masking condition.")
    if (not (primary["mechanism"] == "mcar").all()):
        raise ValueError("Concentration requires the random masking condition.")
    total_width = saved.sum(axis = 1)
    total_width_difference = maximum_difference(total_width, primary["genus_width"])
    require_close(total_width_difference, 1e-10, "Genus widths do not sum to total width.")
    all_spans = (support["orig_high"] - support["orig_low"])
    spans = all_spans.loc[GENERA]
    expected_spans = {
        "Navicula"   : 4.0,
        "Gomphonema" : 4.0,
        "Nitzschia"  : 3.0,
    }
    if (spans.to_dict() != expected_spans):
        raise ValueError("The original coefficient spans for the fixed genera differ.")
    contributions["masked_mass"] = (contributions["masked_valves"] / contributions["N"])
    mass_means = contributions.groupby(genus_keys, sort = True)["masked_mass"].mean()
    mass_matrix = mass_means.unstack("genus", fill_value = 0)
    masses = hierarchy(mass_matrix).loc[saved.index]
    selected_mass = masses[GENERA].sum(axis = 1)
    algebra_mass = saved[GENERA].div((25 * spans), axis = 1).sum(axis = 1)
    mass_difference = maximum_difference(selected_mass, algebra_mass)
    require_close(mass_difference, 1e-12, "Masked mass disagrees with the width algebra.")
    total_mass = masses.sum(axis = 1)
    unrestricted_mass = (primary["unrestricted_width"] / 100)
    total_mass_difference = maximum_difference(total_mass, unrestricted_mass)
    require_close(
        total_mass_difference, 1e-12, "Total mass disagrees with unrestricted width."
    )
    if ((not (total_width > 0).all()) or (not (total_mass > 0).all())):
        raise ValueError("Width and masked mass must have positive catchment totals.")
    unique_analyses = contributions[keys].drop_duplicates()
    counts = unique_analyses.groupby("catchment").agg(
        n_analyses = ("ANALYSIS_ID", "size"), n_sites = ("SITE_ID", "nunique")
    )
    catchments = counts.loc[saved.index].copy()
    for genus in GENERA:
        catchments[(genus + "_width")] = saved[genus]
        catchments[(genus + "_masked_mass")] = masses[genus]
    catchments["three_genus_width"] = saved[GENERA].sum(axis = 1)
    catchments["total_width"] = total_width
    catchments["width_share"] = (catchments["three_genus_width"] / total_width)
    catchments["three_genus_masked_mass"] = selected_mass
    catchments["total_masked_mass"] = total_mass
    catchments["masked_mass_share"] = (selected_mass / total_mass)
    catchments["partition"] = partition
    catchments = catchments.reset_index()
    catchments.to_csv((destination / f"{partition}_catchment_shares.csv"), index = False)
    ratio_results = paired_ratios(catchments, partition)
    ratio_results[0].to_csv((destination / f"{partition}_bootstrap_ratios.csv"), index = False)
    ratio_results[1].to_csv(
        (destination / f"{partition}_leave_one_catchment_out.csv"), index = False
    )
    overall = ratio_results[2]
    observed = overall["width_share"]["ratio_of_means"]
    genus_contributions = genus_summary["catchment_equal_width_contribution"]
    expected = float(genus_contributions.loc[GENERA].sum() / genus_contributions.sum())
    share_difference = abs(observed - expected)
    require_close(share_difference, 1e-12, "The selected genera disagree with the summary.")
    distributions = {}
    extremes = {}
    for name in overall:
        distributions[name] = distribution(catchments[name])
        extremes[name] = {
            "min" : catchments.loc[catchments[name].idxmin(), "catchment"],
            "max" : catchments.loc[catchments[name].idxmax(), "catchment"],
        }
    result = {
        "partition"          : partition,
        "catchments"         : len(catchments),
        "analyses"           : int(contributions["ANALYSIS_ID"].nunique()),
        "sites"              : int(contributions["SITE_ID"].nunique()),
        "genus_score_spans"  : spans.to_dict(),
        "overall"            : overall,
        "distributions"      : distributions,
        "min_max_catchments" : extremes,
        "checks"             : {
            "raw_contribution_matrix_max_abs_diff"                 : reconstruction_difference,
            "width_sum_vs_primary_max_abs_diff"                    : total_width_difference,
            "direct_masked_mass_vs_width_algebra_max_abs_diff"     : mass_difference,
            "total_masked_mass_vs_unrestricted_width_max_abs_diff" : total_mass_difference,
            "share_vs_original_genus_summary_abs_diff"             : share_difference,
            "five_masks_every_analysis_genus"                      : True,
            "absent_genera_zero_filled_before_site_averaging"      : True,
        },
    }
    return (catchments, result)

def write_combined(destination, results, summaries):
    development = set(results[0]["catchment"])
    confirmation = set(results[1]["catchment"])
    if (not development.isdisjoint(confirmation)):
        raise ValueError("Development and confirmation catchments must be disjoint.")
    combined = pd.concat(results, ignore_index = True)
    combined.to_csv((destination / "catchment_shares_all.csv"), index = False)
    distribution_rows = []
    overall_rows = []
    for partition in PARTITIONS:
        for metric in summaries[partition]["distributions"]:
            row = {
                "partition" : partition,
                "metric"    : metric,
            }
            row.update(summaries[partition]["distributions"][metric])
            distribution_rows.append(row)
        for metric in summaries[partition]["overall"]:
            values = summaries[partition]["overall"][metric]
            row = {
                "partition" : partition,
                "metric"    : metric,
            }
            for name in values:
                if (name != "ci95_percentile"):
                    row[name] = values[name]
            row["ci95_lower"] = values["ci95_percentile"][0]
            row["ci95_upper"] = values["ci95_percentile"][1]
            overall_rows.append(row)
    pd.DataFrame(distribution_rows).to_csv(
        (destination / "distribution_summary.csv"), index = False
    )
    pd.DataFrame(overall_rows).to_csv(
        (destination / "overall_ratios_and_influence.csv"), index = False
    )
    report = {
        "status"                         : "COMPLETE",
        "fixed_genera"                   : GENERA,
        "condition"                      : "20% requested MCAR; original coefficients; complete counts; correct genus known",
        "aggregation"                    : "Five-mask means within analyses; analyses within sites; sites within catchments",
        "ratio_definition"               : "Sum of catchment numerators divided by sum of catchment denominators",
        "bootstrap"                      : {
            "resamples"                             : RESAMPLES,
            "seed"                                  : SEED,
            "unit"                                  : "catchment",
            "paired_numerator_denominator"          : True,
            "interval"                              : "95% percentile",
            "same_seed_reinitialized_per_partition" : True,
        },
        "post_selection"                 : "Fixed genera were selected from the earlier descriptive ranking; both partitions are exploratory for concentration.",
        "no_p_values_or_enrichment_test" : True,
        "partitions"                     : summaries,
    }
    write_json((destination / "summary.json"), report)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--experiments", type = Path, required = True)
    parser.add_argument("--summaries", type = Path, required = True)
    parser.add_argument("--output", type = Path, required = True)
    arguments = parser.parse_args()
    experiments = arguments.experiments.resolve()
    summary_directory = arguments.summaries.resolve()
    destination = check_output_directory(arguments.output)
    for partition in PARTITIONS:
        check_run_receipt(experiments / partition)
        if (not (summary_directory / partition).is_dir()):
            raise FileNotFoundError(f"Missing summary partition: {partition}")
    destination.mkdir(parents = True, exist_ok = True)
    results = []
    summaries = {}
    for partition in PARTITIONS:
        result = process(experiments, summary_directory, destination, partition)
        results.append(result[0])
        summaries[partition] = result[1]
    write_combined(destination, results, summaries)
    outputs = {}
    for path in sorted(destination.iterdir()):
        if (path.is_file()):
            outputs[path.name] = file_hash(path)
    write_json((destination / "run_receipt.json"), {
        "status"  : "COMPLETE",
        "outputs" : outputs,
    })

if (__name__ == "__main__"):
    main()
