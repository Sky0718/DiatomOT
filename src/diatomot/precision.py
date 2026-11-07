import argparse
from pathlib import Path
import numpy as np
import pandas as pd

TOLERANCE = 1e-9
THRESHOLDS = [5.0, 10.0, 15.0, 20.0]
QUANTILES = [0.10, 0.25, 0.50, 0.75, 0.90, 0.95]
WIDTH_COLUMNS = {
    "genus"            : "orig_width",
    "unrestricted"     : "orig_blind_width",
    "definition_union" : "definition_union_hull_width",
}

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

def row_weights(records):
    record_keys = ["catchment", "SITE_ID", "ANALYSIS_ID"]
    masks_per_record = records.groupby(record_keys, sort = True)["ANALYSIS_ID"].transform("size")
    records_per_site = records.groupby(["catchment", "SITE_ID"], sort = True)["ANALYSIS_ID"].transform("nunique")
    sites_per_catchment = records.groupby("catchment", sort = True)["SITE_ID"].transform("nunique")
    weights = (1.0 / masks_per_record.to_numpy())
    weights = (weights / records_per_site.to_numpy())
    weights = (weights / sites_per_catchment.to_numpy())
    totals = pd.Series(weights, index = records.index).groupby(records["catchment"]).sum()
    if (np.max(np.abs(totals.to_numpy() - 1.0)) >= 1e-12):
        raise ValueError("Within-catchment weights do not sum to one.")
    return (weights)

def threshold_probabilities(values, weights):
    order = np.argsort(values, kind = "stable")
    sorted_values = values[order]
    sorted_weights = weights[order]
    cumulative_weights = np.concatenate(([0.0], np.cumsum(sorted_weights)))
    positions = np.searchsorted(sorted_values, (np.asarray(THRESHOLDS) + TOLERANCE), side = "right")
    probabilities = cumulative_weights[positions]
    return (np.clip(probabilities, 0.0, 1.0))

def inverse_quantiles(values, weights):
    order = np.argsort(values, kind = "stable")
    sorted_values = values[order]
    cumulative_weights = np.cumsum(weights[order])
    results = {}
    for quantile in QUANTILES:
        position = int(np.searchsorted(cumulative_weights, quantile, side = "left"))
        position = min(position, (len(sorted_values) - 1))
        name = "Q" + str(int(quantile * 100))
        results[name] = float(sorted_values[position])
    return (results)

def condition_rows(records, split, fraction, mechanism, views, bootstrap_indices):
    records = records.reset_index(drop = True)
    catchments = sorted(records["catchment"].unique())
    if (len(catchments) != bootstrap_indices.shape[1]):
        raise ValueError("Catchment membership changed between conditions.")
    weights = row_weights(records)
    overall_weights = (weights / len(catchments))
    if (abs(overall_weights.sum() - 1.0) >= 1e-12):
        raise ValueError("Overall weights do not sum to one.")
    identifiers = {
        "split"      : split,
        "fraction"   : fraction,
        "mechanism"  : mechanism,
        "analyses"   : int(records["ANALYSIS_ID"].nunique()),
        "sites"      : int(records["SITE_ID"].nunique()),
        "catchments" : len(catchments),
        "masks"      : len(records),
    }
    thresholds = []
    quantiles = []
    catchment_rows = []
    for view in views:
        widths = records[WIDTH_COLUMNS[view]].to_numpy(dtype = float)
        if (not np.isfinite(widths).all()):
            raise ValueError("Width values must be finite.")
        if ((widths.min() < -TOLERANCE) or (widths.max() > (100.0 + TOLERANCE))):
            raise ValueError("A width is outside the TDI range.")
        probabilities = np.zeros((len(catchments), len(THRESHOLDS)))
        for index in range(len(catchments)):
            selected = (records["catchment"].to_numpy() == catchments[index])
            probabilities[index] = threshold_probabilities(widths[selected], weights[selected])
        direct = threshold_probabilities(widths, overall_weights)
        if (np.max(np.abs(probabilities.mean(axis = 0) - direct)) >= 2e-11):
            raise ValueError("Hierarchical and directly weighted thresholds disagree.")
        quantile_record = identifiers.copy()
        quantile_record["view"] = view
        quantile_record.update(inverse_quantiles(widths, overall_weights))
        quantile_record["minimum"] = float(widths.min())
        quantile_record["maximum"] = float(widths.max())
        quantiles.append(quantile_record)
        bootstrap_means = probabilities[bootstrap_indices].mean(axis = 1)
        limits = np.quantile(bootstrap_means, [0.025, 0.975], axis = 0, method = "linear")
        for threshold_index in range(len(THRESHOLDS)):
            record = identifiers.copy()
            record["view"] = view
            record["threshold"] = THRESHOLDS[threshold_index]
            record["estimate"] = float(probabilities[:, threshold_index].mean())
            record["ci95_lower"] = float(limits[0, threshold_index])
            record["ci95_upper"] = float(limits[1, threshold_index])
            thresholds.append(record)
            for catchment_index in range(len(catchments)):
                catchment_record = identifiers.copy()
                catchment_record["view"] = view
                catchment_record["threshold"] = THRESHOLDS[threshold_index]
                catchment_record["catchment"] = catchments[catchment_index]
                catchment_record["cdf"] = float(probabilities[catchment_index, threshold_index])
                catchment_rows.append(catchment_record)
    return (thresholds, quantiles, catchment_rows)

def analyse_partition(source, split):
    columns = [
        "ANALYSIS_ID", "SITE_ID", "catchment", "fully_scored",
        "fraction_requested", "mechanism", "replicate_seed",
        "orig_width", "orig_blind_width", "definition_union_hull_width",
    ]
    raw = pd.read_csv(source / split / "mask_results.csv.gz", usecols = columns)
    records = raw.loc[raw["fully_scored"]].copy()
    if (records.empty):
        raise ValueError("No fully scored records are available for " + split)
    if (records["catchment"].isna().any()):
        raise ValueError("Catchment identifiers are missing.")
    if (records["catchment"].eq("MISSING_CATCHMENT").any()):
        raise ValueError("Missing catchments cannot enter these summaries.")
    keys = ["ANALYSIS_ID", "fraction_requested", "mechanism", "replicate_seed"]
    if (records.duplicated(keys).any()):
        raise ValueError("Duplicate masking results.")
    count = int(records["catchment"].nunique())
    generator = np.random.default_rng(20261106)
    bootstrap_indices = generator.integers(0, count, size = (5000, count))
    thresholds = []
    quantiles = []
    catchment_rows = []
    conditions = []
    for fraction in [0.05, 0.10, 0.20, 0.30, 0.50]:
        views = ["genus"]
        if (fraction == 0.20):
            views = list(WIDTH_COLUMNS)
        conditions.append((fraction, "mcar", views, 5))
    for mechanism in ["rarest_first", "low_indicator_first", "high_indicator_first"]:
        conditions.append((0.20, mechanism, ["genus"], 1))
    for condition in conditions:
        fraction = condition[0]
        mechanism = condition[1]
        selected = ((records["fraction_requested"] == fraction) & (records["mechanism"] == mechanism))
        subset = records.loc[selected]
        if (not subset.groupby("ANALYSIS_ID").size().eq(condition[3]).all()):
            raise ValueError("Unexpected replicate count.")
        result = condition_rows(subset, split, fraction, mechanism, condition[2], bootstrap_indices)
        thresholds.extend(result[0])
        quantiles.extend(result[1])
        catchment_rows.extend(result[2])
    return (thresholds, quantiles, catchment_rows)

def main():
    parser = argparse.ArgumentParser(description = "Summarize attainable TDI widths and precision thresholds.")
    parser.add_argument("--experiments", type = Path, required = True)
    parser.add_argument("--output", type = Path, required = True)
    arguments = parser.parse_args()
    source = arguments.experiments.expanduser().resolve()
    destination = output_directory(arguments.output)
    thresholds = []
    quantiles = []
    catchment_rows = []
    for split in ["development", "confirmation"]:
        result = analyse_partition(source, split)
        thresholds.extend(result[0])
        quantiles.extend(result[1])
        catchment_rows.extend(result[2])
    pd.DataFrame(thresholds).to_csv(destination / "threshold_readouts.csv", index = False)
    pd.DataFrame(quantiles).to_csv(destination / "weighted_quantiles.csv", index = False)
    pd.DataFrame(catchment_rows).to_csv(destination / "catchment_threshold_cdf.csv", index = False)

if (__name__ == "__main__"):
    main()
