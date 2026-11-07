import argparse
import hashlib
import json
import re
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
import pyreadr

DIGEST_METHODS = ("ANCA", "ANCL", "ANCO", "ANHO", "ANNI", "ANSU")
UNNAMED_TAXA = {
    "Unknown",
    "Unidentified",
    "Diatom",
    "Diatoms",
    "Bacillariophyta",
    "Bacillariophyceae",
}

def file_digest(path):
    digest = hashlib.sha256()
    with (path.open("rb")) as stream:
        while (True):
            block = stream.read(1048576)
            if (not block):
                break
            digest.update(block)
    return (digest.hexdigest())

def read_text_table(path):
    table = pd.read_csv(path, dtype = str, keep_default_na = False)
    return (table)

def require_columns(table, columns, name):
    missing = []
    for column in columns:
        if (column not in table.columns):
            missing.append(column)
    if (missing):
        raise ValueError(f"Missing columns in {name}: {', '.join(missing)}")

def require_unique(table, column, name):
    if (table[column].duplicated().any()):
        raise ValueError(f"Repeated {column} in {name}")

def read_coefficients(path):
    contents = pyreadr.read_r(str(path))
    if ("darleq3_taxa" not in contents):
        raise ValueError("The coefficient file must contain the darleq3_taxa table")
    dictionary = contents["darleq3_taxa"]
    required = ["NBSCode", "TaxonName", "ValidName", "TDI5LM", "TDI5LM.orig"]
    require_columns(dictionary, required, "coefficient table")
    keyed = dictionary.loc[dictionary["NBSCode"].notna()].copy()
    require_unique(keyed, "NBSCode", "coefficient table")
    for column in ["TDI5LM", "TDI5LM.orig"]:
        scores = pd.to_numeric(keyed[column], errors = "raise")
        if (not scores.dropna().between(1, 5).all()):
            raise ValueError(f"Scores outside 1 to 5 in {column}")
    return (dictionary)

def read_counts(path, metrics, names, dictionary):
    columns = [
        "SITE_ID",
        "SAMPLE_ID",
        "ANALYSIS_ID",
        "ANALYSIS_TYPE",
        "ANALYSIS_METHOD",
        "TAXON_LIST_ITEM_KEY",
        "UNITS_FOUND_OR_SEQUENCE_READS",
    ]
    counts = pd.read_csv(path, dtype = str, keep_default_na = False, usecols = columns)
    row_numbers = np.arange(1, (len(counts) + 1))
    counts.insert(0, "taxa_data_record_1based", row_numbers)
    counts["count"] = pd.to_numeric(counts["UNITS_FOUND_OR_SEQUENCE_READS"], errors = "coerce")
    nonfinite = (~np.isfinite(counts["count"]))
    negative = (counts["count"] < 0)
    fractional = (counts["count"] != counts["count"].round())
    counts["invalid_count"] = (nonfinite | negative | fractional)
    counts["duplicate_analysis_taxon"] = counts.duplicated(
        ["ANALYSIS_ID", "TAXON_LIST_ITEM_KEY"], keep = False
    )
    identity_columns = ["SITE_ID", "SAMPLE_ID", "ANALYSIS_TYPE", "ANALYSIS_METHOD"]
    metric_columns = ["ANALYSIS_ID"]
    metric_columns.extend(identity_columns)
    identities = counts.merge(
        metrics[metric_columns],
        on = "ANALYSIS_ID",
        how = "left",
        validate = "many_to_one",
        suffixes = ("", "_metrics"),
        indicator = True,
        sort = False,
    )
    mismatch = identities["_merge"].ne("both")
    for column in identity_columns:
        metric_column = (column + "_metrics")
        different = identities[column].ne(identities[metric_column])
        mismatch = (mismatch | different)
    counts["identity_mismatch"] = mismatch.to_numpy()
    coefficients = dictionary.dropna(subset = ["NBSCode"]).set_index("NBSCode")
    counts["tdi5lm_indicator"] = counts["TAXON_LIST_ITEM_KEY"].map(coefficients["TDI5LM"])
    counts["tdi5lm_orig_indicator"] = counts["TAXON_LIST_ITEM_KEY"].map(coefficients["TDI5LM.orig"])
    counts["in_dictionary"] = counts["TAXON_LIST_ITEM_KEY"].isin(coefficients.index)
    counts["in_taxon_info"] = counts["TAXON_LIST_ITEM_KEY"].isin(names["TAXON_LIST_ITEM_KEY"])
    counts["scoring_count"] = counts["count"].where(counts["tdi5lm_indicator"].notna(), 0)
    counts["indicator_count"] = (counts["scoring_count"] * counts["tdi5lm_indicator"].fillna(0))
    counts["orig_scoring_count"] = counts["count"].where(counts["tdi5lm_orig_indicator"].notna(), 0)
    counts["orig_indicator_count"] = (
        counts["orig_scoring_count"] * counts["tdi5lm_orig_indicator"].fillna(0)
    )
    counts["unmapped_count"] = counts["count"].where((~counts["in_dictionary"]), 0)
    counts["non_scoring_count"] = counts["count"].where(counts["tdi5lm_indicator"].isna(), 0)
    return (counts)

def make_inventory(metrics, sites, counts):
    totals = counts.groupby("ANALYSIS_ID", sort = False).agg(
        taxa_rows = ("count", "size"),
        raw_count_sum = ("count", "sum"),
        scoring_count_sum = ("scoring_count", "sum"),
        weighted_indicator_sum = ("indicator_count", "sum"),
        orig_scoring_count_sum = ("orig_scoring_count", "sum"),
        orig_weighted_indicator_sum = ("orig_indicator_count", "sum"),
        dictionary_unmapped_count = ("unmapped_count", "sum"),
        non_scoring_count = ("non_scoring_count", "sum"),
        invalid_count_rows = ("invalid_count", "sum"),
        duplicate_taxon_rows = ("duplicate_analysis_taxon", "sum"),
        identity_mismatch_rows = ("identity_mismatch", "sum"),
    ).reset_index()
    inventory = metrics.merge(totals, how = "left", on = "ANALYSIS_ID", validate = "one_to_one")
    site_columns = [
        "SITE_ID",
        "ALKALINITY",
        "BASE_DATA_DATE",
        "CATCHMENT",
        "WFD_WATERBODY_ID",
        "FULL_EASTING",
        "FULL_NORTHING",
    ]
    inventory = inventory.merge(
        sites[site_columns], on = "SITE_ID", how = "left", validate = "many_to_one"
    )
    raw_total = pd.to_numeric(inventory["TOTAL_UNITS_FOUND_LM"], errors = "coerce")
    scoring_total = pd.to_numeric(inventory["TOTAL_TDI5LM_CELLS"], errors = "coerce")
    archive_index = pd.to_numeric(inventory["TDI5LM"], errors = "coerce")
    inventory["raw_total_matches"] = np.isclose(
        inventory["raw_count_sum"], raw_total, rtol = 0, atol = 1e-9
    )
    inventory["scoring_total_matches"] = np.isclose(
        inventory["scoring_count_sum"], scoring_total, rtol = 0, atol = 1e-9
    )
    inventory["recomputed_TDI5LM"] = (
        ((25 * inventory["weighted_indicator_sum"]) / inventory["scoring_count_sum"]) - 25
    )
    inventory["TDI5LM_abs_error"] = (inventory["recomputed_TDI5LM"] - archive_index).abs()
    inventory["TDI5LM_rounding_agreement"] = inventory["TDI5LM_abs_error"].le(0.005000001)
    inventory["recomputed_TDI5LM_orig"] = (
        ((25 * inventory["orig_weighted_indicator_sum"]) / inventory["orig_scoring_count_sum"]) - 25
    )
    inventory["TDI5LM_orig_abs_error"] = (inventory["recomputed_TDI5LM_orig"] - archive_index).abs()
    inventory["TDI5LM_orig_rounding_agreement"] = inventory["TDI5LM_orig_abs_error"].le(0.005000001)
    inventory["orig_scoring_total_matches"] = np.isclose(
        inventory["orig_scoring_count_sum"], scoring_total, rtol = 0, atol = 1e-9
    )
    inventory["scoring_fraction_all_valves"] = (
        inventory["scoring_count_sum"] / inventory["raw_count_sum"]
    )
    inventory["sample_date_parsed"] = pd.to_datetime(
        inventory["SAMPLE_DATE"], dayfirst = True, errors = "coerce"
    )
    inventory["analysis_date_parsed"] = pd.to_datetime(
        inventory["DATE_OF_ANALYSIS"], dayfirst = True, errors = "coerce"
    )
    inventory["analysis_precedes_sample"] = inventory["analysis_date_parsed"].lt(
        inventory["sample_date_parsed"]
    )
    inventory["site_alkalinity_numeric"] = pd.to_numeric(inventory["ALKALINITY"], errors = "coerce")
    inventory["site_alkalinity_present_positive"] = inventory["site_alkalinity_numeric"].gt(0)
    inventory["sample_record_id"] = (inventory["SAMPLE_ID"] + ":v" + inventory["SAMPLE_VERSION"])
    inventory["analysis_id"] = ("EA_BIOSYS:" + inventory["ANALYSIS_ID"])
    inventory["sample_id"] = ("EA_BIOSYS:" + inventory["sample_record_id"])
    inventory["site_cluster_id"] = ("EA_BIOSYS_SITE:" + inventory["SITE_ID"])
    inventory["microscopy_digest_method"] = (
        inventory["ANALYSIS_TYPE"].eq("LABP")
        & inventory["ANALYSIS_METHOD"].isin(DIGEST_METHODS)
    )
    inventory["data_integrity_ok"] = (
        inventory["raw_total_matches"]
        & inventory["taxa_rows"].gt(0)
        & inventory["invalid_count_rows"].eq(0)
        & inventory["duplicate_taxon_rows"].eq(0)
        & inventory["identity_mismatch_rows"].eq(0)
    )
    inventory["candidate_standard_counts"] = (
        inventory["microscopy_digest_method"]
        & inventory["SAMPLE_METHOD"].eq("SDIA")
        & inventory["data_integrity_ok"]
        & inventory["sample_date_parsed"].ge(pd.Timestamp("2006-01-01"))
        & inventory["scoring_count_sum"].ge(300)
    )
    return (inventory)

def lexical_genus(value):
    if (not isinstance(value, str)):
        return (None)
    match = re.match(r"^([A-Z][a-z]+)(?:\s|$)", value.strip())
    if (match is None):
        return (None)
    genus = match.group(1)
    if (genus in UNNAMED_TAXA):
        return (None)
    return (genus)

def make_selected_counts(counts, candidates, names, dictionary):
    columns = [
        "taxa_data_record_1based",
        "ANALYSIS_ID",
        "SAMPLE_ID",
        "SITE_ID",
        "TAXON_LIST_ITEM_KEY",
        "count",
        "tdi5lm_indicator",
        "tdi5lm_orig_indicator",
        "in_dictionary",
        "in_taxon_info",
    ]
    selected = counts.loc[counts["ANALYSIS_ID"].isin(candidates["ANALYSIS_ID"]), columns].copy()
    coefficients = dictionary.dropna(subset = ["NBSCode"]).set_index("NBSCode")
    taxon_names = names.set_index("TAXON_LIST_ITEM_KEY")
    selected["taxon_name_biosys"] = selected["TAXON_LIST_ITEM_KEY"].map(taxon_names["TAXON_NAME"])
    selected["darleq_taxon_name"] = selected["TAXON_LIST_ITEM_KEY"].map(coefficients["TaxonName"])
    selected["darleq_valid_name"] = selected["TAXON_LIST_ITEM_KEY"].map(coefficients["ValidName"])
    selected["genus_darleq_lexical"] = selected["darleq_valid_name"].map(lexical_genus)
    selected["genus_rule"] = "first_capitalized_Latin_token_of_DARLEQ_ValidName_no_synonym_merging"
    selected["v_current"] = np.where(selected["tdi5lm_indicator"].notna(), 1.0, np.nan)
    selected["v_orig"] = np.where(selected["tdi5lm_orig_indicator"].notna(), 1.0, np.nan)
    selected["analysis_id"] = ("EA_BIOSYS:" + selected["ANALYSIS_ID"])
    samples = candidates.set_index("ANALYSIS_ID")["sample_id"]
    selected["sample_id"] = selected["ANALYSIS_ID"].map(samples)
    selected["site_cluster_id"] = ("EA_BIOSYS_SITE:" + selected["SITE_ID"])
    return (selected)

def summarize_eligibility(inventory):
    conditions = {
        "LABP"                        : inventory["ANALYSIS_TYPE"].eq("LABP"),
        "specified_digest_methods"    : inventory["microscopy_digest_method"],
        "SDIA_sampling"               : inventory["SAMPLE_METHOD"].eq("SDIA"),
        "count_identity_integrity"    : inventory["data_integrity_ok"],
        "sample_date_2006_onwards"    : inventory["sample_date_parsed"].ge(pd.Timestamp("2006-01-01")),
        "at_least_300_scoring_valves" : inventory["scoring_count_sum"].ge(300),
    }
    keep = pd.Series(True, index = inventory.index)
    rows = []
    for rule in conditions:
        keep = (keep & conditions[rule])
        subset = inventory.loc[keep]
        catchments = subset["CATCHMENT"].replace("", np.nan)
        row = {
            "rule_cumulative" : rule,
            "analyses"        : len(subset),
            "samples"         : int(subset["SAMPLE_ID"].nunique()),
            "sites"           : int(subset["SITE_ID"].nunique()),
            "catchments"      : int(catchments.nunique()),
        }
        rows.append(row)
    return (pd.DataFrame(rows))

def prepare_output_directory(path):
    destination = path.expanduser().resolve()
    repository = Path(__file__).resolve().parents[2]
    if (destination.is_relative_to(repository)):
        raise ValueError("Choose an output directory outside this repository")
    if (destination.exists()):
        if (not destination.is_dir()):
            raise NotADirectoryError(destination)
        if (any(destination.iterdir())):
            raise FileExistsError(f"Output directory is not empty: {destination}")
    return (destination)

def main():
    parser = argparse.ArgumentParser(description = "Prepare BIOSYS counts and DARLEQ3 coefficients for analysis")
    parser.add_argument("--data-dir", type = Path, required = True)
    parser.add_argument("--coefficients", type = Path, required = True)
    parser.add_argument("--output", type = Path, required = True)
    arguments = parser.parse_args()
    destination = prepare_output_directory(arguments.output)
    data_directory = arguments.data_dir.expanduser().resolve()
    coefficient_path = arguments.coefficients.expanduser().resolve()
    inputs = {
        "metrics"      : (data_directory / "DIAT_OPEN_DATA_METRICS.csv"),
        "sites"        : (data_directory / "DIAT_OPEN_DATA_SITE.csv"),
        "counts"       : (data_directory / "DIAT_OPEN_DATA_TAXA.csv"),
        "taxon_info"   : (data_directory / "OPEN_DATA_TAXON_INFO.csv"),
        "coefficients" : coefficient_path,
    }
    input_hashes = {}
    for label in inputs:
        path = inputs[label]
        if (not path.is_file()):
            raise FileNotFoundError(path)
        input_hashes[label] = file_digest(path)
    metrics = read_text_table(inputs["metrics"])
    sites = read_text_table(inputs["sites"])
    names = read_text_table(inputs["taxon_info"])
    require_unique(metrics, "ANALYSIS_ID", "metrics")
    require_unique(sites, "SITE_ID", "sites")
    require_unique(names, "TAXON_LIST_ITEM_KEY", "taxon information")
    metrics.insert(0, "metrics_data_record_1based", np.arange(1, (len(metrics) + 1)))
    dictionary = read_coefficients(coefficient_path)
    counts = read_counts(inputs["counts"], metrics, names, dictionary)
    inventory = make_inventory(metrics, sites, counts)
    candidates = inventory.loc[inventory["candidate_standard_counts"]].copy()
    if (candidates.empty):
        raise ValueError("No analyses meet the sampling, count and date criteria")
    selected = make_selected_counts(counts, candidates, names, dictionary)
    stages = summarize_eligibility(inventory)
    destination.mkdir(parents = True, exist_ok = True)
    dictionary.to_csv((destination / "darleq3_rda_coefficients.csv"), index = False)
    candidates.to_csv((destination / "candidate_sample_metadata.csv"), index = False)
    selected.to_csv(
        (destination / "candidate_sample_taxon_long.csv.gz"),
        index = False,
        compression = {
            "method" : "gzip",
            "mtime"  : 0,
        },
    )
    inventory.to_csv((destination / "analysis_inventory.csv"), index = False)
    stages.to_csv((destination / "eligibility_summary.csv"), index = False)
    shutil.copyfile(inputs["taxon_info"], (destination / "OPEN_DATA_TAXON_INFO.csv"))
    output_hashes = {}
    for path in sorted(destination.iterdir()):
        output_hashes[path.name] = file_digest(path)
    fully_scored = np.isclose(
        candidates["raw_count_sum"], candidates["scoring_count_sum"], rtol = 0, atol = 0
    )
    summary = {
        "analyses"              : len(candidates),
        "fully_scored_analyses" : int(fully_scored.sum()),
        "sites"                 : int(candidates["SITE_ID"].nunique()),
        "catchments"            : int(candidates["CATCHMENT"].nunique()),
        "taxon_rows"            : len(selected),
        "numpy"                 : np.__version__,
        "pandas"                : pd.__version__,
        "pyreadr"               : pyreadr.__version__,
        "input_sha256"          : input_hashes,
        "output_sha256"         : output_hashes,
    }
    result_path = (destination / "preparation.json")
    result_path.write_text(json.dumps(summary, indent = 2), encoding = "utf-8")

if (__name__ == "__main__"):
    main()
