#!/usr/bin/env python3
"""Standardize retained measurements and calculate record-bin IQRs and EV."""
from __future__ import annotations
import importlib.util
import json
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent


def load_library():
    path = ROOT / "lib" / "integrate_new_proxy_data_v3.py"
    spec = importlib.util.spec_from_file_location("ev_numerics", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def compare_numeric(actual, reference, keys, columns, tolerance=5e-10):
    joined = actual.merge(reference, on=keys, how="outer", suffixes=("_new", "_reference"), indicator=True, validate="one_to_one")
    if not joined._merge.eq("both").all():
        raise AssertionError(f"Reference keys differ: {joined._merge.value_counts().to_dict()}")
    errors = {}
    for column in columns:
        a = joined[column + "_new"].to_numpy(float)
        b = joined[column + "_reference"].to_numpy(float)
        if not np.array_equal(np.isnan(a), np.isnan(b)):
            raise AssertionError(f"Missing-value pattern differs for {column}")
        error = float(np.nanmax(np.abs(a - b)))
        errors[column] = error
        if error > tolerance:
            raise AssertionError(f"{column}: {error} exceeds {tolerance}")
    return errors


def main():
    numerical = load_library()
    source = pd.read_csv(ROOT / "inputs" / "recomputed_measurements_with_zscores_v3.csv", low_memory=False)
    wang = (source.category.eq("pCO2") & source.site.eq("Nanxiong Basin")
            & source.source_DOI.eq("https://doi.org/10.1016/j.palaeo.2024.112177")
            & source.source_status.eq("new_approved_proxy"))
    assert int(wang.sum()) == 64
    source.loc[wang, "dedup_excluded"] = True
    source.loc[wang, "source_status"] = "excluded_derivative_overlap"
    source.loc[wang, "exclusion_reason"] = "pCO2 derived from retained Nanxiong carbonate d13C"
    retained = numerical.recompute_zscores(source.loc[~source.dedup_excluded].copy())
    assert len(retained) == 10881 and retained.record_id.nunique() == 70
    assert retained.loc[retained.category.eq("pCO2"), "record_id"].unique().tolist() == ["pCO2::all data"]
    assert int(retained.category.eq("pCO2").sum()) == 67
    iqrs, binned = numerical.calculate_iqrs(retained, min_n=2)
    n4_sensitivity, _ = numerical.calculate_iqrs(retained, min_n=4, schemes=["50"])
    legacy_global, legacy_proxy = numerical.aggregate_iqrs(iqrs)
    # Calculate EV from positive record-bin IQRs.
    final_global, final_proxy = numerical.aggregate_iqrs(iqrs.loc[iqrs.IQR.gt(0)].copy())
    mean_z = numerical.aggregate_zscore(binned)
    references = ROOT / "reference_results"
    ref_iqr = pd.read_csv(references / "record_level_IQR_all_bin_widths_revised_v4_n2.csv")
    for frame in (iqrs, ref_iqr):
        frame["scheme"] = frame.scheme.astype(str)
    checks = {"record_IQR": compare_numeric(iqrs, ref_iqr,
               ["scheme", "record_id", "bin_center_Ma"], ["IQR", "n_measurements"])}
    reference_retained = pd.read_csv(references / "retained_measurements_with_zscores_v4.csv", low_memory=False)
    checks["standardization"] = compare_numeric(retained, reference_retained,
            ["record_id", "sample_id", "source_row", "age_Ma", "source_file"], ["z_score_recomputed", "bin_center_Ma"])
    ref_global = pd.read_csv(references / "global_EV_revised_v4_all_bin_widths.csv")
    for frame in (legacy_global, ref_global):
        frame["scheme"] = frame.scheme.astype(str)
    checks["v4_including_zero_mean"] = compare_numeric(legacy_global, ref_global,
            ["scheme", "bin_center_Ma"], ["revised_EV", "contributing_record_count"])
    raw = pd.read_csv(references / "Global_EV_50kyr_raw_values.csv")
    fifty = final_global.loc[final_global.scheme.astype(str).eq("50")].rename(
            columns={"bin_center_Ma": "Age_Ma", "revised_EV": "Global_EV"})
    checks["final_50kyr_nonzero_mean"] = compare_numeric(fifty, raw, ["Age_Ma"], ["Global_EV"])
    pre = fifty.loc[fifty.Age_Ma.between(66.95 - 1e-8, 67.2 + 1e-8), "Global_EV"]
    post = fifty.loc[fifty.Age_Ma.between(64.8 - 1e-8, 65.45 + 1e-8), "Global_EV"]
    assert len(pre) == 6 and len(post) == 14
    background = float((pre.mean() + post.mean()) / 2)
    assert abs(background - 0.34503257945945265) < 1e-11
    out = ROOT / "outputs" / "generated"
    out.mkdir(parents=True, exist_ok=True)
    tables = {"retained_measurements_recomputed.csv": retained,
        "record_level_IQR_all_schemes_n2.csv": iqrs,
        "record_level_IQR_50kyr_n4_sensitivity.csv": n4_sensitivity,
        "global_EV_v4_including_zero_reference_equivalent.csv": legacy_global,
        "proxy_EV_v4_including_zero_reference_equivalent.csv": legacy_proxy,
        "global_EV_final_nonzero_mean.csv": final_global,
        "proxy_EV_final_nonzero_mean.csv": final_proxy,
        "global_mean_zscore.csv": mean_z}
    for filename, table in tables.items():
        table.to_csv(out / filename, index=False)
    metadata = {"output_status": "newly_recomputed_deterministic",
        "input_handoff": "compiled sample-level v3 inventory; raw values are restandardized",
        "retained_measurements": len(retained), "operational_records": 70,
        "minimum_measurements_per_record_bin_final": 2,
        "n4_status": "sensitivity only; does not replace final n2 baseline",
        "sample_SD_ddof": 1, "IQR_quantile": "linear; R type 7",
        "pCO2_pooled_measurements": 67, "Wang_derived_pCO2_excluded": 64,
        "zero_IQR_final_mean_rule": "exclude IQR=0 before means",
        "background_equal_window_mean": background,
        "background_pre_window_Ma": [66.95, 67.2], "background_post_window_Ma": [64.8, 65.45],
        "checks": checks}
    (out / "verification.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
