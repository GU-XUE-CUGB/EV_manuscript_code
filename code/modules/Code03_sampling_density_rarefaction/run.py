#!/usr/bin/env python3
"""Run the final stage3 algorithms with packaged inputs and isolated outputs."""
from __future__ import annotations
import argparse
import importlib.util
import json
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent


def load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "lib" / (name + ".py"))
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def numeric_check(actual, reference, keys, values, tolerance=5e-10):
    joined = actual.merge(reference, on=keys, how="outer", suffixes=("_a", "_r"), indicator=True)
    assert joined._merge.eq("both").all(), "Reference key mismatch"
    errors = {}
    for value in values:
        a, r = joined[value + "_a"].to_numpy(float), joined[value + "_r"].to_numpy(float)
        assert np.array_equal(np.isnan(a), np.isnan(r)), value
        error = float(np.nanmax(np.abs(a - r)))
        assert error <= tolerance, (value, error)
        errors[value] = error
    return errors


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--validate-only", action="store_true", help="Reconstruct inventories and coefficients, without simulation")
    parser.add_argument("--iterations", type=int, default=10000, help="Formal=10000; a smaller count is explicitly marked smoke")
    args = parser.parse_args()
    if args.iterations < 2:
        parser.error("iterations must be at least 2")
    stage3 = load("stage3_sampling_density_robustness")
    depth = load("stage3_A2_measurement_depth_sensitivity")
    stage3.ITERATIONS = depth.ITERATIONS = args.iterations
    all_iqrs, contributing, zeros, measurements = stage3.load_current_baseline()
    inventory = stage3.build_inventory(contributing, measurements)
    reference = ROOT / "reference_results"
    check_values = ["contributing_record_count", "raw_measurement_count_contributing_record_bins",
        "raw_measurement_count_all_retained", "contributing_canonical_site_count", "observed_global_EV"]
    checks = {"inventory": numeric_check(inventory,
        pd.read_csv(reference / "audited_50kyr_sampling_composition_table.csv"), ["bin_center_Ma"], check_values)}
    # Validate all 16 point estimates without any bootstrap.
    density_ref = pd.read_csv(reference / "density_associations_moving_block_bootstrap.csv")
    chronological = inventory.sort_values("bin_center_Ma", ascending=False)
    errors = []
    for row in density_ref.drop_duplicates(["predictor", "scale", "method"]).itertuples(index=False):
        x, y = chronological[row.source_column].to_numpy(float), chronological.observed_global_EV.to_numpy(float)
        if row.scale == "first_difference_older_to_younger":
            x, y = np.diff(x), np.diff(y)
        errors.append(abs(stage3.correlation(x, y, row.method) - row.estimate))
    assert max(errors) < 1e-12
    checks["density_point_estimate_max_error"] = max(errors)
    m5_eligible = contributing.loc[contributing.n_measurements.ge(5)].groupby("bin_center_Ma").size()
    m5_ages = inventory.bin_center_Ma[inventory.bin_center_Ma.map(m5_eligible).fillna(0).ge(5)]
    assert len(m5_ages) == 47
    assert int(m5_ages.between(66.95 - 1e-9, 67.2 + 1e-9).sum()) == 4
    label = "validation" if args.validate_only else ("formal_10000" if args.iterations == 10000 else f"smoke_{args.iterations}")
    out = ROOT / "outputs" / label
    out.mkdir(parents=True, exist_ok=True)
    inventory.to_csv(out / "audited_50kyr_sampling_composition_table.csv", index=False)
    zeros.to_csv(out / "excluded_zero_IQR_record_bins.csv", index=False)
    phase, adjacent, composition = stage3.composition_outputs(inventory)
    for filename, frame in [("proxy_composition_phase_summary.csv", phase),
            ("proxy_composition_adjacent_changes.csv", adjacent), ("proxy_composition_EV_associations.csv", composition)]:
        frame.to_csv(out / filename, index=False)
    if not args.validate_only:
        a1, features, convergence, iterations, matrices = stage3.record_rarefaction(contributing, inventory)
        a1.to_csv(out / "A1_fixed_record_rarefaction_bin_summary.csv", index=False)
        features.to_csv(out / "A1_feature_summary.csv", index=False)
        iterations.to_csv(out / "A1_feature_iterations.csv", index=False)
        stage3.density_associations(inventory).to_csv(out / "density_associations_moving_block_bootstrap.csv", index=False)
        depth_summaries, contrast_summaries, eligibility_summaries = [], [], []
        for m in depth.M_VALUES:
            summary, contrasts, eligibility, convergence, matrix = depth.run_design(stage3, m, contributing, measurements, inventory)
            depth_summaries.append(summary)
            contrast_summaries.append(contrasts)
            eligibility_summaries.append(eligibility)
            if m == 4:
                approved = stage3.hierarchical_rarefaction(contributing, measurements, inventory)[-1]
                assert np.allclose(matrix, approved, rtol=0, atol=1e-12, equal_nan=True)
        joint = pd.concat(depth_summaries, ignore_index=True)
        joint.to_csv(out / "A2_measurement_depth_sensitivity_bin_summary.csv", index=False)
        pd.concat(contrast_summaries, ignore_index=True).to_csv(out / "A2_measurement_depth_sensitivity_contrast_summary.csv", index=False)
        pd.concat(eligibility_summaries, ignore_index=True).to_csv(out / "A2_measurement_depth_sensitivity_phase_eligibility.csv", index=False)
        checks["m4_original_stage3_algorithm_equals_depth_algorithm"] = True
        if args.iterations == 10000:
            checks["A1_formal_reference"] = numeric_check(a1,
                pd.read_csv(reference / "A1_fixed_record_rarefaction_bin_summary.csv"),
                ["N_records", "bin_center_Ma"], ["median_rarefied_EV", "p2_5", "p97_5"])
            checks["A2_formal_reference"] = numeric_check(joint,
                pd.read_csv(reference / "A2_measurement_depth_sensitivity_bin_summary.csv"),
                ["m_measurements", "bin_center_Ma"], ["median_rarefied_EV", "p2_5", "p97_5"])
    metadata = {"output_status": label, "iterations": 0 if args.validate_only else args.iterations,
        "full_formal_run": not args.validate_only and args.iterations == 10000,
        "baseline_minimum_measurements": 2, "baseline_nonzero_record_bins": len(contributing),
        "record_rarefaction_N": [3, 5, 7], "record_rarefaction_seeds": stage3.SEEDS_A1,
        "joint_rarefaction_m": [3, 4, 5], "joint_rarefaction_N": 5, "joint_seed_reset_each_m": depth.SEED,
        "block_lengths": list(stage3.BLOCK_LENGTHS), "primary_block": stage3.PRIMARY_BLOCK,
        "block_seed_formula": "20260900 + predictor_index*100 + scale_offset(0/40) + method_index*10 + block_length",
        "density_methods": ["Pearson", "Spearman"], "density_scales": ["raw_bins", "first_difference_older_to_younger"],
        "sampling": "without replacement at each rarefaction level; circular moving-block pairs for density CI",
        "subsample_zscores": "original record-wide sample-SD z-scores retained; zero subsample IQR outcomes retained",
        "m5_eligible_bins": 47, "m5_pre_escalation_eligible_bins": "4/6", "checks": checks}
    (out / "verification.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
