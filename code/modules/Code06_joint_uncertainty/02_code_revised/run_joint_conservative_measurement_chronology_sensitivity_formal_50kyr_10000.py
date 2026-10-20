#!/usr/bin/env python3
"""Formal joint conservative measurement + chronology sensitivity analysis.

This wrapper deliberately reuses the validated 1,000-simulation joint pilot
engine without changing any scientific parameter or execution rule. It expands
the locked analysis to 10,000 simulations, verifies that its first 1,000 rows
reproduce the pilot pointwise, and adds convergence and formal reporting tables.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd


MASTER_SEED = 20260824
ITERATIONS = 10000
PREFIXES = (1000, 5000, 10000)

SCRIPT = Path(__file__).resolve()
WORK_ROOT = SCRIPT.parent.parent
INVENTORY_DIR = WORK_ROOT / "03_results" / "EV_uncertainty_inventory"
PILOT_SCRIPT = WORK_ROOT / "02_code_revised" / "run_joint_conservative_measurement_chronology_sensitivity_pilot_50kyr_1000.py"
PILOT_DIR = INVENTORY_DIR / "Joint_conservative_measurement_chronology_sensitivity_pilot_50kyr_1000"
PILOT_NPZ = PILOT_DIR / "joint_complete_realizations_50kyr_1000.npz"
PILOT_WORKBOOK = PILOT_DIR / "Joint_conservative_measurement_chronology_sensitivity_pilot_50kyr_1000.xlsx"
RESULT_DIR = INVENTORY_DIR / "Joint_conservative_measurement_chronology_sensitivity_formal_50kyr_10000"
RESULT_DIR.mkdir(parents=True, exist_ok=True)

PHASE_NAMES = [
    "Phase 1 - pre-escalation",
    "Phase 2 - pre-escalation",
    "Recovery - Phase 2",
]


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def bin_summary(ev: np.ndarray, ages: np.ndarray, observed: np.ndarray) -> pd.DataFrame:
    q = np.quantile(ev, [0.025, 0.5, 0.975], axis=0)
    result = pd.DataFrame({
        "midpoint_Ma": ages,
        "observed_EV": observed,
        "simulation_median": q[1],
        "simulation_2.5th": q[0],
        "simulation_97.5th": q[2],
        "interval_width": q[2] - q[0],
        "observed_minus_simulation_median": observed - q[1],
        "valid_simulations": np.sum(np.isfinite(ev), axis=0),
    })
    return result.sort_values("midpoint_Ma", ascending=False).reset_index(drop=True)


def contrast_summary(contrasts: np.ndarray) -> pd.DataFrame:
    q = np.quantile(contrasts, [0.025, 0.5, 0.975], axis=0)
    return pd.DataFrame({
        "contrast": PHASE_NAMES,
        "simulation_median": q[1],
        "simulation_2.5th": q[0],
        "simulation_97.5th": q[2],
        "interval_width": q[2] - q[0],
        "valid_simulations": np.sum(np.isfinite(contrasts), axis=0),
        "interval_entirely_above_zero": q[0] > 0,
        "interval_entirely_below_zero": q[2] < 0,
    })


def peak_frequency(peaks: np.ndarray, prefix: int) -> pd.DataFrame:
    counts = pd.Series(peaks[:prefix]).value_counts().sort_index()
    return pd.DataFrame({
        "simulations": prefix,
        "peak_midpoint_Ma": counts.index.astype(float),
        "frequency": counts.values.astype(int),
        "proportion": counts.values.astype(float) / prefix,
    })


def eligibility_by_bin(
    counts: np.ndarray,
    ages: np.ndarray,
    observed_counts: np.ndarray,
    prefix: int,
) -> pd.DataFrame:
    values = counts[:prefix]
    return pd.DataFrame({
        "simulations": prefix,
        "midpoint_Ma": ages,
        "observed_eligible_record_bins": observed_counts,
        "simulation_min": np.min(values, axis=0),
        "simulation_2.5th": np.quantile(values, 0.025, axis=0),
        "simulation_median": np.median(values, axis=0),
        "simulation_97.5th": np.quantile(values, 0.975, axis=0),
        "simulation_max": np.max(values, axis=0),
        "simulations_with_empty_bin": np.sum(values == 0, axis=0),
    }).sort_values("midpoint_Ma", ascending=False).reset_index(drop=True)


def total_eligibility_summary(counts: np.ndarray, observed_counts: np.ndarray, prefix: int) -> dict:
    total = counts[:prefix].sum(axis=1)
    return {
        "simulations": prefix,
        "metric": "Total eligible record-bins per simulation",
        "observed": int(observed_counts.sum()),
        "minimum": int(total.min()),
        "2.5th": float(np.quantile(total, 0.025)),
        "median": float(np.median(total)),
        "97.5th": float(np.quantile(total, 0.975)),
        "maximum": int(total.max()),
    }


def prefix_metrics(
    ev: np.ndarray,
    contrasts: np.ndarray,
    peaks: np.ndarray,
    ages: np.ndarray,
    observed: np.ndarray,
    prefix: int,
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    bins = bin_summary(ev[:prefix], ages, observed)
    phase = contrast_summary(contrasts[:prefix])
    widths = bins.interval_width.to_numpy(float)
    difference = np.abs(bins.observed_minus_simulation_median.to_numpy(float))
    peak_counts = pd.Series(peaks[:prefix]).value_counts()
    modal_peak = float(peak_counts.index[0])
    row = {
        "simulations": prefix,
        "median_interval_width": float(np.median(widths)),
        "maximum_interval_width": float(np.max(widths)),
        "maximum_width_midpoint_Ma": float(bins.iloc[int(np.argmax(widths))].midpoint_Ma),
        "maximum_abs_observed_minus_median": float(np.max(difference)),
        "maximum_difference_midpoint_Ma": float(bins.iloc[int(np.argmax(difference))].midpoint_Ma),
        "modal_peak_midpoint_Ma": modal_peak,
        "modal_peak_frequency": int(peak_counts.loc[modal_peak]),
        "modal_peak_proportion": float(peak_counts.loc[modal_peak] / prefix),
        "all_49_bins_finite_in_every_simulation": bool(np.isfinite(ev[:prefix]).all()),
    }
    return bins, phase, row


def pilot_reproduction_audit(formal: dict, pilot: dict) -> pd.DataFrame:
    rows = []
    numeric_keys = [
        "EV_50kyr",
        "phase_contrasts",
        "peak_midpoint_Ma",
        "eligible_record_bin_counts",
        "measurement_delta_signature",
    ]
    numeric_keys += sorted(key for key in pilot if key.startswith("shift_"))
    for key in numeric_keys:
        reference = pilot[key]
        candidate = formal[key][: len(reference)]
        same_shape = candidate.shape == reference.shape
        exact = bool(same_shape and np.array_equal(candidate, reference, equal_nan=True))
        if same_shape:
            delta = np.abs(candidate.astype(float) - reference.astype(float))
            max_abs = float(np.nanmax(delta)) if delta.size else 0.0
        else:
            max_abs = math.nan
        rows.append({
            "array": key,
            "pilot_shape": str(reference.shape),
            "formal_prefix_shape": str(candidate.shape),
            "exact_pointwise_match": exact,
            "maximum_absolute_difference": max_abs,
        })
    result = pd.DataFrame(rows)
    if not bool(result.exact_pointwise_match.all()):
        raise RuntimeError("Formal first 1,000 simulations do not reproduce the validated pilot pointwise")
    return result


def comparisons_to_single_components(formal_bins: pd.DataFrame, formal_contrasts: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    measurement_dir = INVENTORY_DIR / "Layer2_conservative_measurement_sensitivity_formal_50kyr_10000"
    chronology_dir = INVENTORY_DIR / "Chronology_LayerB_conservative_sensitivity_formal_50kyr_10000"
    scenarios = [
        (
            "Measurement B formal 10,000",
            pd.read_csv(measurement_dir / "Layer2_formal_50kyr_Global_EV_summary.csv"),
            "MC_median", "MC_2.5th", "MC_97.5th",
        ),
        (
            "Chronology B formal 10,000",
            pd.read_csv(chronology_dir / "Chronology_LayerB_formal_50kyr_Global_EV_summary.csv"),
            "simulation_median", "simulation_2.5th", "simulation_97.5th",
        ),
        (
            "Joint B+B formal 10,000",
            formal_bins,
            "simulation_median", "simulation_2.5th", "simulation_97.5th",
        ),
    ]
    envelope_rows = []
    for label, table, median_col, low_col, high_col in scenarios:
        width = table[high_col].to_numpy(float) - table[low_col].to_numpy(float)
        midpoint = table.midpoint_Ma.to_numpy(float)
        observed = table.observed_EV.to_numpy(float)
        median = table[median_col].to_numpy(float)
        envelope_rows.append({
            "scenario": label,
            "median_interval_width": float(np.median(width)),
            "maximum_interval_width": float(np.max(width)),
            "maximum_width_midpoint_Ma": float(midpoint[int(np.argmax(width))]),
            "maximum_abs_observed_minus_median": float(np.max(np.abs(observed - median))),
            "maximum_difference_midpoint_Ma": float(midpoint[int(np.argmax(np.abs(observed - median)))]),
        })

    measurement_contrast = pd.read_csv(measurement_dir / "Layer2_formal_phase_contrasts_50kyr.csv")
    chronology_contrast = pd.read_csv(chronology_dir / "Chronology_LayerB_formal_phase_contrasts.csv").rename(columns={
        "simulation_median": "MC_median",
        "simulation_2.5th": "MC_2.5th",
        "simulation_97.5th": "MC_97.5th",
    })
    joint_contrast = formal_contrasts.rename(columns={
        "simulation_median": "MC_median",
        "simulation_2.5th": "MC_2.5th",
        "simulation_97.5th": "MC_97.5th",
    })
    contrast_rows = []
    for label, table in [
        ("Measurement B formal 10,000", measurement_contrast),
        ("Chronology B formal 10,000", chronology_contrast),
        ("Joint B+B formal 10,000", joint_contrast),
    ]:
        for _, row in table.iterrows():
            contrast_rows.append({
                "scenario": label,
                "contrast": row.contrast,
                "median": float(row.MC_median),
                "2.5th": float(row["MC_2.5th"]),
                "97.5th": float(row["MC_97.5th"]),
                "interval_width": float(row["MC_97.5th"] - row["MC_2.5th"]),
            })
    return pd.DataFrame(envelope_rows), pd.DataFrame(contrast_rows)


def main() -> None:
    if not PILOT_NPZ.exists():
        raise FileNotFoundError(PILOT_NPZ)

    pilot_engine = load_module(PILOT_SCRIPT, "joint_bb_formal_locked_pilot_engine")
    pilot_engine.ITERATIONS = ITERATIONS
    pilot_engine.RESULT_DIR = RESULT_DIR

    # This creates the 10,000-realization arrays with the exact validated pilot
    # parameter mapping and RNG stream construction.
    pilot_engine.main()

    generated_npz = RESULT_DIR / "joint_complete_realizations_50kyr_1000.npz"
    if not generated_npz.exists():
        raise RuntimeError("Locked engine did not create the expected complete realization archive")

    formal_npz_path = RESULT_DIR / "Joint_BB_complete_realizations_50kyr_10000.npz"
    generated_npz.replace(formal_npz_path)
    formal_npz = dict(np.load(formal_npz_path))
    pilot_npz = dict(np.load(PILOT_NPZ))
    reproduction = pilot_reproduction_audit(formal_npz, pilot_npz)

    ev = formal_npz["EV_50kyr"]
    ages = formal_npz["age_50kyr"]
    observed = formal_npz["observed_EV_50kyr"]
    contrasts = formal_npz["phase_contrasts"]
    peaks = formal_npz["peak_midpoint_Ma"]
    eligible = formal_npz["eligible_record_bin_counts"]
    if ev.shape != (ITERATIONS, 49) or contrasts.shape != (ITERATIONS, 3):
        raise RuntimeError(f"Unexpected formal realization shapes: EV {ev.shape}, contrasts {contrasts.shape}")

    final_bins = bin_summary(ev, ages, observed)
    final_contrasts = contrast_summary(contrasts)
    final_peaks = peak_frequency(peaks, ITERATIONS)
    final_eligibility = eligibility_by_bin(eligible, ages, formal_npz["eligible_record_bin_counts"][0] * 0 + pilot_npz["eligible_record_bin_counts"][0] * 0, ITERATIONS)

    # Recover the independently validated observed eligible count vector from
    # the locked engine's generated 10,000 summary table.
    engine_eligibility_path = RESULT_DIR / "joint_eligibility_by_bin.csv"
    engine_eligibility = pd.read_csv(engine_eligibility_path)
    observed_counts_lookup = dict(zip(engine_eligibility.midpoint_Ma, engine_eligibility.observed_eligible_record_bins))
    observed_counts = np.array([observed_counts_lookup[float(age)] for age in ages], dtype=int)
    final_eligibility = eligibility_by_bin(eligible, ages, observed_counts, ITERATIONS)

    convergence_bins = []
    convergence_phases = []
    convergence_metrics = []
    prefix_bin_tables = {}
    for prefix in PREFIXES:
        bins, phase, metrics = prefix_metrics(ev, contrasts, peaks, ages, observed, prefix)
        bins.insert(0, "simulations", prefix)
        phase.insert(0, "simulations", prefix)
        convergence_bins.append(bins)
        convergence_phases.append(phase)
        convergence_metrics.append(metrics)
        prefix_bin_tables[prefix] = bins.set_index("midpoint_Ma")

    convergence_bins = pd.concat(convergence_bins, ignore_index=True)
    convergence_phases = pd.concat(convergence_phases, ignore_index=True)
    convergence_metrics = pd.DataFrame(convergence_metrics)
    final_index = prefix_bin_tables[ITERATIONS]
    difference_rows = []
    for prefix in PREFIXES[:-1]:
        current = prefix_bin_tables[prefix]
        for metric in ["simulation_median", "simulation_2.5th", "simulation_97.5th", "interval_width"]:
            delta = np.abs(current[metric].to_numpy(float) - final_index[metric].to_numpy(float))
            difference_rows.append({
                "simulations": prefix,
                "comparison_target": ITERATIONS,
                "metric": metric,
                "median_absolute_difference": float(np.median(delta)),
                "maximum_absolute_difference": float(np.max(delta)),
                "maximum_difference_midpoint_Ma": float(current.index[int(np.argmax(delta))]),
            })
    convergence_differences = pd.DataFrame(difference_rows)

    peak_convergence = pd.concat([peak_frequency(peaks, prefix) for prefix in PREFIXES], ignore_index=True)
    eligibility_convergence = pd.concat([
        eligibility_by_bin(eligible, ages, observed_counts, prefix) for prefix in PREFIXES
    ], ignore_index=True)
    eligibility_distribution = pd.DataFrame([
        total_eligibility_summary(eligible, observed_counts, prefix) for prefix in PREFIXES
    ])

    envelope_comparison, contrast_comparison = comparisons_to_single_components(final_bins, final_contrasts)

    formal_prefix = "Joint_BB_formal"
    outputs = {
        f"{formal_prefix}_50kyr_Global_EV_summary.csv": final_bins,
        f"{formal_prefix}_phase_contrasts.csv": final_contrasts,
        f"{formal_prefix}_peak_midpoint_frequency.csv": final_peaks,
        f"{formal_prefix}_eligibility_by_bin.csv": final_eligibility,
        f"{formal_prefix}_eligibility_distribution.csv": eligibility_distribution,
        f"{formal_prefix}_convergence_bin_results.csv": convergence_bins,
        f"{formal_prefix}_convergence_phase_contrasts.csv": convergence_phases,
        f"{formal_prefix}_convergence_metrics.csv": convergence_metrics,
        f"{formal_prefix}_convergence_differences_vs_10000.csv": convergence_differences,
        f"{formal_prefix}_peak_convergence.csv": peak_convergence,
        f"{formal_prefix}_eligibility_convergence_by_bin.csv": eligibility_convergence,
        f"{formal_prefix}_pilot_reproduction_audit.csv": reproduction,
        f"{formal_prefix}_single_component_envelope_comparison.csv": envelope_comparison,
        f"{formal_prefix}_single_component_phase_contrast_comparison.csv": contrast_comparison,
    }
    for filename, table in outputs.items():
        table.to_csv(RESULT_DIR / filename, index=False)

    # Rename the locked engine's complete CSV matrices to unambiguous formal names.
    rename_map = {
        "joint_Global_EV_realizations_50kyr_1000.csv": f"{formal_prefix}_Global_EV_realizations_50kyr_10000.csv",
        "joint_phase_contrast_realizations_1000.csv": f"{formal_prefix}_phase_contrast_realizations_10000.csv",
        "joint_eligible_record_bin_counts_1000.csv": f"{formal_prefix}_eligible_record_bin_counts_10000.csv",
        "joint_chronology_shift_draws_1000.csv": f"{formal_prefix}_chronology_shift_draws_10000.csv",
    }
    for source_name, target_name in rename_map.items():
        source = RESULT_DIR / source_name
        target = RESULT_DIR / target_name
        if source.exists():
            source.replace(target)

    old_metadata_path = RESULT_DIR / "joint_pilot_metadata.json"
    engine_metadata = json.loads(old_metadata_path.read_text(encoding="utf-8"))
    widths = final_bins.interval_width.to_numpy(float)
    abs_difference = np.abs(final_bins.observed_minus_simulation_median.to_numpy(float))
    phase_lookup = final_contrasts.set_index("contrast")
    formal_metadata = {
        **engine_metadata,
        "analysis": "Formal joint conservative measurement and chronological sensitivity - standard 50-kyr Global EV",
        "iterations": ITERATIONS,
        "master_random_seed": MASTER_SEED,
        "formal_first_1000_exactly_reproduce_validated_pilot": bool(reproduction.exact_pointwise_match.all()),
        "formal_first_1000_maximum_absolute_difference": float(reproduction.maximum_absolute_difference.max()),
        "interval_width_median": float(np.median(widths)),
        "interval_width_maximum": float(np.max(widths)),
        "maximum_width_midpoint_Ma": float(final_bins.iloc[int(np.argmax(widths))].midpoint_Ma),
        "observed_vs_simulation_median_max_absolute_difference": float(np.max(abs_difference)),
        "maximum_difference_midpoint_Ma": float(final_bins.iloc[int(np.argmax(abs_difference))].midpoint_Ma),
        "phase_2_interval_entirely_above_zero": bool(phase_lookup.loc["Phase 2 - pre-escalation", "simulation_2.5th"] > 0),
        "recovery_interval_entirely_below_zero": bool(phase_lookup.loc["Recovery - Phase 2", "simulation_97.5th"] < 0),
        "phase_1_interval_entirely_above_zero": bool(phase_lookup.loc["Phase 1 - pre-escalation", "simulation_2.5th"] > 0),
        "convergence_prefixes": list(PREFIXES),
        "formal_output_npz": str(formal_npz_path),
    }
    formal_metadata_path = RESULT_DIR / f"{formal_prefix}_metadata.json"
    formal_metadata_path.write_text(json.dumps(formal_metadata, indent=2, ensure_ascii=False), encoding="utf-8")
    old_metadata_path.unlink(missing_ok=True)

    source_manifest_path = RESULT_DIR / "input_source_hash_manifest.csv"
    source_manifest = pd.read_csv(source_manifest_path)
    source_manifest = pd.concat([
        source_manifest,
        pd.DataFrame([{
            "file": SCRIPT.name,
            "path": str(SCRIPT),
            "sha256": sha256(SCRIPT),
        }, {
            "file": PILOT_NPZ.name,
            "path": str(PILOT_NPZ),
            "sha256": sha256(PILOT_NPZ),
        }, {
            "file": PILOT_WORKBOOK.name,
            "path": str(PILOT_WORKBOOK),
            "sha256": sha256(PILOT_WORKBOOK),
        }])
    ], ignore_index=True).drop_duplicates(["path", "sha256"])
    source_manifest.to_csv(source_manifest_path, index=False)

    # Remove pilot-label-only comparison tables generated by the locked engine;
    # their formal equivalents are written above.
    for obsolete in [
        "single_component_envelope_comparison.csv",
        "single_component_phase_contrast_comparison.csv",
    ]:
        (RESULT_DIR / obsolete).unlink(missing_ok=True)

    print(json.dumps(formal_metadata, indent=2, ensure_ascii=False))
    print("\nFORMAL PHASE CONTRASTS")
    print(final_contrasts.to_string(index=False))
    print("\nCONVERGENCE METRICS")
    print(convergence_metrics.to_string(index=False))
    print("\nPILOT REPRODUCTION")
    print(reproduction.to_string(index=False))


if __name__ == "__main__":
    main()
