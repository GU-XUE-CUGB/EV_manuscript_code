#!/usr/bin/env python3
"""Formal 10,000-iteration Chronology Layer B sensitivity analysis.

This executable preserves the validated 1,000-realization pilot design:
standard 50-kyr Global EV only, chronology-only perturbation, six independent
chronology-group random streams, dynamically recomputed bin membership/IQR/
eligibility, fixed measurement values, and no common absolute-timescale
systematic.  The first 1,000 realizations must reproduce the validated pilot
bit-for-bit before formal outputs are accepted.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd


MASTER_SEED = 20260824
ITERATIONS = 10000
PREFIXES = (1000, 5000, 10000)
TOTAL_RETAINED_MEASUREMENTS = 10881
TOTAL_OPERATIONAL_RECORDS = 70

SCRIPT = Path(__file__).resolve()
WORK_ROOT = SCRIPT.parent.parent
DATA_DIR = WORK_ROOT / "01_data_revised"
FORMAL_DIR = DATA_DIR / "formal_datasets"
INVENTORY_DIR = WORK_ROOT / "03_results" / "EV_uncertainty_inventory"
PILOT_DIR = INVENTORY_DIR / "Chronology_LayerB_conservative_sensitivity_pilot_50kyr_1000"
RESULT_DIR = INVENTORY_DIR / "Chronology_LayerB_conservative_sensitivity_formal_50kyr_10000"
RESULT_DIR.mkdir(parents=True, exist_ok=True)

PILOT_SCRIPT = WORK_ROOT / "02_code_revised" / "run_chronology_layerB_conservative_sensitivity_pilot_50kyr_1000.py"
PILOT_NPZ = PILOT_DIR / "Chronology_LayerB_pilot_complete_realizations_50kyr_1000.npz"
PRIMARY_FORMAL_NPZ = INVENTORY_DIR / "Chronology_primary_sensitivity_formal_50kyr_10000" / "Chronology_primary_formal_Global_EV_realizations_50kyr_10000.npz"
RAW_PATH = DATA_DIR / "retained_measurements_with_zscores_v4.csv"
IQR_PATH = DATA_DIR / "record_level_IQR_all_bin_widths_revised_v4_n2.csv"
S1_PATH = FORMAL_DIR / "Dataset_S1_revised_v6.xlsx"
PARAMETER_PATH = INVENTORY_DIR / "Chronology_LayerB_conservative_method_based_sensitivity" / "Chronology_LayerB_conservative_method_based_parameter_table.xlsx"

GROUP_ORDER = [
    "CHR-P1-ODP1209-SHATSKY",
    "CHR-P2-U1403",
    "CHR-P3-SEYMOUR",
    "CHR-B1-ODP1262",
    "CHR-B2-ODP1267",
    "CHR-B3-ODP690-POSTKPG",
]
SEGMENT_GROUP = "CHR-B3-ODP690-POSTKPG"
PRIMARY_SHIFT_KEYS = {
    "CHR-P1-ODP1209-SHATSKY": "ODP1209_shift_kyr",
    "CHR-P2-U1403": "U1403_shift_kyr",
    "CHR-P3-SEYMOUR": "Seymour_shift_kyr",
}
CONTRAST_NAMES = [
    "Phase 1 - pre-escalation",
    "Phase 2 - pre-escalation",
    "Recovery - Phase 2",
]

PROTECTED_PATHS = [
    PARAMETER_PATH,
    PILOT_NPZ,
    PRIMARY_FORMAL_NPZ,
    RAW_PATH,
    IQR_PATH,
    S1_PATH,
    FORMAL_DIR / "Dataset_S2_revised.xlsx",
    FORMAL_DIR / "Dataset_S3_revised_v3.xlsx",
    FORMAL_DIR / "Dataset_S4_revised_v3.xlsx",
    FORMAL_DIR / "Dataset_S5_revised_v4.xlsx",
    FORMAL_DIR / "Dataset_S6_revised_v4.xlsx",
    WORK_ROOT / "05_manuscript_revised" / "Manuscript_working_stage3_sampling_density.docx",
    WORK_ROOT / "06_supplement_revised" / "Supplementary_materials_working_stage3_sampling_density.docx",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_pilot_core():
    spec = importlib.util.spec_from_file_location("chronology_layerB_pilot_core", PILOT_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError("Cannot import validated Chronology Layer B pilot")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    if module.MASTER_SEED != MASTER_SEED or module.ITERATIONS != 1000:
        raise RuntimeError("Validated pilot seed or iteration count changed")
    module.ITERATIONS = ITERATIONS
    return module


def summarize_bins(ev: np.ndarray, ages: np.ndarray, observed: np.ndarray, prefix: int) -> pd.DataFrame:
    quantiles = np.quantile(ev[:prefix], [0.025, 0.5, 0.975], axis=0)
    return pd.DataFrame({
        "iterations": prefix,
        "midpoint_Ma": ages,
        "observed_EV": observed,
        "simulation_median": quantiles[1],
        "simulation_2.5th": quantiles[0],
        "simulation_97.5th": quantiles[2],
        "envelope_width": quantiles[2] - quantiles[0],
        "observed_minus_simulation_median": observed - quantiles[1],
        "valid_realizations": np.sum(np.isfinite(ev[:prefix]), axis=0),
    })


def summarize_contrasts(contrasts: np.ndarray, prefix: int) -> pd.DataFrame:
    quantiles = np.quantile(contrasts[:prefix], [0.025, 0.5, 0.975], axis=0)
    return pd.DataFrame([
        {
            "iterations": prefix,
            "contrast": name,
            "simulation_median": float(quantiles[1, index]),
            "simulation_2.5th": float(quantiles[0, index]),
            "simulation_97.5th": float(quantiles[2, index]),
            "interval_width": float(quantiles[2, index] - quantiles[0, index]),
            "valid_realizations": int(np.isfinite(contrasts[:prefix, index]).sum()),
        }
        for index, name in enumerate(CONTRAST_NAMES)
    ])


def peak_frequency(peaks: np.ndarray, prefix: int, ages: np.ndarray) -> pd.DataFrame:
    counts = Counter(peaks[:prefix])
    return pd.DataFrame([
        {
            "iterations": prefix,
            "peak_midpoint_Ma": float(age),
            "frequency": int(counts.get(float(age), 0)),
            "proportion": float(counts.get(float(age), 0) / prefix),
        }
        for age in ages
        if counts.get(float(age), 0) > 0
    ])


def convergence_tables(ev, contrasts, peaks, ages, observed):
    bin_prefixes = pd.concat(
        [summarize_bins(ev, ages, observed, prefix) for prefix in PREFIXES],
        ignore_index=True,
    )
    contrast_prefixes = pd.concat(
        [summarize_contrasts(contrasts, prefix) for prefix in PREFIXES],
        ignore_index=True,
    )
    peak_prefixes = pd.concat(
        [peak_frequency(peaks, prefix, ages) for prefix in PREFIXES],
        ignore_index=True,
    )

    final_bins = bin_prefixes.loc[bin_prefixes.iterations.eq(10000)].set_index("midpoint_Ma")
    metric_rows = []
    for prefix in PREFIXES:
        current = bin_prefixes.loc[bin_prefixes.iterations.eq(prefix)].set_index("midpoint_Ma")
        for statistic in ["simulation_median", "simulation_2.5th", "simulation_97.5th", "envelope_width"]:
            difference = (current[statistic] - final_bins[statistic]).abs()
            metric_rows.append({
                "iterations": prefix,
                "statistic": statistic,
                "median_absolute_difference_vs_10000": float(difference.median()),
                "maximum_absolute_difference_vs_10000": float(difference.max()),
                "midpoint_of_maximum_difference_Ma": float(difference.idxmax()),
            })
    convergence_metrics = pd.DataFrame(metric_rows)

    final_contrasts = contrast_prefixes.loc[contrast_prefixes.iterations.eq(10000)].set_index("contrast")
    contrast_rows = []
    for prefix in PREFIXES:
        current = contrast_prefixes.loc[contrast_prefixes.iterations.eq(prefix)].set_index("contrast")
        for contrast in CONTRAST_NAMES:
            for statistic in ["simulation_median", "simulation_2.5th", "simulation_97.5th", "interval_width"]:
                contrast_rows.append({
                    "iterations": prefix,
                    "contrast": contrast,
                    "statistic": statistic,
                    "absolute_difference_vs_10000": float(abs(current.loc[contrast, statistic] - final_contrasts.loc[contrast, statistic])),
                })
    contrast_convergence = pd.DataFrame(contrast_rows)

    final_counts = Counter(peaks)
    final_probability = np.array([final_counts.get(float(age), 0) / ITERATIONS for age in ages])
    peak_rows = []
    for prefix in PREFIXES:
        counts = Counter(peaks[:prefix])
        probability = np.array([counts.get(float(age), 0) / prefix for age in ages])
        modal_peak = max(counts, key=counts.get)
        peak_rows.append({
            "iterations": prefix,
            "modal_peak_midpoint_Ma": float(modal_peak),
            "modal_peak_frequency": int(counts[modal_peak]),
            "modal_peak_proportion": float(counts[modal_peak] / prefix),
            "total_variation_distance_vs_10000": float(0.5 * np.sum(np.abs(probability - final_probability))),
            "maximum_absolute_probability_difference_vs_10000": float(np.max(np.abs(probability - final_probability))),
        })
    peak_convergence = pd.DataFrame(peak_rows)
    return bin_prefixes, convergence_metrics, contrast_prefixes, contrast_convergence, peak_prefixes, peak_convergence


def primary_comparison(layer_ev, layer_contrasts, layer_peaks, primary_npz, ages, observed):
    primary_ev = np.asarray(primary_npz["primary_EV_50kyr"], dtype=float)
    primary_contrasts = np.asarray(primary_npz["primary_phase_contrasts"], dtype=float)
    primary_peaks = np.asarray(primary_npz["primary_peak_midpoint_Ma"], dtype=float)
    if primary_ev.shape != layer_ev.shape or primary_contrasts.shape != layer_contrasts.shape:
        raise RuntimeError("Primary and Layer B formal realization dimensions differ")
    if not np.array_equal(primary_npz["age_50kyr"], ages):
        raise RuntimeError("Primary formal ages differ from Layer B ages")
    if np.max(np.abs(primary_npz["observed_EV_50kyr"] - observed)) > 1e-12:
        raise RuntimeError("Primary formal and Layer B observed baselines differ")

    scenario_arrays = {
        "Primary formal 10,000": primary_ev,
        "Layer B formal 10,000": layer_ev,
    }
    envelope_rows, widths = [], {}
    for label, values in scenario_arrays.items():
        quantiles = np.quantile(values, [0.025, 0.975], axis=0)
        widths[label] = quantiles[1] - quantiles[0]
        envelope_rows.append({
            "scenario": label,
            "median_envelope_width": float(np.median(widths[label])),
            "maximum_envelope_width": float(np.max(widths[label])),
            "maximum_width_midpoint_Ma": float(ages[int(np.argmax(widths[label]))]),
        })
    envelope_comparison = pd.DataFrame(envelope_rows)
    width_change = widths["Layer B formal 10,000"] - widths["Primary formal 10,000"]
    envelope_by_bin = pd.DataFrame({
        "midpoint_Ma": ages,
        "primary_formal_10000_width": widths["Primary formal 10,000"],
        "LayerB_formal_10000_width": widths["Layer B formal 10,000"],
        "LayerB_minus_primary_width": width_change,
        "absolute_width_change": np.abs(width_change),
    })

    primary_q = np.quantile(primary_contrasts, [0.025, 0.5, 0.975], axis=0)
    layer_q = np.quantile(layer_contrasts, [0.025, 0.5, 0.975], axis=0)
    difference_q = np.quantile(layer_contrasts - primary_contrasts, [0.025, 0.5, 0.975], axis=0)
    contrast_comparison = pd.DataFrame([
        {
            "contrast": name,
            "primary_median": float(primary_q[1, index]),
            "primary_2.5th": float(primary_q[0, index]),
            "primary_97.5th": float(primary_q[2, index]),
            "LayerB_median": float(layer_q[1, index]),
            "LayerB_2.5th": float(layer_q[0, index]),
            "LayerB_97.5th": float(layer_q[2, index]),
            "paired_LayerB_minus_primary_median": float(difference_q[1, index]),
            "paired_difference_2.5th": float(difference_q[0, index]),
            "paired_difference_97.5th": float(difference_q[2, index]),
        }
        for index, name in enumerate(CONTRAST_NAMES)
    ])

    primary_counts, layer_counts = Counter(primary_peaks), Counter(layer_peaks)
    peak_comparison = pd.DataFrame([
        {
            "peak_midpoint_Ma": float(age),
            "primary_frequency": int(primary_counts.get(float(age), 0)),
            "primary_proportion": float(primary_counts.get(float(age), 0) / ITERATIONS),
            "LayerB_frequency": int(layer_counts.get(float(age), 0)),
            "LayerB_proportion": float(layer_counts.get(float(age), 0) / ITERATIONS),
            "LayerB_minus_primary_percentage_points": float(100 * (layer_counts.get(float(age), 0) - primary_counts.get(float(age), 0)) / ITERATIONS),
        }
        for age in ages
        if primary_counts.get(float(age), 0) or layer_counts.get(float(age), 0)
    ])
    metrics = {
        "LayerB_minus_primary_median_envelope_width": float(np.median(widths["Layer B formal 10,000"]) - np.median(widths["Primary formal 10,000"])),
        "LayerB_minus_primary_maximum_envelope_width": float(np.max(widths["Layer B formal 10,000"]) - np.max(widths["Primary formal 10,000"])),
        "median_bin_specific_envelope_width_change": float(np.median(width_change)),
        "maximum_absolute_bin_specific_width_change": float(np.max(np.abs(width_change))),
        "maximum_absolute_width_change_midpoint_Ma": float(ages[int(np.argmax(np.abs(width_change)))]),
        "primary_modal_peak_Ma": float(pd.Series(primary_peaks).mode().iloc[0]),
        "LayerB_modal_peak_Ma": float(pd.Series(layer_peaks).mode().iloc[0]),
        "paired_realizations_with_changed_peak_fraction": float(np.mean(primary_peaks != layer_peaks)),
        "peak_total_variation_distance": float(0.5 * sum(abs(layer_counts.get(float(age), 0) - primary_counts.get(float(age), 0)) / ITERATIONS for age in ages)),
    }
    return envelope_comparison, envelope_by_bin, contrast_comparison, peak_comparison, metrics


def main() -> None:
    core = load_pilot_core()
    protected_before = {str(path): sha256(path) for path in PROTECTED_PATHS if path.exists()}

    raw = pd.read_csv(RAW_PATH).reset_index(drop=True)
    raw["record_id"] = raw.record_id.map(core.clean_id)
    if len(raw) != TOTAL_RETAINED_MEASUREMENTS or raw.record_id.nunique() != TOTAL_OPERATIONAL_RECORDS:
        raise RuntimeError("Final retained baseline must contain 10,881 measurements / 70 records")
    parameter_records, parameter_groups, group_members = core.load_parameter_table(raw)
    record_arrays = core.load_module(core.PRIMARY_PILOT_SCRIPT, "primary_chronology_pilot_formal_core").prepare_record_arrays(raw)
    chronology_core = sys.modules["primary_chronology_pilot_formal_core"]
    ages, observed, observed_counts, s1_difference, observed_zero_bins = chronology_core.validate_observed_baseline(raw, record_arrays)

    active_records = set(parameter_records.loc[parameter_records.active_flag.eq(1), "record_id"])
    fixed_records = [record_id for record_id in record_arrays if record_id not in active_records]
    if len(active_records) != 38 or len(fixed_records) != 32:
        raise RuntimeError("Expected 38 active and 32 fixed operational records")
    fixed_sums, fixed_counts, fixed_diagnostics = chronology_core.record_contributions(record_arrays, fixed_records, 0.0)

    baseline_groups = {}
    for group_id in GROUP_ORDER:
        if group_id == SEGMENT_GROUP:
            baseline_groups[group_id] = core.segment_record_contributions(chronology_core, record_arrays, group_members[group_id], 0.0)
        else:
            baseline_groups[group_id] = chronology_core.record_contributions(record_arrays, group_members[group_id], 0.0)
    baseline_sums = fixed_sums + sum(baseline_groups[group_id][0] for group_id in GROUP_ORDER)
    baseline_counts = fixed_counts + sum(baseline_groups[group_id][1] for group_id in GROUP_ORDER)
    if not np.array_equal(baseline_counts, observed_counts):
        raise RuntimeError("Formal groups do not reproduce observed eligibility")
    if np.max(np.abs(baseline_sums / baseline_counts - observed)) > 1e-12:
        raise RuntimeError("Formal groups do not reproduce observed Global EV")

    primary_npz = np.load(PRIMARY_FORMAL_NPZ)
    shifts, seed_info = core.draw_group_shifts(parameter_groups, primary_npz)

    group_sums, group_counts, cache_diagnostics = {}, {}, {}
    for group_id in GROUP_ORDER:
        sums, counts, diagnostics = core.cached_group_contributions(
            chronology_core,
            record_arrays,
            group_members[group_id],
            shifts[group_id],
            group_id,
            segment_only=group_id == SEGMENT_GROUP,
        )
        group_sums[group_id], group_counts[group_id], cache_diagnostics[group_id] = sums, counts, diagnostics

    total_sums = np.broadcast_to(fixed_sums, (ITERATIONS, chronology_core.BIN_SPEC["count"])).copy()
    total_counts = np.broadcast_to(fixed_counts, (ITERATIONS, chronology_core.BIN_SPEC["count"])).copy()
    for group_id in GROUP_ORDER:
        total_sums += group_sums[group_id]
        total_counts += group_counts[group_id]
    ev = np.divide(total_sums, total_counts, out=np.full_like(total_sums, np.nan), where=total_counts > 0)
    contrasts = chronology_core.phase_contrasts(ev, ages)
    peaks = core.peak_midpoints(ev, ages)

    if not PILOT_NPZ.exists():
        raise RuntimeError("Validated Layer B pilot realization archive is missing")
    pilot_npz = np.load(PILOT_NPZ)
    prefix_checks = {
        "EV_max_absolute_difference": float(np.max(np.abs(ev[:1000] - pilot_npz["LayerB_EV_50kyr"]))),
        "eligible_count_max_absolute_difference": int(np.max(np.abs(total_counts[:1000] - pilot_npz["LayerB_eligible_record_counts"]))),
        "phase_contrast_max_absolute_difference": float(np.max(np.abs(contrasts[:1000] - pilot_npz["LayerB_phase_contrasts"]))),
        "peak_midpoint_max_absolute_difference": float(np.max(np.abs(peaks[:1000] - pilot_npz["LayerB_peak_midpoint_Ma"]))),
    }
    for group_id in GROUP_ORDER:
        prefix_checks[f"{group_id}_shift_max_absolute_difference"] = float(
            np.max(np.abs(shifts[group_id][:1000] - pilot_npz[f"shift_{group_id}"]))
        )
    if any(value != 0 for value in prefix_checks.values()):
        raise RuntimeError(f"Formal first 1,000 realizations do not exactly reproduce pilot: {prefix_checks}")

    direct_recalculation_checks = {}
    for iteration in [0, 999, 4999, 9999]:
        direct_sums, direct_counts = fixed_sums.copy(), fixed_counts.copy()
        for group_id in GROUP_ORDER:
            shift_ma = float(shifts[group_id][iteration] / 1000.0)
            if group_id == SEGMENT_GROUP:
                sums, counts, _ = core.segment_record_contributions(chronology_core, record_arrays, group_members[group_id], shift_ma)
            else:
                sums, counts, _ = chronology_core.record_contributions(record_arrays, group_members[group_id], shift_ma)
            direct_sums += sums
            direct_counts += counts
        direct_ev = np.divide(direct_sums, direct_counts, out=np.full_like(direct_sums, np.nan), where=direct_counts > 0)
        check = {
            "max_EV_difference": float(np.nanmax(np.abs(direct_ev - ev[iteration]))),
            "max_record_count_difference": int(np.max(np.abs(direct_counts - total_counts[iteration]))),
        }
        direct_recalculation_checks[str(iteration + 1)] = check
        if check["max_EV_difference"] > 1e-12 or check["max_record_count_difference"] != 0:
            raise RuntimeError(f"Cached contribution verification failed at iteration {iteration + 1}")

    bin_prefixes, convergence_metrics, contrast_prefixes, contrast_convergence, peak_prefixes, peak_convergence = convergence_tables(
        ev, contrasts, peaks, ages, observed
    )
    final_ev = bin_prefixes.loc[bin_prefixes.iterations.eq(10000)].copy()
    final_contrasts = contrast_prefixes.loc[contrast_prefixes.iterations.eq(10000)].copy()
    final_peaks = peak_prefixes.loc[peak_prefixes.iterations.eq(10000)].copy()

    eligibility = pd.DataFrame({
        "midpoint_Ma": ages,
        "observed_eligible_record_bins": observed_counts,
        "simulation_record_bins_min": np.min(total_counts, axis=0),
        "simulation_record_bins_median": np.median(total_counts, axis=0),
        "simulation_record_bins_max": np.max(total_counts, axis=0),
        "simulation_record_bins_2.5th": np.quantile(total_counts, 0.025, axis=0),
        "simulation_record_bins_97.5th": np.quantile(total_counts, 0.975, axis=0),
        "iterations_with_empty_bin": np.sum(total_counts == 0, axis=0),
    })
    total_eligible_record_bins = total_counts.sum(axis=1)
    eligibility_distribution = pd.DataFrame([{
        "metric": "Total eligible record-bins per realization",
        "minimum": int(total_eligible_record_bins.min()),
        "2.5th": float(np.quantile(total_eligible_record_bins, 0.025)),
        "median": float(np.median(total_eligible_record_bins)),
        "97.5th": float(np.quantile(total_eligible_record_bins, 0.975)),
        "maximum": int(total_eligible_record_bins.max()),
        "observed": int(observed_counts.sum()),
    }])

    envelope_comparison, envelope_by_bin, contrast_comparison, peak_comparison, comparison_metrics = primary_comparison(
        ev, contrasts, peaks, primary_npz, ages, observed
    )

    segment_diagnostics = cache_diagnostics[SEGMENT_GROUP]["diagnostic_totals"]
    invalid_iqr_events = int(fixed_diagnostics["invalid_iqr_record_bins"])
    for diagnostics in cache_diagnostics.values():
        invalid_iqr_events += int(diagnostics["diagnostic_totals"].get("invalid_iqr_record_bins", 0))
    execution_diagnostics = {
        "age_inversions": int(segment_diagnostics.get("age_inversions", 0)),
        "segment_boundary_order_violations": int(segment_diagnostics.get("segment_boundary_order_violations", 0)),
        "ODP690_post_KPg_measurements_shifted_per_realization": int(
            sum(np.sum(record_arrays[record_id][0] < core.KPG_BOUNDARY_MA) for record_id in group_members[SEGMENT_GROUP])
        ),
        "ODP690_shift_draw_min_kyr": float(shifts[SEGMENT_GROUP].min()),
        "ODP690_shift_draw_max_kyr": float(shifts[SEGMENT_GROUP].max()),
        "iterations_with_at_least_one_empty_bin": int(np.sum(np.any(total_counts == 0, axis=1))),
        "total_empty_bin_events": int(np.sum(total_counts == 0)),
        "invalid_or_missing_EV_trajectory_iterations": int(np.sum(~np.all(np.isfinite(ev), axis=1))),
        "nonfinite_EV_values": int(np.sum(~np.isfinite(ev))),
        "invalid_IQR_events": invalid_iqr_events,
        "observed_total_eligible_record_bins": int(observed_counts.sum()),
        "simulated_total_eligible_record_bins_min": int(total_eligible_record_bins.min()),
        "simulated_total_eligible_record_bins_median": float(np.median(total_eligible_record_bins)),
        "simulated_total_eligible_record_bins_max": int(total_eligible_record_bins.max()),
        "maximum_absolute_per_bin_record_count_change": int(np.max(np.abs(total_counts - observed_counts[None, :]))),
        "observed_zero_IQR_record_bins_excluded": int(observed_zero_bins),
    }
    if any(execution_diagnostics[key] != 0 for key in [
        "age_inversions", "segment_boundary_order_violations", "total_empty_bin_events",
        "invalid_or_missing_EV_trajectory_iterations", "nonfinite_EV_values", "invalid_IQR_events",
    ]):
        raise RuntimeError(f"Formal Layer B execution diagnostics failed: {execution_diagnostics}")

    final_ev.to_csv(RESULT_DIR / "Chronology_LayerB_formal_50kyr_Global_EV_summary.csv", index=False)
    final_contrasts.to_csv(RESULT_DIR / "Chronology_LayerB_formal_phase_contrasts.csv", index=False)
    final_peaks.to_csv(RESULT_DIR / "Chronology_LayerB_formal_peak_midpoint_frequency.csv", index=False)
    eligibility.to_csv(RESULT_DIR / "Chronology_LayerB_formal_eligibility_by_bin.csv", index=False)
    eligibility_distribution.to_csv(RESULT_DIR / "Chronology_LayerB_formal_eligibility_distribution.csv", index=False)
    bin_prefixes.to_csv(RESULT_DIR / "Chronology_LayerB_formal_convergence_bin_prefixes.csv", index=False)
    convergence_metrics.to_csv(RESULT_DIR / "Chronology_LayerB_formal_convergence_metrics.csv", index=False)
    contrast_prefixes.to_csv(RESULT_DIR / "Chronology_LayerB_formal_convergence_phase_contrasts.csv", index=False)
    contrast_convergence.to_csv(RESULT_DIR / "Chronology_LayerB_formal_convergence_phase_metrics.csv", index=False)
    peak_prefixes.to_csv(RESULT_DIR / "Chronology_LayerB_formal_convergence_peak_distributions.csv", index=False)
    peak_convergence.to_csv(RESULT_DIR / "Chronology_LayerB_formal_convergence_peak_metrics.csv", index=False)
    envelope_comparison.to_csv(RESULT_DIR / "Chronology_LayerB_formal_vs_primary_envelope_comparison.csv", index=False)
    envelope_by_bin.to_csv(RESULT_DIR / "Chronology_LayerB_formal_vs_primary_envelope_change_by_bin.csv", index=False)
    contrast_comparison.to_csv(RESULT_DIR / "Chronology_LayerB_formal_vs_primary_phase_contrast_comparison.csv", index=False)
    peak_comparison.to_csv(RESULT_DIR / "Chronology_LayerB_formal_vs_primary_peak_distribution_comparison.csv", index=False)

    shift_table = pd.DataFrame({"iteration": np.arange(1, ITERATIONS + 1)})
    for group_id in GROUP_ORDER:
        shift_table[f"{group_id}_shift_kyr"] = shifts[group_id]
    shift_table.to_csv(RESULT_DIR / "Chronology_LayerB_formal_group_shift_draws.csv", index=False)

    group_parameters = parameter_groups.copy()
    raw_counts = raw.groupby("record_id").size()
    group_parameters["actual_record_count"] = group_parameters.perturbation_group_id.map(lambda group_id: len(group_members[group_id]))
    group_parameters["actual_measurement_count"] = group_parameters.perturbation_group_id.map(lambda group_id: int(raw_counts.loc[group_members[group_id]].sum()))
    group_parameters["random_seed"] = MASTER_SEED
    group_parameters.to_csv(RESULT_DIR / "Chronology_LayerB_formal_group_parameters.csv", index=False)

    active_membership = parameter_records.loc[parameter_records.active_flag.eq(1), [
        "record_id", "category", "sheet", "n_measurements", "chronology_method_class", "actual_chronology_method",
        "scenario_role", "perturbation_group_id", "recommended_perturbation_structure", "bound_min_kyr",
        "bound_max_kyr", "distribution", "applicable_interval", "correlation_rule", "range_basis",
    ]].copy()
    active_membership.to_csv(RESULT_DIR / "Chronology_LayerB_formal_active_record_membership.csv", index=False)

    realizations_wide = pd.DataFrame(ev, columns=[f"EV_{age:.2f}_Ma" for age in ages])
    realizations_wide.insert(0, "total_eligible_record_bins", total_eligible_record_bins)
    realizations_wide.insert(0, "peak_midpoint_Ma", peaks)
    realizations_wide.insert(0, "iteration", np.arange(1, ITERATIONS + 1))
    realizations_wide.to_csv(RESULT_DIR / "Chronology_LayerB_formal_Global_EV_realizations_50kyr_10000.csv", index=False)

    phase_realizations = pd.DataFrame(contrasts, columns=CONTRAST_NAMES)
    phase_realizations.insert(0, "iteration", np.arange(1, ITERATIONS + 1))
    phase_realizations.to_csv(RESULT_DIR / "Chronology_LayerB_formal_phase_contrast_realizations_10000.csv", index=False)

    np.savez_compressed(
        RESULT_DIR / "Chronology_LayerB_formal_complete_realizations_50kyr_10000.npz",
        age_50kyr=ages,
        observed_EV_50kyr=observed,
        LayerB_EV_50kyr=ev,
        LayerB_eligible_record_counts=total_counts,
        LayerB_phase_contrasts=contrasts,
        LayerB_peak_midpoint_Ma=peaks,
        total_eligible_record_bins=total_eligible_record_bins,
        **{f"shift_{group_id}": shifts[group_id] for group_id in GROUP_ORDER},
    )

    widths = final_ev.envelope_width.to_numpy(float)
    observed_differences = np.abs(final_ev.observed_EV - final_ev.simulation_median).to_numpy(float)
    phase2 = final_contrasts.loc[final_contrasts.contrast.eq("Phase 2 - pre-escalation")].iloc[0]
    metadata = {
        "analysis": "Formal Chronology Layer B conservative method-based sensitivity — standard 50-kyr Global EV",
        "iterations": ITERATIONS,
        "master_random_seed": MASTER_SEED,
        "active_chronology_groups": 6,
        "active_operational_records": 38,
        "active_measurements": 8843,
        "active_measurement_coverage_fraction": 8843 / TOTAL_RETAINED_MEASUREMENTS,
        "active_measurement_coverage_percent": 100 * 8843 / TOTAL_RETAINED_MEASUREMENTS,
        "fixed_operational_records": 32,
        "measurement_values_perturbed": False,
        "common_absolute_timescale_systematic_included": False,
        "common_astronomical_target_systematic_included": False,
        "common_GPTS_systematic_included": False,
        "observed_record_bin_mask_reused": False,
        "dynamic_rebinning_and_eligibility": True,
        "KPg_boundary_Ma_for_ODP690_segment_definition": core.KPG_BOUNDARY_MA,
        "Dataset_S1_observed_EV_max_absolute_validation_difference": s1_difference,
        "formal_prefix_vs_validated_pilot": prefix_checks,
        "formal_prefix_exactly_reproduces_pilot": True,
        "cache_diagnostics": cache_diagnostics,
        "direct_recalculation_checks": direct_recalculation_checks,
        "observed_EV_vs_simulation_median_max_absolute_difference": float(np.max(observed_differences)),
        "maximum_difference_midpoint_Ma": float(final_ev.iloc[int(np.argmax(observed_differences))].midpoint_Ma),
        "envelope_median_width": float(np.median(widths)),
        "envelope_maximum_width": float(np.max(widths)),
        "maximum_width_midpoint_Ma": float(final_ev.iloc[int(np.argmax(widths))].midpoint_Ma),
        "observed_peak_midpoint_Ma": float(ages[int(np.argmax(observed))]),
        "Phase2_minus_pre_median": float(phase2.simulation_median),
        "Phase2_minus_pre_2.5th": float(phase2["simulation_2.5th"]),
        "Phase2_minus_pre_97.5th": float(phase2["simulation_97.5th"]),
        "Phase2_stable_positive_interval": bool(phase2["simulation_2.5th"] > 0),
        "execution_diagnostics": execution_diagnostics,
        "primary_comparison": comparison_metrics,
    }
    (RESULT_DIR / "Chronology_LayerB_formal_metadata.json").write_text(json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8")
    (RESULT_DIR / "Chronology_LayerB_formal_group_cache_diagnostics.json").write_text(json.dumps(cache_diagnostics, indent=2, ensure_ascii=False), encoding="utf-8")

    protected_after = {str(path): sha256(path) for path in PROTECTED_PATHS if path.exists()}
    protected_check = pd.DataFrame([
        {
            "path": path,
            "sha256_before": digest,
            "sha256_after": protected_after.get(path),
            "unchanged": digest == protected_after.get(path),
        }
        for path, digest in protected_before.items()
    ])
    protected_check.to_csv(RESULT_DIR / "protected_files_hash_check.csv", index=False)
    if not protected_check.unchanged.all():
        raise RuntimeError("A protected Dataset, source result, figure, or manuscript changed")

    manifest_paths = [RAW_PATH, IQR_PATH, S1_PATH, PARAMETER_PATH, PRIMARY_FORMAL_NPZ, PILOT_NPZ, PILOT_SCRIPT, SCRIPT]
    pd.DataFrame([
        {"file": path.name, "path": str(path), "sha256": sha256(path)} for path in manifest_paths
    ]).to_csv(RESULT_DIR / "input_source_hash_manifest.csv", index=False)

    print(json.dumps(metadata, indent=2, ensure_ascii=False), flush=True)
    print("\nFORMAL PHASE CONTRASTS", flush=True)
    print(final_contrasts.to_string(index=False), flush=True)
    print("\nFORMAL PEAK DISTRIBUTION", flush=True)
    print(final_peaks.to_string(index=False), flush=True)
    print("\nCONVERGENCE METRICS", flush=True)
    print(convergence_metrics.to_string(index=False), flush=True)
    print("\nPRIMARY COMPARISON", flush=True)
    print(contrast_comparison.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
