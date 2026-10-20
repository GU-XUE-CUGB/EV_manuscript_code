#!/usr/bin/env python3
"""Joint conservative measurement + chronological sensitivity pilot.

Locked scope
------------
* Standard 50-kyr Global EV only; 1,000 simulations.
* Measurement perturbations exactly reuse the validated Fig. S6B mappings and
  perturbation structures (11 documented records plus Tier 2 Gaussian and
  bounded observation-level records; non-executable/fixed records unchanged).
* Chronology perturbations exactly reuse the validated Fig. S7B six-group
  parameter table, including shared group shifts and the ODP 690 post-K-Pg
  segment-only perturbation.
* Measurement and chronology draws are generated from independent RNG streams
  derived from one master SeedSequence (entropy 20260824).
* Every simulation starts from retained observed measurements and ages,
  recalculates record-wide z-scores, dynamically re-bins perturbed ages,
  recalculates record-bin IQRs, reapplies n >= 2 and IQR > 0 eligibility, and
  computes record-weighted Global EV.
* No Dataset, manuscript, Supplement, or existing figure is modified.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import math
import sys
import unicodedata
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd


MASTER_SEED = 20260824
ITERATIONS = 1000
TOTAL_MEASUREMENTS = 10881
TOTAL_RECORDS = 70
KPG_BOUNDARY_MA = 66.016

SCRIPT = Path(__file__).resolve()
WORK_ROOT = SCRIPT.parent.parent
DATA_DIR = WORK_ROOT / "01_data_revised"
FORMAL_DIR = DATA_DIR / "formal_datasets"
INVENTORY_DIR = WORK_ROOT / "03_results" / "EV_uncertainty_inventory"
RESULT_DIR = INVENTORY_DIR / "Joint_conservative_measurement_chronology_sensitivity_pilot_50kyr_1000"
RESULT_DIR.mkdir(parents=True, exist_ok=True)

RAW_PATH = DATA_DIR / "retained_measurements_with_zscores_v4.csv"
IQR_PATH = DATA_DIR / "record_level_IQR_all_bin_widths_revised_v4_n2.csv"
S1_PATH = FORMAL_DIR / "Dataset_S1_revised_v6.xlsx"
MEASUREMENT_PARAMETER_PATH = INVENTORY_DIR / "Executable_uncertainty_parameter_table.xlsx"
CHRONOLOGY_PARAMETER_PATH = (
    INVENTORY_DIR
    / "Chronology_LayerB_conservative_method_based_sensitivity"
    / "Chronology_LayerB_conservative_method_based_parameter_table.xlsx"
)

MEASUREMENT_PILOT_SCRIPT = WORK_ROOT / "02_code_revised" / "run_layer2_conservative_measurement_sensitivity_pilot.py"
MEASUREMENT_FORMAL_SCRIPT = WORK_ROOT / "02_code_revised" / "run_layer2_conservative_measurement_sensitivity_formal_50kyr_10000.py"
CHRONOLOGY_PILOT_SCRIPT = WORK_ROOT / "02_code_revised" / "run_chronology_layerB_conservative_sensitivity_pilot_50kyr_1000.py"
PRIMARY_CHRONOLOGY_SCRIPT = WORK_ROOT / "02_code_revised" / "run_primary_chronology_sensitivity_pilot_50kyr_1000.py"

MEASUREMENT_FORMAL_DIR = INVENTORY_DIR / "Layer2_conservative_measurement_sensitivity_formal_50kyr_10000"
CHRONOLOGY_FORMAL_DIR = INVENTORY_DIR / "Chronology_LayerB_conservative_sensitivity_formal_50kyr_10000"

MEASUREMENT_FORMAL_EV = MEASUREMENT_FORMAL_DIR / "Layer2_formal_50kyr_Global_EV_summary.csv"
MEASUREMENT_FORMAL_CONTRAST = MEASUREMENT_FORMAL_DIR / "Layer2_formal_phase_contrasts_50kyr.csv"
CHRONOLOGY_FORMAL_EV = CHRONOLOGY_FORMAL_DIR / "Chronology_LayerB_formal_50kyr_Global_EV_summary.csv"
CHRONOLOGY_FORMAL_CONTRAST = CHRONOLOGY_FORMAL_DIR / "Chronology_LayerB_formal_phase_contrasts.csv"

PHASES = {
    "pre-escalation": (66.95, 67.20),
    "Phase 1": (66.45, 66.90),
    "Phase 2": (66.05, 66.40),
    "Recovery": (64.80, 65.95),
}
CONTRAST_NAMES = [
    "Phase 1 - pre-escalation",
    "Phase 2 - pre-escalation",
    "Recovery - Phase 2",
]

PROTECTED_PATHS = [
    RAW_PATH,
    IQR_PATH,
    S1_PATH,
    MEASUREMENT_PARAMETER_PATH,
    CHRONOLOGY_PARAMETER_PATH,
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


def clean_id(value) -> str:
    return unicodedata.normalize("NFC", str(value)).strip()


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def measurement_execution_audit(raw, params, layer1_maps, layer1_notes, tier2_audit, formal_module):
    audit = formal_module.build_execution_audit(raw, params, layer1_maps, layer1_notes, tier2_audit)
    if int(audit.actually_active.sum()) != 53:
        raise RuntimeError("Fig. S6B execution scope must contain 53 active records")
    if int(audit.actually_perturbed_measurements.sum()) != 10015:
        raise RuntimeError("Fig. S6B execution scope must contain 10,015 perturbed measurements")

    documented = audit.primary_structure.str.startswith("Tier 1") & audit.actually_active
    gaussian = audit.primary_structure.str.contains("Conservative Gaussian") & audit.actually_active
    bounded = audit.primary_structure.str.contains("Bounded uniform") & audit.actually_active
    checks = {
        "documented": (int(documented.sum()), int(audit.loc[documented, "actually_perturbed_measurements"].sum())),
        "tier2_gaussian": (int(gaussian.sum()), int(audit.loc[gaussian, "actually_perturbed_measurements"].sum())),
        "tier2_bounded": (int(bounded.sum()), int(audit.loc[bounded, "actually_perturbed_measurements"].sum())),
    }
    expected = {
        "documented": (11, 375),
        "tier2_gaussian": (19, 4414),
        "tier2_bounded": (23, 5226),
    }
    if checks != expected:
        raise RuntimeError(f"Fig. S6B structure counts changed: {checks}")
    return audit, formal_module.structure_participation(audit), checks


def prepare_measurement_draw_arrays(mappings):
    normal = [item for item in mappings if item.distribution == "normal_direct"]
    uniform = [item for item in mappings if item.distribution == "uniform_direct"]
    isotope = [item for item in mappings if item.distribution == "odp1049_isotope_normal"]
    return {
        "normal_ix": np.array([item.row_index for item in normal], dtype=int),
        "normal_sigma": np.array([item.sigma for item in normal], dtype=float),
        "uniform_ix": np.array([item.row_index for item in uniform], dtype=int),
        "uniform_lo": np.array([item.lower for item in uniform], dtype=float),
        "uniform_hi": np.array([item.upper for item in uniform], dtype=float),
        "isotope_ix": np.array([item.row_index for item in isotope], dtype=int),
        "isotope_sigma": np.array([item.sigma for item in isotope], dtype=float),
        "isotope_base": np.array([item.isotope_value for item in isotope], dtype=float),
        "isotope_species": np.array([item.species for item in isotope]),
    }


def perturb_measurements(base, draw_arrays, rng):
    values = base.copy()
    normal_ix = draw_arrays["normal_ix"]
    uniform_ix = draw_arrays["uniform_ix"]
    isotope_ix = draw_arrays["isotope_ix"]
    if len(normal_ix):
        values[normal_ix] += rng.normal(0.0, draw_arrays["normal_sigma"])
    if len(uniform_ix):
        values[uniform_ix] += rng.uniform(draw_arrays["uniform_lo"], draw_arrays["uniform_hi"])
    if len(isotope_ix):
        isotope_values = draw_arrays["isotope_base"] + rng.normal(0.0, draw_arrays["isotope_sigma"])
        species_s = draw_arrays["isotope_species"] == "S"
        values[isotope_ix[species_s]] = (
            0.03 * isotope_values[species_s] ** 2
            - 4.43 * isotope_values[species_s]
            + 10.2855
        )
        values[isotope_ix[~species_s]] = (
            0.03 * isotope_values[~species_s] ** 2
            - 4.427 * isotope_values[~species_s]
            + 10.064075
        )
    return values


def draw_chronology_shifts(groups, group_order, chronology_seed_sequence):
    children = chronology_seed_sequence.spawn(len(group_order))
    shifts = {}
    metadata = []
    for group_id, child in zip(group_order, children):
        row = groups.loc[groups.perturbation_group_id.eq(group_id)].iloc[0]
        rng = np.random.default_rng(child)
        draws = rng.uniform(float(row.bound_min_kyr), float(row.bound_max_kyr), size=ITERATIONS)
        shifts[group_id] = draws
        metadata.append({
            "component": "chronology",
            "stream": group_id,
            "master_entropy": int(child.entropy),
            "spawn_key": ".".join(str(value) for value in child.spawn_key),
            "bit_generator": rng.bit_generator.__class__.__name__,
            "distribution": "bounded uniform",
            "bound_min_kyr": float(row.bound_min_kyr),
            "bound_max_kyr": float(row.bound_max_kyr),
            "draw_min_kyr": float(draws.min()),
            "draw_max_kyr": float(draws.max()),
        })
    return shifts, pd.DataFrame(metadata)


def chronology_row_masks(raw, group_members, group_order, segment_group):
    record_values = raw.record_id.to_numpy(str)
    masks = {}
    for group_id in group_order:
        mask = np.isin(record_values, group_members[group_id])
        if group_id == segment_group:
            mask &= raw.age_Ma.to_numpy(float) < KPG_BOUNDARY_MA
        masks[group_id] = mask
    return masks


def perturb_ages(base_ages, shifts, iteration, group_order, row_masks):
    ages = base_ages.copy()
    for group_id in group_order:
        ages[row_masks[group_id]] += shifts[group_id][iteration] / 1000.0
    return ages


def compute_recordwise_z(values, record_indices):
    z = np.empty(len(values), dtype=float)
    invalid_sd = 0
    for indices in record_indices:
        part = values[indices]
        sd = float(part.std(ddof=1))
        if not math.isfinite(sd) or sd <= 0:
            z[indices] = np.nan
            invalid_sd += 1
        else:
            z[indices] = (part - float(part.mean())) / sd
    return z, invalid_sd


def dynamic_ev(primary_module, z, ages, record_indices):
    assigned = primary_module.assign_bins(ages)
    n_bins = primary_module.BIN_SPEC["count"]
    sums = np.zeros(n_bins, dtype=float)
    counts = np.zeros(n_bins, dtype=np.int16)
    diagnostics = Counter()
    for indices in record_indices:
        record_bins = assigned[indices]
        record_z = z[indices]
        diagnostics["measurements_assigned_to_domain"] += int(np.sum(record_bins >= 0))
        for bin_index in np.unique(record_bins[record_bins >= 0]):
            values = record_z[record_bins == bin_index]
            if len(values) < 2:
                continue
            diagnostics["record_bins_with_n_ge_2"] += 1
            iqr = primary_module.iqr_linear(values)
            if not math.isfinite(iqr) or iqr < 0:
                diagnostics["invalid_iqr_record_bins"] += 1
            elif iqr == 0:
                diagnostics["zero_iqr_record_bins_excluded"] += 1
            else:
                sums[int(bin_index)] += iqr
                counts[int(bin_index)] += 1
                diagnostics["eligible_nonzero_record_bins"] += 1
    ev = np.divide(sums, counts, out=np.full(n_bins, np.nan), where=counts > 0)
    return ev, counts, diagnostics


def phase_contrasts(ev, ages):
    phase_means = {}
    for name, (lo, hi) in PHASES.items():
        selected = (ages >= lo - 1e-10) & (ages <= hi + 1e-10)
        phase_means[name] = float(np.mean(ev[selected]))
    return np.array([
        phase_means["Phase 1"] - phase_means["pre-escalation"],
        phase_means["Phase 2"] - phase_means["pre-escalation"],
        phase_means["Recovery"] - phase_means["Phase 2"],
    ])


def observed_coverage(raw, observed_mask, measurement_active_rows, chronology_active_records):
    chronology_active_rows = set(raw.index[raw.record_id.isin(chronology_active_records)].astype(int))
    rows = []
    per_bin = {index: set() for index in range(len(observed_mask["centers"]))}
    for bin_index, indices, _record_id in observed_mask["groups"]:
        per_bin[int(bin_index)].update(int(index) for index in indices)
    for bin_index, midpoint in enumerate(observed_mask["centers"]):
        contributing = per_bin[bin_index]
        measurement = contributing.intersection(measurement_active_rows)
        chronology = contributing.intersection(chronology_active_rows)
        both = measurement.intersection(chronology)
        either = measurement.union(chronology)
        denominator = len(contributing)
        rows.append({
            "midpoint_Ma": float(midpoint),
            "contributing_measurements": denominator,
            "measurement_active_measurements": len(measurement),
            "measurement_coverage_proportion": len(measurement) / denominator,
            "chronology_active_measurements": len(chronology),
            "chronology_coverage_proportion": len(chronology) / denominator,
            "active_in_both_components": len(both),
            "both_components_coverage_proportion": len(both) / denominator,
            "active_in_either_component": len(either),
            "either_component_coverage_proportion": len(either) / denominator,
        })
    return pd.DataFrame(rows).sort_values("midpoint_Ma", ascending=False)


def summarize_bins(ev, ages, observed):
    quantiles = np.quantile(ev, [0.025, 0.5, 0.975], axis=0)
    return pd.DataFrame({
        "midpoint_Ma": ages,
        "observed_EV": observed,
        "simulation_median": quantiles[1],
        "simulation_2.5th": quantiles[0],
        "simulation_97.5th": quantiles[2],
        "interval_width": quantiles[2] - quantiles[0],
        "observed_minus_simulation_median": observed - quantiles[1],
        "valid_simulations": np.sum(np.isfinite(ev), axis=0),
    }).sort_values("midpoint_Ma", ascending=False)


def summarize_contrasts(contrasts):
    quantiles = np.quantile(contrasts, [0.025, 0.5, 0.975], axis=0)
    return pd.DataFrame({
        "contrast": CONTRAST_NAMES,
        "simulation_median": quantiles[1],
        "simulation_2.5th": quantiles[0],
        "simulation_97.5th": quantiles[2],
        "interval_width": quantiles[2] - quantiles[0],
        "valid_simulations": np.sum(np.isfinite(contrasts), axis=0),
    })


def compare_single_components(joint_bins, joint_contrasts):
    measurement_ev = pd.read_csv(MEASUREMENT_FORMAL_EV)
    chronology_ev = pd.read_csv(CHRONOLOGY_FORMAL_EV)
    measurement_contrast = pd.read_csv(MEASUREMENT_FORMAL_CONTRAST)
    chronology_contrast = pd.read_csv(CHRONOLOGY_FORMAL_CONTRAST)

    chronology_contrast = chronology_contrast.rename(columns={
        "simulation_median": "MC_median",
        "simulation_2.5th": "MC_2.5th",
        "simulation_97.5th": "MC_97.5th",
    })
    scenarios = [
        ("Measurement B formal 10,000", measurement_ev, "MC_median", "MC_2.5th", "MC_97.5th"),
        ("Chronology B formal 10,000", chronology_ev, "simulation_median", "simulation_2.5th", "simulation_97.5th"),
        ("Joint pilot 1,000", joint_bins, "simulation_median", "simulation_2.5th", "simulation_97.5th"),
    ]
    envelope_rows = []
    for label, table, median_column, low_column, high_column in scenarios:
        width = table[high_column].to_numpy(float) - table[low_column].to_numpy(float)
        midpoint = table.midpoint_Ma.to_numpy(float)
        observed = table.observed_EV.to_numpy(float)
        median = table[median_column].to_numpy(float)
        envelope_rows.append({
            "scenario": label,
            "median_interval_width": float(np.median(width)),
            "maximum_interval_width": float(np.max(width)),
            "maximum_width_midpoint_Ma": float(midpoint[int(np.argmax(width))]),
            "maximum_abs_observed_minus_median": float(np.max(np.abs(observed - median))),
            "maximum_difference_midpoint_Ma": float(midpoint[int(np.argmax(np.abs(observed - median)))]),
        })

    joint_c = joint_contrasts.rename(columns={
        "simulation_median": "MC_median",
        "simulation_2.5th": "MC_2.5th",
        "simulation_97.5th": "MC_97.5th",
    })
    contrast_rows = []
    for label, table in [
        ("Measurement B formal 10,000", measurement_contrast),
        ("Chronology B formal 10,000", chronology_contrast),
        ("Joint pilot 1,000", joint_c),
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


def main():
    measurement_pilot = load_module(MEASUREMENT_PILOT_SCRIPT, "joint_measurement_pilot_core")
    measurement_formal = load_module(MEASUREMENT_FORMAL_SCRIPT, "joint_measurement_formal_core")
    chronology_pilot = load_module(CHRONOLOGY_PILOT_SCRIPT, "joint_chronology_layerb_core")
    primary_chronology = load_module(PRIMARY_CHRONOLOGY_SCRIPT, "joint_primary_chronology_core")

    protected_before = {str(path): sha256(path) for path in PROTECTED_PATHS if path.exists()}
    raw = pd.read_csv(RAW_PATH).reset_index(drop=True)
    raw["record_id"] = raw.record_id.map(clean_id)
    if len(raw) != TOTAL_MEASUREMENTS or raw.record_id.nunique() != TOTAL_RECORDS:
        raise RuntimeError("Retained baseline must contain 10,881 measurements / 70 records")

    # Exact Fig. S6B measurement mappings.
    measurement_params = measurement_pilot.parameter_table()
    layer1_maps, layer1_notes, layer1_sources = measurement_pilot.build_layer1_maps(raw)
    tier2_maps, tier2_audit = measurement_pilot.tier2_maps(raw, measurement_params)
    mappings = layer1_maps + tier2_maps
    mapped_rows = [item.row_index for item in mappings]
    if len(mapped_rows) != 10015 or len(set(mapped_rows)) != 10015:
        raise RuntimeError("Fig. S6B mappings must uniquely cover 10,015 retained observations")
    measurement_audit, measurement_participation, measurement_structure_counts = measurement_execution_audit(
        raw, measurement_params, layer1_maps, layer1_notes, tier2_audit, measurement_formal
    )
    draw_arrays = prepare_measurement_draw_arrays(mappings)

    # Exact Fig. S7B chronology membership, grouping, and bounds.
    chronology_records, chronology_groups, chronology_members = chronology_pilot.load_parameter_table(raw)
    chronology_active = chronology_records[chronology_records.active_flag.eq(1)]
    if len(chronology_active) != 38 or int(chronology_active.n_measurements.sum()) != 8843:
        raise RuntimeError("Fig. S7B chronology scope must contain 38 records / 8,843 measurements")
    group_order = list(chronology_pilot.GROUP_ORDER)
    segment_group = chronology_pilot.SEGMENT_GROUP
    age_row_masks = chronology_row_masks(raw, chronology_members, group_order, segment_group)

    # Independent component RNG streams derived from one master SeedSequence.
    master_sequence = np.random.SeedSequence(MASTER_SEED)
    measurement_sequence, chronology_sequence = master_sequence.spawn(2)
    measurement_rng = np.random.default_rng(measurement_sequence)
    chronology_shifts, chronology_rng_metadata = draw_chronology_shifts(
        chronology_groups, group_order, chronology_sequence
    )
    rng_metadata = pd.concat([
        pd.DataFrame([{
            "component": "measurement",
            "stream": "all locked Fig. S6B measurement draws",
            "master_entropy": int(measurement_sequence.entropy),
            "spawn_key": ".".join(str(value) for value in measurement_sequence.spawn_key),
            "bit_generator": measurement_rng.bit_generator.__class__.__name__,
            "distribution": "mixed: documented Gaussian/bounded; Tier 2 Gaussian/bounded; ODP1049 exact reconstruction",
            "bound_min_kyr": np.nan,
            "bound_max_kyr": np.nan,
            "draw_min_kyr": np.nan,
            "draw_max_kyr": np.nan,
        }]),
        chronology_rng_metadata,
    ], ignore_index=True)

    base_values = raw.value.to_numpy(float)
    base_ages = raw.age_Ma.to_numpy(float)
    categories = raw.category.to_numpy(str)
    record_ids = pd.unique(raw.record_id)
    record_lookup = {record_id: index for index, record_id in enumerate(record_ids)}
    record_code = raw.record_id.map(record_lookup).to_numpy(int)
    record_indices = [np.flatnonzero(record_code == index) for index in range(len(record_ids))]

    # Observed baseline is independently validated through both existing engines.
    iqr_table = pd.read_csv(IQR_PATH)
    iqr_table["record_id"] = iqr_table.record_id.map(clean_id)
    observed_masks = measurement_pilot.build_masks(raw, iqr_table, record_code)
    if set(observed_masks) != {"50"}:
        raise RuntimeError("Joint pilot must remain restricted to the standard 50-kyr reconstruction")
    s1_difference = measurement_pilot.validate_s1_50(observed_masks, measurement_pilot.load_layer1_module())
    record_arrays = primary_chronology.prepare_record_arrays(raw)
    ages, observed_ev, observed_counts, chronology_s1_difference, observed_zero_bins = (
        primary_chronology.validate_observed_baseline(raw, record_arrays)
    )
    if max(s1_difference, chronology_s1_difference) > 1e-11:
        raise RuntimeError("Observed bottom-up baseline no longer reproduces Dataset S1")
    if np.max(np.abs(observed_masks["50"]["observed_ev"] - observed_ev)) > 1e-11:
        raise RuntimeError("Measurement and chronology engines do not share the same observed baseline")

    baseline_z, invalid_baseline_sd = compute_recordwise_z(base_values, record_indices)
    baseline_dynamic_ev, baseline_dynamic_counts, baseline_diagnostics = dynamic_ev(
        primary_chronology, baseline_z, base_ages, record_indices
    )
    if invalid_baseline_sd or np.max(np.abs(baseline_dynamic_ev - observed_ev)) > 1e-11:
        raise RuntimeError("Joint dynamic engine does not reproduce the observed Global EV baseline")
    if not np.array_equal(baseline_dynamic_counts, observed_counts):
        raise RuntimeError("Joint dynamic engine does not reproduce observed eligible record-bin counts")

    measurement_active_rows = set(int(value) for value in mapped_rows)
    chronology_active_records = set(chronology_active.record_id)
    coverage = observed_coverage(
        raw, observed_masks["50"], measurement_active_rows, chronology_active_records
    )

    ev_realizations = np.empty((ITERATIONS, len(ages)), dtype=float)
    contrast_realizations = np.empty((ITERATIONS, len(CONTRAST_NAMES)), dtype=float)
    eligible_counts = np.empty((ITERATIONS, len(ages)), dtype=np.int16)
    peak_midpoints = np.empty(ITERATIONS, dtype=float)
    measurement_delta_signature = np.empty(ITERATIONS, dtype=float)
    diagnostics = Counter()
    max_abs_z = 0.0
    min_ev, max_ev = math.inf, -math.inf
    category_min = {category: math.inf for category in np.unique(categories)}
    category_max = {category: -math.inf for category in np.unique(categories)}

    segment_record_indices = [
        np.flatnonzero(raw.record_id.to_numpy(str) == record_id)
        for record_id in chronology_members[segment_group]
    ]

    for iteration in range(ITERATIONS):
        values = perturb_measurements(base_values, draw_arrays, measurement_rng)
        measurement_delta_signature[iteration] = float(np.mean(values[list(measurement_active_rows)] - base_values[list(measurement_active_rows)]))
        if not np.all(np.isfinite(values)):
            diagnostics["nonfinite_measurement_value_iterations"] += 1
        diagnostics["temperature_below_absolute_zero_events"] += int(
            np.sum(values[categories == "temperature"] < -273.15)
        )
        for category in category_min:
            part = values[categories == category]
            category_min[category] = min(category_min[category], float(np.min(part)))
            category_max[category] = max(category_max[category], float(np.max(part)))

        z, invalid_sd = compute_recordwise_z(values, record_indices)
        diagnostics["zero_or_invalid_record_SD_events"] += invalid_sd
        if not np.all(np.isfinite(z)):
            diagnostics["nonfinite_z_iterations"] += 1
        max_abs_z = max(max_abs_z, float(np.nanmax(np.abs(z))))

        perturbed_ages = perturb_ages(
            base_ages, chronology_shifts, iteration, group_order, age_row_masks
        )
        for indices in segment_record_indices:
            original_order = np.argsort(base_ages[indices], kind="mergesort")
            if np.any(np.diff(perturbed_ages[indices][original_order]) < -1e-12):
                diagnostics["age_inversion_events"] += 1
            post = base_ages[indices] < KPG_BOUNDARY_MA
            if np.any(post) and np.any(~post):
                if float(np.max(perturbed_ages[indices][post])) >= float(np.min(perturbed_ages[indices][~post])) + 1e-12:
                    diagnostics["segment_boundary_order_violations"] += 1

        ev, counts, iteration_diagnostics = dynamic_ev(
            primary_chronology, z, perturbed_ages, record_indices
        )
        diagnostics.update(iteration_diagnostics)
        eligible_counts[iteration] = counts
        ev_realizations[iteration] = ev
        contrast_realizations[iteration] = phase_contrasts(ev, ages)
        if not np.all(np.isfinite(ev)):
            diagnostics["simulations_with_nonfinite_or_missing_EV"] += 1
        diagnostics["empty_bin_events"] += int(np.sum(counts == 0))
        min_ev = min(min_ev, float(np.nanmin(ev)))
        max_ev = max(max_ev, float(np.nanmax(ev)))
        peak_midpoints[iteration] = ages[int(np.nanargmax(ev))]

        if (iteration + 1) % 100 == 0:
            print(f"Joint pilot: {iteration + 1}/{ITERATIONS} simulations", flush=True)

    bin_summary = summarize_bins(ev_realizations, ages, observed_ev)
    contrast_summary = summarize_contrasts(contrast_realizations)
    peak_counts = pd.Series(peak_midpoints).value_counts().sort_index()
    peak_frequency = pd.DataFrame({
        "peak_midpoint_Ma": peak_counts.index.astype(float),
        "frequency": peak_counts.values.astype(int),
        "proportion": peak_counts.values / ITERATIONS,
    })
    eligibility_by_bin = pd.DataFrame({
        "midpoint_Ma": ages,
        "observed_eligible_record_bins": observed_counts,
        "simulation_min": np.min(eligible_counts, axis=0),
        "simulation_2.5th": np.quantile(eligible_counts, 0.025, axis=0),
        "simulation_median": np.median(eligible_counts, axis=0),
        "simulation_97.5th": np.quantile(eligible_counts, 0.975, axis=0),
        "simulation_max": np.max(eligible_counts, axis=0),
        "simulations_with_empty_bin": np.sum(eligible_counts == 0, axis=0),
    }).sort_values("midpoint_Ma", ascending=False)
    total_eligible = eligible_counts.sum(axis=1)
    eligibility_distribution = pd.DataFrame([{
        "metric": "Total eligible record-bins per simulation",
        "observed": int(observed_counts.sum()),
        "minimum": int(total_eligible.min()),
        "2.5th": float(np.quantile(total_eligible, 0.025)),
        "median": float(np.median(total_eligible)),
        "97.5th": float(np.quantile(total_eligible, 0.975)),
        "maximum": int(total_eligible.max()),
    }])

    envelope_comparison, contrast_comparison = compare_single_components(
        bin_summary, contrast_summary
    )

    # RNG-stream audit: construction is the proof of independence; empirical
    # correlations are retained only as a finite-sample diagnostic.
    rng_correlation_rows = []
    for group_id in group_order:
        value = float(np.corrcoef(measurement_delta_signature, chronology_shifts[group_id])[0, 1])
        rng_correlation_rows.append({
            "measurement_signature_vs_chronology_stream": group_id,
            "Pearson_correlation": value,
            "absolute_correlation": abs(value),
            "interpretation": "Finite-sample diagnostic only; independence is guaranteed by distinct SeedSequence spawn keys.",
        })
    rng_correlation = pd.DataFrame(rng_correlation_rows)

    joint_record_audit = measurement_audit[[
        "record_id", "n_measurements", "actually_active", "actually_perturbed_measurements",
        "primary_structure", "table_distribution", "execution_implementation"
    ]].rename(columns={
        "actually_active": "measurement_active",
        "actually_perturbed_measurements": "measurement_perturbed_measurements",
        "primary_structure": "measurement_primary_structure",
        "table_distribution": "measurement_distribution",
        "execution_implementation": "measurement_implementation",
    }).merge(
        chronology_records[[
            "record_id", "active_flag", "perturbation_group_id", "recommended_perturbation_structure",
            "bound_min_kyr", "bound_max_kyr", "execution_status"
        ]],
        on="record_id", how="left", validate="one_to_one"
    ).rename(columns={
        "active_flag": "chronology_active",
        "perturbation_group_id": "chronology_group_id",
        "recommended_perturbation_structure": "chronology_structure",
        "bound_min_kyr": "chronology_bound_min_kyr",
        "bound_max_kyr": "chronology_bound_max_kyr",
        "execution_status": "chronology_execution_status",
    })

    # Write all machine-readable pilot outputs before workbook construction.
    tables = {
        "joint_50kyr_Global_EV_summary.csv": bin_summary,
        "joint_phase_contrasts.csv": contrast_summary,
        "joint_peak_midpoint_frequency.csv": peak_frequency,
        "joint_eligibility_by_bin.csv": eligibility_by_bin,
        "joint_eligibility_distribution.csv": eligibility_distribution,
        "joint_observed_perturbation_coverage_by_bin.csv": coverage,
        "joint_record_execution_audit.csv": joint_record_audit,
        "measurement_structure_participation.csv": measurement_participation,
        "chronology_group_parameters.csv": chronology_groups,
        "rng_streams.csv": rng_metadata,
        "rng_independence_diagnostic.csv": rng_correlation,
        "single_component_envelope_comparison.csv": envelope_comparison,
        "single_component_phase_contrast_comparison.csv": contrast_comparison,
    }
    for filename, table in tables.items():
        table.to_csv(RESULT_DIR / filename, index=False)

    pd.DataFrame(
        ev_realizations,
        columns=[f"EV_{age:.2f}_Ma" for age in ages],
    ).assign(simulation=np.arange(1, ITERATIONS + 1)).loc[:, ["simulation"] + [f"EV_{age:.2f}_Ma" for age in ages]].to_csv(
        RESULT_DIR / "joint_Global_EV_realizations_50kyr_1000.csv", index=False
    )
    pd.DataFrame(
        contrast_realizations, columns=CONTRAST_NAMES
    ).assign(simulation=np.arange(1, ITERATIONS + 1)).loc[:, ["simulation"] + CONTRAST_NAMES].to_csv(
        RESULT_DIR / "joint_phase_contrast_realizations_1000.csv", index=False
    )
    pd.DataFrame(
        eligible_counts,
        columns=[f"eligible_records_{age:.2f}_Ma" for age in ages],
    ).assign(simulation=np.arange(1, ITERATIONS + 1)).loc[:, ["simulation"] + [f"eligible_records_{age:.2f}_Ma" for age in ages]].to_csv(
        RESULT_DIR / "joint_eligible_record_bin_counts_1000.csv", index=False
    )
    shift_table = pd.DataFrame({"simulation": np.arange(1, ITERATIONS + 1)})
    for group_id in group_order:
        shift_table[f"{group_id}_shift_kyr"] = chronology_shifts[group_id]
    shift_table.to_csv(RESULT_DIR / "joint_chronology_shift_draws_1000.csv", index=False)

    np.savez_compressed(
        RESULT_DIR / "joint_complete_realizations_50kyr_1000.npz",
        EV_50kyr=ev_realizations,
        age_50kyr=ages,
        observed_EV_50kyr=observed_ev,
        phase_contrasts=contrast_realizations,
        peak_midpoint_Ma=peak_midpoints,
        eligible_record_bin_counts=eligible_counts,
        measurement_delta_signature=measurement_delta_signature,
        **{f"shift_{group_id}": chronology_shifts[group_id] for group_id in group_order},
    )

    width = bin_summary.interval_width.to_numpy(float)
    observed_difference = np.abs(
        bin_summary.observed_EV.to_numpy(float) - bin_summary.simulation_median.to_numpy(float)
    )
    measurement_active_records = set(joint_record_audit.loc[joint_record_audit.measurement_active, "record_id"])
    chronology_active_record_set = set(joint_record_audit.loc[joint_record_audit.chronology_active.eq(1), "record_id"])
    union_records = measurement_active_records | chronology_active_record_set
    intersection_records = measurement_active_records & chronology_active_record_set
    raw_counts = raw.groupby("record_id").size()

    protected_after = {str(path): sha256(path) for path in PROTECTED_PATHS if path.exists()}
    protected_check = pd.DataFrame([
        {
            "path": path,
            "sha256_before": before,
            "sha256_after": protected_after.get(path),
            "unchanged": before == protected_after.get(path),
        }
        for path, before in protected_before.items()
    ])
    protected_check.to_csv(RESULT_DIR / "protected_files_hash_check.csv", index=False)

    diagnostics.update({
        "maximum_absolute_simulated_z_score": max_abs_z,
        "minimum_simulated_Global_EV": min_ev,
        "maximum_simulated_Global_EV": max_ev,
    })
    metadata = {
        "analysis": "Joint conservative measurement and chronological sensitivity pilot - standard 50-kyr Global EV",
        "iterations": ITERATIONS,
        "master_random_seed": MASTER_SEED,
        "rng_implementation": "numpy.random.SeedSequence(master_seed).spawn(2); child 0 -> one PCG64 measurement stream; child 1 -> spawn(6) independent PCG64 chronology-group streams",
        "measurement_component_spawn_key": list(measurement_sequence.spawn_key),
        "chronology_component_spawn_key": list(chronology_sequence.spawn_key),
        "measurement_and_chronology_draws_independent": True,
        "retained_records": TOTAL_RECORDS,
        "retained_measurements": TOTAL_MEASUREMENTS,
        "measurement_active_records": 53,
        "measurement_perturbed_measurements": 10015,
        "measurement_coverage_percent": 100 * 10015 / TOTAL_MEASUREMENTS,
        "measurement_structure_counts": measurement_structure_counts,
        "chronology_active_records": 38,
        "chronology_active_measurements": 8843,
        "chronology_coverage_percent": 100 * 8843 / TOTAL_MEASUREMENTS,
        "chronology_independent_groups": len(group_order),
        "active_record_union_count": len(union_records),
        "active_record_intersection_count": len(intersection_records),
        "active_measurement_entries_in_record_union": int(raw_counts.loc[list(union_records)].sum()),
        "active_measurement_entries_in_record_intersection": int(raw_counts.loc[list(intersection_records)].sum()),
        "chronology_rebinned_each_simulation": True,
        "observed_record_bin_mask_fixed": False,
        "eligibility_rule": "dynamic n >= 2 and nonzero IQR",
        "record_wide_z_standardization_recomputed_each_simulation": True,
        "common_absolute_timescale_systematics_included": False,
        "other_bin_widths_analyzed": False,
        "s1_max_absolute_validation_difference": float(max(s1_difference, chronology_s1_difference)),
        "observed_zero_IQR_record_bins_excluded": int(observed_zero_bins),
        "observed_vs_simulation_median_max_absolute_difference": float(np.max(observed_difference)),
        "maximum_difference_midpoint_Ma": float(bin_summary.iloc[int(np.argmax(observed_difference))].midpoint_Ma),
        "interval_width_median": float(np.median(width)),
        "interval_width_maximum": float(np.max(width)),
        "maximum_width_midpoint_Ma": float(bin_summary.iloc[int(np.argmax(width))].midpoint_Ma),
        "modal_peak_midpoint_Ma": float(pd.Series(peak_midpoints).mode().iloc[0]),
        "protected_input_and_document_hashes_unchanged": bool(protected_check.unchanged.all()),
        "execution_diagnostics": dict(diagnostics),
        "simulated_measurement_min_by_category": category_min,
        "simulated_measurement_max_by_category": category_max,
        "implementation_notes": [
            "This is a joint bottom-up simulation; Fig. S6B and Fig. S7B interval limits were not added or otherwise combined.",
            "Measurement and chronology perturbations are sampled independently within each simulation realization.",
            "The seven non-executable record-wide calibration/reconstruction measurement structures and all other fixed measurements follow the locked Fig. S6B rules.",
            "Chronology perturbations use the exact six Fig. S7B group memberships and bounds; ODP 690 is perturbed only in the specified post-K-Pg segment.",
            "The K-Pg boundary remains 66.016 Ma for the ODP 690 segment definition; 66.05 Ma remains a normal 50-kyr bin midpoint and the last Phase 2 midpoint.",
        ],
    }
    (RESULT_DIR / "joint_pilot_metadata.json").write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    source_manifest = []
    for path in [
        RAW_PATH, IQR_PATH, S1_PATH, MEASUREMENT_PARAMETER_PATH,
        CHRONOLOGY_PARAMETER_PATH, MEASUREMENT_PILOT_SCRIPT,
        MEASUREMENT_FORMAL_SCRIPT, CHRONOLOGY_PILOT_SCRIPT,
        PRIMARY_CHRONOLOGY_SCRIPT, SCRIPT,
    ]:
        source_manifest.append({"file": path.name, "path": str(path), "sha256": sha256(path)})
    for name, path in layer1_sources.items():
        source_manifest.append({"file": str(name), "path": str(path), "sha256": sha256(path)})
    pd.DataFrame(source_manifest).drop_duplicates("sha256").to_csv(
        RESULT_DIR / "input_source_hash_manifest.csv", index=False
    )

    print(json.dumps(metadata, indent=2, ensure_ascii=False))
    print("\nJOINT PHASE CONTRASTS")
    print(contrast_summary.to_string(index=False))
    print("\nSINGLE-COMPONENT COMPARISON")
    print(envelope_comparison.to_string(index=False))


if __name__ == "__main__":
    main()
