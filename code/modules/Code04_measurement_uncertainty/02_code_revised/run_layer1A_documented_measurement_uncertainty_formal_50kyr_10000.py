#!/usr/bin/env python3
"""Formal Layer 1A documented measurement-uncertainty propagation.

Locked scope: 10,000 iterations, standard 50-kyr Global EV only, fixed
chronology, the validated observation-level mappings only, and the observed
eligible/non-zero record-bin mask.  The RNG call order is identical to the
validated 1,000-iteration pilot so the first 1,000 realizations can be checked
directly.
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


SEED = 20260823
ITERATIONS = 10000
BATCH_SIZE = 250
PREFIXES = (1000, 5000, 10000)

SCRIPT = Path(__file__).resolve()
WORK_ROOT = SCRIPT.parent.parent
DATA_DIR = WORK_ROOT / "01_data_revised"
FORMAL_DATA_DIR = DATA_DIR / "formal_datasets"
INVENTORY_DIR = WORK_ROOT / "03_results" / "EV_uncertainty_inventory"
RESULT_DIR = INVENTORY_DIR / "Layer1A_documented_measurement_uncertainty_formal_50kyr_10000"
RESULT_DIR.mkdir(parents=True, exist_ok=True)

PILOT_SCRIPT = WORK_ROOT / "02_code_revised" / "run_layer1_measurement_uncertainty_pilot.py"
PILOT_DIR = INVENTORY_DIR / "Layer1_measurement_uncertainty_pilot_1000"
PILOT_NPZ = PILOT_DIR / "Layer1_MC_Global_EV_realizations_all_widths.npz"

RAW_PATH = DATA_DIR / "retained_measurements_with_zscores_v4.csv"
IQR_PATH = DATA_DIR / "record_level_IQR_all_bin_widths_revised_v4_n2.csv"
S1_PATH = FORMAL_DATA_DIR / "Dataset_S1_revised_v6.xlsx"

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
    FORMAL_DATA_DIR / "Dataset_S1_revised_v6.xlsx",
    FORMAL_DATA_DIR / "Dataset_S2_revised.xlsx",
    FORMAL_DATA_DIR / "Dataset_S3_revised_v3.xlsx",
    FORMAL_DATA_DIR / "Dataset_S4_revised_v3.xlsx",
    FORMAL_DATA_DIR / "Dataset_S5_revised_v4.xlsx",
    FORMAL_DATA_DIR / "Dataset_S6_revised_v4.xlsx",
    WORK_ROOT / "05_manuscript_revised" / "Manuscript_working_stage3_sampling_density.docx",
    WORK_ROOT / "06_supplement_revised" / "Supplementary_materials_working_stage3_sampling_density.docx",
]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def load_pilot_module():
    spec = importlib.util.spec_from_file_location("layer1_documented_pilot", PILOT_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError("Unable to load the validated Layer 1 pilot")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def build_50kyr_masks(pilot, raw: pd.DataFrame, iqr: pd.DataFrame, record_code: np.ndarray):
    locked_specs = dict(pilot.BIN_SPECS)
    try:
        pilot.BIN_SPECS = {"50": locked_specs["50"]}
        return pilot.build_record_bin_masks(raw, iqr, record_code)
    finally:
        pilot.BIN_SPECS = locked_specs


def coverage_table(raw: pd.DataFrame, info: dict, active_rows: set[int]) -> pd.DataFrame:
    rows = []
    for b, age in enumerate(info["centers"]):
        groups = [ix for bb, ix in info["groups"] if bb == b]
        contributing = np.concatenate(groups)
        active = sum(int(i) in active_rows for i in contributing)
        rows.append({
            "midpoint_Ma": float(age),
            "contributing_measurements": int(len(contributing)),
            "perturbed_contributing_measurements": int(active),
            "perturbation_coverage_proportion": float(active / len(contributing)),
            "perturbation_coverage_percent": float(100 * active / len(contributing)),
            "contributing_record_bins": int(info["bin_counts"][b]),
        })
    return pd.DataFrame(rows)


def simulate_batched(pilot, raw: pd.DataFrame, mappings: list, masks: dict):
    if set(masks) != {"50"}:
        raise RuntimeError(f"Formal Layer 1A run is 50-kyr only; received {list(masks)}")

    rng = np.random.default_rng(SEED)
    record_ids = pd.unique(raw.record_id)
    record_lookup = {rid: i for i, rid in enumerate(record_ids)}
    record_code = raw.record_id.map(record_lookup).to_numpy(int)
    record_indices = [np.flatnonzero(record_code == code) for code in range(len(record_ids))]
    base = raw.value.to_numpy(float)
    categories = raw.category.to_numpy(str)

    normal = [m for m in mappings if m.distribution == "normal"]
    bounded = [m for m in mappings if m.distribution == "bounded_uniform"]
    isotope = [m for m in mappings if m.distribution == "odp1049_isotope_normal"]
    normal_ix = np.array([m.row_index for m in normal], dtype=int)
    normal_sigma = np.array([m.sigma for m in normal], dtype=float)
    bounded_ix = np.array([m.row_index for m in bounded], dtype=int)
    bounded_u = np.array([m.half_width for m in bounded], dtype=float)
    isotope_ix = np.array([m.row_index for m in isotope], dtype=int)
    isotope_sigma = np.array([m.sigma for m in isotope], dtype=float)
    isotope_base = np.array([m.isotope_value for m in isotope], dtype=float)
    isotope_species = np.array([m.species for m in isotope])

    info = masks["50"]
    outputs = np.empty((ITERATIONS, len(info["centers"])), dtype=float)
    contrasts = np.empty((ITERATIONS, 3), dtype=float)
    anomaly = {
        "nonfinite_value_iterations": 0,
        "nonfinite_z_iterations": 0,
        "nonfinite_EV_iterations": 0,
        "zero_or_invalid_record_SD_events": 0,
        "missing_EV_trajectory_iterations": 0,
        "temperature_below_absolute_zero_events": 0,
    }
    max_abs_z = 0.0
    min_ev, max_ev = math.inf, -math.inf
    category_min = {c: math.inf for c in np.unique(categories)}
    category_max = {c: -math.inf for c in np.unique(categories)}

    for start in range(0, ITERATIONS, BATCH_SIZE):
        stop = min(start + BATCH_SIZE, ITERATIONS)
        n_batch = stop - start
        values = np.tile(base, (n_batch, 1))

        # Keep the pilot's RNG order iteration-by-iteration.
        for j in range(n_batch):
            if len(normal_ix):
                values[j, normal_ix] += rng.normal(0.0, normal_sigma)
            if len(bounded_ix):
                values[j, bounded_ix] += rng.uniform(-bounded_u, bounded_u)
            if len(isotope_ix):
                iso = isotope_base + rng.normal(0.0, isotope_sigma)
                s = isotope_species == "S"
                values[j, isotope_ix[s]] = 0.03 * iso[s] ** 2 - 4.43 * iso[s] + 10.2855
                values[j, isotope_ix[~s]] = 0.03 * iso[~s] ** 2 - 4.427 * iso[~s] + 10.064075

        anomaly["nonfinite_value_iterations"] += int((~np.all(np.isfinite(values), axis=1)).sum())
        anomaly["temperature_below_absolute_zero_events"] += int(np.sum(values[:, categories == "temperature"] < -273.15))
        for cat in category_min:
            vv = values[:, categories == cat]
            category_min[cat] = min(category_min[cat], float(np.min(vv)))
            category_max[cat] = max(category_max[cat], float(np.max(vv)))

        z = np.empty_like(values)
        for ix in record_indices:
            vv = values[:, ix]
            means = vv.mean(axis=1)
            sds = vv.std(axis=1, ddof=1)
            bad = (~np.isfinite(sds)) | (sds <= 0)
            anomaly["zero_or_invalid_record_SD_events"] += int(bad.sum())
            z[:, ix] = (vv - means[:, None]) / sds[:, None]
        anomaly["nonfinite_z_iterations"] += int((~np.all(np.isfinite(z), axis=1)).sum())
        max_abs_z = max(max_abs_z, float(np.nanmax(np.abs(z))))

        sums = np.zeros((n_batch, len(info["centers"])), dtype=float)
        for b, ix in info["groups"]:
            q = np.quantile(z[:, ix], [0.25, 0.75], axis=1, method="linear")
            sums[:, b] += q[1] - q[0]
        ev = sums / info["bin_counts"][None, :]
        outputs[start:stop] = ev
        nonfinite_ev = ~np.all(np.isfinite(ev), axis=1)
        anomaly["nonfinite_EV_iterations"] += int(nonfinite_ev.sum())
        anomaly["missing_EV_trajectory_iterations"] += int(nonfinite_ev.sum())
        min_ev = min(min_ev, float(np.nanmin(ev)))
        max_ev = max(max_ev, float(np.nanmax(ev)))

        phase_means = {}
        ages = info["centers"]
        for name, (lo, hi) in PHASES.items():
            sel = (ages >= lo - 1e-10) & (ages <= hi + 1e-10)
            phase_means[name] = ev[:, sel].mean(axis=1)
        contrasts[start:stop] = np.column_stack([
            phase_means["Phase 1"] - phase_means["pre-escalation"],
            phase_means["Phase 2"] - phase_means["pre-escalation"],
            phase_means["Recovery"] - phase_means["Phase 2"],
        ])

    anomaly.update({
        "maximum_absolute_simulated_z_score": max_abs_z,
        "minimum_simulated_global_EV_50kyr": min_ev,
        "maximum_simulated_global_EV_50kyr": max_ev,
        "simulated_measurement_min_by_category": category_min,
        "simulated_measurement_max_by_category": category_max,
    })
    return outputs, contrasts, anomaly


def summarize_prefixes(outputs: np.ndarray, contrasts: np.ndarray, masks: dict):
    ages = masks["50"]["centers"]
    observed = masks["50"]["observed_ev"]
    bin_rows, metric_rows, contrast_rows = [], [], []
    stats = {}
    for prefix in PREFIXES:
        q = np.quantile(outputs[:prefix], [0.025, 0.5, 0.975], axis=0)
        stats[prefix] = q
        for age, obs, lo, med, hi in zip(ages, observed, q[0], q[1], q[2]):
            bin_rows.append({
                "iterations": prefix,
                "midpoint_Ma": float(age),
                "observed_EV": float(obs),
                "MC_median": float(med),
                "MC_2.5th": float(lo),
                "MC_97.5th": float(hi),
                "envelope_width": float(hi - lo),
                "observed_minus_MC_median": float(obs - med),
            })
        cq = np.quantile(contrasts[:prefix], [0.025, 0.5, 0.975], axis=0)
        for j, name in enumerate(CONTRAST_NAMES):
            contrast_rows.append({
                "iterations": prefix,
                "contrast": name,
                "MC_median": float(cq[1, j]),
                "MC_2.5th": float(cq[0, j]),
                "MC_97.5th": float(cq[2, j]),
                "interval_width": float(cq[2, j] - cq[0, j]),
            })

    final = stats[10000]
    for prefix in PREFIXES:
        q = stats[prefix]
        for label, x, ref in [
            ("MC_median", q[1], final[1]),
            ("MC_2.5th", q[0], final[0]),
            ("MC_97.5th", q[2], final[2]),
            ("envelope_width", q[2] - q[0], final[2] - final[0]),
        ]:
            diff = np.abs(x - ref)
            metric_rows.append({
                "iterations": prefix,
                "statistic": label,
                "median_absolute_difference_vs_10000": float(np.median(diff)),
                "maximum_absolute_difference_vs_10000": float(np.max(diff)),
                "midpoint_of_maximum_difference_Ma": float(ages[int(np.argmax(diff))]),
            })
    return pd.DataFrame(bin_rows), pd.DataFrame(metric_rows), pd.DataFrame(contrast_rows)


def mapping_audit(pilot, raw: pd.DataFrame, mappings: list, notes: dict) -> pd.DataFrame:
    mapped_by_record = pd.Series([raw.loc[m.row_index, "record_id"] for m in mappings]).value_counts().to_dict()
    distribution_by_record, source_by_record = {}, {}
    for m in mappings:
        rid = raw.loc[m.row_index, "record_id"]
        distribution_by_record.setdefault(rid, set()).add(m.distribution)
        source_by_record.setdefault(rid, set()).add(m.source)
    rows = []
    for rid in pilot.LAYER1_RECORDS:
        retained_n = int((raw.record_id == rid).sum())
        mapped_n = int(mapped_by_record.get(rid, 0))
        rows.append({
            "record_id": rid,
            "retained_measurements": retained_n,
            "perturbed_measurements": mapped_n,
            "coverage_percent": 100 * mapped_n / retained_n if retained_n else 0.0,
            "formal_status": "perturbed" if mapped_n else "fixed_unresolved_observation_mapping",
            "distribution": "; ".join(sorted(distribution_by_record.get(rid, {"fixed"}))),
            "source_uncertainty": "; ".join(sorted(source_by_record.get(rid, {"not executable at observation level"}))),
            "mapping_note": notes.get(rid, ""),
        })
    return pd.DataFrame(rows)


def main():
    pilot = load_pilot_module()
    if pilot.SEED != SEED or pilot.ITERATIONS != 1000:
        raise RuntimeError("Validated Layer 1A pilot seed/iteration settings changed")
    pilot.validate_layer1_parameter_membership()

    protected_before = {str(p): sha256(p) for p in PROTECTED_PATHS if p.exists()}
    sources = pilot.ensure_sources()
    raw = pd.read_csv(RAW_PATH).reset_index(drop=True)
    if len(raw) != 10881 or raw.record_id.nunique() != 70:
        raise RuntimeError("Retained baseline must contain 10,881 measurements and 70 records")
    iqr = pd.read_csv(IQR_PATH)

    mappings, notes = [], {}
    for mapper in (pilot.map_temperature_sources, pilot.map_jouini_sources, pilot.map_os_sources):
        mm, nn = mapper(raw, sources)
        mappings.extend(mm)
        notes.update(nn)
    row_ids = [m.row_index for m in mappings]
    if len(row_ids) != len(set(row_ids)):
        raise RuntimeError("An observation received multiple Layer 1A mappings")

    audit = mapping_audit(pilot, raw, mappings, notes)
    if int((audit.perturbed_measurements > 0).sum()) != 11 or int(audit.perturbed_measurements.sum()) != 375:
        raise RuntimeError("Locked Layer 1A scope must be 11 active records / 375 measurements")

    record_ids = pd.unique(raw.record_id)
    record_lookup = {rid: i for i, rid in enumerate(record_ids)}
    record_code = raw.record_id.map(record_lookup).to_numpy(int)
    masks = build_50kyr_masks(pilot, raw, iqr, record_code)
    s1_diff = pilot.validate_s1_50(masks)
    if s1_diff > 1e-12:
        raise RuntimeError(f"Observed EV differs from Dataset S1: {s1_diff}")

    coverage = coverage_table(raw, masks["50"], set(row_ids))
    outputs, contrasts, anomaly = simulate_batched(pilot, raw, mappings, masks)

    if not PILOT_NPZ.exists():
        raise RuntimeError("Validated Layer 1A pilot realizations are missing")
    pilot_npz = np.load(PILOT_NPZ)
    pilot_ev_diff = float(np.max(np.abs(outputs[:1000] - pilot_npz["EV_50kyr"])))
    pilot_contrast_diff = float(np.max(np.abs(contrasts[:1000] - pilot_npz["phase_contrasts_50kyr"])))
    if pilot_ev_diff > 1e-12 or pilot_contrast_diff > 1e-12:
        raise RuntimeError(f"Formal prefix fails pilot reproduction: EV={pilot_ev_diff}; contrasts={pilot_contrast_diff}")

    bin_prefix, convergence_metrics, contrast_prefix = summarize_prefixes(outputs, contrasts, masks)
    final_ev = bin_prefix[bin_prefix.iterations.eq(10000)].copy()
    final_contrasts = contrast_prefix[contrast_prefix.iterations.eq(10000)].copy()

    final_ev.to_csv(RESULT_DIR / "Layer1A_formal_50kyr_Global_EV_summary.csv", index=False)
    final_contrasts.to_csv(RESULT_DIR / "Layer1A_formal_phase_contrasts_50kyr.csv", index=False)
    coverage.to_csv(RESULT_DIR / "Layer1A_formal_perturbation_coverage_50kyr.csv", index=False)
    bin_prefix.to_csv(RESULT_DIR / "Layer1A_formal_convergence_bin_prefixes.csv", index=False)
    convergence_metrics.to_csv(RESULT_DIR / "Layer1A_formal_convergence_metrics.csv", index=False)
    contrast_prefix.to_csv(RESULT_DIR / "Layer1A_formal_convergence_phase_contrasts.csv", index=False)
    audit.to_csv(RESULT_DIR / "Layer1A_formal_observation_mapping_audit.csv", index=False)

    protected_after = {str(p): sha256(p) for p in PROTECTED_PATHS if p.exists()}
    protected = pd.DataFrame([
        {"path": p, "sha256_before": h, "sha256_after": protected_after.get(p), "unchanged": h == protected_after.get(p)}
        for p, h in protected_before.items()
    ])
    protected.to_csv(RESULT_DIR / "protected_files_hash_check.csv", index=False)

    widths = final_ev.envelope_width.to_numpy(float)
    diffs = np.abs(final_ev.observed_EV - final_ev.MC_median)
    metadata = {
        "analysis": "Formal Layer 1A documented measurement uncertainty propagation — standard 50-kyr Global EV only",
        "iterations": ITERATIONS,
        "random_seed": SEED,
        "batch_size": BATCH_SIZE,
        "chronology_perturbed": False,
        "tier2_method_based_uncertainty_included": False,
        "other_bin_widths_analyzed": False,
        "retained_operational_records": int(raw.record_id.nunique()),
        "retained_measurements": int(len(raw)),
        "layer1_parameter_records": int(len(audit)),
        "actually_active_records": int((audit.perturbed_measurements > 0).sum()),
        "actually_perturbed_measurements": int(audit.perturbed_measurements.sum()),
        "perturbed_measurement_percent": float(100 * audit.perturbed_measurements.sum() / len(raw)),
        "fixed_unresolved_records": audit.loc[audit.perturbed_measurements.eq(0), "record_id"].tolist(),
        "observed_nonzero_record_bin_mask_fixed": True,
        "s1_50kyr_max_absolute_validation_difference": float(s1_diff),
        "formal_prefix_vs_pilot_max_EV_difference": pilot_ev_diff,
        "formal_prefix_vs_pilot_max_contrast_difference": pilot_contrast_diff,
        "observed_EV_vs_MC_median_max_absolute_difference": float(np.max(diffs)),
        "max_difference_midpoint_Ma": float(final_ev.iloc[int(np.argmax(diffs))].midpoint_Ma),
        "50kyr_envelope_width_median": float(np.median(widths)),
        "50kyr_envelope_width_maximum": float(np.max(widths)),
        "50kyr_maximum_width_midpoint_Ma": float(final_ev.iloc[int(np.argmax(widths))].midpoint_Ma),
        "perturbation_coverage_minimum": float(coverage.perturbation_coverage_proportion.min()),
        "perturbation_coverage_median": float(coverage.perturbation_coverage_proportion.median()),
        "perturbation_coverage_maximum": float(coverage.perturbation_coverage_proportion.max()),
        "protected_dataset_and_document_hashes_unchanged": bool(protected.unchanged.all()),
        "anomaly_diagnostics": anomaly,
        "implementation_notes": [
            "Only the 11 records / 375 observations with reliable one-to-one documented uncertainty mapping were perturbed.",
            "ODP525A Os, ODP577 Os, and IODP U1403 Os remained fixed.",
            "Every realization restarted from retained raw measurements and recomputed record-wide sample-SD standardization, z-scores, within-record IQR, and record-weighted Global EV.",
            "Ages, bin assignments, and the observed eligible/non-zero record-bin mask were fixed.",
            "The first 1,000 realizations reproduce the validated pilot exactly within 1e-12.",
        ],
    }
    (RESULT_DIR / "formal_metadata.json").write_text(json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8")

    manifest = []
    for p in [RAW_PATH, IQR_PATH, S1_PATH, PILOT_SCRIPT, PILOT_NPZ, SCRIPT]:
        manifest.append({"file": p.name, "sha256": sha256(p), "source": str(p)})
    for name, p in sources.items():
        manifest.append({"file": name, "sha256": sha256(p), "source": str(p)})
    pd.DataFrame(manifest).drop_duplicates("sha256").to_csv(RESULT_DIR / "input_source_hash_manifest.csv", index=False)

    np.savez_compressed(
        RESULT_DIR / "Layer1A_formal_Global_EV_realizations_50kyr_10000.npz",
        EV_50kyr=outputs,
        age_50kyr=masks["50"]["centers"],
        phase_contrasts_50kyr=contrasts,
    )

    print(json.dumps(metadata, indent=2, ensure_ascii=False))
    print("\nFORMAL PHASE CONTRASTS")
    print(final_contrasts.to_string(index=False))
    print("\nCONVERGENCE METRICS")
    print(convergence_metrics.to_string(index=False))


if __name__ == "__main__":
    main()
