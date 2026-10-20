#!/usr/bin/env python3
"""Formal 10,000-iteration Layer 2 measurement-sensitivity analysis.

Scope is restricted to the standard 50-kyr Global EV reconstruction. The
validated pilot's seed, uncertainty mappings, perturbation structures, fixed
records, chronology, bin assignments, and observed non-zero record-bin mask
are retained. Prefixes of the same realization sequence provide the 1,000 vs
5,000 vs 10,000 convergence check.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import math
import sys
import unicodedata
from pathlib import Path

import numpy as np
import pandas as pd


SEED = 20260824
ITERATIONS = 10000
BATCH_SIZE = 250
PREFIXES = (1000, 5000, 10000)

SCRIPT = Path(__file__).resolve()
WORK_ROOT = SCRIPT.parent.parent
DATA_DIR = WORK_ROOT / "01_data_revised"
FORMAL_DIR = DATA_DIR / "formal_datasets"
INVENTORY_DIR = WORK_ROOT / "03_results" / "EV_uncertainty_inventory"
RESULT_DIR = INVENTORY_DIR / "Layer2_conservative_measurement_sensitivity_formal_50kyr_10000"
RESULT_DIR.mkdir(parents=True, exist_ok=True)

PILOT_SCRIPT = WORK_ROOT / "02_code_revised" / "run_layer2_conservative_measurement_sensitivity_pilot.py"
PILOT_RESULT_DIR = INVENTORY_DIR / "Layer2_conservative_measurement_sensitivity_pilot_50kyr_1000"
PILOT_NPZ = PILOT_RESULT_DIR / "Layer2_Global_EV_realizations_50kyr.npz"

RAW_PATH = DATA_DIR / "retained_measurements_with_zscores_v4.csv"
IQR_PATH = DATA_DIR / "record_level_IQR_all_bin_widths_revised_v4_n2.csv"
S1_PATH = FORMAL_DIR / "Dataset_S1_revised_v6.xlsx"
PARAM_PATH = INVENTORY_DIR / "Executable_uncertainty_parameter_table.xlsx"

PHASES = {
    "pre-escalation": (66.95, 67.20),
    "Phase 1": (66.45, 66.90),
    "Phase 2": (66.05, 66.40),
    "Recovery": (64.80, 65.95),
}
CONTRAST_NAMES = ["Phase 1 - pre-escalation", "Phase 2 - pre-escalation", "Recovery - Phase 2"]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def load_pilot_module():
    spec = importlib.util.spec_from_file_location("layer2_pilot", PILOT_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError("Unable to load the validated Layer 2 pilot implementation")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def build_execution_audit(raw, params, layer1_maps, layer1_notes, tier2_audit):
    layer1_counts = pd.Series([raw.loc[m.row_index, "record_id"] for m in layer1_maps]).value_counts().to_dict()
    rows = []
    layer1_params = params[params.measurement_analysis_layer.eq("documented measurement propagation")]
    for _, r in layer1_params.iterrows():
        mapped = int(layer1_counts.get(r.record_id, 0))
        rows.append({
            "record_id": r.record_id,
            "n_measurements": int(r.n_measurements),
            "primary_structure": r.measurement_primary_structure,
            "table_distribution": r.measurement_distribution,
            "actually_active": mapped > 0,
            "actually_perturbed_measurements": mapped,
            "random_analytical": bool(r.random_analytical),
            "calibration_reconstruction": bool(r.calibration_reconstruction),
            "systematic_upper_bound": bool(r.systematic_upper_bound),
            "execution_implementation": layer1_notes.get(r.record_id, ""),
            "table_execution_status": r.measurement_execution_status,
        })
    audit = pd.concat([pd.DataFrame(rows), tier2_audit], ignore_index=True)
    if len(audit) != 70 or audit.record_id.nunique() != 70:
        raise RuntimeError("Formal execution audit must contain 70 unique records")
    return audit


def structure_participation(audit):
    rows = []
    primary = [
        ("Documented Layer 1A source-mapped", audit.primary_structure.str.startswith("Tier 1")),
        ("Tier 2 conservative Gaussian observation-level", audit.primary_structure.str.contains("Conservative Gaussian")),
        ("Tier 2 bounded uniform/truncated observation-level", audit.primary_structure.str.contains("Bounded uniform")),
        ("Tier 2 record-wide correlated calibration/reconstruction (fixed)", audit.primary_structure.str.contains("Record-wide correlated")),
        ("Fixed: no defensible quantitative basis", audit.primary_structure.str.contains("No defensible")),
    ]
    for label, selector in primary:
        part = audit[selector]
        rows.append({
            "summary_type": "mutually exclusive primary execution structure",
            "structure": label,
            "records_in_parameter_table": len(part),
            "measurements_in_parameter_table": int(part.n_measurements.sum()),
            "actually_active_records": int(part.actually_active.sum()),
            "actually_perturbed_measurements": int(part.actually_perturbed_measurements.sum()),
        })
    for flag, label in [
        ("random_analytical", "Random analytical flag"),
        ("calibration_reconstruction", "Calibration/reconstruction flag"),
        ("systematic_upper_bound", "Systematic/upper-bound flag"),
    ]:
        active = audit[audit[flag] & audit.actually_active]
        rows.append({
            "summary_type": "non-mutually-exclusive active uncertainty flag",
            "structure": label,
            "records_in_parameter_table": int(audit[flag].sum()),
            "measurements_in_parameter_table": int(audit.loc[audit[flag], "n_measurements"].sum()),
            "actually_active_records": len(active),
            "actually_perturbed_measurements": int(active.actually_perturbed_measurements.sum()),
        })
    return pd.DataFrame(rows)


def simulate_batched(raw, mappings, masks):
    if set(masks) != {"50"}:
        raise RuntimeError(f"Formal run must be 50-kyr only; received masks {list(masks)}")
    rng = np.random.default_rng(SEED)
    base = raw.value.to_numpy(float)
    categories = raw.category.to_numpy(str)
    record_ids = pd.unique(raw.record_id)
    record_lookup = {rid: i for i, rid in enumerate(record_ids)}
    record_code = raw.record_id.map(record_lookup).to_numpy(int)
    record_indices = [np.flatnonzero(record_code == code) for code in range(len(record_ids))]

    normal = [m for m in mappings if m.distribution == "normal_direct"]
    uniform = [m for m in mappings if m.distribution == "uniform_direct"]
    isotope = [m for m in mappings if m.distribution == "odp1049_isotope_normal"]
    normal_ix = np.array([m.row_index for m in normal], dtype=int)
    normal_sigma = np.array([m.sigma for m in normal], dtype=float)
    uniform_ix = np.array([m.row_index for m in uniform], dtype=int)
    uniform_lo = np.array([m.lower for m in uniform], dtype=float)
    uniform_hi = np.array([m.upper for m in uniform], dtype=float)
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

        # Preserve the pilot's exact RNG call order iteration by iteration.
        for j in range(n_batch):
            if len(normal_ix):
                values[j, normal_ix] += rng.normal(0.0, normal_sigma)
            if len(uniform_ix):
                values[j, uniform_ix] += rng.uniform(uniform_lo, uniform_hi)
            if len(isotope_ix):
                iso = isotope_base + rng.normal(0.0, isotope_sigma)
                s = isotope_species == "S"
                values[j, isotope_ix[s]] = 0.03 * iso[s] ** 2 - 4.43 * iso[s] + 10.2855
                values[j, isotope_ix[~s]] = 0.03 * iso[~s] ** 2 - 4.427 * iso[~s] + 10.064075

        nonfinite_value_rows = ~np.all(np.isfinite(values), axis=1)
        anomaly["nonfinite_value_iterations"] += int(nonfinite_value_rows.sum())
        anomaly["temperature_below_absolute_zero_events"] += int(np.sum(values[:, categories == "temperature"] < -273.15))
        for cat in category_min:
            vv = values[:, categories == cat]
            category_min[cat] = min(category_min[cat], float(np.min(vv)))
            category_max[cat] = max(category_max[cat], float(np.max(vv)))

        z = np.empty_like(values)
        invalid_row = np.zeros(n_batch, dtype=bool)
        for ix in record_indices:
            vv = values[:, ix]
            means = vv.mean(axis=1)
            sds = vv.std(axis=1, ddof=1)
            bad = (~np.isfinite(sds)) | (sds <= 0)
            anomaly["zero_or_invalid_record_SD_events"] += int(bad.sum())
            invalid_row |= bad
            z[:, ix] = (vv - means[:, None]) / sds[:, None]
        anomaly["nonfinite_z_iterations"] += int((~np.all(np.isfinite(z), axis=1)).sum())
        max_abs_z = max(max_abs_z, float(np.nanmax(np.abs(z))))

        sums = np.zeros((n_batch, len(info["centers"])), dtype=float)
        for b, ix, _rid in info["groups"]:
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


def summarize_prefixes(outputs, contrasts, masks):
    ages = masks["50"]["centers"]
    observed = masks["50"]["observed_ev"]
    bin_rows = []
    metric_rows = []
    contrast_rows = []
    stats_by_prefix = {}
    for prefix in PREFIXES:
        q = np.quantile(outputs[:prefix], [0.025, 0.5, 0.975], axis=0)
        stats_by_prefix[prefix] = q
        for age, obs, lo, med, hi in zip(ages, observed, q[0], q[1], q[2]):
            bin_rows.append({
                "iterations": prefix,
                "midpoint_Ma": float(age),
                "observed_EV": float(obs),
                "MC_median": float(med),
                "MC_2.5th": float(lo),
                "MC_97.5th": float(hi),
                "envelope_width": float(hi - lo),
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

    final = stats_by_prefix[10000]
    for prefix in PREFIXES:
        q = stats_by_prefix[prefix]
        labels = [("MC_median", q[1], final[1]), ("MC_2.5th", q[0], final[0]), ("MC_97.5th", q[2], final[2]), ("envelope_width", q[2]-q[0], final[2]-final[0])]
        for stat, x, ref in labels:
            diff = np.abs(x - ref)
            metric_rows.append({
                "iterations": prefix,
                "statistic": stat,
                "median_absolute_difference_vs_10000": float(np.median(diff)),
                "maximum_absolute_difference_vs_10000": float(np.max(diff)),
                "midpoint_of_maximum_difference_Ma": float(ages[int(np.argmax(diff))]),
            })
    return pd.DataFrame(bin_rows), pd.DataFrame(metric_rows), pd.DataFrame(contrast_rows)


def main():
    pilot = load_pilot_module()
    if pilot.SEED != SEED or set(pilot.BIN_SPECS) != {"50"}:
        raise RuntimeError("Validated pilot seed/scope no longer matches the locked formal analysis")

    protected_before = {str(p): sha256(p) for p in pilot.PROTECTED_PATHS if p.exists()}
    raw = pd.read_csv(RAW_PATH).reset_index(drop=True)
    raw["record_id"] = raw.record_id.map(lambda x: unicodedata.normalize("NFC", str(x)).strip())
    if len(raw) != 10881 or raw.record_id.nunique() != 70:
        raise RuntimeError("Retained baseline must contain 10,881 measurements and 70 records")
    params = pilot.parameter_table()
    layer1_maps, layer1_notes, sources = pilot.build_layer1_maps(raw)
    tier2_maps, tier2_audit = pilot.tier2_maps(raw, params)
    mappings = layer1_maps + tier2_maps
    row_ids = [m.row_index for m in mappings]
    if len(row_ids) != len(set(row_ids)):
        raise RuntimeError("An observation received more than one uncertainty mapping")

    record_ids = pd.unique(raw.record_id)
    record_lookup = {rid: i for i, rid in enumerate(record_ids)}
    record_code = raw.record_id.map(record_lookup).to_numpy(int)
    iqr = pd.read_csv(IQR_PATH)
    iqr["record_id"] = iqr.record_id.map(lambda x: unicodedata.normalize("NFC", str(x)).strip())
    masks = pilot.build_masks(raw, iqr, record_code)
    s1_diff = pilot.validate_s1_50(masks, pilot.load_layer1_module())
    if s1_diff > 1e-12:
        raise RuntimeError(f"Observed EV does not reproduce Dataset S1: {s1_diff}")

    active_rows = set(row_ids)
    coverage = pilot.coverage_table(raw, masks, active_rows)
    outputs, contrasts, anomaly = simulate_batched(raw, mappings, masks)

    if not PILOT_NPZ.exists():
        raise RuntimeError("Validated 1,000-iteration pilot realization file is missing")
    pilot_npz = np.load(PILOT_NPZ)
    pilot_ev_diff = float(np.max(np.abs(outputs[:1000] - pilot_npz["EV_50kyr"])))
    pilot_contrast_diff = float(np.max(np.abs(contrasts[:1000] - pilot_npz["phase_contrasts_50kyr"])))
    if pilot_ev_diff > 1e-12 or pilot_contrast_diff > 1e-12:
        raise RuntimeError(f"Formal prefix does not reproduce validated pilot: EV={pilot_ev_diff}; contrasts={pilot_contrast_diff}")

    bin_prefix, convergence_metrics, contrast_prefix = summarize_prefixes(outputs, contrasts, masks)
    final_ev = bin_prefix[bin_prefix.iterations.eq(10000)].copy()
    final_contrasts = contrast_prefix[contrast_prefix.iterations.eq(10000)].copy()
    final_ev.to_csv(RESULT_DIR / "Layer2_formal_50kyr_Global_EV_summary.csv", index=False)
    final_contrasts.to_csv(RESULT_DIR / "Layer2_formal_phase_contrasts_50kyr.csv", index=False)
    coverage.to_csv(RESULT_DIR / "Layer2_formal_perturbation_coverage_50kyr.csv", index=False)
    bin_prefix.to_csv(RESULT_DIR / "Layer2_formal_convergence_bin_prefixes.csv", index=False)
    convergence_metrics.to_csv(RESULT_DIR / "Layer2_formal_convergence_metrics.csv", index=False)
    contrast_prefix.to_csv(RESULT_DIR / "Layer2_formal_convergence_phase_contrasts.csv", index=False)

    audit = build_execution_audit(raw, params, layer1_maps, layer1_notes, tier2_audit)
    participation = structure_participation(audit)
    audit.to_csv(RESULT_DIR / "Layer2_formal_record_execution_audit.csv", index=False)
    participation.to_csv(RESULT_DIR / "Layer2_formal_structure_participation_summary.csv", index=False)

    protected_after = {str(p): sha256(p) for p in pilot.PROTECTED_PATHS if p.exists()}
    protected_check = pd.DataFrame([
        {"path": p, "sha256_before": h, "sha256_after": protected_after.get(p), "unchanged": h == protected_after.get(p)}
        for p, h in protected_before.items()
    ])
    protected_check.to_csv(RESULT_DIR / "protected_files_hash_check.csv", index=False)

    width = final_ev.envelope_width.to_numpy(float)
    diff = np.abs(final_ev.observed_EV - final_ev.MC_median)
    metadata = {
        "analysis": "Formal Layer 2 conservative method-based measurement sensitivity — standard 50-kyr Global EV only",
        "iterations": ITERATIONS,
        "random_seed": SEED,
        "batch_size": BATCH_SIZE,
        "chronology_perturbed": False,
        "other_bin_widths_analyzed": False,
        "retained_operational_records": int(raw.record_id.nunique()),
        "retained_measurements": int(len(raw)),
        "actually_active_records": int(audit.actually_active.sum()),
        "actually_perturbed_measurements": len(active_rows),
        "perturbed_measurement_percent": 100 * len(active_rows) / len(raw),
        "documented_active_records": int((audit.primary_structure.str.startswith("Tier 1") & audit.actually_active).sum()),
        "documented_perturbed_measurements": int(audit.loc[audit.primary_structure.str.startswith("Tier 1"), "actually_perturbed_measurements"].sum()),
        "tier2_gaussian_active_records": int((audit.primary_structure.str.contains("Conservative Gaussian") & audit.actually_active).sum()),
        "tier2_gaussian_perturbed_measurements": int(audit.loc[audit.primary_structure.str.contains("Conservative Gaussian"), "actually_perturbed_measurements"].sum()),
        "tier2_bounded_active_records": int((audit.primary_structure.str.contains("Bounded uniform") & audit.actually_active).sum()),
        "tier2_bounded_perturbed_measurements": int(audit.loc[audit.primary_structure.str.contains("Bounded uniform"), "actually_perturbed_measurements"].sum()),
        "tier2_correlated_active_records": int((audit.primary_structure.str.contains("Record-wide correlated") & audit.actually_active).sum()),
        "tier2_correlated_perturbed_measurements": int(audit.loc[audit.primary_structure.str.contains("Record-wide correlated"), "actually_perturbed_measurements"].sum()),
        "s1_50kyr_max_absolute_validation_difference": s1_diff,
        "formal_prefix_vs_validated_pilot_max_EV_difference": pilot_ev_diff,
        "formal_prefix_vs_validated_pilot_max_contrast_difference": pilot_contrast_diff,
        "observed_EV_vs_MC_median_max_absolute_difference": float(np.max(diff)),
        "max_difference_midpoint_Ma": float(final_ev.iloc[int(np.argmax(diff))].midpoint_Ma),
        "50kyr_envelope_width_median": float(np.median(width)),
        "50kyr_envelope_width_maximum": float(np.max(width)),
        "50kyr_maximum_width_midpoint_Ma": float(final_ev.iloc[int(np.argmax(width))].midpoint_Ma),
        "perturbation_coverage_minimum": float(coverage.perturbation_coverage_proportion.min()),
        "perturbation_coverage_median": float(coverage.perturbation_coverage_proportion.median()),
        "perturbation_coverage_maximum": float(coverage.perturbation_coverage_proportion.max()),
        "protected_dataset_and_document_hashes_unchanged": bool(protected_check.unchanged.all()),
        "anomaly_diagnostics": anomaly,
        "implementation_notes": [
            "The formal run uses the same seed and RNG call order as the validated pilot; the first 1,000 EV trajectories and phase contrasts were checked directly against the pilot realization file.",
            "Only the standard 50-kyr Global EV reconstruction was analyzed. No 100-, 150-, or 200-kyr uncertainty propagation was run.",
            "Layer 1A documented mappings, Tier 2 Gaussian and bounded observation-level structures, fixed chronology, and the observed non-zero eligible record-bin mask were unchanged.",
            "The seven record-wide correlated calibration/reconstruction records remained fixed because the executable table does not provide the raw calibration variables/equations required for a non-cancelling realization.",
        ],
    }
    (RESULT_DIR / "formal_metadata.json").write_text(json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8")

    manifest = []
    for p in [RAW_PATH, IQR_PATH, S1_PATH, PARAM_PATH, PILOT_SCRIPT, PILOT_NPZ, SCRIPT]:
        manifest.append({"file": p.name, "sha256": sha256(p), "source": str(p)})
    for name, p in sources.items():
        manifest.append({"file": name, "sha256": sha256(p), "source": str(p)})
    pd.DataFrame(manifest).drop_duplicates("sha256").to_csv(RESULT_DIR / "input_source_hash_manifest.csv", index=False)

    np.savez_compressed(
        RESULT_DIR / "Layer2_formal_Global_EV_realizations_50kyr_10000.npz",
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
