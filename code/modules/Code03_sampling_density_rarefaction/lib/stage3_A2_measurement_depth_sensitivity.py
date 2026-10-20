#!/usr/bin/env python3
"""Measurement-depth sensitivity for the existing Stage 3 A2 design.

This script adds m = 3, 4 and 5 measurements per record-bin while keeping
N = 5 operational records per 50-kyr bin, 10,000 iterations, seed 20260821,
the formal non-zero observed-IQR eligibility rule, phase definitions and
quantile method used by the existing Stage 3 A2 analysis.

It does not modify the formal Stage 3 A2 output, manuscript or figures.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd


WORK = Path(__file__).resolve().parents[1]
BASE_SCRIPT = Path(__file__).with_name("stage3_sampling_density_robustness.py")
RESULTS = WORK / "outputs" / "generated" / "A2_measurement_depth_sensitivity"
FIGURES = WORK / "outputs" / "generated" / "figures"

ITERATIONS = 10_000
SEED = 20260821
N_RECORDS = 5
M_VALUES = (3, 4, 5)
CONTRASTS = ("phase1_minus_pre", "phase2_minus_pre", "recovery_minus_phase2")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def load_stage3_module():
    spec = importlib.util.spec_from_file_location("stage3_baseline", BASE_SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def run_design(stage3, m_depth: int, contributing: pd.DataFrame,
               measurements: pd.DataFrame, inventory: pd.DataFrame):
    ages = inventory["bin_center_Ma"].to_numpy(float)
    meas = measurements.copy()
    meas["bin_center_Ma"] = meas["bin_center_Ma_audited"]

    eligible_rows = contributing.loc[contributing["n_measurements"].ge(m_depth)].copy()
    values: dict[tuple[float, str], np.ndarray] = {}
    for (age, record), group in meas.groupby(["bin_center_Ma", "record_id"]):
        values[(float(age), record)] = group["z_score_recomputed"].to_numpy(float)

    eligible_by_bin = {
        float(age): group["record_id"].tolist()
        for age, group in eligible_rows.groupby("bin_center_Ma")
    }
    eligible_counts = {
        float(age): len(eligible_by_bin.get(float(age), [])) for age in ages
    }

    # The established A2 seed is reset for each m design so all settings use
    # exactly the same fixed-seed rule and m=4 reproduces the approved result.
    rng = np.random.default_rng(SEED)
    matrix = np.full((ITERATIONS, len(ages)), np.nan, dtype=float)
    iteration_index = np.arange(ITERATIONS)[:, None]

    for j, age in enumerate(ages):
        records = eligible_by_bin.get(float(age), [])
        if len(records) < N_RECORDS:
            continue

        record_subsample_iqrs = np.empty((ITERATIONS, len(records)), dtype=float)
        for r, record in enumerate(records):
            z = values[(float(age), record)]
            if len(z) < m_depth:
                raise AssertionError("Eligibility/value mismatch")
            priorities = rng.random((ITERATIONS, len(z)))
            measurement_indices = np.argpartition(
                priorities, m_depth - 1, axis=1
            )[:, :m_depth]
            chosen = z[measurement_indices]
            q1 = np.quantile(chosen, 0.25, axis=1, method="linear")
            q3 = np.quantile(chosen, 0.75, axis=1, method="linear")
            record_subsample_iqrs[:, r] = q3 - q1

        record_priorities = rng.random((ITERATIONS, len(records)))
        selected_record_indices = np.argpartition(
            record_priorities, N_RECORDS - 1, axis=1
        )[:, :N_RECORDS]
        matrix[:, j] = record_subsample_iqrs[
            iteration_index, selected_record_indices
        ].mean(axis=1)

    eligible = np.array([eligible_counts[float(age)] >= N_RECORDS for age in ages])
    summary = pd.DataFrame({
        "m_measurements": m_depth,
        "N_records": N_RECORDS,
        "seed": SEED,
        "iterations": ITERATIONS,
        "bin_center_Ma": ages,
        "eligible_record_count_n_ge_m": [eligible_counts[float(age)] for age in ages],
        "eligible": eligible,
        "observed_global_EV": inventory["observed_global_EV"].to_numpy(float),
        "median_rarefied_EV": np.nanmedian(matrix, axis=0),
        "p2_5": np.nanquantile(matrix, 0.025, axis=0),
        "p97_5": np.nanquantile(matrix, 0.975, axis=0),
    })

    features = stage3.phase_features(matrix, ages)
    contrast_rows = []
    for contrast in CONTRASTS:
        values_contrast = features[contrast].dropna().to_numpy(float)
        contrast_rows.append({
            "m_measurements": m_depth,
            "N_records": N_RECORDS,
            "seed": SEED,
            "iterations": ITERATIONS,
            "contrast": contrast,
            "median": float(np.median(values_contrast)),
            "p2_5": float(np.quantile(values_contrast, 0.025)),
            "p97_5": float(np.quantile(values_contrast, 0.975)),
            "probability_gt_zero": float(np.mean(values_contrast > 0)),
            "valid_iterations": len(values_contrast),
        })

    phase_rows = []
    for phase, (lo, hi) in stage3.WINDOWS.items():
        in_phase = (ages >= lo - 1e-9) & (ages <= hi + 1e-9)
        phase_rows.append({
            "m_measurements": m_depth,
            "N_records": N_RECORDS,
            "phase": phase,
            "phase_start_Ma": hi,
            "phase_end_Ma": lo,
            "available_bins": int(in_phase.sum()),
            "eligible_bins": int((eligible & in_phase).sum()),
        })

    convergence_rows = []
    for prefix in (1_000, 5_000):
        sub = matrix[:prefix]
        final = matrix
        convergence_rows.append({
            "m_measurements": m_depth,
            "N_records": N_RECORDS,
            "prefix_iterations": prefix,
            "max_abs_median_change_vs_10000": float(np.nanmax(
                np.abs(np.nanmedian(sub, axis=0) - np.nanmedian(final, axis=0))
            )),
            "max_abs_p2_5_change_vs_10000": float(np.nanmax(
                np.abs(np.nanquantile(sub, 0.025, axis=0) - np.nanquantile(final, 0.025, axis=0))
            )),
            "max_abs_p97_5_change_vs_10000": float(np.nanmax(
                np.abs(np.nanquantile(sub, 0.975, axis=0) - np.nanquantile(final, 0.975, axis=0))
            )),
        })

    return (
        summary,
        pd.DataFrame(contrast_rows),
        pd.DataFrame(phase_rows),
        pd.DataFrame(convergence_rows),
        matrix,
    )


def main():
    stage3 = load_stage3_module()
    standard, contributing, zeros, measurements = stage3.load_current_baseline()
    inventory = stage3.build_inventory(contributing, measurements)

    bin_outputs = []
    contrast_outputs = []
    phase_outputs = []
    convergence_outputs = []
    matrices = {}
    for m_depth in M_VALUES:
        bin_summary, contrasts, phase, convergence, matrix = run_design(
            stage3, m_depth, contributing, measurements, inventory
        )
        bin_outputs.append(bin_summary)
        contrast_outputs.append(contrasts)
        phase_outputs.append(phase)
        convergence_outputs.append(convergence)
        matrices[m_depth] = matrix

    bin_summary = pd.concat(bin_outputs, ignore_index=True)
    contrast_summary = pd.concat(contrast_outputs, ignore_index=True)
    phase_eligibility = pd.concat(phase_outputs, ignore_index=True)
    convergence_summary = pd.concat(convergence_outputs, ignore_index=True)

    # Mandatory checks: m=4 reproduces the approved A2 result exactly within
    # numerical tolerance; m=5 retains 47/49 total and 4/6 pre-escalation bins.
    approved = pd.read_csv(
        WORK / "03_results" / "stage3_sampling_density" /
        "A2_joint_measurement_record_rarefaction_bin_summary.csv"
    ).sort_values("bin_center_Ma")
    current_m4 = bin_summary.loc[bin_summary["m_measurements"].eq(4)].sort_values("bin_center_Ma")
    for column in ("median_rarefied_EV", "p2_5", "p97_5"):
        error = np.nanmax(np.abs(
            current_m4[column].to_numpy(float) - approved[column].to_numpy(float)
        ))
        if error > 1e-12:
            raise AssertionError(f"m=4 does not reproduce approved A2 {column}: {error}")

    m5 = bin_summary.loc[bin_summary["m_measurements"].eq(5)]
    if int(m5["eligible"].sum()) != 47:
        raise AssertionError("m=5 does not retain 47/49 bins")
    m5_pre = phase_eligibility.loc[
        phase_eligibility["m_measurements"].eq(5)
        & phase_eligibility["phase"].eq("pre_escalation")
    ].iloc[0]
    if int(m5_pre["eligible_bins"]) != 4 or int(m5_pre["available_bins"]) != 6:
        raise AssertionError("m=5 pre-escalation eligibility is not 4/6")

    bin_summary.to_csv(RESULTS / "A2_measurement_depth_sensitivity_bin_summary.csv", index=False)
    contrast_summary.to_csv(RESULTS / "A2_measurement_depth_sensitivity_contrast_summary.csv", index=False)
    phase_eligibility.to_csv(RESULTS / "A2_measurement_depth_sensitivity_phase_eligibility.csv", index=False)
    convergence_summary.to_csv(RESULTS / "A2_measurement_depth_sensitivity_convergence.csv", index=False)

    provenance = {
        "scope": "standard 50-kyr bins, 67.2-64.8 Ma",
        "design": "joint measurement-depth and record-number rarefaction",
        "m_measurements": list(M_VALUES),
        "N_records": N_RECORDS,
        "iterations": ITERATIONS,
        "seed_rule": "reset NumPy default_rng(20260821) independently for each m design",
        "seed": SEED,
        "record_sampling": "without replacement within each temporal bin",
        "measurement_sampling": "without replacement within each record-bin",
        "standardization": "existing record-wide z-scores; no restandardization of subsamples",
        "IQR_quantile": "linear interpolation; equivalent to R type 7",
        "eligibility": "observed record-bin IQR > 0, n_measurements >= m, at least 5 eligible records",
        "phase_definitions": stage3.WINDOWS,
        "source_script": str(BASE_SCRIPT),
        "source_script_sha256": sha256(BASE_SCRIPT),
        "record_iqr_source": str(stage3.IQR_SOURCE),
        "record_iqr_source_sha256": sha256(stage3.IQR_SOURCE),
        "measurement_source": str(stage3.MEAS_SOURCE),
        "measurement_source_sha256": sha256(stage3.MEAS_SOURCE),
        "m4_reproduces_approved_A2": True,
        "m5_total_eligible_bins": 47,
        "m5_pre_escalation_eligible_bins": "4/6",
        "manuscript_modified": False,
        "formal_Fig_A2_modified": False,
    }
    (RESULTS / "A2_measurement_depth_sensitivity_provenance.json").write_text(
        json.dumps(provenance, indent=2, ensure_ascii=False) + "\n"
    )

    print(contrast_summary.to_string(index=False))
    print("\nEligibility")
    print(phase_eligibility.to_string(index=False))


if __name__ == "__main__":
    raise SystemExit("Preserved numerical library: run ../run.py for the packaged final workflow.")
