#!/usr/bin/env python3
"""Stage 3 sampling-density robustness and proxy-composition diagnostics.

Scope is limited to the standard 50-kyr, 67.2--64.8 Ma EV reconstruction.
The formal observed baseline is the arithmetic mean of non-zero contributing
record-level IQRs. Wider bins, EV--biotic analyses, Monte Carlo propagation,
and composition-controlled estimators are intentionally outside this script.
"""

from __future__ import annotations

import hashlib
import json
import math
import unicodedata
from pathlib import Path

import numpy as np
import pandas as pd


WORK = Path(__file__).resolve().parents[1]
DATA = WORK / "inputs"
RESULTS = WORK / "outputs" / "generated"
FIGURES = WORK / "outputs" / "generated" / "figures"
LOGS = WORK / "outputs" / "generated"

IQR_SOURCE = DATA / "record_level_IQR_all_bin_widths_revised_v4_n2.csv"
MEAS_SOURCE = DATA / "retained_measurements_with_zscores_v4.csv"
SITE_MAP_SOURCE = DATA / "canonical_site_mapping_v5.csv"

ITERATIONS = 10_000
SEEDS_A1 = {3: 20260824, 5: 20260826, 7: 20260828}
SEED_A2 = 20260821
SEED_B_BASE = 20260900
BLOCK_LENGTHS = (3, 4, 5)
PRIMARY_BLOCK = 4
CATEGORIES = ["temperature", "d18O", "d13C", "weathering", "pCO2"]
CATEGORY_LABELS = {
    "temperature": "Temperature",
    "d18O": "δ18O",
    "d13C": "δ13C",
    "weathering": "Weathering",
    "pCO2": "pCO₂",
}
WINDOWS = {
    "pre_escalation": (66.95, 67.20),
    "phase_1": (66.45, 66.90),
    "phase_2": (66.05, 66.40),
    "recovery": (64.80, 65.95),
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def center50(age: float) -> float | None:
    for i in range(49):
        lo = 64.775 + i * 0.05
        hi = lo + 0.05
        if lo <= age < hi:
            return round(lo + 0.025, 6)
    return None


def quantile_iqr(values: np.ndarray) -> float:
    q1, q3 = np.quantile(values, [0.25, 0.75], method="linear")
    return float(q3 - q1)


def load_current_baseline() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    iqrs_all = pd.read_csv(IQR_SOURCE)
    standard = iqrs_all.loc[
        (iqrs_all["scheme"].astype(str) == "50")
        & (iqrs_all["alignment"] == "standard")
    ].copy()
    standard["bin_center_Ma"] = standard["bin_center_Ma"].round(6)

    zeros = standard.loc[standard["IQR"].eq(0)].copy()
    contributing = standard.loc[standard["IQR"].gt(0)].copy()
    if len(standard) != 912 or len(zeros) != 3 or len(contributing) != 909:
        raise AssertionError(
            f"Unexpected 50-kyr inventory: all={len(standard)}, zero={len(zeros)}, "
            f"contributing={len(contributing)}"
        )

    measurements = pd.read_csv(MEAS_SOURCE, low_memory=False)
    # Use the bin assignments supplied in the retained-measurement table.
    measurements["bin_center_Ma_audited"] = pd.to_numeric(measurements["bin_center_Ma"], errors="coerce")

    mapping = pd.read_csv(SITE_MAP_SOURCE)
    mapping = mapping.rename(columns={"raw_site": "site", "canonical_site_v5": "canonical_site_audited"})
    mapping = mapping[["record_id", "site", "canonical_site_audited"]].drop_duplicates()
    for frame in (measurements, mapping):
        frame["record_id_key"] = frame["record_id"].astype(str).map(lambda x: unicodedata.normalize("NFC", x).strip())
        frame["site_key"] = frame["site"].astype(str).map(lambda x: unicodedata.normalize("NFC", x).strip())
    measurements = measurements.merge(
        mapping[["record_id_key", "site_key", "canonical_site_audited"]].drop_duplicates(),
        on=["record_id_key", "site_key"], how="left", validate="many_to_one"
    )
    missing = measurements["canonical_site_audited"].isna()
    if missing.any():
        examples = measurements.loc[missing, ["record_id", "site"]].drop_duplicates().head(20)
        raise AssertionError(f"Canonical-site mapping incomplete:\n{examples.to_string(index=False)}")

    domain_measurements = measurements.loc[measurements["bin_center_Ma_audited"].notna()].copy()
    domain_measurements["bin_center_Ma_audited"] = domain_measurements["bin_center_Ma_audited"].astype(float)

    # Recompute every record-bin directly from retained sample-level z-scores.
    recalculated = (
        domain_measurements.groupby(["record_id", "category", "bin_center_Ma_audited"], as_index=False)
        .agg(
            n_recomputed=("z_score_recomputed", "size"),
            IQR_recomputed=("z_score_recomputed", lambda x: quantile_iqr(x.to_numpy(float))),
        )
        .rename(columns={"bin_center_Ma_audited": "bin_center_Ma"})
    )
    matched = standard.merge(
        recalculated,
        on=["record_id", "category", "bin_center_Ma"],
        how="left",
        validate="one_to_one",
    )
    if matched["IQR_recomputed"].isna().any():
        raise AssertionError("At least one current record-bin cannot be traced to sample-level z-scores")
    max_iqr_error = float((matched["IQR"] - matched["IQR_recomputed"]).abs().max())
    max_n_error = int((matched["n_measurements"] - matched["n_recomputed"]).abs().max())
    if max_iqr_error > 5e-10 or max_n_error != 0:
        raise AssertionError(f"Bottom-up IQR mismatch: max IQR error={max_iqr_error}, max n error={max_n_error}")

    return standard, contributing, zeros, domain_measurements


def build_inventory(contributing: pd.DataFrame, measurements: pd.DataFrame) -> pd.DataFrame:
    bins = pd.DataFrame({"bin_center_Ma": np.round(np.arange(64.8, 67.2001, 0.05), 6)})
    pairs = contributing[["record_id", "category", "bin_center_Ma"]].drop_duplicates()
    meas = measurements.copy()
    meas["bin_center_Ma"] = meas["bin_center_Ma_audited"]

    contributing_meas = meas.merge(pairs, on=["record_id", "category", "bin_center_Ma"], how="inner")

    inventory = bins.copy()
    inventory = inventory.merge(
        contributing.groupby("bin_center_Ma").size().rename("contributing_record_count").reset_index(),
        on="bin_center_Ma", how="left"
    )
    inventory = inventory.merge(
        contributing_meas.groupby("bin_center_Ma").size().rename("raw_measurement_count_contributing_record_bins").reset_index(),
        on="bin_center_Ma", how="left"
    )
    inventory = inventory.merge(
        meas.groupby("bin_center_Ma").size().rename("raw_measurement_count_all_retained").reset_index(),
        on="bin_center_Ma", how="left"
    )
    inventory = inventory.merge(
        contributing_meas.groupby("bin_center_Ma")["canonical_site_audited"].nunique().rename("contributing_canonical_site_count").reset_index(),
        on="bin_center_Ma", how="left"
    )
    inventory = inventory.merge(
        contributing.groupby("bin_center_Ma")["category"].nunique().rename("proxy_category_count").reset_index(),
        on="bin_center_Ma", how="left"
    )
    inventory = inventory.merge(
        contributing.groupby("bin_center_Ma")["IQR"].mean().rename("observed_global_EV").reset_index(),
        on="bin_center_Ma", how="left"
    )
    for category in CATEGORIES:
        counts = (
            contributing.loc[contributing["category"].eq(category)]
            .groupby("bin_center_Ma").size().rename(f"records_{category}").reset_index()
        )
        inventory = inventory.merge(counts, on="bin_center_Ma", how="left")
        inventory[f"records_{category}"] = inventory[f"records_{category}"].fillna(0).astype(int)
        inventory[f"proportion_{category}"] = (
            inventory[f"records_{category}"] / inventory["contributing_record_count"]
        )
        inventory[f"absent_{category}"] = inventory[f"records_{category}"].eq(0)

    count_cols = [
        "contributing_record_count", "raw_measurement_count_contributing_record_bins",
        "raw_measurement_count_all_retained", "contributing_canonical_site_count", "proxy_category_count",
    ]
    inventory[count_cols] = inventory[count_cols].fillna(0).astype(int)
    if inventory["observed_global_EV"].isna().any():
        raise AssertionError("Observed Global EV is missing in at least one of 49 bins")
    return inventory


def phase_features(matrix: np.ndarray, ages: np.ndarray) -> pd.DataFrame:
    out: dict[str, np.ndarray] = {}
    for name, (lo, hi) in WINDOWS.items():
        mask = (ages >= lo - 1e-9) & (ages <= hi + 1e-9)
        if not mask.any():
            raise AssertionError(f"No bins in feature window {name}")
        out[f"{name}_mean"] = np.nanmean(matrix[:, mask], axis=1)
        out[f"{name}_max"] = np.nanmax(matrix[:, mask], axis=1)
    out["phase1_minus_pre"] = out["phase_1_mean"] - out["pre_escalation_mean"]
    out["phase2_minus_pre"] = out["phase_2_mean"] - out["pre_escalation_mean"]
    out["recovery_minus_phase2"] = out["recovery_mean"] - out["phase_2_mean"]
    safe = np.where(np.isnan(matrix), -np.inf, matrix)
    max_index = np.argmax(safe, axis=1)
    out["peak_magnitude"] = matrix[np.arange(matrix.shape[0]), max_index]
    out["peak_age_Ma"] = ages[max_index]
    return pd.DataFrame(out)


def summarise_feature_distribution(features: pd.DataFrame, test: str, design: str) -> pd.DataFrame:
    rows = []
    for column in features.columns:
        x = features[column].dropna().to_numpy(float)
        rows.append({
            "test": test,
            "design": design,
            "feature": column,
            "median": float(np.median(x)),
            "p2_5": float(np.quantile(x, 0.025)),
            "p97_5": float(np.quantile(x, 0.975)),
            "n_iterations_valid": len(x),
            "probability_gt_zero": float(np.mean(x > 0)) if "minus" in column else np.nan,
        })
    return pd.DataFrame(rows)


def convergence_rows(matrix: np.ndarray, ages: np.ndarray, test: str, design: str) -> pd.DataFrame:
    finals = {}
    for n in (1_000, 5_000, 10_000):
        sub = matrix[:n]
        finals[n] = {
            "median": np.nanmedian(sub, axis=0),
            "lo": np.nanquantile(sub, 0.025, axis=0),
            "hi": np.nanquantile(sub, 0.975, axis=0),
        }
    rows = []
    for n in (1_000, 5_000):
        rows.append({
            "test": test,
            "design": design,
            "prefix_iterations": n,
            "max_abs_median_change_vs_10000": float(np.nanmax(np.abs(finals[n]["median"] - finals[10_000]["median"]))),
            "max_abs_p2_5_change_vs_10000": float(np.nanmax(np.abs(finals[n]["lo"] - finals[10_000]["lo"]))),
            "max_abs_p97_5_change_vs_10000": float(np.nanmax(np.abs(finals[n]["hi"] - finals[10_000]["hi"]))),
            "median_abs_median_change_vs_10000": float(np.nanmedian(np.abs(finals[n]["median"] - finals[10_000]["median"]))),
        })
    return pd.DataFrame(rows)


def record_rarefaction(contributing: pd.DataFrame, inventory: pd.DataFrame):
    ages = inventory["bin_center_Ma"].to_numpy(float)
    summaries, feature_summaries, convergence, feature_iterations = [], [], [], []
    matrices: dict[int, np.ndarray] = {}
    by_bin = {age: g["IQR"].to_numpy(float) for age, g in contributing.groupby("bin_center_Ma")}
    for n_draw in (3, 5, 7):
        rng = np.random.default_rng(SEEDS_A1[n_draw])
        matrix = np.full((ITERATIONS, len(ages)), np.nan)
        for j, age in enumerate(ages):
            values = by_bin[age]
            if len(values) < n_draw:
                continue
            # Draw without replacement by selecting the N smallest random priorities.
            priorities = rng.random((ITERATIONS, len(values)))
            indices = np.argpartition(priorities, n_draw - 1, axis=1)[:, :n_draw]
            matrix[:, j] = values[indices].mean(axis=1)
        matrices[n_draw] = matrix

        med = np.nanmedian(matrix, axis=0)
        lo = np.nanquantile(matrix, 0.025, axis=0)
        hi = np.nanquantile(matrix, 0.975, axis=0)
        summaries.append(pd.DataFrame({
            "test": "A1_fixed_record_number",
            "N_records": n_draw,
            "seed": SEEDS_A1[n_draw],
            "iterations": ITERATIONS,
            "bin_center_Ma": ages,
            "available_record_count": inventory["contributing_record_count"].to_numpy(int),
            "eligible": inventory["contributing_record_count"].to_numpy(int) >= n_draw,
            "observed_global_EV": inventory["observed_global_EV"].to_numpy(float),
            "median_rarefied_EV": med,
            "p2_5": lo,
            "p97_5": hi,
        }))
        features = phase_features(matrix, ages)
        features.insert(0, "iteration", np.arange(1, ITERATIONS + 1))
        features.insert(0, "design", f"N={n_draw}")
        features.insert(0, "test", "A1_fixed_record_number")
        feature_iterations.append(features)
        feature_summaries.append(summarise_feature_distribution(features.drop(columns=["test", "design", "iteration"]), "A1_fixed_record_number", f"N={n_draw}"))
        convergence.append(convergence_rows(matrix, ages, "A1_fixed_record_number", f"N={n_draw}"))
    return (
        pd.concat(summaries, ignore_index=True),
        pd.concat(feature_summaries, ignore_index=True),
        pd.concat(convergence, ignore_index=True),
        pd.concat(feature_iterations, ignore_index=True),
        matrices,
    )


def hierarchical_rarefaction(contributing: pd.DataFrame, measurements: pd.DataFrame, inventory: pd.DataFrame):
    ages = inventory["bin_center_Ma"].to_numpy(float)
    meas = measurements.copy()
    meas["bin_center_Ma"] = meas["bin_center_Ma_audited"]
    eligible_rows = contributing.loc[contributing["n_measurements"].ge(4)].copy()
    values: dict[tuple[float, str], np.ndarray] = {}
    for (age, record), g in meas.groupby(["bin_center_Ma", "record_id"]):
        values[(float(age), record)] = g["z_score_recomputed"].to_numpy(float)
    eligible_by_bin = {
        float(age): g["record_id"].tolist() for age, g in eligible_rows.groupby("bin_center_Ma")
    }
    counts = {age: len(eligible_by_bin.get(float(age), [])) for age in ages}
    if min(counts.values()) < 5:
        raise AssertionError(f"A2 has an ineligible bin; minimum eligible records={min(counts.values())}")

    rng = np.random.default_rng(SEED_A2)
    matrix = np.empty((ITERATIONS, len(ages)), dtype=float)
    iteration_index = np.arange(ITERATIONS)[:, None]
    for j, age in enumerate(ages):
        records = eligible_by_bin[float(age)]
        record_subsample_iqrs = np.empty((ITERATIONS, len(records)), dtype=float)
        for r, record in enumerate(records):
            z = values[(float(age), record)]
            priorities = rng.random((ITERATIONS, len(z)))
            measurement_indices = np.argpartition(priorities, 3, axis=1)[:, :4]
            chosen = z[measurement_indices]
            q1 = np.quantile(chosen, 0.25, axis=1, method="linear")
            q3 = np.quantile(chosen, 0.75, axis=1, method="linear")
            record_subsample_iqrs[:, r] = q3 - q1
        record_priorities = rng.random((ITERATIONS, len(records)))
        selected_record_indices = np.argpartition(record_priorities, 4, axis=1)[:, :5]
        matrix[:, j] = record_subsample_iqrs[iteration_index, selected_record_indices].mean(axis=1)

    summary = pd.DataFrame({
        "test": "A2_measurement_and_record_rarefaction",
        "m_measurements": 4,
        "N_records": 5,
        "seed": SEED_A2,
        "iterations": ITERATIONS,
        "bin_center_Ma": ages,
        "eligible_record_count_n_ge_4": [counts[float(age)] for age in ages],
        "eligible": True,
        "observed_global_EV": inventory["observed_global_EV"].to_numpy(float),
        "median_rarefied_EV": np.median(matrix, axis=0),
        "p2_5": np.quantile(matrix, 0.025, axis=0),
        "p97_5": np.quantile(matrix, 0.975, axis=0),
        "zero_EV_outcome_count": np.sum(matrix == 0, axis=0),
    })
    features = phase_features(matrix, ages)
    feature_iterations = features.copy()
    feature_iterations.insert(0, "iteration", np.arange(1, ITERATIONS + 1))
    feature_iterations.insert(0, "design", "m=4,N=5")
    feature_iterations.insert(0, "test", "A2_measurement_and_record_rarefaction")
    feature_summary = summarise_feature_distribution(features, "A2_measurement_and_record_rarefaction", "m=4,N=5")
    convergence = convergence_rows(matrix, ages, "A2_measurement_and_record_rarefaction", "m=4,N=5")
    return summary, feature_summary, convergence, feature_iterations, matrix


def correlation(x: np.ndarray, y: np.ndarray, method: str) -> float:
    if method == "Pearson":
        return float(np.corrcoef(x, y)[0, 1])
    rx = pd.Series(x).rank(method="average").to_numpy(float)
    ry = pd.Series(y).rank(method="average").to_numpy(float)
    return float(np.corrcoef(rx, ry)[0, 1])


def circular_block_indices(n: int, block: int, rng: np.random.Generator) -> np.ndarray:
    n_blocks = math.ceil(n / block)
    starts = rng.integers(0, n, size=n_blocks)
    indices = np.concatenate([(start + np.arange(block)) % n for start in starts])
    return indices[:n]


def moving_block_ci(x: np.ndarray, y: np.ndarray, method: str, block: int, seed: int) -> tuple[float, float, int]:
    rng = np.random.default_rng(seed)
    n = len(x)
    n_blocks = math.ceil(n / block)
    starts = rng.integers(0, n, size=(ITERATIONS, n_blocks))
    idx = ((starts[:, :, None] + np.arange(block)[None, None, :]) % n).reshape(ITERATIONS, -1)[:, :n]
    bx = x[idx]
    by = y[idx]
    if method == "Spearman":
        bx = pd.DataFrame(bx).rank(axis=1, method="average").to_numpy(float)
        by = pd.DataFrame(by).rank(axis=1, method="average").to_numpy(float)
    bx = bx - bx.mean(axis=1, keepdims=True)
    by = by - by.mean(axis=1, keepdims=True)
    denom = np.sqrt(np.sum(bx * bx, axis=1) * np.sum(by * by, axis=1))
    estimates = np.sum(bx * by, axis=1) / denom
    estimates = estimates[np.isfinite(estimates)]
    valid = len(estimates)
    if valid < int(ITERATIONS * 0.99):
        raise AssertionError(f"Too few valid block-bootstrap correlations: {valid}/{ITERATIONS}")
    return float(np.quantile(estimates, 0.025)), float(np.quantile(estimates, 0.975)), valid


def density_associations(inventory: pd.DataFrame) -> pd.DataFrame:
    chronological = inventory.sort_values("bin_center_Ma", ascending=False).reset_index(drop=True)
    predictors = {
        "contributing_operational_record_count": "contributing_record_count",
        "raw_measurement_count_contributing_record_bins_primary": "raw_measurement_count_contributing_record_bins",
        "raw_measurement_count_all_retained_secondary": "raw_measurement_count_all_retained",
        "contributing_canonical_site_count": "contributing_canonical_site_count",
    }
    rows = []
    for pi, (label, column) in enumerate(predictors.items()):
        raw_x = chronological[column].to_numpy(float)
        raw_y = chronological["observed_global_EV"].to_numpy(float)
        delta_x = np.diff(raw_x)
        delta_y = np.diff(raw_y)
        for scale, x, y in (("raw_bins", raw_x, raw_y), ("first_difference_older_to_younger", delta_x, delta_y)):
            for mi, method in enumerate(("Pearson", "Spearman")):
                estimate = correlation(x, y, method)
                for block in BLOCK_LENGTHS:
                    seed = SEED_B_BASE + pi * 100 + (0 if scale == "raw_bins" else 40) + mi * 10 + block
                    lo, hi, valid = moving_block_ci(x, y, method, block, seed)
                    rows.append({
                        "predictor": label,
                        "source_column": column,
                        "scale": scale,
                        "method": method,
                        "estimate": estimate,
                        "block_length_bins": block,
                        "primary_block_length": block == PRIMARY_BLOCK,
                        "CI_2_5": lo,
                        "CI_97_5": hi,
                        "bootstrap_iterations_requested": ITERATIONS,
                        "bootstrap_iterations_valid": valid,
                        "seed": seed,
                        "n_pairs": len(x),
                    })
    return pd.DataFrame(rows)


def composition_outputs(inventory: pd.DataFrame):
    chronological = inventory.sort_values("bin_center_Ma", ascending=False).reset_index(drop=True)
    phase_rows = []
    for phase, (lo, hi) in WINDOWS.items():
        sub = chronological.loc[chronological["bin_center_Ma"].between(lo, hi)]
        for category in CATEGORIES:
            phase_rows.append({
                "phase": phase,
                "category": category,
                "mean_record_count": sub[f"records_{category}"].mean(),
                "mean_record_proportion": sub[f"proportion_{category}"].mean(),
                "bins_absent": int(sub[f"absent_{category}"].sum()),
                "n_bins": len(sub),
            })
    phase = pd.DataFrame(phase_rows)

    adjacent = chronological[["bin_center_Ma"]].copy()
    adjacent["younger_bin_center_Ma"] = chronological["bin_center_Ma"].shift(-1)
    total_variation = np.zeros(len(chronological))
    for category in CATEGORIES:
        delta = chronological[f"proportion_{category}"].shift(-1) - chronological[f"proportion_{category}"]
        adjacent[f"delta_proportion_{category}_older_to_younger"] = delta
        total_variation += delta.fillna(0).abs().to_numpy()
    adjacent["composition_total_variation_distance"] = total_variation / 2
    adjacent = adjacent.iloc[:-1].copy()

    correlations = []
    for category in CATEGORIES:
        p = chronological[f"proportion_{category}"].to_numpy(float)
        ev = chronological["observed_global_EV"].to_numpy(float)
        dp = np.diff(p)
        dev = np.diff(ev)
        for scale, x, y in (("raw_bins", p, ev), ("first_difference_older_to_younger", dp, dev)):
            for method in ("Pearson", "Spearman"):
                correlations.append({
                    "category": category,
                    "scale": scale,
                    "method": method,
                    "estimate": correlation(x, y, method),
                    "n_pairs": len(x),
                })
    return phase, adjacent, pd.DataFrame(correlations)


def observed_feature_table(inventory: pd.DataFrame) -> pd.DataFrame:
    ages = inventory["bin_center_Ma"].to_numpy(float)
    matrix = inventory["observed_global_EV"].to_numpy(float)[None, :]
    f = phase_features(matrix, ages).iloc[0]
    return pd.DataFrame([{"test": "Observed baseline", "design": "all non-zero contributing record IQRs", **f.to_dict()}])


def plot_outputs(inventory, a1, a2, density, composition_correlations):
    plt.rcParams.update({"font.family": "Arial", "font.size": 9, "axes.linewidth": 0.8})
    boundary = 66.0

    fig, axes = plt.subplots(3, 1, figsize=(7.0, 8.4), sharex=True, sharey=True)
    for ax, n in zip(axes, (3, 5, 7)):
        d = a1.loc[a1["N_records"].eq(n)].sort_values("bin_center_Ma")
        ax.fill_between(d["bin_center_Ma"], d["p2_5"], d["p97_5"], color="#9ecae1", alpha=0.55, linewidth=0)
        ax.plot(d["bin_center_Ma"], d["median_rarefied_EV"], color="#2171b5", lw=1.3, label="Rarefied median")
        ax.plot(d["bin_center_Ma"], d["observed_global_EV"], color="#333333", lw=0.9, label="Observed EV")
        ax.axvline(boundary, color="#777777", ls="--", lw=0.8)
        ax.set_title(f"N = {n} records per bin", loc="left", fontweight="bold")
        ax.set_ylabel("Global EV")
    axes[-1].set_xlabel("Age (Ma)")
    axes[-1].invert_xaxis()
    axes[0].legend(frameon=False, ncol=2, loc="upper left")
    fig.tight_layout()
    save_figure(fig, "Fig_stage3_A1_fixed_record_rarefaction")

    fig, ax = plt.subplots(figsize=(7.0, 3.7))
    d = a2.sort_values("bin_center_Ma")
    ax.fill_between(d["bin_center_Ma"], d["p2_5"], d["p97_5"], color="#fdae6b", alpha=0.5, linewidth=0)
    ax.plot(d["bin_center_Ma"], d["median_rarefied_EV"], color="#d94801", lw=1.4, label="m = 4, N = 5 median")
    ax.plot(d["bin_center_Ma"], d["observed_global_EV"], color="#333333", lw=0.9, label="Observed EV")
    ax.axvline(boundary, color="#777777", ls="--", lw=0.8)
    ax.set(xlabel="Age (Ma)", ylabel="Global EV")
    ax.invert_xaxis()
    ax.legend(frameon=False)
    fig.tight_layout()
    save_figure(fig, "Fig_stage3_A2_joint_measurement_record_rarefaction")

    predictors = [
        ("contributing_record_count", "Contributing records"),
        ("raw_measurement_count_contributing_record_bins", "Measurements in contributing record-bins"),
        ("raw_measurement_count_all_retained", "All retained measurements"),
        ("contributing_canonical_site_count", "Contributing canonical sites"),
    ]
    chronological = inventory.sort_values("bin_center_Ma", ascending=False).reset_index(drop=True)
    fig, axes = plt.subplots(2, 4, figsize=(12.0, 6.0))
    for col_idx, (column, label) in enumerate(predictors):
        for row_idx, scale in enumerate(("raw_bins", "first_difference_older_to_younger")):
            ax = axes[row_idx, col_idx]
            if scale == "raw_bins":
                x = chronological[column].to_numpy(float)
                y = chronological["observed_global_EV"].to_numpy(float)
                xlab = label
            else:
                x = np.diff(chronological[column].to_numpy(float))
                y = np.diff(chronological["observed_global_EV"].to_numpy(float))
                xlab = "Δ " + label
                ax.axhline(0, color="#bbbbbb", lw=0.6)
                ax.axvline(0, color="#bbbbbb", lw=0.6)
            ax.scatter(x, y, s=16, color="#4c78a8", alpha=0.8)
            if np.ptp(x) > 0:
                slope, intercept = np.polyfit(x, y, 1)
                xx = np.linspace(x.min(), x.max(), 100)
                ax.plot(xx, intercept + slope * xx, color="#b44b4b", lw=1.0)
            label_key = {
                "contributing_record_count": "contributing_operational_record_count",
                "raw_measurement_count_contributing_record_bins": "raw_measurement_count_contributing_record_bins_primary",
                "raw_measurement_count_all_retained": "raw_measurement_count_all_retained_secondary",
                "contributing_canonical_site_count": "contributing_canonical_site_count",
            }[column]
            ann = density.loc[
                density["predictor"].eq(label_key) & density["scale"].eq(scale)
                & density["method"].eq("Spearman") & density["block_length_bins"].eq(4)
            ].iloc[0]
            ax.text(0.03, 0.97, f"ρ={ann.estimate:.2f}\n95% CI [{ann.CI_2_5:.2f}, {ann.CI_97_5:.2f}]",
                    transform=ax.transAxes, ha="left", va="top", fontsize=7.5)
            ax.set_xlabel(xlab, fontsize=8)
            if col_idx == 0:
                ax.set_ylabel("Global EV" if row_idx == 0 else "Δ Global EV")
    fig.tight_layout()
    save_figure(fig, "Fig_stage3_B_density_associations")

    chronological = inventory.sort_values("bin_center_Ma", ascending=False)
    colors = {"temperature": "#d95f02", "d18O": "#1b9e77", "d13C": "#7570b3", "weathering": "#e6ab02", "pCO2": "#666666"}
    fig, axes = plt.subplots(2, 1, figsize=(8.0, 6.4), sharex=True)
    x = chronological["bin_center_Ma"].to_numpy(float)
    bottom = np.zeros(len(x))
    for category in CATEGORIES:
        y = chronological[f"records_{category}"].to_numpy(float)
        axes[0].bar(x, y, bottom=bottom, width=0.045, color=colors[category], label=CATEGORY_LABELS[category])
        bottom += y
    bottom = np.zeros(len(x))
    for category in CATEGORIES:
        y = chronological[f"proportion_{category}"].to_numpy(float)
        axes[1].bar(x, y, bottom=bottom, width=0.045, color=colors[category])
        bottom += y
    for ax in axes:
        ax.axvline(boundary, color="#333333", ls="--", lw=0.8)
        ax.invert_xaxis()
    axes[0].set_ylabel("Contributing records")
    axes[1].set_ylabel("Record proportion")
    axes[1].set_xlabel("Age (Ma)")
    axes[0].legend(frameon=False, ncol=5, loc="upper center", bbox_to_anchor=(0.5, 1.18))
    fig.tight_layout()
    save_figure(fig, "Fig_stage3_proxy_composition_diagnostic")


def save_figure(fig, stem: str):
    for suffix in ("pdf", "png", "svg"):
        kwargs = {"dpi": 450} if suffix == "png" else {}
        fig.savefig(FIGURES / f"{stem}.{suffix}", bbox_inches="tight", **kwargs)
    plt.close(fig)


def simple_markdown(frame: pd.DataFrame) -> str:
    shown = frame.copy()
    for column in shown.columns:
        if pd.api.types.is_float_dtype(shown[column]):
            shown[column] = shown[column].map(lambda x: "" if pd.isna(x) else f"{x:.4f}")
    headers = [str(c) for c in shown.columns]
    rows = [[str(v) for v in row] for row in shown.itertuples(index=False, name=None)]
    return "\n".join([
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
        *["| " + " | ".join(row) + " |" for row in rows],
    ])


def main():
    standard, contributing, zeros, measurements = load_current_baseline()
    inventory = build_inventory(contributing, measurements)

    a1, a1_features, a1_convergence, a1_feature_iterations, a1_matrices = record_rarefaction(contributing, inventory)
    a2, a2_features, a2_convergence, a2_feature_iterations, a2_matrix = hierarchical_rarefaction(contributing, measurements, inventory)
    density = density_associations(inventory)
    phase_composition, adjacent_composition, composition_correlations = composition_outputs(inventory)
    observed_features = observed_feature_table(inventory)

    feature_summary = pd.concat([a1_features, a2_features], ignore_index=True)
    convergence = pd.concat([a1_convergence, a2_convergence], ignore_index=True)
    feature_iterations = pd.concat([a1_feature_iterations, a2_feature_iterations], ignore_index=True)

    # Verification: observed Global EV must equal the non-zero record-IQR mean in every bin.
    independently_observed = contributing.groupby("bin_center_Ma")["IQR"].mean().sort_index()
    compared = inventory.set_index("bin_center_Ma")["observed_global_EV"].sort_index()
    max_ev_error = float((independently_observed - compared).abs().max())
    if max_ev_error > 1e-12:
        raise AssertionError(f"Observed EV mean identity failed: {max_ev_error}")

    outputs = {
        "audited_50kyr_sampling_composition_table.csv": inventory,
        "excluded_zero_IQR_record_bins.csv": zeros,
        "A1_fixed_record_rarefaction_bin_summary.csv": a1,
        "A1_fixed_record_rarefaction_feature_iterations.csv": a1_feature_iterations,
        "A2_joint_measurement_record_rarefaction_bin_summary.csv": a2,
        "A2_joint_measurement_record_rarefaction_feature_iterations.csv": a2_feature_iterations,
        "rarefaction_feature_distribution_summary.csv": feature_summary,
        "rarefaction_convergence_summary.csv": convergence,
        "density_associations_moving_block_bootstrap.csv": density,
        "proxy_composition_phase_summary.csv": phase_composition,
        "proxy_composition_adjacent_changes.csv": adjacent_composition,
        "proxy_composition_EV_associations.csv": composition_correlations,
        "observed_EV_phase_feature_summary.csv": observed_features,
    }
    for filename, frame in outputs.items():
        frame.to_csv(RESULTS / filename, index=False)

    primary_density = density.loc[density["block_length_bins"].eq(PRIMARY_BLOCK)].copy()
    phase_key = feature_summary.loc[feature_summary["feature"].isin(["phase1_minus_pre", "phase2_minus_pre", "recovery_minus_phase2"])]
    comp_ranges = pd.DataFrame([
        {
            "category": category,
            "minimum_proportion": inventory[f"proportion_{category}"].min(),
            "maximum_proportion": inventory[f"proportion_{category}"].max(),
            "range": inventory[f"proportion_{category}"].max() - inventory[f"proportion_{category}"].min(),
            "bins_absent": int(inventory[f"absent_{category}"].sum()),
        }
        for category in CATEGORIES
    ])
    comp_ranges.to_csv(RESULTS / "proxy_composition_range_and_absence.csv", index=False)

    metadata = {
        "analysis_scope": "standard 50-kyr bins, 67.2-64.8 Ma",
        "baseline_definition": "arithmetic mean of non-zero contributing record-level IQRs",
        "record_iqr_source": str(IQR_SOURCE),
        "measurement_source": str(MEAS_SOURCE),
        "record_iqr_source_sha256": sha256(IQR_SOURCE),
        "measurement_source_sha256": sha256(MEAS_SOURCE),
        "all_record_bins_before_zero_exclusion": len(standard),
        "excluded_zero_IQR_record_bins": len(zeros),
        "contributing_record_bins": len(contributing),
        "operational_records_contributing": int(contributing["record_id"].nunique()),
        "bins": len(inventory),
        "A1_iterations": ITERATIONS,
        "A1_seeds": SEEDS_A1,
        "A2_iterations": ITERATIONS,
        "A2_seed": SEED_A2,
        "A2_m_measurements": 4,
        "A2_N_records": 5,
        "TestB_iterations": ITERATIONS,
        "TestB_block_lengths": BLOCK_LENGTHS,
        "TestB_primary_block_length": PRIMARY_BLOCK,
        "maximum_observed_EV_identity_error": max_ev_error,
        "quantile_definition": "linear interpolation (NumPy method='linear'; equivalent to R type 7)",
        "composition_control_performed": False,
        "EV_biotic_analysis_performed": False,
        "monte_carlo_performed": False,
    }
    (RESULTS / "stage3_analysis_provenance.json").write_text(json.dumps(metadata, indent=2, ensure_ascii=False) + "\n")

    report = [
        "# Stage 3 sampling-density robustness and proxy-composition diagnosis",
        "",
        "## Audited baseline",
        "",
        f"- Standard 50-kyr bins: {len(inventory)} (67.2-64.8 Ma).",
        f"- Current source contains {len(standard)} record-bin IQRs; {len(zeros)} zero-IQR rows are excluded from the formal mean, leaving {len(contributing)} contributing record-bins from {contributing['record_id'].nunique()} operational records.",
        f"- Per-bin contributing records: {inventory.contributing_record_count.min()}-{inventory.contributing_record_count.max()} (median {inventory.contributing_record_count.median():.0f}).",
        f"- Per-bin primary raw measurements: {inventory.raw_measurement_count_contributing_record_bins.min()}-{inventory.raw_measurement_count_contributing_record_bins.max()} (median {inventory.raw_measurement_count_contributing_record_bins.median():.0f}).",
        f"- Per-bin secondary all-retained measurements: {inventory.raw_measurement_count_all_retained.min()}-{inventory.raw_measurement_count_all_retained.max()} (median {inventory.raw_measurement_count_all_retained.median():.0f}).",
        f"- Per-bin contributing canonical sites: {inventory.contributing_canonical_site_count.min()}-{inventory.contributing_canonical_site_count.max()} (median {inventory.contributing_canonical_site_count.median():.0f}).",
        "",
        "## Rarefaction results",
        "",
        simple_markdown(phase_key),
        "",
        "## Density associations (moving-block bootstrap, primary block length = 4 bins)",
        "",
        simple_markdown(primary_density),
        "",
        "## Proxy-composition diagnosis",
        "",
        simple_markdown(comp_ranges),
        "",
        "Composition varies through time, so a composition-controlled estimator may be scientifically informative, but it would change the estimand. No proxy-equal weighting, stratified rarefaction, leave-one-proxy-out analysis, or fixed-composition reweighting was run in this stage.",
        "",
        "Candidate control options for author review are: (1) proxy-stratified rarefaction with predeclared per-category quotas, restricted to bins satisfying all quotas; or (2) direct standardization to fixed proxy weights over a predeclared common-coverage interval. Both require explicit confirmation before implementation.",
        "",
        "## Scope",
        "",
        "No Monte Carlo uncertainty propagation or EV-biotic relationship analysis was performed.",
    ]
    (RESULTS / "Stage3_sampling_density_robustness_report.md").write_text("\n".join(report) + "\n")

    print(json.dumps({
        "inventory_rows": len(inventory),
        "contributing_record_bins": len(contributing),
        "A1_rows": len(a1),
        "A2_rows": len(a2),
        "density_rows": len(density),
        "figure_dir": str(FIGURES),
        "result_dir": str(RESULTS),
    }, indent=2))


if __name__ == "__main__":
    raise SystemExit("Preserved numerical library: run ../run.py for the packaged final workflow.")
