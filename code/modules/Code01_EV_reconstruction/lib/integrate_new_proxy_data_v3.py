#!/usr/bin/env python3
"""Integrate approved 2024-2025 proxy records into the revised K-Pg EV baseline.

This script only updates measurement- and EV-level numerical products. It does
not update manuscripts, figures, EV-biotic analyses, composition-controlled
resampling, or Monte Carlo uncertainty propagation.
"""

from __future__ import annotations

import hashlib
import json
import math
from collections import OrderedDict
from pathlib import Path

import numpy as np
import pandas as pd


PACKAGE = Path(__file__).resolve().parents[1]
DATA = PACKAGE / "inputs"
RESULTS = PACKAGE / "outputs" / "generated"
LOGS = PACKAGE / "outputs" / "generated"
SOURCE = DATA / "source_data_2024_2025"

SEED = 20260817
ITERATIONS = 10_000

BIN_SPECS = OrderedDict(
    [
        ("200", {"bin_width_kyr": 200, "alignment": "standard", "end_age_Ma": 63.1, "n_bins": 29}),
        ("150", {"bin_width_kyr": 150, "alignment": "standard", "end_age_Ma": 63.975, "n_bins": 27}),
        ("100", {"bin_width_kyr": 100, "alignment": "standard", "end_age_Ma": 64.75, "n_bins": 25}),
        ("50", {"bin_width_kyr": 50, "alignment": "standard", "end_age_Ma": 64.775, "n_bins": 49}),
        ("50_offset25", {"bin_width_kyr": 50, "alignment": "offset_25_kyr", "end_age_Ma": 64.8, "n_bins": 48}),
    ]
)

CATEGORY_ORDER = {"temperature": 0, "d18O": 1, "d13C": 2, "pCO2": 3, "weathering": 4}
DEDUP_REASONS = {
    "temperature::ODP1209 2": "Same Hull et al. bulk-foraminiferal d18O measurements as retained d18O record",
    "temperature::ODP1209 3": "Same Hull et al. planktonic-foraminiferal d18O measurements as retained d18O record",
    "temperature::ODP1262": "Same benthic-foraminiferal d18O measurements as retained d18O record",
    "temperature::DSDP465": "Same benthic-foraminiferal d18O measurements as retained d18O record",
    "temperature::ODP1210": "Same benthic-foraminiferal d18O measurements as retained d18O record",
}


def save_csv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False, float_format="%.12g")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_pangaea(path: Path) -> pd.DataFrame:
    lines = path.read_text(encoding="utf-8").splitlines()
    header = next(i for i, line in enumerate(lines) if line.startswith("Event\t"))
    return pd.read_csv(path, sep="\t", skiprows=header, dtype=str)


def numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def rss_mean_error(values: pd.Series) -> float:
    x = numeric(values).dropna().to_numpy(float)
    return float(np.sqrt(np.sum(x * x)) / len(x)) if len(x) else np.nan


def new_row(
    *, record_id: str, category: str, sheet: str, source_row: int, study: str,
    site: str, canonical_site: str, proxy_material: str, sample_id: str,
    age_ma: float, value: float, analytical_error: float, replicate_count: int,
    source_doi: str, source_file: str, source_age_model: str,
) -> dict:
    return {
        "record_id": record_id,
        "category": category,
        "sheet": sheet,
        "source_row": source_row,
        "study": study,
        "site": site,
        "canonical_site": canonical_site,
        "proxy_material": proxy_material,
        "sample_id": sample_id,
        "age_Ma": float(age_ma),
        "age_ka_as_compiled": float(age_ma) * 1000,
        "value": float(value),
        "analytical_error": float(analytical_error) if np.isfinite(analytical_error) else np.nan,
        "replicate_count": int(replicate_count),
        "source_DOI": source_doi,
        "source_file": source_file,
        "source_age_model": source_age_model,
        "source_status": "new_approved_proxy",
        "dedup_excluded": False,
    }


def parse_westerhold_stable() -> list[dict]:
    path = SOURCE / "PANGAEA_964956.tab"
    raw = read_pangaea(path)
    age_col = "Age [Ma] (CycleCount time)"
    sample_col = "Sample label (DSDP/ODP/IODP sample designation)"
    c_col = "Foram benth δ13C [‰ PDB] (Species Corrected, Corrected)"
    o_col = "Foram benth δ18O [‰ PDB] (Species Corrected, Corrected)"
    raw["_age"] = numeric(raw[age_col])
    raw["_c"] = numeric(raw[c_col])
    raw["_o"] = numeric(raw[o_col])
    valid = raw.dropna(subset=["_age", "_c", "_o"]).copy()
    grouped = (
        valid.groupby([sample_col, "_age"], sort=True, dropna=False)
        .agg(event=("Event", "first"), d13C=("_c", "mean"), d18O=("_o", "mean"), replicate_count=("_c", "size"))
        .reset_index()
    )
    assert len(valid) == 299 and len(grouped) == 292
    rows: list[dict] = []
    for i, r in grouped.iterrows():
        common = dict(
            source_row=int(i + 1), study="Westerhold et al. (2025)", site="ODP1209",
            canonical_site="ODP 1209", proxy_material="Benthic foraminifera",
            sample_id=str(r[sample_col]), age_ma=float(r["_age"]), analytical_error=np.nan,
            replicate_count=int(r["replicate_count"]), source_doi="https://doi.org/10.1594/PANGAEA.964956",
            source_file=path.name, source_age_model="CycleCount time",
        )
        rows.append(new_row(record_id="d13C::ODP1209 6", category="d13C", sheet="ODP1209 6", value=float(r.d13C), **common))
        rows.append(new_row(record_id="d18O::ODP1209 6", category="d18O", sheet="ODP1209 6", value=float(r.d18O), **common))
    return rows


def parse_westerhold_os(filename: str, record_id: str, sheet: str, site: str, canonical_site: str, doi: str) -> list[dict]:
    path = SOURCE / filename
    raw = read_pangaea(path)
    age_col = next(c for c in raw.columns if c.startswith("Age [Ma] (CycleCount"))
    sample_col = "Sample label"
    value_col = "187Os/188Os"
    error_col = "187Os/188Os e [±]"
    raw["_age"] = numeric(raw[age_col])
    raw["_value"] = numeric(raw[value_col])
    raw["_error"] = numeric(raw[error_col])
    valid = raw.dropna(subset=["_age", "_value"]).copy()
    grouped = (
        valid.groupby([sample_col, "_age"], sort=True, dropna=False)
        .agg(value=("_value", "mean"), analytical_error=("_error", rss_mean_error), replicate_count=("_value", "size"))
        .reset_index()
    )
    expected = 16 if "964772" in filename else 22
    assert len(grouped) == expected
    return [
        new_row(
            record_id=record_id, category="weathering", sheet=sheet, source_row=int(i + 1),
            study="Westerhold et al. (2025)", site=site, canonical_site=canonical_site,
            proxy_material="Sediment", sample_id=str(r[sample_col]), age_ma=float(r["_age"]),
            value=float(r.value), analytical_error=float(r.analytical_error),
            replicate_count=int(r.replicate_count), source_doi=doi, source_file=path.name,
            source_age_model="CycleCount time",
        )
        for i, r in grouped.iterrows()
    ]


def parse_nanxiong() -> list[dict]:
    path = SOURCE / "nanxiong_pco2_mmc2.xlsx"
    raw = pd.read_excel(path, sheet_name="data", header=None).iloc[5:].copy()
    raw["_sample"] = raw.iloc[:, 3].astype(str)
    raw["_age"] = numeric(raw.iloc[:, 8])
    raw["_d18O"] = numeric(raw.iloc[:, 9])
    raw["_d13C"] = numeric(raw.iloc[:, 10])
    raw["_pCO2"] = numeric(raw.iloc[:, 21])
    raw["_d13C_error"] = numeric(raw.iloc[:, 23])
    valid = raw[raw["_age"].between(63, 69)].dropna(subset=["_age", "_d18O", "_d13C", "_pCO2"]).copy()
    grouped = (
        valid.groupby("_age", sort=True)
        .agg(
            sample_id=("_sample", lambda x: " | ".join(x.astype(str))),
            d18O=("_d18O", "mean"), d13C=("_d13C", "mean"), pCO2=("_pCO2", "mean"),
            d13C_error=("_d13C_error", rss_mean_error), replicate_count=("_pCO2", "size"),
        )
        .reset_index()
    )
    assert len(valid) == 120 and len(grouped) == 64
    rows: list[dict] = []
    for i, r in grouped.iterrows():
        common = dict(
            source_row=int(i + 1), study="Wang et al. (2024)", site="Nanxiong Basin",
            canonical_site="Nanxiong Basin", sample_id=str(r.sample_id), age_ma=float(r["_age"]),
            replicate_count=int(r.replicate_count), source_doi="https://doi.org/10.1016/j.palaeo.2024.112177",
            source_file=path.name, source_age_model="Calibrated age model (second Age column)",
        )
        rows.append(new_row(record_id="d13C::Nanxiong Basin", category="d13C", sheet="Nanxiong Basin", proxy_material="Pedogenic carbonate", value=float(r.d13C), analytical_error=float(r.d13C_error), **common))
        rows.append(new_row(record_id="d18O::Nanxiong Basin", category="d18O", sheet="Nanxiong Basin", proxy_material="Pedogenic carbonate", value=float(r.d18O), analytical_error=np.nan, **common))
        rows.append(new_row(record_id="pCO2::all data", category="pCO2", sheet="all data", proxy_material="Paleosol carbonate", value=float(r.pCO2), analytical_error=np.nan, **common))
    return rows


def canonical_site_for_legacy(record_id: str, site: str) -> str:
    sheet = record_id.split("::", 1)[1]
    rules = [
        ("1209", "ODP 1209"), ("1262", "ODP 1262"), ("1267", "ODP 1267"),
        ("1210", "ODP 1210"), ("U1403", "IODP U1403"), ("690", "ODP 690"),
        ("525", "DSDP/ODP 525"), ("577", "DSDP/ODP 577"), ("465", "DSDP 465"),
        ("Seymour", "Seymour Island"), ("Songliao", "Songliao Basin"),
        ("SK-In", "Songliao Basin"), ("Bajada", "Bajada del Jagüel"),
    ]
    for token, canonical in rules:
        if token.lower() in sheet.lower() or token.lower() in str(site).lower():
            return canonical
    return str(site).strip() or sheet


def center_for_age(age: float, spec: dict) -> float | None:
    if not np.isfinite(age):
        return None
    width = spec["bin_width_kyr"] / 1000
    for i in range(spec["n_bins"]):
        lower = spec["end_age_Ma"] + i * width
        upper = spec["end_age_Ma"] + (i + 1) * width
        if lower <= age < upper:
            return round(float(lower + width / 2), 6)
    return None


def recompute_zscores(measurements: pd.DataFrame) -> pd.DataFrame:
    df = measurements.copy()
    means = df.groupby("record_id")["value"].transform("mean")
    sds = df.groupby("record_id")["value"].transform(lambda x: x.std(ddof=1))
    assert not (sds <= 0).any()
    df["z_score_recomputed"] = (df["value"] - means) / sds
    df["bin_center_Ma"] = [center_for_age(a, BIN_SPECS["50"]) for a in df.age_Ma]
    return df


def calculate_iqrs(measurements: pd.DataFrame, min_n: int = 2, schemes: list[str] | None = None) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows, binned = [], []
    use = measurements[~measurements.dedup_excluded].copy()
    for scheme in schemes or list(BIN_SPECS):
        spec = BIN_SPECS[scheme]
        part = use.copy()
        part["scheme"] = scheme
        part["bin_width_kyr"] = spec["bin_width_kyr"]
        part["alignment"] = spec["alignment"]
        part["bin_center_Ma"] = [center_for_age(a, spec) for a in part.age_Ma]
        part = part[part.bin_center_Ma.notna()].copy()
        binned.append(part)
        for (record_id, category, sheet, center), g in part.groupby(["record_id", "category", "sheet", "bin_center_Ma"], sort=False):
            if len(g) < min_n:
                continue
            z = g.z_score_recomputed.to_numpy(float)
            rows.append({
                "scheme": scheme, "bin_width_kyr": spec["bin_width_kyr"], "alignment": spec["alignment"],
                "record_id": record_id, "category": category, "sheet": sheet,
                "bin_center_Ma": float(center), "n_measurements": len(g),
                "IQR": float(np.quantile(z, 0.75) - np.quantile(z, 0.25)),
                "sites_in_record_bin": " | ".join(sorted(set(g.site.astype(str)))),
            })
    out = pd.DataFrame(rows)
    order = pd.read_csv(DATA / "record_sheet_order_v3.csv").set_index("record_id").sheet_order.to_dict()
    out["category_order"] = out.category.map(CATEGORY_ORDER)
    out["sheet_order"] = out.record_id.map(order)
    out = out.sort_values(["scheme", "category_order", "sheet_order", "bin_center_Ma"]).drop(columns="category_order")
    return out.reset_index(drop=True), pd.concat(binned, ignore_index=True)


def aggregate_iqrs(iqrs: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    global_ev = (
        iqrs.groupby(["scheme", "bin_width_kyr", "alignment", "bin_center_Ma"], as_index=False)
        .agg(revised_EV=("IQR", "mean"), contributing_record_count=("IQR", "size"))
    )
    proxy_ev = (
        iqrs.groupby(["scheme", "bin_width_kyr", "alignment", "bin_center_Ma", "category"], as_index=False)
        .agg(proxy_EV=("IQR", "mean"), contributing_record_count=("IQR", "size"))
    )
    return global_ev, proxy_ev


def aggregate_zscore(binned: pd.DataFrame) -> pd.DataFrame:
    return (
        binned.groupby(["scheme", "bin_width_kyr", "alignment", "bin_center_Ma"], as_index=False)
        .agg(mean_z_score=("z_score_recomputed", "mean"), measurement_count=("z_score_recomputed", "size"))
    )


def correlation_estimate(x, y, method: str) -> float:
    x, y = pd.Series(x, dtype=float), pd.Series(y, dtype=float)
    good = x.notna() & y.notna()
    x, y = x[good], y[good]
    if method == "Spearman":
        x, y = x.rank(method="average"), y.rank(method="average")
    return float(np.corrcoef(x.to_numpy(float), y.to_numpy(float))[0, 1])


def bootstrap_ci(x, y, method: str, seed: int) -> tuple[float, float]:
    x, y = np.asarray(x, float), np.asarray(y, float)
    good = np.isfinite(x) & np.isfinite(y)
    x, y = x[good], y[good]
    rng, values = np.random.default_rng(seed), []
    for _ in range(ITERATIONS):
        idx = rng.integers(0, len(x), len(x))
        if np.std(x[idx]) and np.std(y[idx]):
            values.append(correlation_estimate(x[idx], y[idx], method))
    return tuple(np.quantile(values, [0.025, 0.975]))


def sampling_inventory(measurements: pd.DataFrame, iqrs50: pd.DataFrame) -> pd.DataFrame:
    spec = BIN_SPECS["50"]
    retained = measurements[~measurements.dedup_excluded].copy()
    retained["center"] = [center_for_age(a, spec) for a in retained.age_Ma]
    rows = []
    for i in range(spec["n_bins"]):
        center = round(spec["end_age_Ma"] + (i + .5) * .05, 6)
        raw = retained[retained.center == center]
        contrib = iqrs50[iqrs50.bin_center_Ma == center]
        ids = set(contrib.record_id)
        used = raw[raw.record_id.isin(ids)]
        row = {
            "bin_center_Ma": center,
            "raw_measurement_count_all_revised_records": len(raw),
            "raw_measurement_count_contributing_records": len(used),
            "contributing_record_count": len(contrib),
            "contributing_site_count": used.site.nunique(),
            "proxy_category_count": contrib.category.nunique(),
        }
        for cat in CATEGORY_ORDER:
            n = int((contrib.category == cat).sum())
            row[f"records_{cat}"] = n
            row[f"proportion_{cat}"] = n / len(contrib) if len(contrib) else np.nan
        rows.append(row)
    return pd.DataFrame(rows)


def rarefaction(iqrs50: pd.DataFrame) -> pd.DataFrame:
    rows = []
    groups = {float(k): g.IQR.to_numpy(float) for k, g in iqrs50.groupby("bin_center_Ma")}
    centers = [round(BIN_SPECS["50"]["end_age_Ma"] + (i + .5) * .05, 6) for i in range(49)]
    for n_draw in [3, 5, 7]:
        rng = np.random.default_rng(SEED + n_draw)
        for center in centers:
            x = groups.get(center, np.array([]))
            row = {"N": n_draw, "bin_center_Ma": center, "available_record_count": len(x), "eligible": len(x) >= n_draw, "iterations": ITERATIONS, "replacement": False, "seed": SEED + n_draw}
            if len(x) >= n_draw:
                priorities = rng.random((ITERATIONS, len(x)))
                selected = np.argpartition(priorities, n_draw - 1, axis=1)[:, :n_draw]
                means = x[selected].mean(axis=1)
                row.update({"median_resampled_EV": np.median(means), "p2_5": np.quantile(means, .025), "p97_5": np.quantile(means, .975)})
            else:
                row.update({"median_resampled_EV": np.nan, "p2_5": np.nan, "p97_5": np.nan})
            rows.append(row)
    return pd.DataFrame(rows)


def density_correlations(global50: pd.DataFrame, inventory: pd.DataFrame) -> pd.DataFrame:
    merged = global50[["bin_center_Ma", "revised_EV"]].merge(inventory, on="bin_center_Ma")
    variables = ["contributing_record_count", "raw_measurement_count_all_revised_records", "raw_measurement_count_contributing_records", "contributing_site_count"]
    rows = []
    for vi, variable in enumerate(variables):
        for mi, method in enumerate(["Pearson", "Spearman"]):
            estimate = correlation_estimate(merged[variable], merged.revised_EV, method)
            lo, hi = bootstrap_ci(merged[variable], merged.revised_EV, method, SEED + 500 + 20 * vi + mi)
            rows.append({"response": "revised_EV_v3", "predictor": variable, "method": method, "n_bins": len(merged), "estimate": estimate, "p_value": np.nan, "bootstrap_95CI_low": lo, "bootstrap_95CI_high": hi, "bootstrap_iterations": ITERATIONS, "seed": SEED + 500 + 20 * vi + mi})
    return pd.DataFrame(rows)


def proxy_correlations(global50: pd.DataFrame, proxy50: pd.DataFrame) -> pd.DataFrame:
    wide = proxy50.pivot(index="bin_center_Ma", columns="category", values="proxy_EV").reset_index()
    merged = global50[["bin_center_Ma", "revised_EV"]].merge(wide, on="bin_center_Ma", how="left")
    phases = {"full_67.20_to_64.80_Ma": (64.80, 67.20), "phase_1_66.90_to_66.45_Ma": (66.45, 66.90), "phase_2_66.40_to_66.05_Ma": (66.05, 66.40)}
    rows = []
    for phase, (lo, hi) in phases.items():
        part = merged[merged.bin_center_Ma.between(lo, hi)]
        for category in CATEGORY_ORDER:
            if category not in part:
                continue
            xy = part[["revised_EV", category]].dropna()
            if len(xy) < 3 or xy[category].nunique() < 2:
                continue
            for method in ["Pearson", "Spearman"]:
                rows.append({"phase": phase, "proxy": category, "method": method, "n_bins": len(xy), "estimate": correlation_estimate(xy.revised_EV, xy[category], method), "p_value": np.nan})
    return pd.DataFrame(rows)


def comparison_outputs(global_v3: pd.DataFrame, proxy_v3: pd.DataFrame, iqrs_v3: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    old_global = pd.read_csv(DATA / "global_EV_revised_all_bin_widths.csv").rename(columns={"revised_EV": "v2_EV", "contributing_record_count": "v2_record_count"})
    global_cmp = old_global.merge(global_v3.rename(columns={"revised_EV": "v3_EV", "contributing_record_count": "v3_record_count"}), on=["scheme", "bin_width_kyr", "alignment", "bin_center_Ma"], how="outer")
    global_cmp["delta_v3_minus_v2"] = global_cmp.v3_EV - global_cmp.v2_EV

    old_proxy = pd.read_csv(DATA / "single_proxy_EV_revised_all_bin_widths.csv").rename(columns={"proxy_EV": "v2_proxy_EV", "contributing_record_count": "v2_record_count"})
    proxy_cmp = old_proxy.merge(proxy_v3.rename(columns={"proxy_EV": "v3_proxy_EV", "contributing_record_count": "v3_record_count"}), on=["scheme", "bin_width_kyr", "alignment", "bin_center_Ma", "category"], how="outer")
    proxy_cmp["delta_v3_minus_v2"] = proxy_cmp.v3_proxy_EV - proxy_cmp.v2_proxy_EV

    metrics = []
    features = []
    for scheme, g in global_cmp.groupby("scheme"):
        pair = g.dropna(subset=["v2_EV", "v3_EV"]).copy()
        delta = pair.v3_EV - pair.v2_EV
        metrics.append({
            "scheme": scheme, "bin_width_kyr": pair.bin_width_kyr.iloc[0], "alignment": pair.alignment.iloc[0],
            "n_common_bins": len(pair), "pearson_r": correlation_estimate(pair.v2_EV, pair.v3_EV, "Pearson"),
            "mean_absolute_difference": delta.abs().mean(), "root_mean_square_difference": np.sqrt(np.mean(delta ** 2)),
            "maximum_absolute_difference": delta.abs().max(), "maximum_change_age_Ma": pair.loc[delta.abs().idxmax(), "bin_center_Ma"],
        })
        for curve in ["v2_EV", "v3_EV"]:
            series = g.dropna(subset=[curve]).set_index("bin_center_Ma")[curve].sort_index()
            phase1 = series[(series.index >= 66.45) & (series.index <= 66.90)]
            phase2 = series[(series.index >= 66.05) & (series.index <= 66.40)]
            recovery = series[(series.index >= 64.80) & (series.index <= 65.95)]
            common = series[(series.index >= 64.80) & (series.index <= 67.20)]
            slope = np.polyfit(recovery.index.to_numpy(float), recovery.to_numpy(float), 1)[0] if len(recovery) >= 2 else np.nan
            features.append({
                "scheme": scheme, "curve": curve, "full_peak_EV": series.max(), "full_peak_age_Ma": series.idxmax(),
                "common_window_peak_EV": common.max() if len(common) else np.nan, "common_window_peak_age_Ma": common.idxmax() if len(common) else np.nan,
                "phase1_mean": phase1.mean(), "phase2_mean": phase2.mean(), "phase2_max": phase2.max(), "phase2_max_age_Ma": phase2.idxmax() if len(phase2) else np.nan,
                "recovery_mean": recovery.mean(), "recovery_min": recovery.min(), "recovery_slope_vs_age": slope,
            })

    old_iqr = pd.read_csv(DATA / "record_level_IQR_all_bin_widths_revised_deduplicated_n2.csv")
    keys = ["scheme", "record_id", "category", "sheet", "bin_center_Ma"]
    iqr_cmp = old_iqr[keys + ["IQR", "n_measurements"]].rename(columns={"IQR": "v2_IQR", "n_measurements": "v2_n"}).merge(
        iqrs_v3[keys + ["IQR", "n_measurements"]].rename(columns={"IQR": "v3_IQR", "n_measurements": "v3_n"}), on=keys, how="outer", indicator=True
    )
    iqr_cmp["delta_v3_minus_v2"] = iqr_cmp.v3_IQR - iqr_cmp.v2_IQR
    iqr_cmp["status"] = np.select(
        [iqr_cmp._merge.eq("right_only"), iqr_cmp._merge.eq("left_only"), iqr_cmp.delta_v3_minus_v2.abs().fillna(0).gt(1e-10)],
        ["new_record_bin", "removed_record_bin", "changed_record_bin"], default="unchanged_record_bin"
    )
    return global_cmp, proxy_cmp, pd.DataFrame(metrics), pd.DataFrame(features), iqr_cmp


def build_order_and_latitude(records: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    old = pd.read_csv(DATA / "record_sheet_order.csv")
    rows = []
    for category in ["temperature", "d13C", "d18O", "weathering", "pCO2"]:
        names = old[old.category == category].sort_values("sheet_order").sheet.tolist()
        if category in {"d13C", "d18O"}:
            insert_at = names.index("ODP1209 5") + 1
            names.insert(insert_at, "ODP1209 6")
            names.append("Nanxiong Basin")
        for i, sheet in enumerate(names):
            rid = f"{category}::{sheet}"
            if rid in records:
                rows.append({"category": category, "sheet": sheet, "record_id": rid, "sheet_order": i})
    order = pd.DataFrame(rows)
    latitude = pd.read_csv(DATA / "latitude_record_mapping_from_submitted_S1.csv")
    latitude = pd.concat([latitude, pd.DataFrame([
        {"record_id": "d13C::ODP1209 6", "category": "d13C", "latitude_group": "low", "source_sheet": "IQR (low latitudes_50 bin)"},
        {"record_id": "d18O::ODP1209 6", "category": "d18O", "latitude_group": "low", "source_sheet": "IQR (low latitudes_50 bin)"},
        {"record_id": "d13C::Nanxiong Basin", "category": "d13C", "latitude_group": "low", "source_sheet": "IQR (low latitudes_50 bin)"},
        {"record_id": "d18O::Nanxiong Basin", "category": "d18O", "latitude_group": "low", "source_sheet": "IQR (low latitudes_50 bin)"},
    ])], ignore_index=True).drop_duplicates("record_id", keep="last")
    return order, latitude


def main() -> None:
    for path in [DATA, RESULTS, LOGS]:
        path.mkdir(parents=True, exist_ok=True)

    legacy = pd.read_csv(DATA / "recomputed_measurements_with_zscores.csv")
    legacy["canonical_site"] = [canonical_site_for_legacy(rid, site) for rid, site in zip(legacy.record_id, legacy.site)]
    legacy["sample_id"] = [f"legacy::{rid}::{row}" for rid, row in zip(legacy.record_id, legacy.source_row)]
    legacy["analytical_error"] = np.nan
    legacy["replicate_count"] = 1
    legacy["source_DOI"] = ""
    legacy["source_file"] = "submitted_Datasets_S2-S6"
    legacy["source_age_model"] = "as compiled in submitted workbook"
    legacy["source_status"] = "legacy_revised_baseline"

    new_rows = []
    new_rows.extend(parse_westerhold_stable())
    new_rows.extend(parse_westerhold_os("PANGAEA_964772.tab", "weathering::ODP1209C Os", "ODP1209C Os", "ODP1209", "ODP 1209", "https://doi.org/10.1594/PANGAEA.964772"))
    new_rows.extend(parse_westerhold_os("PANGAEA_964777.tab", "weathering::ODP1262B Os", "ODP1262B Os", "ODP1262B", "ODP 1262", "https://doi.org/10.1594/PANGAEA.964777"))
    new_rows.extend(parse_nanxiong())
    new = pd.DataFrame(new_rows)

    expected_counts = {
        "d13C::ODP1209 6": 292, "d18O::ODP1209 6": 292,
        "weathering::ODP1209C Os": 16, "weathering::ODP1262B Os": 22,
        "d13C::Nanxiong Basin": 64, "d18O::Nanxiong Basin": 64,
    }
    assert len(new) == 814
    for rid, n in expected_counts.items():
        count = len(new[new.record_id == rid]) if "weathering" not in rid else len(new[(new.record_id == rid)])
        assert count == n, (rid, count)
    assert len(new[(new.category == "pCO2") & (new.site == "Nanxiong Basin")]) == 64
    assert not new.duplicated(["record_id", "sample_id", "age_Ma"]).any()

    columns = list(dict.fromkeys(list(legacy.columns) + list(new.columns)))
    combined = pd.concat([legacy.reindex(columns=columns), new.reindex(columns=columns)], ignore_index=True)
    combined = recompute_zscores(combined)
    assert int((~combined.dedup_excluded).sum()) == 10945
    assert combined.record_id.nunique() == 75
    assert combined.loc[~combined.dedup_excluded, "record_id"].nunique() == 70

    order, latitude = build_order_and_latitude(combined.record_id.unique().tolist())
    save_csv(order, DATA / "record_sheet_order_v3.csv")
    save_csv(latitude, DATA / "latitude_record_mapping_v3.csv")

    records = []
    for rid, g in combined.groupby("record_id", sort=False):
        records.append({
            "record_id": rid, "category": g.category.iloc[0], "sheet": g.sheet.iloc[0],
            "n_measurements_total": len(g), "age_min_Ma": g.age_Ma.min(), "age_max_Ma": g.age_Ma.max(),
            "value_header": {"d13C": "δ13C", "d18O": "δ18O", "weathering": "187Os/188Os or retained weathering proxy", "pCO2": "pCO2", "temperature": "Temperature"}[g.category.iloc[0]],
            "studies": " | ".join(sorted(set(g.study.dropna().astype(str)))),
            "sites": " | ".join(sorted(set(g.site.dropna().astype(str)))),
            "canonical_sites": " | ".join(sorted(set(g.canonical_site.dropna().astype(str)))),
            "proxy_material": " | ".join(sorted(set(g.proxy_material.dropna().astype(str)))),
            "record_mean": g.value.mean(), "record_sd_sample": g.value.std(ddof=1),
            "dedup_excluded": bool(g.dedup_excluded.iloc[0]), "dedup_reason": DEDUP_REASONS.get(rid, ""),
            "source_DOIs": " | ".join(sorted(set(x for x in g.source_DOI.dropna().astype(str) if x))),
        })
    record_inventory = pd.DataFrame(records).merge(order[["record_id", "sheet_order"]], on="record_id", how="left")

    save_csv(new.sort_values(["category", "record_id", "age_Ma"]), DATA / "new_proxy_processed_measurements_v3.csv")
    save_csv(combined, DATA / "recomputed_measurements_with_zscores_v3.csv")
    save_csv(record_inventory, DATA / "record_inventory_v3.csv")

    iqrs, binned = calculate_iqrs(combined, min_n=2)
    iqrs_n4, _ = calculate_iqrs(combined, min_n=4, schemes=["50"])
    global_ev, proxy_ev = aggregate_iqrs(iqrs)
    mean_z = aggregate_zscore(binned)
    save_csv(iqrs, DATA / "record_level_IQR_all_bin_widths_revised_v3_n2.csv")
    save_csv(iqrs_n4, DATA / "record_level_IQR_revised_v3_50kyr_n4.csv")
    save_csv(binned, DATA / "measurements_binned_all_schemes_revised_v3.csv")
    save_csv(global_ev, DATA / "global_EV_revised_v3_all_bin_widths.csv")
    save_csv(proxy_ev, DATA / "single_proxy_EV_revised_v3_all_bin_widths.csv")
    save_csv(mean_z, DATA / "global_mean_zscore_revised_v3_all_bin_widths.csv")

    iqr50 = iqrs[iqrs.scheme.astype(str) == "50"].copy()
    global50 = global_ev[global_ev.scheme.astype(str) == "50"].copy()
    proxy50 = proxy_ev[proxy_ev.scheme.astype(str) == "50"].copy()
    inventory = sampling_inventory(combined, iqr50)
    rare = rarefaction(iqr50)
    density = density_correlations(global50, inventory)
    proxy_corr = proxy_correlations(global50, proxy50)
    save_csv(inventory, DATA / "per_bin_sampling_and_proxy_composition_v3.csv")
    save_csv(rare, RESULTS / "fixed_number_rarefaction_N3_N5_N7_v3.csv")
    save_csv(density, RESULTS / "EV_sampling_density_correlations_v3.csv")
    save_csv(proxy_corr, RESULTS / "EV_single_proxy_correlations_v3.csv")

    global_cmp, proxy_cmp, metrics, features, iqr_cmp = comparison_outputs(global_ev, proxy_ev, iqrs)
    save_csv(global_cmp, DATA / "global_EV_v2_vs_v3_all_bin_widths.csv")
    save_csv(proxy_cmp, DATA / "single_proxy_EV_v2_vs_v3_all_bin_widths.csv")
    save_csv(metrics, RESULTS / "EV_v2_v3_curve_comparison_metrics.csv")
    save_csv(features, RESULTS / "EV_v2_v3_feature_comparison.csv")
    save_csv(iqr_cmp, RESULTS / "new_proxy_record_bin_change_detail.csv")

    series_summary = pd.DataFrame([
        {"series": "ODP1209 benthic d13C", "record_id": "d13C::ODP1209 6", "processed_measurements": 292, "action": "new operational record", "source_DOI": "https://doi.org/10.1594/PANGAEA.964956"},
        {"series": "ODP1209 benthic d18O", "record_id": "d18O::ODP1209 6", "processed_measurements": 292, "action": "new operational record", "source_DOI": "https://doi.org/10.1594/PANGAEA.964956"},
        {"series": "ODP1209 measured 187Os/188Os", "record_id": "weathering::ODP1209C Os", "processed_measurements": 16, "action": "appended to existing record", "source_DOI": "https://doi.org/10.1594/PANGAEA.964772"},
        {"series": "ODP1262 measured 187Os/188Os", "record_id": "weathering::ODP1262B Os", "processed_measurements": 22, "action": "appended to existing record", "source_DOI": "https://doi.org/10.1594/PANGAEA.964777"},
        {"series": "Nanxiong pedogenic-carbonate d13C", "record_id": "d13C::Nanxiong Basin", "processed_measurements": 64, "action": "new operational record", "source_DOI": "https://doi.org/10.1016/j.palaeo.2024.112177"},
        {"series": "Nanxiong pedogenic-carbonate d18O", "record_id": "d18O::Nanxiong Basin", "processed_measurements": 64, "action": "new operational record", "source_DOI": "https://doi.org/10.1016/j.palaeo.2024.112177"},
        {"series": "Nanxiong nominal pCO2", "record_id": "pCO2::all data", "processed_measurements": 64, "action": "appended to pooled record", "source_DOI": "https://doi.org/10.1016/j.palaeo.2024.112177"},
    ])
    save_csv(series_summary, RESULTS / "new_proxy_integration_series_summary.csv")
    site_summary = pd.DataFrame([
        {"definition": "canonical geographic sites", "v2_count": 31, "new_count": 1, "v3_count": 32, "note": "Only Nanxiong Basin is geographically new; ODP1209 and ODP1262 already existed."},
        {"definition": "retained operational records", "v2_count": 66, "new_count": 4, "v3_count": 70, "note": "Two ODP1209 isotope records and two Nanxiong isotope records added; Os and pCO2 appended to existing records."},
        {"definition": "retained proxy measurement entries", "v2_count": 10131, "new_count": 814, "v3_count": 10945, "note": "Counts are replicate-resolved measurements entering retained records."},
    ])
    save_csv(site_summary, DATA / "site_record_measurement_count_summary_v3.csv")

    source_files = [SOURCE / "PANGAEA_964956.tab", SOURCE / "PANGAEA_964772.tab", SOURCE / "PANGAEA_964777.tab", SOURCE / "nanxiong_pco2_mmc2.xlsx"]
    source_hashes = pd.DataFrame([{"file": p.name, "path": str(p), "sha256": sha256(p)} for p in source_files])
    save_csv(source_hashes, LOGS / "new_proxy_source_hashes_v3.csv")

    zcheck = combined.groupby("record_id").z_score_recomputed.agg(["count", "mean", "std"]).reset_index()
    assert zcheck["mean"].abs().max() < 1e-10
    assert (zcheck["std"] - 1).abs().max() < 1e-10
    for _, row in global_ev.iterrows():
        subset = iqrs[(iqrs.scheme.astype(str) == str(row.scheme)) & (iqrs.bin_center_Ma == row.bin_center_Ma)]
        assert math.isclose(row.revised_EV, subset.IQR.mean(), rel_tol=0, abs_tol=1e-12)
    for _, row in proxy_ev.iterrows():
        subset = iqrs[(iqrs.scheme.astype(str) == str(row.scheme)) & (iqrs.bin_center_Ma == row.bin_center_Ma) & (iqrs.category == row.category)]
        assert math.isclose(row.proxy_EV, subset.IQR.mean(), rel_tol=0, abs_tol=1e-12)

    peak_rows = global_ev.loc[global_ev.groupby("scheme").revised_EV.idxmax(), ["scheme", "bin_width_kyr", "alignment", "bin_center_Ma", "revised_EV"]]
    summary = {
        "analysis_date": "2026-08-17", "canonical_geographic_sites": 32,
        "records_total_including_excluded": 75, "retained_operational_records": 70,
        "retained_measurements": 10945, "new_measurements": 814,
        "new_series_counts": series_summary[["series", "processed_measurements"]].to_dict("records"),
        "global_EV_peaks": peak_rows.to_dict("records"),
        "standardization": "record-wide z score using sample SD (ddof=1)",
        "record_bin_IQR": "Q75-Q25 using numpy linear quantiles; n>=2",
        "global_EV": "arithmetic mean of retained contributing record-level IQRs",
        "technical_replicates": "arithmetic mean by sample/age; Os mean error propagated as RSS/n",
        "composition_control_run": False, "biotic_analysis_run": False, "monte_carlo_run": False,
    }
    (RESULTS / "new_proxy_integration_v3_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    (LOGS / "new_proxy_integration_v3_provenance.json").write_text(json.dumps({**summary, "source_hashes": source_hashes.to_dict("records")}, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    raise SystemExit("Preserved numerical library: run ../run.py for the packaged final workflow.")
