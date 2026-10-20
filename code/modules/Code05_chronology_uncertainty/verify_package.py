#!/usr/bin/env python3
"""Recompute a tiny raw-input prefix using the retained executed kernels."""
from __future__ import annotations

import argparse
import ast
import dataclasses
import importlib.util
import json
import sys
import unicodedata
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
CODE = ROOT / "02_code_revised"
REFERENCE = ROOT / "reference_outputs"
OUT = ROOT / "verification"


def load(name):
    path = CODE / (name + ".py")
    spec = importlib.util.spec_from_file_location("package_" + name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def raw_data():
    raw = pd.read_csv(ROOT / "01_data_revised" / "retained_measurements_with_zscores_v4.csv").reset_index(drop=True)
    raw["record_id"] = raw.record_id.map(lambda value: unicodedata.normalize("NFC", str(value)).strip())
    return raw


def row_codes(raw):
    lookup = {rid: i for i, rid in enumerate(pd.unique(raw.record_id))}
    return raw.record_id.map(lookup).to_numpy(int)


def output_check(case, ev, ages, observed, archive_key, q_columns, contrasts, contrast_key):
    archive = np.load(REFERENCE / case / "realizations.npz")
    reference_ev = archive[archive_key]
    difference = float(np.max(np.abs(ev - reference_ev[:len(ev)])))
    if difference > 1e-12:
        raise RuntimeError(f"{case} raw recomputation differs from archived prefix: {difference}")
    summary = pd.read_csv(REFERENCE / case / "summary.csv").sort_values("midpoint_Ma")
    q = np.quantile(reference_ev, [0.025, 0.5, 0.975], axis=0)
    sorted_age = np.argsort(ages)
    summary_difference = float(np.max(np.abs(summary[list(q_columns)].to_numpy(float).T - q[:, sorted_age])))
    observed_difference = float(np.max(np.abs(summary.observed_EV.to_numpy(float) - observed[sorted_age])))
    if max(summary_difference, observed_difference) > 1e-12:
        raise RuntimeError(f"{case} archived summary is inconsistent with realization archive/baseline")
    contrast_difference = float(np.max(np.abs(contrasts - archive[contrast_key][:len(ev)])))
    if contrast_difference > 1e-12:
        raise RuntimeError(f"{case} raw recomputation phase contrasts differ from archived prefix")
    pd.DataFrame(ev, columns=[f"{age:.8f}" for age in ages]).to_csv(OUT / f"{case}_smoke_EV.csv", index=False)
    return {
        "case": case, "smoke_iterations": len(ev), "formal_archive_iterations": len(reference_ev),
        "raw_prefix_max_absolute_EV_difference": difference,
        "raw_prefix_max_absolute_phase_contrast_difference": contrast_difference,
        "archived_10000_quantile_summary_max_absolute_difference": summary_difference,
        "observed_baseline_summary_max_absolute_difference": observed_difference,
        "all_smoke_EV_finite": bool(np.isfinite(ev).all()),
    }


def export_maps(raw, maps, name):
    frame = pd.DataFrame([dataclasses.asdict(item) for item in maps])
    frame["record_id"] = [raw.loc[item.row_index, "record_id"] for item in maps]
    frame["retained_age_Ma"] = [raw.loc[item.row_index, "age_Ma"] for item in maps]
    frame["retained_value"] = [raw.loc[item.row_index, "value"] for item in maps]
    frame.to_csv(OUT / f"{name}_exact_observation_parameters.csv", index=False)
    return {"active_records": int(frame.record_id.nunique()), "active_measurements": len(frame)}


def measurement(n):
    raw = raw_data()
    iqr = pd.read_csv(ROOT / "01_data_revised" / "record_level_IQR_all_bin_widths_revised_v4_n2.csv")
    iqr["record_id"] = iqr.record_id.map(lambda value: unicodedata.normalize("NFC", str(value)).strip())
    codes = row_codes(raw)
    l1 = load("run_layer1_measurement_uncertainty_pilot")
    l1.validate_layer1_parameter_membership()
    sources = l1.ensure_sources()
    mappings = []
    notes = {}
    for mapper in (l1.map_temperature_sources, l1.map_jouini_sources, l1.map_os_sources):
        mm, nn = mapper(raw, sources)
        mappings.extend(mm)
        notes.update(nn)
    m1 = load("run_layer1A_documented_measurement_uncertainty_formal_50kyr_10000")
    masks1 = m1.build_50kyr_masks(l1, raw, iqr, codes)
    m1.ITERATIONS = n
    ev1, contrasts1, diagnostics1 = m1.simulate_batched(l1, raw, mappings, masks1)
    result1 = output_check("documented", ev1, masks1["50"]["centers"], masks1["50"]["observed_ev"], "EV_50kyr", ("MC_2.5th", "MC_median", "MC_97.5th"), contrasts1, "phase_contrasts_50kyr")
    result1.update(export_maps(raw, mappings, "documented"))
    result1["sample_SD_ddof"] = 1
    result1["observed_record_bin_mask_fixed"] = True
    result1["baseline_vs_Dataset_S1_max_difference"] = l1.validate_s1_50(masks1)
    result1["diagnostics"] = diagnostics1
    assert (result1["active_records"], result1["active_measurements"]) == (11, 375)

    l2 = load("run_layer2_conservative_measurement_sensitivity_pilot")
    params = l2.parameter_table()
    maps1, notes1, _ = l2.build_layer1_maps(raw)
    maps2, audit2 = l2.tier2_maps(raw, params)
    maps = maps1 + maps2
    masks2 = l2.build_masks(raw, iqr, codes)
    m2 = load("run_layer2_conservative_measurement_sensitivity_formal_50kyr_10000")
    m2.ITERATIONS = n
    ev2, contrasts2, diagnostics2 = m2.simulate_batched(raw, maps, masks2)
    result2 = output_check("conservative", ev2, masks2["50"]["centers"], masks2["50"]["observed_ev"], "EV_50kyr", ("MC_2.5th", "MC_median", "MC_97.5th"), contrasts2, "phase_contrasts_50kyr")
    result2.update(export_maps(raw, maps, "conservative"))
    result2["diagnostics"] = diagnostics2
    result2["observed_record_bin_mask_fixed"] = True
    audit = m2.build_execution_audit(raw, params, maps1, notes1, audit2)
    audit.to_csv(OUT / "measurement_record_execution_audit.csv", index=False)
    assert (result2["active_records"], result2["active_measurements"]) == (53, 10015)
    return [result1, result2]


def chronology(n):
    raw = raw_data()
    primary = load("run_primary_chronology_sensitivity_pilot_50kyr_1000")
    membership, groups = primary.load_grouping_audit()
    members = primary.validate_grouping(raw, membership, groups)
    arrays = primary.prepare_record_arrays(raw)
    ages, observed, observed_counts, s1_difference, zero_bins = primary.validate_observed_baseline(raw, arrays)
    primary.ITERATIONS = n
    shifts, stream_info = primary.draw_group_shifts()
    active = set(membership.operational_record)
    sums, counts, _ = primary.record_contributions(arrays, [rid for rid in arrays if rid not in active], 0.0)
    ev = []
    for i in range(n):
        ss, cc = sums.copy(), counts.copy()
        for gid in primary.GROUP_ORDER:
            s, c, _ = primary.record_contributions(arrays, members[gid], float(shifts[gid][i] / 1000.0))
            ss += s
            cc += c
        ev.append(ss / cc)
    result1 = output_check("primary", np.array(ev), ages, observed, "primary_EV_50kyr", ("simulation_2.5th", "simulation_median", "simulation_97.5th"), primary.phase_contrasts(np.array(ev), ages), "primary_phase_contrasts")
    result1.update({"active_records": 24, "active_measurements": 4230, "groups": primary.EXPECTED_GROUPS, "RNG_streams": stream_info, "dynamic_n_ge_2_and_nonzero_IQR": True, "baseline_vs_Dataset_S1_max_difference": s1_difference})
    expanded_core = load("run_primary_chronology_sensitivity_formal_50kyr_10000")
    expanded = expanded_core.read_expanded_membership(primary, raw)
    result1["expanded_ODP1210_records"] = expanded.operational_record.tolist()
    result1["expanded_ODP1210_measurements"] = int(expanded.n_measurements.sum())
    result1["expanded_ODP1210_rng"] = "shared CHR-P1-ODP1209 shift; no additional stream"

    layerb = load("run_chronology_layerB_conservative_sensitivity_pilot_50kyr_1000")
    records, bgroups, bmembers = layerb.load_parameter_table(raw)
    layerb.ITERATIONS = n
    paired = {"ODP1209_shift_kyr": shifts[primary.GROUP_ORDER[0]], "U1403_shift_kyr": shifts[primary.GROUP_ORDER[1]], "Seymour_shift_kyr": shifts[primary.GROUP_ORDER[2]]}
    bshifts, bstream_info = layerb.draw_group_shifts(bgroups, paired)
    active_b = set(records.loc[records.active_flag.eq(1), "record_id"])
    sums, counts, _ = primary.record_contributions(arrays, [rid for rid in arrays if rid not in active_b], 0.0)
    bev = []
    for i in range(n):
        ss, cc = sums.copy(), counts.copy()
        for gid in layerb.GROUP_ORDER:
            if gid == layerb.SEGMENT_GROUP:
                s, c, _ = layerb.segment_record_contributions(primary, arrays, bmembers[gid], float(bshifts[gid][i] / 1000.0))
            else:
                s, c, _ = primary.record_contributions(arrays, bmembers[gid], float(bshifts[gid][i] / 1000.0))
            ss += s
            cc += c
        bev.append(ss / cc)
    result2 = output_check("layerB", np.array(bev), ages, observed, "LayerB_EV_50kyr", ("simulation_2.5th", "simulation_median", "simulation_97.5th"), primary.phase_contrasts(np.array(bev), ages), "LayerB_phase_contrasts")
    result2.update({"active_records": 38, "active_measurements_in_active_records": 8843, "groups": bgroups.to_dict("records"), "RNG_streams": bstream_info, "dynamic_n_ge_2_and_nonzero_IQR": True, "ODP690_post_boundary_Ma": layerb.KPG_BOUNDARY_MA, "ODP690_shifted_rows_per_realization": int(sum(np.sum(arrays[rid][0] < layerb.KPG_BOUNDARY_MA) for rid in bmembers[layerb.SEGMENT_GROUP]))})
    membership.to_csv(OUT / "primary_chronology_exact_membership.csv", index=False)
    records.to_csv(OUT / "layerB_chronology_exact_record_parameters.csv", index=False)
    bgroups.to_csv(OUT / "layerB_chronology_exact_group_parameters.csv", index=False)
    return [result1, result2]


def joint(n):
    raw = raw_data()
    engine = load("run_joint_conservative_measurement_chronology_sensitivity_pilot_50kyr_1000")
    m = engine.load_module(engine.MEASUREMENT_PILOT_SCRIPT, "verify_joint_measurement")
    mf = engine.load_module(engine.MEASUREMENT_FORMAL_SCRIPT, "verify_joint_measurement_formal")
    c = engine.load_module(engine.CHRONOLOGY_PILOT_SCRIPT, "verify_joint_layerB")
    p = engine.load_module(engine.PRIMARY_CHRONOLOGY_SCRIPT, "verify_joint_primary")
    params = m.parameter_table()
    maps1, notes1, _ = m.build_layer1_maps(raw)
    maps2, audit2 = m.tier2_maps(raw, params)
    maps = maps1 + maps2
    audit, participation, structure = engine.measurement_execution_audit(raw, params, maps1, notes1, audit2, mf)
    draw_arrays = engine.prepare_measurement_draw_arrays(maps)
    records, groups, members = c.load_parameter_table(raw)
    row_masks = engine.chronology_row_masks(raw, members, c.GROUP_ORDER, c.SEGMENT_GROUP)
    master = np.random.SeedSequence(engine.MASTER_SEED)
    measurement_sequence, chronology_sequence = master.spawn(2)
    rng = np.random.default_rng(measurement_sequence)
    engine.ITERATIONS = n
    shifts, stream_info = engine.draw_chronology_shifts(groups, c.GROUP_ORDER, chronology_sequence)
    codes = row_codes(raw)
    indices = [np.flatnonzero(codes == code) for code in np.unique(codes)]
    arrays = p.prepare_record_arrays(raw)
    ages, observed, observed_counts, baseline_difference, _ = p.validate_observed_baseline(raw, arrays)
    ev, eligibility = [], []
    for i in range(n):
        values = engine.perturb_measurements(raw.value.to_numpy(float), draw_arrays, rng)
        z, invalid_sd = engine.compute_recordwise_z(values, indices)
        assert invalid_sd == 0
        shifted_ages = engine.perturb_ages(raw.age_Ma.to_numpy(float), shifts, i, c.GROUP_ORDER, row_masks)
        trajectory, counts, _ = engine.dynamic_ev(p, z, shifted_ages, indices)
        ev.append(trajectory)
        eligibility.append(counts)
    result = output_check("jointBB", np.array(ev), ages, observed, "EV_50kyr", ("simulation_2.5th", "simulation_median", "simulation_97.5th"), np.array([engine.phase_contrasts(trajectory, ages) for trajectory in ev]), "phase_contrasts")
    reference = np.load(REFERENCE / "jointBB" / "realizations.npz")
    eligibility_difference = int(np.max(np.abs(np.array(eligibility) - reference["eligible_record_bin_counts"][:n])))
    assert eligibility_difference == 0
    result.update({"measurement_scope": structure, "measurement_active_records": 53, "measurement_active_measurements": 10015, "chronology_active_records": 38, "chronology_active_record_measurements": 8843, "chronology_groups": groups.to_dict("records"), "master_seed": engine.MASTER_SEED, "measurement_spawn_key": list(measurement_sequence.spawn_key), "chronology_RNG_streams": stream_info.to_dict("records"), "smoke_eligibility_max_difference": eligibility_difference, "baseline_vs_Dataset_S1_max_difference": baseline_difference, "dynamic_n_ge_2_and_nonzero_IQR": True})
    export_maps(raw, maps, "joint_measurement")
    audit.to_csv(OUT / "joint_measurement_record_execution_audit.csv", index=False)
    records.to_csv(OUT / "joint_chronology_record_parameters.csv", index=False)
    groups.to_csv(OUT / "joint_chronology_group_parameters.csv", index=False)
    return [result]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--iterations", type=int, default=3)
    args = parser.parse_args()
    if not 1 <= args.iterations <= 20:
        parser.error("This verification command allows only 1–20 smoke iterations; use run_analysis.py for formal runs.")
    OUT.mkdir(exist_ok=True)
    syntax = []
    for path in sorted(CODE.glob("*.py")):
        ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        syntax.append(path.name)
    code = ROOT.name[:6]
    results = {"Code04": measurement, "Code05": chronology, "Code06": joint}[code](args.iterations)
    report = {"code": code, "verification_level": "syntax plus tiny raw-input recomputation and archived formal output consistency", "full_10000_simulations_rerun": False, "Python_sources_parsed": syntax, "results": results}
    (OUT / "verification_report.json").write_text(json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(report, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
