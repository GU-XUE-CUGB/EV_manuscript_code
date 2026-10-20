#!/usr/bin/env python3
"""Run the packaged formal analysis, or a tiny verification prefix."""
from __future__ import annotations
import argparse
import importlib.util
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--smoke", type=int, metavar="N", help="Verify 1–20 raw-input realizations instead of running the formal analysis")
    parser.add_argument("--scenario", default="all", help="Run one configured formal scenario, or all in dependency order")
    args = parser.parse_args()
    if args.smoke is not None:
        subprocess.run([sys.executable, str(ROOT / "verify_package.py"), "--iterations", str(args.smoke)], check=True)
        return
    config = json.loads((ROOT / "package_config.json").read_text(encoding="utf-8"))
    scenarios = config["formal_scenarios"]
    if args.scenario != "all":
        scenarios = [item for item in scenarios if item["name"] == args.scenario]
        if not scenarios:
            parser.error("Unknown scenario. Choose all or " + ", ".join(item["name"] for item in config["formal_scenarios"]))
    for item in scenarios:
        script = ROOT / config["code_directory"] / item["source"]
        print(f"Running {item['name']}: formal 10,000 realizations", flush=True)
        subprocess.run([sys.executable, str(script)], cwd=ROOT, check=True)


if __name__ == "__main__":
    main()
