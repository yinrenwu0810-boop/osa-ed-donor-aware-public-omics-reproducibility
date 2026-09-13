#!/usr/bin/env python3
"""Read perturbation labels from an H5AD without loading the expression matrix."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import h5py
import numpy as np


def decode(value: object) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    return str(value)


def values(group: h5py.Group, key: str) -> list[str]:
    node = group[key]
    if isinstance(node, h5py.Dataset):
        return [decode(item) for item in node[()]]
    if "categories" in node and "codes" in node:
        categories = [decode(item) for item in node["categories"][()]]
        codes = node["codes"][()]
        return [categories[int(code)] if int(code) >= 0 else "NOT_REPORTED" for code in codes]
    raise ValueError(f"Unsupported H5AD encoding for obs/{key}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("h5ad", type=Path)
    parser.add_argument("--targets", nargs="+", default=["EFNB2", "LRRC17", "TYMS"])
    args = parser.parse_args()

    with h5py.File(args.h5ad, "r") as handle:
        obs = handle["obs"]
        keys = sorted(obs.keys())
        condition_key = next((key for key in ["condition", "perturbation", "gene", "target"] if key in obs), None)
        if condition_key is None:
            raise SystemExit(f"No supported perturbation column; obs keys={keys}")
        labels = values(obs, condition_key)
        unique, counts = np.unique(np.asarray(labels, dtype=object), return_counts=True)
        matches_by_column = {}
        for candidate_key in ["condition", "condition_name", "perturbation", "gene", "target"]:
            if candidate_key not in obs:
                continue
            candidate_values = values(obs, candidate_key)
            candidate_unique, candidate_counts = np.unique(np.asarray(candidate_values, dtype=object), return_counts=True)
            matches = {}
            for target in args.targets:
                exact = []
                for label, count in zip(candidate_unique.tolist(), candidate_counts.tolist()):
                    normalized = str(label).replace("_", "+").replace("-", "+")
                    parts = [part.strip() for part in normalized.split("+")]
                    if target in parts or target in str(label):
                        exact.append({"label": str(label), "cells": int(count)})
                matches[target] = exact
            matches_by_column[candidate_key] = matches
        selected_obs = {}
        for key in ["cell_type", "condition", "condition_name", "perturbation"]:
            if key in obs:
                vals = values(obs, key)
                vals_unique, vals_counts = np.unique(np.asarray(vals, dtype=object), return_counts=True)
                if len(vals_unique) <= 25:
                    selected_obs[key] = [
                        {"value": str(value), "cells": int(count)}
                        for value, count in zip(vals_unique.tolist(), vals_counts.tolist())
                    ]
                else:
                    selected_obs[key] = {"unique_values": int(len(vals_unique))}
        result = {
            "path": str(args.h5ad),
            "obs_keys": keys,
            "condition_key": condition_key,
            "cells": len(labels),
            "unique_conditions": int(len(unique)),
            "target_matches_by_column": matches_by_column,
            "condition_examples": [str(item) for item in unique[:20].tolist()],
            "selected_obs_summary": selected_obs,
        }
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
