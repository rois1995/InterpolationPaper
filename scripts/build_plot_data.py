#!/usr/bin/env python3
"""Build plot-ready CSV files and convergence-order summaries.

The script reads the raw comparison files produced by the interpolation test
suite and writes compact, wide tables consumed directly by pgfplots.  It also
recomputes convergence orders from the error data, so the paper does not rely
on manually transcribed values.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

import numpy as np
import pandas as pd

METHOD_MAP = {
    "NN": "NN",
    "Barycentric": "BAR",
    "ConsGalerkinProj": "GP",
}
ORDER_MAP = {
    "FirstOrder": "1",
    "SecondOrder": "2",
    "ThirdOrder": "3",
}
SCHEME_ORDER = ["NN1", "NN2", "NN3", "BAR2", "BAR3", "GP2", "GP3"]
NOMINAL_ORDER = {
    "NN1": 1,
    "NN2": 2,
    "NN3": 3,
    "BAR2": 2,
    "BAR3": 3,
    "GP2": 2,
    "GP3": 3,
}


def scheme_name(method: str, interp_order: str) -> str:
    try:
        return METHOD_MAP[method] + ORDER_MAP[interp_order]
    except KeyError as exc:
        raise ValueError(f"Unsupported method/order combination: {method}/{interp_order}") from exc


def require_columns(df: pd.DataFrame, columns: list[str], filename: Path) -> None:
    missing = [column for column in columns if column not in df.columns]
    if missing:
        raise ValueError(f"{filename} is missing columns: {', '.join(missing)}")


def fit_order(h: np.ndarray, error: np.ndarray) -> float:
    mask = np.isfinite(h) & np.isfinite(error) & (h > 0.0) & (error > 0.0)
    if mask.sum() < 2:
        return float("nan")
    return float(np.polyfit(np.log(h[mask]), np.log(error[mask]), 1)[0])


def fine_order(h: np.ndarray, error: np.ndarray) -> float:
    mask = np.isfinite(h) & np.isfinite(error) & (h > 0.0) & (error > 0.0)
    h = h[mask]
    error = error[mask]
    if len(h) < 2:
        return float("nan")
    idx = np.argsort(h)[::-1]  # coarse -> fine
    h = h[idx]
    error = error[idx]
    return float(np.log(error[-2] / error[-1]) / np.log(h[-2] / h[-1]))


def make_wide_from_long(
    df: pd.DataFrame,
    *,
    solution: str,
    field: str,
    value_column: str,
) -> pd.DataFrame:
    sub = df[(df["solution"] == solution) & (df["field"] == field)].copy()
    if sub.empty:
        raise ValueError(f"No data for solution={solution}, field={field}")
    sub["scheme"] = [scheme_name(m, o) for m, o in zip(sub["method"], sub["interp_order"])]

    resolution_sets = {
        grid: set(group["resolution"].astype(int).tolist()) for grid, group in sub.groupby("grid")
    }
    if len({tuple(sorted(values)) for values in resolution_sets.values()}) != 1:
        raise ValueError(f"Structured/mixed resolution sets differ: {resolution_sets}")
    resolutions = sorted(next(iter(resolution_sets.values())))

    out = pd.DataFrame({"N": resolutions})
    for grid in ["Struct", "Mixed"]:
        g = sub[sub["grid"] == grid]
        h = g.groupby("resolution", as_index=True)["h"].first().reindex(resolutions)
        out[f"h_{grid}"] = h.to_numpy()
        for scheme in SCHEME_ORDER:
            values = (
                g[g["scheme"] == scheme]
                .set_index("resolution")[value_column]
                .reindex(resolutions)
            )
            if values.isna().any():
                raise ValueError(f"Missing {solution}/{field}/{grid}/{scheme}/{value_column} data")
            out[f"{grid}_{scheme}"] = values.to_numpy()
    return out


def make_order_summary(wide: pd.DataFrame, expected: dict[str, float]) -> pd.DataFrame:
    rows: list[dict[str, float | str]] = []
    for scheme in SCHEME_ORDER:
        row: dict[str, float | str] = {"scheme": scheme, "nominal": expected[scheme]}
        for grid in ["Struct", "Mixed"]:
            h = wide[f"h_{grid}"].to_numpy(dtype=float)
            e = wide[f"{grid}_{scheme}"].to_numpy(dtype=float)
            row[f"{grid}_fit"] = fit_order(h, e)
            row[f"{grid}_fine"] = fine_order(h, e)
        rows.append(row)
    return pd.DataFrame(rows)


def make_norm_order_summary(errors: pd.DataFrame, solution: str) -> pd.DataFrame:
    sub = errors[(errors["solution"] == solution) & (errors["field"] == "Density")].copy()
    sub["scheme"] = [scheme_name(m, o) for m, o in zip(sub["method"], sub["interp_order"])]
    rows = []
    for scheme in SCHEME_ORDER:
        row: dict[str, float | str] = {"scheme": scheme}
        for grid in ["Struct", "Mixed"]:
            g = sub[(sub["grid"] == grid) & (sub["scheme"] == scheme)].sort_values("h", ascending=False)
            for norm in ["rel_L1", "rel_L2", "rel_Linf"]:
                row[f"{grid}_{norm}"] = fit_order(
                    g["h"].to_numpy(dtype=float), g[norm].to_numpy(dtype=float)
                )
        rows.append(row)
    return pd.DataFrame(rows)


def make_boundedness_plot_data(df: pd.DataFrame) -> pd.DataFrame:
    order_short = {"First": "1", "Second": "2", "Third": "3"}
    method_short = {"NN": "NN", "Barycentric": "BAR", "ConsGalerkinProj": "GP"}
    labels = []
    for row in df.itertuples(index=False):
        base = method_short[row.method] + order_short[row.order]
        if row.method == "ConsGalerkinProj":
            suffix = {"soft": "S", "hard": "H", "none": "U"}.get(str(row.limiter), str(row.limiter))
        else:
            if row.limiter == "none":
                suffix = "U"
            else:
                suffix = str(row.limiter)
            if str(row.enforce_cons) == "Global":
                suffix += "-GC"
        labels.append(f"{base}-{suffix}")
    out = df.copy()
    out.insert(0, "idx", np.arange(len(out), dtype=int))
    out.insert(1, "plot_label", labels)
    out["abs_cons_defect"] = out["cons_defect"].abs()
    return out[[
        "idx", "plot_label", "method", "order", "limiter", "enforce_cons",
        "overshoot", "undershoot", "abs_cons_defect", "rho_min", "rho_max"
    ]]


def verify_comparison_orders(errors: pd.DataFrame, supplied: pd.DataFrame) -> float:
    keys = ["grid", "dim", "solution", "interp_order", "method", "field", "resolution"]
    supplied_map = supplied.set_index(keys)["order_L2"]
    max_diff = 0.0
    for _, group in errors.groupby(["grid", "dim", "solution", "interp_order", "method", "field"]):
        group = group.sort_values("h", ascending=False)
        previous = None
        for row in group.itertuples(index=False):
            key = tuple(getattr(row, name) for name in keys)
            expected = float("nan")
            if previous is not None:
                expected = np.log(previous.rel_L2 / row.rel_L2) / np.log(previous.h / row.h)
            provided = supplied_map.get(key, np.nan)
            if np.isfinite(expected) and np.isfinite(provided):
                max_diff = max(max_diff, abs(float(expected) - float(provided)))
            previous = row
    return max_diff


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, default=Path("data"))
    parser.add_argument("--out-dir", type=Path, default=Path("generated"))
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    errors_path = args.data_dir / "Comparison_Errors.csv"
    cons_path = args.data_dir / "Comparison_Conservation.csv"
    supplied_orders_path = args.data_dir / "Comparison_Orders.csv"
    disc_l2_path = args.data_dir / "discontinuous_convergence_L2.csv"
    disc_cons_path = args.data_dir / "discontinuous_conservation_defect.csv"
    bounded_path = args.data_dir / "discontinuous_boundedness.csv"

    errors = pd.read_csv(errors_path)
    cons = pd.read_csv(cons_path)
    supplied_orders = pd.read_csv(supplied_orders_path)
    disc_l2 = pd.read_csv(disc_l2_path)
    disc_cons = pd.read_csv(disc_cons_path)
    bounded = pd.read_csv(bounded_path)

    require_columns(errors, ["grid", "solution", "method", "interp_order", "resolution", "h", "field", "rel_L1", "rel_L2", "rel_Linf"], errors_path)
    require_columns(cons, ["grid", "solution", "method", "interp_order", "resolution", "h", "field", "geom_defect"], cons_path)

    asym_l2 = make_wide_from_long(errors, solution="Asym", field="Density", value_column="rel_L2")
    asym_cons = make_wide_from_long(cons, solution="Asym", field="Density", value_column="geom_defect")
    asym_orders = make_order_summary(asym_l2, NOMINAL_ORDER)
    asym_norm_orders = make_norm_order_summary(errors, "Asym")

    # The discontinuous files are already in the desired wide layout.  Preserve
    # the data exactly and compute order summaries from them.
    source_suffix = {
        "NN1": "NN1", "NN2": "NN2", "NN3": "NN3",
        "BAR2": "Bary2", "BAR3": "Bary3",
        "GP2": "Gal2", "GP3": "Gal3",
    }
    required_wide = ["N", "h_Struct", "h_Mixed"] + [
        f"{grid}_{source_suffix[scheme]}"
        for grid in ["Struct", "Mixed"] for scheme in SCHEME_ORDER
    ]
    require_columns(disc_l2, required_wide, disc_l2_path)
    require_columns(disc_cons, required_wide, disc_cons_path)

    # Normalize Bary/Gal column names to the BAR/GP notation used in the paper.
    rename_columns = {
        f"{grid}_{source_suffix[scheme]}": f"{grid}_{scheme}"
        for grid in ["Struct", "Mixed"] for scheme in SCHEME_ORDER
    }
    disc_l2_normalized = disc_l2.rename(columns=rename_columns)
    disc_cons_normalized = disc_cons.rename(columns=rename_columns)
    disc_orders = make_order_summary(disc_l2_normalized, {scheme: 0.5 for scheme in SCHEME_ORDER})
    disc_norm_orders = make_norm_order_summary(errors, "Discontinuous")
    bounded_plot = make_boundedness_plot_data(bounded)

    outputs = {
        "asymmetric_convergence_L2.csv": asym_l2,
        "asymmetric_conservation_defect.csv": asym_cons,
        "asymmetric_orders.csv": asym_orders,
        "asymmetric_norm_orders.csv": asym_norm_orders,
        "discontinuous_convergence_L2.csv": disc_l2_normalized,
        "discontinuous_conservation_defect.csv": disc_cons_normalized,
        "discontinuous_orders.csv": disc_orders,
        "discontinuous_norm_orders.csv": disc_norm_orders,
        "discontinuous_boundedness_plot.csv": bounded_plot,
    }
    for name, frame in outputs.items():
        frame.to_csv(args.out_dir / name, index=False, float_format="%.16g")

    max_order_difference = verify_comparison_orders(errors, supplied_orders)
    report = [
        "Plot-data generation completed.",
        f"Maximum difference between recomputed and supplied pairwise L2 orders: {max_order_difference:.3e}",
        "Generated files:",
    ]
    report.extend(f"  - {name}" for name in outputs)
    (args.out_dir / "generation_report.txt").write_text("\n".join(report) + "\n", encoding="utf-8")
    print("\n".join(report))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
