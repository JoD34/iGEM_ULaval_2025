import sys
from pathlib import Path
import pandas as pd
import numpy as np

DISPLAY_COLS = [
    "model","KO","carbon","carbon_uptake","o2_lb",
    "pH_lb","pH_label","ion","ion_exch",
    "temp_factor","temp_label","growth_frac",
    "WT_growth","growth_max","citrate_FBA",
    "growth_at_citrate","siderophore_FBA"
]

NUMERIC_COLS = [
    "pH_lb","temp_factor","growth_frac","WT_growth","growth_max",
    "citrate_FBA","growth_at_citrate","siderophore_FBA",
    "o2_lb","carbon_uptake"
]

SECTIONS = [
    ("growth_max",        "Top 10 — Croissance maximale"),
    ("citrate_FBA",       "Top 10 — Export de citrate (FBA)"),
    ("siderophore_FBA",   "Top 10 — Production de sidérophore (FBA)"),
]

def coerce_numeric(df: pd.DataFrame) -> pd.DataFrame:
    for col in NUMERIC_COLS:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df

def top_n(df: pd.DataFrame, metric: str, n: int = 10) -> pd.DataFrame:
    if metric not in df.columns:
        return pd.DataFrame()
    sub = df.copy()
    sub = coerce_numeric(sub)
    sub = sub[~sub[metric].isna()]
    if sub.empty:
        return sub
    sort_cols = [metric]
    ascending = [False]
    if "WT_growth" in sub.columns:
        sort_cols += ["WT_growth"]
        ascending += [False]
    return sub.sort_values(sort_cols, ascending=ascending).head(n)

def print_section(title: str, table: pd.DataFrame):
    print("=" * 100)
    print(title)
    print("=" * 100)
    if table.empty:
        print("Aucune donnée disponible pour cette métrique.\n")
        return
    cols = [c for c in DISPLAY_COLS if c in table.columns]
    print(table[cols].to_string(index=False))
    print()

def env_best_summary(df: pd.DataFrame, metric: str, top_k: int = 10):
    """Affiche les meilleurs cas par (carbon, o2_lb, carbon_uptake) pour une métrique donnée."""
    needed = {"carbon", "o2_lb", "carbon_uptake", metric}
    if not needed.issubset(df.columns):
        return
    sub = coerce_numeric(df[list(needed | set(DISPLAY_COLS))].copy())
    sub = sub.dropna(subset=[metric])
    if sub.empty:
        return
    idx = sub.groupby(["carbon", "o2_lb", "carbon_uptake"])[metric].idxmax()
    best = sub.loc[idx]
    best = best.sort_values([metric, "carbon", "o2_lb", "carbon_uptake"], ascending=[False, True, True, True]).head(top_k)

    title = f"Meilleurs {metric} par (carbon, o2_lb, carbon_uptake) — Top {top_k}"
    print_section(title, best)

def main():
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <results.csv> [N_top]")
        sys.exit(1)

    csv_path = Path(sys.argv[1])
    if not csv_path.exists():
        print(f"Erreur : fichier introuvable : {csv_path}")
        sys.exit(1)

    try:
        N = int(sys.argv[2]) if len(sys.argv) >= 3 else 10
    except Exception:
        N = 10

    pd.set_option("display.width", 180)
    pd.set_option("display.max_columns", None)

    df = pd.read_csv(csv_path)

    for metric, title in SECTIONS:
        t = top_n(df, metric, n=N)
        print_section(title, t)

    env_best_summary(df, "citrate_FBA", top_k=N)
    env_best_summary(df, "siderophore_FBA", top_k=N)

if __name__ == "__main__":
    main()
