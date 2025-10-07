import sys
from pathlib import Path
import pandas as pd
import numpy as np

# Columns to display (subset is taken based on availability in the input CSV)
DISPLAY_COLS = [
    "model","KO","carbon","carbon_uptake","o2_lb",
    "pH_lb","pH_label","ion","ion_exch",
    "temp_factor","temp_label","growth_frac",
    "WT_growth","growth_max","citrate_FBA",
    "growth_at_citrate","siderophore_FBA"
]

# Columns that should be numeric
NUMERIC_COLS = [
    "pH_lb","temp_factor","growth_frac","WT_growth","growth_max",
    "citrate_FBA","growth_at_citrate","siderophore_FBA",
    "o2_lb","carbon_uptake"
]

SECTIONS = [
    ("growth_max",        "Top 10 — Maximum growth"),
    ("citrate_FBA",       "Top 10 — Citrate export (FBA)"),
    ("siderophore_FBA",   "Top 10 — Production of siderophore (FBA)"),
]

def coerce_numeric(df: pd.DataFrame) -> pd.DataFrame:
    """
    Ensure selected columns are numeric.

    Purpose
    -------
    Converts values in NUMERIC_COLS to numeric dtype; non-parsable values
    are coerced to NaN (errors='coerce') to avoid sorting/type issues later.

    Parameters
    ----------
    df : pd.DataFrame
        Input dataframe loaded from the simulation CSV.

    Returns
    -------
    pd.DataFrame
        Same dataframe object with columns in NUMERIC_COLS converted to numeric
        where present. Operation is in-place on those columns (copy semantics
        depend on pandas; callers commonly pass a copy).
    """
    for col in NUMERIC_COLS:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df

def top_n(df: pd.DataFrame, metric: str, n: int = 10) -> pd.DataFrame:
    """
    Compute the top-N rows for a given metric, with optional tie-breaker.

    Purpose
    -------
    Returns the N highest rows for `metric` (descending). If the column
    'WT_growth' exists, it is used as a secondary sort key (also descending)
    to break ties in a biologically meaningful way.

    Parameters
    ----------
    df : pd.DataFrame
        Input dataframe.
    metric : str
        Column name to rank by (e.g., 'growth_max', 'citrate_FBA').
    n : int, default=10
        Number of top rows to return.

    Returns
    -------
    pd.DataFrame
        A dataframe with up to N rows, sorted by the specified metric (and
        'WT_growth' if present). Empty dataframe is returned if the metric
        is missing or if all values are NaN.
    """
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
    """
    Pretty-print a titled section to stdout.

    Purpose
    -------
    Prints a formatted header and then a tabular view of the dataframe
    restricted to DISPLAY_COLS (if available), or a friendly message if
    the table is empty.

    Parameters
    ----------
    title : str
        Section title to print before the table.
    table : pd.DataFrame
        Dataframe to render; commonly produced by `top_n` or `env_best_summary`.

    Returns
    -------
    None
    """
    print("=" * 100)
    print(title)
    print("=" * 100)
    if table.empty:
        print("No data available for this metric.\n")
        return
    cols = [c for c in DISPLAY_COLS if c in table.columns]
    print(table[cols].to_string(index=False))
    print()

def env_best_summary(df: pd.DataFrame, metric: str, top_k: int = 10):
    """
    Print best cases by environment tuple (carbon, o2_lb, carbon_uptake).

    Purpose
    -------
    For a given metric, finds the best-performing row within each unique
    environment bucket defined by (carbon, o2_lb, carbon_uptake), then
    sorts those maxima descending by the metric and prints the top_k rows.

    Parameters
    ----------
    df : pd.DataFrame
        Input dataframe containing results across many environments.
    metric : str
        Column name to evaluate (e.g., 'citrate_FBA', 'siderophore_FBA').
    top_k : int, default=10
        Number of best environment tuples to display.

    Returns
    -------
    None
        (Side-effect: prints a titled section to stdout.)
    """
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

    title = f"Best {metric} by (carbon, o2_lb, carbon_uptake) — Top {top_k}"
    print_section(title, best)

def main():
    """
    Entrypoint: read CSV, compute and print ranked sections and env summaries.

    Behavior
    --------
    1) Parse CLI arguments: path to results CSV and optional N_top (default: 10).
    2) Load the CSV into a dataframe.
    3) For each (metric, title) in SECTIONS, print the top-N table.
    4) Print environment-best summaries for 'citrate_FBA' and 'siderophore_FBA'.

    CLI
    ---
    Usage: python summarize_results.py <results.csv> [N_top]

    Returns
    -------
    None
        (Process exits with status 1 on usage error or missing file.)
    """
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <results.csv> [N_top]")
        sys.exit(1)

    csv_path = Path(sys.argv[1])
    if not csv_path.exists():
        print(f"Error : file not found : {csv_path}")
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
