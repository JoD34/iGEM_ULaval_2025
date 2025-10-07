import cobra
from cobra import io as cbio
import pandas as pd
import numpy as np
import warnings
from pathlib import Path
from itertools import product
from multiprocessing import get_context
from multiprocessing import Pool
from fetch_gems import fetch_models_based_on_org
from threadpoolctl import threadpool_limits

warnings.filterwarnings("ignore", category=UserWarning)

MODELS = {}

PH_LEVELS = [
    (-20, "extrême acide"),
    (-10, "acide pH~5"),
    (-1, "neutre pH~7"),
    (0, "basique pH~9"),
    (5, "très basique")
]
CARBON_SOURCES = {
    "glucose": "EX_glc__D_e",
    "acetate": "EX_ac_e",
    "glycerol": "EX_glyc_e",
    "succinate": "EX_succ_e",
    "ethanol": "EX_etoh_e",
}
ION_BUFFERS = {
    "HCO3": ("EX_hco3_e", np.linspace(-20, 0, 9)),
    "Pi":   ("EX_pi_e",   np.linspace(-20, 0, 9)),
    "NH4":  ("EX_nh4_e",  np.linspace(-20, 0, 9)),
    "CO2":  ("EX_co2_e",  np.linspace(-5, 0, 11)),
    "Fe":   ("EX_fe2_e",  np.linspace(-1000, -1, 10))
}

CARBON_UPTAKES = [-2.0, -5.0, -10.0, -20.0]
O2_LBS = [-20.0, -10.0, -5.0, -1.0, 0.0]

TEMP_LEVELS = [
    (0.7, "low (~30 °C)"),
    (0.9, "light (~34 °C)"),
    (1.0, "standard (~37 °C)"),
    (1.1, "moderate (~40 °C)"),
    (1.3, "high (~45 °C)")
]
TEMP_SENSITIVE = ["ACONTa", "ACONTb", "ICDHyr", "SUCDi"]
GROWTH_FRACS = [
    0.0,
    0.05,
    0.1,
    0.2
]
KOS = [None, "ICDHyr", "SUCDi", "AKGDH", "SUCOAS", "ACONTa", "ACONTb", "PPC", "PPCK"]

CSV_OUTPUT = Path('citrate_growth_siderophore_scan.csv')

def flatten_ion_buffers(buffers: dict) -> list:
    """
    Expand ion-buffer configuration into a flat list of tuples.

    Purpose
    -------
    Converts a dict like {"Pi": ("EX_pi_e", np.linspace(...)), ...}
    into a list of (ion_label, exchange_rxn_id, lower_bound_value) tuples,
    one entry per bound value, suitable for cartesian product sweeps.

    Parameters
    ----------
    buffers : dict
        Mapping from ion label (str) to a tuple:
        (exchange_reaction_id : str, iterable_of_lower_bounds : iterable[float])

    Returns
    -------
    list[tuple[str, str, float]]
        Flat list of (ion_label, exchange_rxn_id, lower_bound_value).
    """
    return [
        (ion, rxn_id, buf_lb)
        for ion, (rxn_id, buf_list) in buffers.items()
        for buf_lb in buf_list
    ]


def generate_param_grid(models: dict, carb_sources: dict, ph_levels: list,
                        ion_params: list, temp_levels: list,
                        growth_fracs: list, kos: list,
                        carbon_uptakes: list, o2_lbs: list):
    """
    Build a cartesian product iterator over all simulation parameters.

    Purpose
    -------
    Produces the full combination space of model file entries, media
    composition, physicochemical proxies, genetic perturbations, and
    objective-gating parameters for batch simulation.

    Parameters
    ----------
    models : dict
        Mapping {model_name: PathLike} of model files to evaluate.
    carb_sources : dict
        Mapping {carbon_label: exchange_reaction_id}.
    ph_levels : list[tuple[float, str]]
        [(EX_h_e_lower_bound, human_label), ...] (more negative = more acidic).
    ion_params : list[tuple[str, str, float]]
        Output of `flatten_ion_buffers`: (ion_label, exchange_rxn_id, lower_bound).
    temp_levels : list[tuple[float, str]]
        [(upper_bound_scale_factor, human_label), ...] temperature proxies.
    growth_fracs : list[float]
        Minimum biomass fraction of WT to enforce (0.0 = none).
    kos : list[Optional[str]]
        Reaction IDs to knock out; None means no KO.
    carbon_uptakes : list[float]
        Carbon source uptake lower bounds (mmol·gDW^-1·h^-1, negative for uptake).
    o2_lbs : list[float]
        Oxygen exchange lower bounds (negative uptake → aerobic; 0 → closed).

    Returns
    -------
    itertools.product
        Iterator over tuples in the exact order expected by `run_simulation`.
    """
    return product(
        models.items(),
        carb_sources.items(),
        ph_levels,
        ion_params,
        temp_levels,
        growth_fracs,
        kos,
        carbon_uptakes,
        o2_lbs
    )

def configure_model(m, c_exch, ph_lb, rxn_id, buf_lb,
                    temp_fact, ko, bio_rxn, WT, frac,
                    *, carbon_uptake=-10.0, o2_lb=None):
    """
    Apply environment, physicochemical proxies, and genetic edits to a model.

    Purpose
    -------
    Configures a COBRA model copy before optimization by:
    1) closing all carbon exchanges, then opening the selected carbon source;
    2) setting pH proxy via `EX_h_e` lower bound;
    3) setting an ion/buffer exchange lower bound;
    4) constraining O2 lower bound (aerobic → anoxic);
    5) scaling UBs of temperature-sensitive reactions by `temp_fact`;
    6) applying a reaction KO if requested;
    7) enforcing a minimum biomass flux as a fraction of WT, if provided.

    Parameters
    ----------
    m : cobra.Model
        Model to modify in place.
    c_exch : str
        Carbon source exchange reaction ID to enable.
    ph_lb : float
        Lower bound for `EX_h_e` (more negative = more acidic).
    rxn_id : str
        Exchange reaction ID for the selected ion/buffer.
    buf_lb : float
        Lower bound for `rxn_id` (COBRA: negative allows uptake).
    temp_fact : float
        Upper-bound scaling factor for reactions in `TEMP_SENSITIVE`.
    ko : Optional[str]
        Reaction ID to knock out; None = no KO.
    bio_rxn : Optional[str]
        Biomass reaction ID (used only if `WT` and `frac` are provided).
    WT : Optional[float]
        Wild-type max growth rate (h^-1) for the same environment.
    frac : float
        Fraction of WT growth to enforce as a biomass lower bound (0–1).
    carbon_uptake : float, keyword-only
        Lower bound for the chosen carbon exchange (mmol·gDW^-1·h^-1).
    o2_lb : Optional[float], keyword-only
        Lower bound for `EX_o2_e`.

    Returns
    -------
    None
    """
    for exch in CARBON_SOURCES.values():
        if exch in m.reactions:
            m.reactions.get_by_id(exch).lower_bound = 0.0

    if c_exch in m.reactions:
        m.reactions.get_by_id(c_exch).lower_bound = float(carbon_uptake)

    if "EX_h_e" in m.reactions:
        m.reactions.EX_h_e.lower_bound = ph_lb

    if rxn_id in m.reactions:
        m.reactions.get_by_id(rxn_id).lower_bound = buf_lb

    if o2_lb is not None and "EX_o2_e" in m.reactions:
        m.reactions.EX_o2_e.lower_bound = float(o2_lb)

    for rid in TEMP_SENSITIVE:
        if rid in m.reactions:
            rx = m.reactions.get_by_id(rid)
            ub0 = rx.upper_bound if abs(rx.upper_bound) < 1e5 else 1000.0
            rx.upper_bound = ub0 * temp_fact

    if ko and ko in m.reactions:
        m.reactions.get_by_id(ko).knock_out()

    if WT and bio_rxn:
        lb_min = max(0.0, WT * float(frac))
        m.reactions.get_by_id(bio_rxn).lower_bound = lb_min


def optimize_citrate(m, bio_rxn):
    """
    Maximize citrate export and report the optimal citrate flux and biomass rate.

    Purpose
    -------
    Sets the model objective to `EX_cit_e`, performs FBA, and returns:
      (citrate_flux, pFBA_flux_placeholder, biomass_flux_at_citrate_optimum).

    Parameters
    ----------
    m : cobra.Model
        Configured model (environment + KO already applied).
    bio_rxn : Optional[str]
        Biomass reaction ID to read the biomass flux at the citrate optimum.

    Returns
    -------
    tuple[float, float, float]
        (citrate_FBA, citrate_pFBA_placeholder, growth_at_citrate)
        NaN triplet if the optimization is infeasible or `EX_cit_e` is absent.
    """
    if "EX_cit_e" not in m.reactions:
        return (np.nan, np.nan, np.nan)
    m.objective = "EX_cit_e"
    sol = m.optimize()
    if sol.status != "optimal":
        return (np.nan, np.nan, np.nan)
        
    gr = np.nan
    if bio_rxn and bio_rxn in m.reactions:
        try:
            gr = sol.fluxes[bio_rxn]
        except Exception:
            gr = np.nan
    return sol.fluxes["EX_cit_e"], np.nan, gr

def optimize_growth(m, bio_rxn) -> float:
     """
    Maximize biomass growth rate.

    Purpose
    -------
    Sets the objective to `bio_rxn` and runs FBA.

    Parameters
    ----------
    m : cobra.Model
        Configured model.
    bio_rxn : str
        Biomass reaction ID.

    Returns
    -------
    float or None
        Optimal growth rate (h^-1), or None if the problem is infeasible.

    Raises
    ------
    KeyError
        If `bio_rxn` does not exist in `m.reactions` (not caught here).
    """
    m.objective = bio_rxn
    sol = m.optimize()
    if sol.status != "optimal":
        return None
    return sol.objective_value


def optimize_siderophore(m):
    """
    Maximize siderophore secretion proxy (enterobactin exchange).

    Purpose
    -------
    Sets the objective to `EX_feenter_e` and performs FBA.

    Parameters
    ----------
    m : cobra.Model
        Configured model.

    Returns
    -------
    tuple[float, float]
        (siderophore_FBA, siderophore_pFBA_placeholder)
        NaNs if infeasible or the exchange reaction is absent.
    """
    sider_id = "EX_feenter_e"
    if sider_id not in m.reactions:
        return np.nan, np.nan
    m.objective = sider_id
    sol = m.optimize()
    if sol.status != "optimal":
        return np.nan, np.nan
    return sol.fluxes[sider_id], np.nan

def run_simulation_wrapper(params):
    """
    Wrapper for parallel workers to limit internal thread fan-out.

    Purpose
    -------
    Ensures that each multiprocessing worker runs with restricted BLAS/OpenMP
    threads to avoid oversubscription on multi-core systems, then delegates to
    `run_simulation`.

    Parameters
    ----------
    params : tuple
        A single parameter tuple as produced by `generate_param_grid`.

    Returns
    -------
    dict
        One result row.
    """
    with threadpool_limits(limits=1):
        return run_simulation(params)

def _ensure_minimal_medium(m):
    """
    Enable a minimal medium to prevent infeasible models due to missing nutrients.

    Purpose
    -------
    Sets conservative uptake lower bounds for standard nutrients:
    O2, NH4+, Pi, SO4, H2O, and H+ (protons).

    Parameters
    ----------
    m : cobra.Model
        Model to modify in place.

    Returns
    -------
    None
    """
    needed = [
        ("EX_o2_e",  -20.0),
        ("EX_nh4_e", -10.0),
        ("EX_pi_e",  -10.0),
        ("EX_so4_e", -10.0),
        ("EX_h2o_e", -1000.0),
        ("EX_h_e",   -1000.0),
    ]
    for rid, lb in needed:
        if rid in m.reactions:
            m.reactions.get_by_id(rid).lower_bound = lb

def _init_worker(models_paths):
    """
    Worker initializer to load models once per process.

    Purpose
    -------
    Populates the global MODELS cache with {model_name: cobra.Model} by reading
    JSON model files (e.g., BiGG) from `models_paths`. This minimizes I/O
    overhead during parallel evaluation.

    Parameters
    ----------
    models_paths : dict
        Mapping {model_name: pathlib.Path} to JSON model files.

    Returns
    -------
    None
    """
    import cobra
    global MODELS
    MODELS = {name: cobra.io.load_json_model(str(p)) for name, p in models_paths.items()}


def run_simulation(params) -> dict:
    """
    Execute a single parameterized simulation and return a result row.

    Purpose
    -------
    Given a parameter tuple (model, carbon source, pH, ion buffer, temperature,
    growth fraction, KO, carbon uptake, O2 LB), this function:
      1) computes WT growth in the same environment (no KO, no growth fraction);
      2) configures a second model with KO and optional minimum growth;
      3) optimizes growth, citrate export (and reads biomass at that optimum),
         and siderophore export;
      4) returns a structured dict of metrics and the parameter annotations.

    Parameters
    ----------
    params : tuple
        ((model_name, model_path),
         (carbon_label, carbon_exch_id),
         (ph_lb, ph_label),
         (ion_label, ion_exch_id, ion_lb),
         (temp_factor, temp_label),
         growth_fraction,
         ko_reaction_id_or_None,
         carbon_uptake_lb,
         o2_lower_bound)

    Returns
    -------
    dict
        Result row with fields:
        {
          'model', 'KO', 'carbon', 'carbon_exch', 'pH_lb', 'pH_label',
          'ion', 'ion_exch', 'temp_factor', 'temp_label', 'growth_frac',
          'WT_growth', 'growth_max',
          'citrate_FBA', 'citrate_pFBA', 'growth_at_citrate',
          'siderophore_FBA', 'siderophore_pFBA',
          'o2_lb', 'carbon_uptake'
        }
    """
    import numpy as np

    (model_name, path), (c_name, c_exch), (ph_lb, ph_lbl), \
        (ion, rxn_id, buf_lb), (temp_fact, temp_lbl), frac, ko, carbon_uptake, o2_lb = params

    def nan_tuple(n=1):
        return tuple(np.nan for _ in range(n))

    base = MODELS[model_name]

    bio_rxn = next((r.id for r in base.reactions if "biomass" in r.id.lower()), None)

    m_wt = base.copy()
    _ensure_minimal_medium(m_wt)
    if "EX_cit_e" in m_wt.reactions:
        m_wt.reactions.EX_cit_e.lower_bound = -1000.0
    configure_model(
        m_wt, c_exch, ph_lb, rxn_id, buf_lb, temp_fact,
        ko=None, bio_rxn=bio_rxn, WT=None, frac=0.0,
        carbon_uptake=carbon_uptake, o2_lb=o2_lb
    )
    WT_growth = optimize_growth(m_wt, bio_rxn) if bio_rxn else np.nan

    m = base.copy()
    _ensure_minimal_medium(m)
    if "EX_cit_e" in m.reactions:
        m.reactions.EX_cit_e.lower_bound = -1000.0
    configure_model(
        m, c_exch, ph_lb, rxn_id, buf_lb, temp_fact,
        ko=ko, bio_rxn=bio_rxn, WT=WT_safe, frac=frac,
        carbon_uptake=carbon_uptake, o2_lb=o2_lb
    )
    growth_max = optimize_growth(m.copy(), bio_rxn) if bio_rxn else np.nan
    if growth_max is None:
        growth_max = np.nan

    cit_fba, cit_pfb, gr_cit = optimize_citrate(m.copy(), bio_rxn)

    sid_fba, sid_pfb = optimize_siderophore(m.copy())

    return {
        'model': model_name,
        'KO': ko or 'None',
        'carbon': c_name,
        'carbon_exch': c_exch,
        'pH_lb': ph_lb,
        'pH_label': ph_lbl,
        'ion': ion,
        'ion_exch': rxn_id,
        'temp_factor': temp_fact,
        'temp_label': temp_lbl,
        'growth_frac': float(frac),

        'WT_growth': WT_growth if WT_growth is not None else np.nan,
        'growth_max': growth_max,

        'citrate_FBA': cit_fba,
        'citrate_pFBA': cit_pfb,
        'growth_at_citrate': gr_cit,

        'siderophore_FBA': sid_fba,
        'siderophore_pFBA': sid_pfb,
        'o2_lb': float(o2_lb),
        'carbon_uptake': float(carbon_uptake),


    }

def save_results(results: list, output: Path, write_header: bool):
    """
    Append a batch of result dicts to CSV and print quick rankings.

    Purpose
    -------
    Converts a list of per-condition dicts into a DataFrame, appends them to
    `output` (creating header if requested), and prints top-k summaries for
    sanity checking during long runs.

    Parameters
    ----------
    results : list[dict]
        Batch of result rows as returned by `run_simulation`.
    output : pathlib.Path
        Destination CSV path.
    write_header : bool
        If True, write the CSV header; afterwards use False to append.

    Returns
    -------
    None
    """
    df = pd.DataFrame(results)
    df.to_csv(output, index=False, mode='a', header=write_header)
    print(f"Appending {len(df)} rows -> {output}")
    print("Top 5 proliferation (growth_max):")
    print(df.nlargest(5, 'growth_max'), "\n")
    print("Top 10 citrate export (citrate_FBA):")
    print(df.nlargest(10, 'citrate_FBA'), "\n")
    print("Top 5 production of siderophore (FBA):")
    print(df.nlargest(5, 'siderophore_FBA'), "\n")

def main():
    """
    Orchestrate model fetching, parameter grid construction, and parallel runs.

    Steps
    -----
    1) Ensure `models/` exists. If empty, fetch BiGG models filtered by organism
       (e.g., "coli") using `fetch_models_based_on_org`.
    2) Build the parameter iterator (cartesian product of all factors).
    3) Print a run banner with environment metadata and an SHA1 of this script.
    4) Prepare the CSV (remove if exists) and set append mode with header-once.
    5) Launch a multiprocessing pool with per-worker model initialization.
    6) Consume results with unordered iteration, flush to CSV every N rows.
    7) Final flush; print the absolute path of the aggregated CSV.

    Returns
    -------
    None
    """
    import os, hashlib
    from datetime import datetime

    models_loc = Path('models')
    models_loc.mkdir(parents=True, exist_ok=True)

    if models_loc.is_dir() and not any(models_loc.iterdir()):
        url = "http://bigg.ucsd.edu/api/v2/models"
        fetch_models_based_on_org(url=url, outdir=str(models_loc), organisms=["coli"])

    models = {m.stem: m for m in models_loc.iterdir() if m.is_file()}

    ion_params = flatten_ion_buffers(ION_BUFFERS)
    param_iter = generate_param_grid(
        models,
        CARBON_SOURCES,
        PH_LEVELS,
        ion_params,
        TEMP_LEVELS,
        GROWTH_FRACS,
        KOS,
        CARBON_UPTAKES,
        O2_LBS
    )

    here = Path(__file__).resolve()
    try:
        script_sha1 = hashlib.sha1(here.read_bytes()).hexdigest()[:12]
    except Exception:
        script_sha1 = "n/a"

    total_expected = (
        len(models) *
        len(CARBON_SOURCES) *
        len(PH_LEVELS) *
        len(ion_params) *
        len(TEMP_LEVELS) *
        len(GROWTH_FRACS) *
        len(KOS) *
        len(CARBON_UPTAKES) *
        len(O2_LBS)
    )

    print("=" * 100)
    print("RUN BANNER @", datetime.now().isoformat())
    print("CWD:", os.getcwd())
    print("Script:", str(here))
    print("Script SHA1:", script_sha1)
    print(f"Models: {sorted(models.keys())}")
    print(f"CARBON_SOURCES: {list(CARBON_SOURCES.keys())}")
    print(f"PH labels: {[lbl for _, lbl in PH_LEVELS]}")
    print(f"IONS: {sorted(set([ion for (ion, _, _) in ion_params]))}")
    print(f"TEMP labels: {[lbl for _, lbl in TEMP_LEVELS]}")
    print(f"GROWTH_FRACS: {GROWTH_FRACS}")
    print(f"KOs: {KOS}")
    print("Expected combinations:", total_expected)
    print("=" * 100)
    print(f"O2 LBs: {O2_LBS}")
    print(f"CARBON_UPTAKES: {CARBON_UPTAKES}")

    if CSV_OUTPUT.exists():
        CSV_OUTPUT.unlink()
    CSV_OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    write_header = True

    flush_every = 200
    n_procs = 10
    ctx = get_context("spawn")
    chunk = []

    with ctx.Pool(processes=n_procs, initializer=_init_worker, initargs=(models,)) as pool:
        for res in pool.imap_unordered(run_simulation_wrapper, param_iter, chunksize=10):
            if res is None:
                continue
            chunk.append(res)
            if len(chunk) >= flush_every:
                save_results(chunk, CSV_OUTPUT, write_header)
                write_header = False
                chunk = []
    if chunk:
        save_results(chunk, CSV_OUTPUT, write_header)

    print(f"Completed. Aggregated results in: {CSV_OUTPUT.resolve()}")


if __name__ == '__main__':
    main()






