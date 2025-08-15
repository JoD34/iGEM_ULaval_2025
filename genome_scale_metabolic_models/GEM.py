#from cobra.flux_analysis.parsimonious import pfba as run_pfba
import cobra    # Pour l'analyse du métabolisme
#from cobra.flux_analysis import pfba
from cobra import io as cbio
import pandas as pd
import numpy as np  # Pour les calculs scientifiques
import warnings
from pathlib import Path
from itertools import product
from multiprocessing import get_context
from multiprocessing import Pool
from fetch_gems import fetch_models_based_on_org
from threadpoolctl import threadpool_limits

warnings.filterwarnings("ignore", category=UserWarning)

MODELS = {}

# -- Configuration des constantes et paramètres --
PH_LEVELS = [
#    (-20, "extrême acide"),
    (-10, "acide pH~5"),
    (-1, "neutre pH~7"),
    (0, "basique pH~9"),
#    (5, "très basique")
]
CARBON_SOURCES = {
    "glucose": "EX_glc__D_e"
#    "acetate": "EX_ac_e"
#    "glycerol": "EX_glyc_e"
#    "succinate": "EX_succ_e",
#    "ethanol": "EX_etoh_e",
}
ION_BUFFERS = {
    "HCO3": ("EX_hco3_e", np.linspace(-20, 0, 9)),
    "Pi":   ("EX_pi_e",   np.linspace(-20, 0, 9)),
    "NH4":  ("EX_nh4_e",  np.linspace(-20, 0, 9)),
    "CO2":  ("EX_co2_e",  np.linspace(-5, 0, 11)),
    "Fe":   ("EX_fe2_e",  np.linspace(-1000, -1, 10))
}

CARBON_UPTAKES = [-2.0, -5.0, -10.0, -20.0]   # mmol/gDW/h (borne inférieure des EX_*)
O2_LBS = [-20.0, -10.0, -5.0, -1.0, 0.0]

TEMP_LEVELS = [
#    (0.7, "faible (~30 °C)"),
    (0.9, "légère (~34 °C)"),
    (1.0, "standard (~37 °C)"),
    (1.1, "modérée (~40 °C)")
#    (1.3, "élevée (~45 °C)")
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

# -- Fonctions utilitaires --
def flatten_ion_buffers(buffers: dict) -> list:
    return [
        (ion, rxn_id, buf_lb)
        for ion, (rxn_id, buf_list) in buffers.items()
        for buf_lb in buf_list
    ]


def generate_param_grid(models: dict, carb_sources: dict, ph_levels: list,
                        ion_params: list, temp_levels: list,
                        growth_fracs: list, kos: list,
                        carbon_uptakes: list, o2_lbs: list):
    # NOTE: l'ordre ici DOIT matcher l'unpacking dans run_simulation(...)
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
    # 1) Désactive toutes les sources de carbone
    for exch in CARBON_SOURCES.values():
        if exch in m.reactions:
            m.reactions.get_by_id(exch).lower_bound = 0.0

    # 2) Active la source sélectionnée avec le débit demandé
    if c_exch in m.reactions:
        m.reactions.get_by_id(c_exch).lower_bound = float(carbon_uptake)

    # 3) pH (via H+)
    if "EX_h_e" in m.reactions:
        m.reactions.EX_h_e.lower_bound = ph_lb

    # 4) Tampon ionique / nutriments
    if rxn_id in m.reactions:
        m.reactions.get_by_id(rxn_id).lower_bound = buf_lb

    # 5) O2 (limitation / anaérobie)
    if o2_lb is not None and "EX_o2_e" in m.reactions:
        m.reactions.EX_o2_e.lower_bound = float(o2_lb)

    # 6) Température (borne sup. des réactions sensibles)
    for rid in TEMP_SENSITIVE:
        if rid in m.reactions:
            rx = m.reactions.get_by_id(rid)
            ub0 = rx.upper_bound if abs(rx.upper_bound) < 1e5 else 1000.0
            rx.upper_bound = ub0 * temp_fact

    # 7) Knock-out
    if ko and ko in m.reactions:
        m.reactions.get_by_id(ko).knock_out()

    # 8) Contrainte de croissance minimale (si besoin)
    if WT and bio_rxn:
        lb_min = max(0.0, WT * float(frac))
        m.reactions.get_by_id(bio_rxn).lower_bound = lb_min


def optimize_citrate(m, bio_rxn):
    if "EX_cit_e" not in m.reactions:
        return (np.nan, np.nan, np.nan)
    m.objective = "EX_cit_e"
    sol = m.optimize()
    if sol.status != "optimal":
        return (np.nan, np.nan, np.nan)

#    m2 = m.copy()
#    m2.solver = "glpk"       # nouvelle instance
#    m2.solver.problem = None # <<< purge l’état interne du solveur
#     m2.objective = "EX_cit_e"  # optionnel, déjà défini
#    pf_sol = run_pfba(m2)

#    return sol.fluxes["EX_cit_e"], pf_sol.fluxes["EX_cit_e"], sol.fluxes[bio_rxn]
    gr = np.nan
    if bio_rxn and bio_rxn in m.reactions:
        try:
            gr = sol.fluxes[bio_rxn]
        except Exception:
            gr = np.nan
    return sol.fluxes["EX_cit_e"], np.nan, gr

def optimize_growth(m, bio_rxn) -> float:
    m.objective = bio_rxn
    sol = m.optimize()
    if sol.status != "optimal":
        return None
    return sol.objective_value


def optimize_siderophore(m):
    sider_id = "EX_feenter_e"
    if sider_id not in m.reactions:
        return np.nan, np.nan
    m.objective = sider_id
    sol = m.optimize()
    if sol.status != "optimal":
        return np.nan, np.nan

#    m2 = m.copy()
#    m2.solver = "glpk"
#    m2.solver.problem = None  # <<< purge
#    m2.objective = sider_id
#    pf_sol = run_pfba(m2)

#    return sol.fluxes[sider_id], pf_sol.fluxes[sider_id]
    return sol.fluxes[sider_id], np.nan

def run_simulation_wrapper(params):
    # Chaque worker n'utilise qu'un thread pour BLAS/OpenMP internes
    with threadpool_limits(limits=1):
        return run_simulation(params)

def _ensure_minimal_medium(m):
    """Active un milieu minimal standard pour éviter les infeasibilités arbitraires."""
    needed = [
        ("EX_o2_e",  -20.0),   # oxygène
        ("EX_nh4_e", -10.0),   # ammonium
        ("EX_pi_e",  -10.0),   # phosphate
        ("EX_so4_e", -10.0),   # sulfate
        ("EX_h2o_e", -1000.0), # eau
        ("EX_h_e",   -1000.0), # protons (utile pour certaines balances)
    ]
    for rid, lb in needed:
        if rid in m.reactions:
            m.reactions.get_by_id(rid).lower_bound = lb

def _init_worker(models_paths):
    # Chargé une fois par worker
    import cobra
    global MODELS
    MODELS = {name: cobra.io.load_json_model(str(p)) for name, p in models_paths.items()}


def run_simulation(params) -> dict:
    import numpy as np

    (model_name, path), (c_name, c_exch), (ph_lb, ph_lbl), \
        (ion, rxn_id, buf_lb), (temp_fact, temp_lbl), frac, ko, carbon_uptake, o2_lb = params

    def nan_tuple(n=1):
        return tuple(np.nan for _ in range(n))

    # Charger modèle
    base = MODELS[model_name]

    # Trouver la réaction biomasse
    bio_rxn = next((r.id for r in base.reactions if "biomass" in r.id.lower()), None)

    # ---- WT growth dans l’environnement (sans KO, sans contrainte frac) ----
    # --- WT growth (sans KO, sans contrainte frac) ---
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

    # --- Modèle final (avec KO + contrainte de fraction) ---
    m = base.copy()
    _ensure_minimal_medium(m)
    if "EX_cit_e" in m.reactions:
        m.reactions.EX_cit_e.lower_bound = -1000.0
    configure_model(
        m, c_exch, ph_lb, rxn_id, buf_lb, temp_fact,
        ko=ko, bio_rxn=bio_rxn, WT=WT_safe, frac=frac,
        carbon_uptake=carbon_uptake, o2_lb=o2_lb
    )
    # Growth max
    growth_max = optimize_growth(m.copy(), bio_rxn) if bio_rxn else np.nan
    if growth_max is None:
        growth_max = np.nan

    # Citrate
    cit_fba, cit_pfb, gr_cit = optimize_citrate(m.copy(), bio_rxn)

    # Sidérophore
    sid_fba, sid_pfb = optimize_siderophore(m.copy())

    # Ligne finale (toujours renvoyée)
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
    df = pd.DataFrame(results)
    df.to_csv(output, index=False, mode='a', header=write_header)
    print(f"Appending {len(df)} rows -> {output}")
    print("Top 5 prolifération (growth_max):")
    print(df.nlargest(5, 'growth_max'), "\n")
    print("Top 10 export citrate (citrate_FBA):")
    print(df.nlargest(10, 'citrate_FBA'), "\n")
    print("Top 5 production de sidérophore (FBA):")
    print(df.nlargest(5, 'siderophore_FBA'), "\n")

def main():
    import os, hashlib
    from datetime import datetime

    # 1) Préparer le dossier des modèles
    models_loc = Path('models')
    models_loc.mkdir(parents=True, exist_ok=True)

    # Si vide, on télécharge les modèles (BiGG)
    if models_loc.is_dir() and not any(models_loc.iterdir()):
        url = "http://bigg.ucsd.edu/api/v2/models"
        fetch_models_based_on_org(url=url, outdir=str(models_loc), organisms=["coli"])

    # Charger la liste des modèles (fichiers .json)
    models = {m.stem: m for m in models_loc.iterdir() if m.is_file()}

    # 2) Construire le plan d'expériences
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

    # 2bis) Bannière de debug
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
        len(CARBON_UPTAKES)
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

    # 3) Écriture CSV en append (header une seule fois)
    if CSV_OUTPUT.exists():
        CSV_OUTPUT.unlink()
    CSV_OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    write_header = True

    flush_every = 200  # <-- défini ici
    n_procs = 10
    ctx = get_context("spawn")
    chunk = []

    with ctx.Pool(processes=n_procs, initializer=_init_worker, initargs=(models,)) as pool:
        for res in pool.imap_unordered(run_simulation_wrapper, param_iter, chunksize=10):
            if res is None:
                continue
            chunk.append(res)
            if len(chunk) >= flush_every:
                save_results(chunk, CSV_OUTPUT, write_header)  # append + header si write_header=True
                write_header = False
                chunk = []

    # Flush final
    if chunk:
        save_results(chunk, CSV_OUTPUT, write_header)

    print(f"Terminé. Résultats agrégés dans: {CSV_OUTPUT.resolve()}")


if __name__ == '__main__':
    main()
