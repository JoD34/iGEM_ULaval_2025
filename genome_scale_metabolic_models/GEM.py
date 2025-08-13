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
    (0, "basique pH~9")
#    (5, "très basique")
]
CARBON_SOURCES = {
    "glucose": "EX_glc__D_e",
    "acetate": "EX_ac_e"
#    "glycerol": "EX_glyc_e",
#    "succinate": "EX_succ_e",
#    "ethanol": "EX_etoh_e",
}
ION_BUFFERS = {
    "HCO3": ("EX_hco3_e", np.linspace(-20, 0, 9)),
#    "Pi":   ("EX_pi_e",   np.linspace(-20, 0, 9)),
#    "NH4":  ("EX_nh4_e",  np.linspace(-20, 0, 9)),
#    "CO2":  ("EX_co2_e",  np.linspace(-5, 0, 11)),
    "Fe":   ("EX_fe2_e",  np.linspace(-1000, -1, 10)),
}
TEMP_LEVELS = [
    #(0.7, "faible (~30 °C)"),
    (0.9, "légère (~34 °C)"),
    (1.0, "standard (~37 °C)"),
    (1.1, "modérée (~40 °C)"),
    #(1.3, "élevée (~45 °C)")
]
TEMP_SENSITIVE = ["ACONTa", "ACONTb", "ICDHyr", "SUCDi"]
GROWTH_FRACS = [
    0.0,
    #0.05,
    0.1,
    #0.2
]
KOS = [
    None
#    "ICDHyr"
]
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
                        growth_fracs: list, kos: list):
    return product(
        models.items(),
        carb_sources.items(),
        ph_levels,
        ion_params,
        temp_levels,
        growth_fracs,
        kos
    )


def configure_model(m, c_exch, ph_lb, rxn_id, buf_lb,
                    temp_fact, ko, bio_rxn, WT, frac):
    # Désactive toutes les sources de carbone
    for exch in CARBON_SOURCES.values():
        if exch in m.reactions:
            m.reactions.get_by_id(exch).lower_bound = 0.0
    # Active la source sélectionnée
    if c_exch in m.reactions:
        m.reactions.get_by_id(c_exch).lower_bound = -10.0
    # pH
    if "EX_h_e" in m.reactions:
        m.reactions.EX_h_e.lower_bound = ph_lb
    # tampon ionique
    if rxn_id in m.reactions:
        m.reactions.get_by_id(rxn_id).lower_bound = buf_lb
    # température
    for rid in TEMP_SENSITIVE:
        if rid in m.reactions:
            rx = m.reactions.get_by_id(rid)
            ub0 = rx.upper_bound if abs(rx.upper_bound) < 1e5 else 1000
            rx.upper_bound = ub0 * temp_fact
    # knock-out
    if ko and ko in m.reactions:
        m.reactions.get_by_id(ko).knock_out()
 # croissance minimale (si besoin)
    if WT and bio_rxn:
        lb_min = max(0.0, WT * float(frac))
        m.reactions.get_by_id(bio_rxn).lower_bound = lb_min

def optimize_citrate(m, bio_rxn):
    if "EX_cit_e" not in m.reactions:
        return None
    m.objective = "EX_cit_e"
    sol = m.optimize()
    if sol.status != "optimal":
        return None

#    m2 = m.copy()
#    m2.solver = "glpk"       # nouvelle instance
#    m2.solver.problem = None # <<< purge l’état interne du solveur
#     m2.objective = "EX_cit_e"  # optionnel, déjà défini
#    pf_sol = run_pfba(m2)

#    return sol.fluxes["EX_cit_e"], pf_sol.fluxes["EX_cit_e"], sol.fluxes[bio_rxn]
    return sol.fluxes["EX_cit_e"], np.nan, sol.fluxes[bio_rxn]


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

def _init_worker(models_paths):
    # Chargé une fois par worker
    import cobra
    global MODELS
    MODELS = {name: cobra.io.load_json_model(str(p)) for name, p in models_paths.items()}


def run_simulation(params) -> dict:
    (model_name, path), (c_name, c_exch), (ph_lb, ph_lbl), \
        (ion, rxn_id, buf_lb), (temp_fact, temp_lbl), frac, ko = params

    base = MODELS[model_name]
    bio_rxn = next(r.id for r in base.reactions if "biomass" in r.id.lower()) # Cétecte la réaction de biomass
    WT = base.optimize().objective_value # Calcul la croissance de référence
    if "EX_cit_e" in base.reactions:
        base.reactions.EX_cit_e.lower_bound = -1000 # Permet exportation de citrate libre

    m = base.copy() # Cloné pour modifier sans altérer base
    m.solver = "glpk" 
    configure_model(m, c_exch, ph_lb, rxn_id, buf_lb, temp_fact, ko, bio_rxn, WT, frac)
    cit = optimize_citrate(m, bio_rxn)
    if cit is None:
        return None
    cit_fba, cit_pfb, gr_cit = cit

    growth_max = optimize_growth(m, bio_rxn)
    if growth_max is None:
        return None

    sid_fba, sid_pfb = optimize_siderophore(m)

    return {
        'model': model_name,
        'C_source': c_name,
        'pH': ph_lbl,
        'ion': ion,
        'ion_lb': round(buf_lb, 2),
        'température': temp_lbl,
        'growth_%': int(frac * 100),
        'KO': ko or 'none',
        'citrate_FBA': round(cit_fba, 4),
        'citrate_pFBA': round(cit_pfb, 4),
        'growth_under_cit': round(gr_cit, 4),
        'growth_max': round(growth_max, 4),
        'siderophore_FBA': round(sid_fba, 4),
        'siderophore_pFBA': round(sid_pfb, 4),
    }


def save_results(results: list, output: Path):
    df = pd.DataFrame(results)
    df.to_csv(output, index=False)
    print(f"Scénarios simulés : {len(df)}\n")
    print("Top 5 prolifération (growth_max):")
    print(df.nlargest(5, 'growth_max'), "\n")
    print("Top 10 export citrate (citrate_FBA):")
    print(df.nlargest(10, 'citrate_FBA'), "\n")
    print("Top 5 production de sidérophore (FBA):")
    print(df.nlargest(5, 'siderophore_FBA'), "\n")


def main():
    # Step 1: Assert models présents
    models_loc = Path('models')
    models_loc.mkdir(exist_ok=True)

    if models_loc.is_dir() and not any(models_loc.iterdir()):
        # Si vide, on fetch tous les modèles depuis BIGG
        url = "http://bigg.ucsd.edu/api/v2/models"
        fetch_models_based_on_org(url=url, outdir=str(models_loc), organisms=["coli"])
    models = {m.stem: m for m in models_loc.iterdir()}

    # Step 2: Prépare le plan d'expériences
    ion_params = flatten_ion_buffers(ION_BUFFERS)
    param_iter = generate_param_grid(models,
                                     CARBON_SOURCES,
                                     PH_LEVELS,
                                     ion_params,
                                     TEMP_LEVELS,
                                     GROWTH_FRACS,
                                     KOS)

    # Step 3: Exécution parallèle des simulations
    n_procs = 10 # Exploration parallèle avec 10 processus
    ctx = get_context("spawn")

    with ctx.Pool(processes=n_procs, initializer=_init_worker, initargs=(models,)) as pool:
        results = []
        for res in pool.imap_unordered(run_simulation_wrapper, param_iter, chunksize=10):
            if res:
                results.append(res)
                if len(results) % 200 == 0:
                    save_results(results, CSV_OUTPUT)  # checkpoint
                    results = []
    # à la fin, sauve le reste
    if results:
        save_results(results, CSV_OUTPUT)


if __name__ == '__main__':
    main()
