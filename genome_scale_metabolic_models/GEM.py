import cobra    # Pour l'analyse du métabolisme
from cobra.flux_analysis import pfba
import pandas as pd
import numpy as np  # Pour les calculs scientifiques
import warnings
from pathlib import Path
from itertools import product
from multiprocessing import Pool
from fetch_gems import fetch_models_based_on_org

warnings.filterwarnings("ignore", category=UserWarning)

# -- Configuration des constantes et paramètres --
PH_LEVELS = [(-20, "extrême acide"), (-10, "acide pH~5"), (-1, "neutre pH~7"),
             (0, "basique pH~9"), (5, "très basique")]
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
    "Fe":   ("EX_fe2_e",  np.linspace(-1000, -1, 10)),
}
TEMP_LEVELS = [(0.7, "faible (~30 °C)"), (0.9, "légère (~34 °C)"),
               (1.0, "standard (~37 °C)"), (1.1, "modérée (~40 °C)"),
               (1.3, "élevée (~45 °C)")]
TEMP_SENSITIVE = ["ACONTa", "ACONTb", "ICDHyr", "SUCDi"]
GROWTH_FRACS = [0.0, 0.05, 0.1, 0.2]
KOS = [None, "ICDHyr"]
CSV_OUTPUT = Path('citrate_growth_siderophore_scan.csv')

# -- Fonctions utilitaires --
def flatten_ion_buffers(buffers: dict) -> list:
    """
    Aplatie la structure ION_BUFFERS en une liste de tuples (ion, rxn_id, buf_lb).

    Args:
        buffers (dict): map ion -> (rxn_id, liste de bornes)
    Returns:
        list of tuples: [(ion, reaction_id, buffer_lower_bound), ...]
    """
    return [
        (ion, rxn_id, buf_lb)
        for ion, (rxn_id, buf_list) in buffers.items()
        for buf_lb in buf_list
    ]


def generate_param_grid(models: dict, carb_sources: dict, ph_levels: list,
                        ion_params: list, temp_levels: list,
                        growth_fracs: list, kos: list):
    """
    Génère un itérateur cartésien de toutes les combinaisons de paramètres.

    Args:
        models (dict): map model_name -> Path du modèle JSON
        carb_sources (dict): map nom_source -> id de réaction
        ph_levels (list): liste de tuples (pH_lb, étiquette)
        ion_params (list): liste de tuples (ion, reaction_id, buffer_lb)
        temp_levels (list): liste de tuples (facteur_temp, étiquette)
        growth_fracs (list): fractions de croissance minimale
        kos (list): knock-outs possibles ou None
    Returns:
        Iterator: itertools.product générant un tuple complet par combinaison
    """
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
    """
    Applique les contraintes expérimentales sur le modèle COBRApy cloné.

    Args:
        m: modèle COBRApy (contexte `with base as m`)
        c_exch (str): id de la réaction de carbone à activer
        ph_lb (float): borne inférieure pour la réaction EX_h_e
        rxn_id (str): id de la réaction tampon ionique
        buf_lb (float): borne inférieure pour la réaction tampon
        temp_fact (float): facteur multiplicatif appliqué aux réactions sensibles
        ko (str|None): id de réaction à KO, ou None
        bio_rxn (str): id de la réaction de biomass
        WT (float): flux de biomass en condition WT
    """
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
        cons = m.solver.interface.Constraint(
            m.reactions.get_by_id(bio_rxn).flux_expression,
            lb=0, ub=1e6, name="min_growth"
        )
        m.add_cons_vars(cons)

def optimize_growth(m, bio_rxn) -> float:
    """
    Optimise la croissance maximale (biomass objective).
    """
    m.objective = bio_rxn
    sol = m.optimize()
    if sol.status != "optimal":
        return None
    return sol.objective_value

def optimize_citrate(m, bio_rxn) -> tuple:
    """
    Optimise l'export de citrate et renvoie ses flux.

    Args:
        m: modèle COBRApy
        bio_rxn: id de la réaction de biomass
    Returns:
        (flux_FBA, flux_pFBA, flux_biomass_durant_export) ou None
    """
    m.objective = "EX_cit_e"
    sol = m.optimize()
    if sol.status != "optimal":
        return None
    pfba_fluxes = pfba(m)
    return sol.fluxes["EX_cit_e"], pfba_fluxes["EX_cit_e"], sol.fluxes[bio_rxn]


def optimize_siderophore(m) -> tuple:
    """
    Optimise l'export du sidérophore et renvoie ses flux.

    Args:
        m: modèle COBRApy
    Returns:
        (flux_FBA, flux_pFBA)
    """
    sider_id = "EX_feenter_e"
    if sider_id not in m.reactions:
        return np.nan, np.nan

    m.objective = sider_id
    sol = m.optimize()
    if sol.status != "optimal":
        return np.nan, np.nan

    pfba_fluxes = pfba(m)
    return sol.fluxes[sider_id], pfba_fluxes[sider_id]


def run_simulation(params) -> dict:
    """
    Exécute une simulation complète pour un jeu de paramètres.

    Args:
        params: tuple fourni par generate_param_grid
    Returns:
        dict des résultats ou None si échec d'une étape
    """
    (model_name, path), (c_name, c_exch), (ph_lb, ph_lbl), \
        (ion, rxn_id, buf_lb), (temp_fact, temp_lbl), frac, ko = params

    base = cobra.io.load_json_model(path)
    bio_rxn = next(r.id for r in base.reactions if "biomass" in r.id.lower())
    WT = base.optimize().objective_value
    base.reactions.EX_cit_e.lower_bound = -1000

    with base as m:
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
        'citrine_pFBA': round(cit_pfb, 4),
        'growth_under_cit': round(gr_cit, 4),
        'growth_max': round(growth_max, 4),
        'siderophore_FBA': round(sid_fba, 4),
        'siderophore_pFBA': round(sid_pfb, 4),
    }


def save_results(results: list, output: Path):
    """
    Sauvegarde les résultats dans un CSV et affiche un résumé.

    Args:
        results (list): liste de dicts de résultats
        output (Path): chemin vers le CSV de sortie
    """
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
    """
    Point d'entrée principal :
    - Vérifie/télécharge les modèles
    - Génère le plan d'expériences
    - Exécute les simulations en parallèle
    - Sauvegarde et affiche les résultats
    """
    # Step 1: Assert models présents
    models_loc = Path('models')
    models_loc.mkdir(exist_ok=True)

    if models_loc.is_dir() and not any(models_loc.iterdir()):
        # Si vide, on fetch tous les modèles depuis BIGG
        url = "http://bigg.ucsd.edu/api/v2/models"
        fetch_models_based_on_org(url=url, outdir=str(models_loc))
    models = {m.stem: m for m in models_loc.iterdir()}

    # Step 2: Prépare le plan d'expériences
    ion_params = flatten_ion_buffers(ION_BUFFERS)
    param_grid = generate_param_grid(models,
                                     CARBON_SOURCES,
                                     PH_LEVELS,
                                     ion_params,
                                     TEMP_LEVELS,
                                     GROWTH_FRACS,
                                     KOS)

    # Step 3: Exécution parallèle des simulations
    with Pool() as pool:
        results = pool.map(run_simulation, param_grid)
    results = [r for r in results if r]

    # Step 4: Sauvegarde et résumé
    save_results(results, CSV_OUTPUT)


if __name__ == '__main__':
    main()
