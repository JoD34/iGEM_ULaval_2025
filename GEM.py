import cobra    # Pour l'analyse du métabolisme
from cobra.flux_analysis import flux_variability_analysis as pfba

import pandas as pd
import numpy as np  # Pour les calculs scientifiques
import warnings
from pathlib import Path



warnings.filterwarnings("ignore", category=UserWarning)

# Step 1: Assert model use are present in the directory
curr_loc = Path(__file__)
models_loc = Path('models')
curr_loc.mkdir(models_loc, exist_ok=True) 

if models_loc.is_dir() and not any(models_loc.iterdir()):
    # Fetch models if needed
    
    pass
models = {m.stem: m for m in models_loc.iterdir()}

"""
models = {
    "iML1515": "C:/Users/aliso/OneDrive - Université Laval/iGEM_E.coli/iML1515.json",
    "iAF1260" : "C:/Users/aliso/OneDrive - Université Laval/iGEM_E.coli/iAF1260.json",
    "iJO1366" : "C:/Users/aliso/OneDrive - Université Laval/iGEM_E.coli/iJO1366.json"
}
"""

# 2) Paramètres environnementaux
ph_levels = [
    (-20, "extrême acide"),
    (-10, "acide pH~5"),
    (-1,  "neutre pH~7"),
    (0,   "basique pH~9"),
    (5,   "très basique")
]

carbon_sources = {
    "glucose":   "EX_glc__D_e",
    "acetate":   "EX_ac_e",
    "glycerol":  "EX_glyc_e",
    "succinate": "EX_succ_e",
    "ethanol":   "EX_etoh_e",
}

ion_buffers = {
    "HCO3": ("EX_hco3_e", np.linspace(-20,   0, 9)),
    "Pi":   ("EX_pi_e",   np.linspace(-20,   0, 9)),
    "NH4":  ("EX_nh4_e",  np.linspace(-20,   0, 9)),
    "CO2":  ("EX_co2_e",  np.linspace(-5,    0,11)),
    "Fe":   ("EX_fe2_e",  np.linspace(-1000, -1,10)),  # tampon fer
}
temp_levels = [
    (0.7, "faible (~30 °C)"),
    (0.9, "légère (~34 °C)"),
    (1.0, "standard (~37 °C)"),
    (1.1, "modérée (~40 °C)"),
    (1.3, "élevée (~45 °C)")
]
temp_sensitive = ["ACONTa", "ACONTb", "ICDHyr", "SUCDi"]
growth_fracs = [0.0, 0.05, 0.1, 0.2]
kos = [None, "ICDHyr"]

# 3) Boucle de simulation
results = []
for model_name, path in models.items():
    base = cobra.io.load_json_model(path)
    bio_rxn = next(r.id for r in base.reactions if "biomass" in r.id.lower())
    WT = base.optimize().objective_value
    base.reactions.EX_cit_e.lower_bound = -1000

    for c_name, c_exch in carbon_sources.items():
        for ph_lb, ph_lbl in ph_levels:
            for ion, (rxn_id, buf_list) in ion_buffers.items():
                for buf_lb in buf_list:
                    for temp_fact, temp_lbl in temp_levels:
                        for frac in growth_fracs:
                            for ko in kos:
                                with base as m:
                                    # 0) Configurer la source de carbone
                                    # désactive toutes
                                    for exch in carbon_sources.values():
                                        if exch in m.reactions:
                                            m.reactions.get_by_id(exch).lower_bound = 0.0
                                    # active la source courante à -10 mmol/gDW/h
                                    if c_exch in m.reactions:
                                        m.reactions.get_by_id(c_exch).lower_bound = -10.0

                                    # pH
                                    if "EX_h_e" in m.reactions:
                                        m.reactions.EX_h_e.lower_bound = ph_lb
                                    # tampon ionique
                                    if rxn_id in m.reactions:
                                        m.reactions.get_by_id(rxn_id).lower_bound = buf_lb
                                    # température
                                    for rid in temp_sensitive:
                                        if rid in m.reactions:
                                            rx = m.reactions.get_by_id(rid)
                                            ub0 = rx.upper_bound if abs(rx.upper_bound)<1e5 else 1000
                                            rx.upper_bound = ub0 * temp_fact
                                    # knock-out
                                    if ko and ko in m.reactions:
                                        m.reactions.get_by_id(ko).knock_out()
                                    # contrainte de croissance min
                                    if frac > 0:
                                        cons = m.solver.interface.Constraint(
                                            m.reactions.get_by_id(bio_rxn).flux_expression,
                                            lb=frac*WT, ub=1e6, name="min_growth"
                                        )
                                        m.add_cons_vars(cons)

                                    # 3.1) export citrate
                                    m.objective = "EX_cit_e"
                                    sol_cit = m.optimize()
                                    if sol_cit.status != "optimal":
                                        continue
                                    cit_fba = sol_cit.fluxes["EX_cit_e"]
                                    cit_pfb = pfba(m).fluxes["EX_cit_e"]
                                    # cit_var = fva(m, ["EX_cit_e"]).loc["EX_cit_e"].to_dict()
                                    gr_cit  = sol_cit.fluxes[bio_rxn]

                                    # 3.2) croissance max
                                    m.objective = bio_rxn
                                    sol_gr = m.optimize()
                                    if sol_gr.status != "optimal":
                                        continue
                                    growth_max = sol_gr.objective_value
                                    # gr_var     = fva(m, [bio_rxn]).loc[bio_rxn].to_dict()

                                    # 3.3) export sidérophore
                                    sider_id = "EX_feenter_e"
                                    m.objective = sider_id
                                    sol_sid = m.optimize()
                                    if sol_sid.status == "optimal":
                                        sid_fba = sol_sid.fluxes[sider_id]
                                        sid_pfb = pfba(m).fluxes[sider_id]
                                        # sid_var = fva(m, [sider_id]).loc[sider_id].to_dict()
                                    else:
                                        sid_fba = sid_pfb = np.nan
                                        # sid_var = {"minimum": np.nan, "maximum": np.nan}

                                    # Stockage
                                    results.append({
                                        "model":            model_name,
                                        "C_source":         c_name,
                                        "pH":               ph_lbl,
                                        "ion":              ion,
                                        "ion_lb":           round(buf_lb,2),
                                        "température":      temp_lbl,
                                        "growth_%":         int(frac*100),
                                        "KO":               ko or "none",
                                        # citrate
                                        "citrate_FBA":      round(cit_fba,    4),
                                        "citrate_pFBA":     round(cit_pfb,    4),
                                        # "cit_min_FVA":      round(cit_var["minimum"],4),
                                        # "cit_max_FVA":      round(cit_var["maximum"],4),
                                        "growth_under_cit": round(gr_cit,      4),
                                        # croissance
                                        "growth_max":       round(growth_max,  4),
                                        # "min_growth":       round(gr_var["minimum"],4),
                                        # "max_growth":       round(gr_var["maximum"],4),
                                        # sidérophore
                                        "siderophore_FBA":  round(sid_fba,    4),
                                        "siderophore_pFBA": round(sid_pfb,    4)
                                        #"sid_min_FVA":      round(sid_var["minimum"],4),
                                        #"sid_max_FVA":      round(sid_var["maximum"],4),
                                    })

# 4) DataFrame et CSV
df = pd.DataFrame(results)
df.to_csv(
    "C:/Users/aliso/OneDrive - Université Laval/iGEM_E.coli/citrate_growth_siderophore_scan.csv",
    index=False
)

# --- Résumé rapide ---
print(f"Scénarios simulés : {len(df)}\n")

print("Top 5 prolifération (growth_max):")
print(df.sort_values('growth_max', ascending=False).head(5), "\n")

print("Top 10 export citrate (citrate_FBA):")
print(df.sort_values('citrate_FBA', ascending=False).head(10), "\n")

print("Top 5 production de sidérophore (FBA):")
print(df.sort_values('siderophore_FBA', ascending=False).head(5), "\n")


# --- Analyses post-simulation ---

# 1) Top 10 conditions pour croissance
top10_growth = df.sort_values('growth_max', ascending=False).head(10)
print("Top 10 scénarios (growth_max):")
print(top10_growth, "\n")

# 2) Top 10 conditions pour export citrate (FBA)
top10_citrate = df.sort_values('citrate_FBA', ascending=False).head(10)
print("Top 10 scénarios (citrate_FBA):")
print(top10_citrate, "\n")

# 3) Taux moyen par KO
mean_by_ko = (
    df.groupby('KO')[['growth_max', 'citrate_FBA']]
      .mean()
      .reset_index()
)
print("Taux moyen par KO:")
print(mean_by_ko, "\n")

# 4) Taux moyen par fraction de croissance
mean_by_growth = (
    df.groupby('growth_%')[['growth_max', 'citrate_FBA']]
      .mean()
      .reset_index()
)
print("Taux moyen par growth_%:")
print(mean_by_growth, "\n")

# 5) Taux moyen par ion
mean_by_ion = (
    df.groupby('ion')[['growth_max', 'citrate_FBA']]
      .mean()
      .reset_index()
      .sort_values('growth_max', ascending=False)
)
print("Taux moyen par ion:")
print(mean_by_ion, "\n")

# 6) Taux moyen par température
mean_by_temp = (
    df.groupby('température')[['growth_max', 'citrate_FBA']]
      .mean()
      .reset_index()
      .sort_values('growth_max', ascending=False)
)
print("Taux moyen par température:")
print(mean_by_temp, "\n")

# Calcul des moyennes par source de carbone (comme défini précédemment)
mean_by_carbon = (
    df.groupby('C_source')[['growth_max','citrate_FBA','siderophore_FBA']]
      .mean()
      .reset_index()
      .sort_values('growth_max', ascending=False)
)

# Affichage complet
print("Taux moyen par source de carbone :")
print(mean_by_carbon, "\n")

# Si tu veux ne voir que les 5 meilleures sources pour la croissance :
print("Top 5 sources de carbone pour growth_max :")
print(mean_by_carbon[['C_source','growth_max']].head(5), "\n")

# Et pour voir les scénarios individuels les plus performants (croissance) :
print("Top 10 scénarios (croissance) avec leur source de carbone :")
print(
    df.sort_values('growth_max', ascending=False)
      .head(10)[['model','C_source','growth_max','pH','ion','température','KO']],
    "\n"
)

# Début : 11h37
# Fin : 

# Début narval: 11h24
# Fin narval: 