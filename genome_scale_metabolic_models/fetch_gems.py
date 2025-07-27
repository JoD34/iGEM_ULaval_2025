import requests
from typing import List, Optional

def fetch_models(url:str, outdir:Optional[str]=None, organisms:Optional[List[str]]=None):
    """ 
    Fetch GEModels from the given database (URL)

    ARGS:
        - url <STRING>: URL of the database from which to request the models
        - outdir <STRING>: Output directory to which save the fetched models
        - organisms: <LIST(STRING)>: Organism names for which to get the appropriate models  
    """
    try:
        response = requests.get(url=url)
        response.raise_for_status() 
        print("Request successful!")
    except requests.exceptions.HTTPError as err:
        print(f"HTTP Error: {err}")
        return
        
    models = response.json()["results"]
    print(f"Fetched {len(models)} models.")

    # Filter for E. coli models
    selected_models = [m for m in models if "coli" in m["organism"].lower()] if organism else list(m)
    for model in selected_models:
        print(model["bigg_id"], "-", model["organism"])

def get_organisms_choice(url:str, simplified:Optional[bool]=True):
    """ 
    Get available organisms name.

    ARGS:
        - url <STRING>: URL of the database from which to request the models 
    """
    try:
        response = requests.get(url=url)
        response.raise_for_status() 
        print("Request successful!")

    except requests.exceptions.HTTPError as err:
        print(f"HTTP Error: {err}")
        return
        
    models = response.json()["results"]
    display_models = set([' '.join(m['organism'].split()[:2]) for m in models]) if simplified else [m["organism"] for m in models]
    for organism in display_models:
        print(f'- {organism}')
    
if __name__=='__main__':
    url = "http://bigg.ucsd.edu/api/v2/models"
    get_organisms_choice(url=url)
    #fetch_models(url=url)