import requests
from typing import List, Optional

def fetch_models(url:str, outdir:Optional[str]=None, organisms:Optional[List[str]]=None) -> None:
    """ 
    Fetch GEModels from the given database (URL)

    ARGS:
        - url <STRING>: URL of the database from which to request the models
        - outdir <STRING>: Output directory to which save the fetched models
        - organisms: <LIST(STRING)>: Organism names for which to get the appropriate models  
    """
    response = make_request(url=url)
    models = response.json()["results"]

    # Filter for specific organisms if any
    selected_models = [m for m in models if "coli" in m["organism"].lower()] if organisms else models
    for model in selected_models:
        print(model["bigg_id"], "-", model["organism"])

def get_organisms_choice(url:str, organisms:Optional[List[str]]=None, simplified:bool=True) -> Optional[List[str]]:
    """ 
    Get available organism names from the given URL.

    Args:
        url (str): URL of the database from which to request the models.
        organisms (List[str], optional): List of organisms to filter by. Defaults to None.
        simplified (bool, optional): Whether to simplify names to Genus + species. Defaults to True.

    Returns:
        Optional[List[str]]: List of selected organism names or None if request failed.
    
    """
    response = make_request(url=url)
    models = response.json().get("results", [])

    if simplified:
        selected_orgm = set([' '.join(m['organism'].split()[:2]) for m in models])  
    else :
        selected_orgm = set([m["organism"] for m in models])

    # Filter for organisms if a list has been provided
    if organisms is not None:
        selected_orgm = [org for org in selected_orgm if any([o in org for o in organisms]) ]
    else:
        # If simplified, convert set to list for consistent output
        if simplified:
            selected_orgm = list(selected_orgm)
    
    print(f'Your request fetch {len(selected_orgm)} different organism names!')
    for orgm in selected_orgm:
        print(f"- {orgm}")

def make_request(url:str) -> Optional[requests.Response]:
    """
    Request information from provided URL.

    ARGS
        url (str): URL of the database from which to request the models.

    Returns
        Optional[requests.Response]: Information from the successful request to provided URL or None if Request failed
    """
    try:
        response = requests.get(url=url)
        response.raise_for_status() 
        return response

    except requests.exceptions.HTTPError as err:
        print(f"HTTP Error: {err}")
        return None

    except requests.exceptions.RequestException as err:
        print(f"Request Error: {err}")
        return None
    
if __name__=='__main__':
    url = "http://bigg.ucsd.edu/api/v2/models"
    #get_organisms_choice(url=url, organisms = ['Escherichia coli'], simplified=False)
    fetch_models(url=url, organisms=['Escherichia coli'])