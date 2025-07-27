import requests
from typing import List, Optional
from pathlib import Path

def fetch_models_based_on_org(url:str, outdir:[Path]=None, organisms:Optional[List[str]]=None) -> None:
    """ 
    Fetch GEModels from the given database (URL)

    ARGS:
        - url <STRING>: URL of the database from which to request the models
        - outdir <STRING>: Output directory to which save the fetched models
        - organisms: <LIST(STRING)>: Organism names for which to get the appropriate models  
    """
    # Template URL for model download given from BIGG MODELS
    DOWNLOAD_URL_TEMPLATE = 'http://bigg.ucsd.edu/api/v2/models/MODEL_ID/download'
     
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    
    response = make_request(url=url)
    if not response:
        return

    models = response.json().get("results", [])

    # Filter for specific organisms if any
    selected_models = [m for m in models if any([org.lower() in m["organism"].lower() for org in organisms])] if organisms else models
    for model in selected_models:
        model_id = model['bigg_id']
        organism = model['organism']
        download_url = DOWNLOAD_URL_TEMPLATE.replace('MODEL_ID', model_id)

        print(f"Downloading model: {model_id} ({organism})")

        try:
            req = requests.get(download_url)
            req.raise_for_status()
            
            if outdir:
                file_path = outdir / f"{model_id}.json"
                with open(file_path, "wb") as f:
                    f.write(req.content)
                print(f"Saved: {file_path}")
        except requests.exceptions.RequestException as e:
            print(f"Failed to download {model_id}: {e}")

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
    fetch_models_based_on_org(url=url, outdir='models', organisms=['Escherichia coli'])