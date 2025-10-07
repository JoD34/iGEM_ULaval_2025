import requests
from typing import List, Optional
from pathlib import Path

def fetch_models_based_on_org(url:str, outdir:Optional[Path]=None, organisms:Optional[List[str]]=None) -> None:
    """ 
    Fetch GEModels from the given database (URL)

    ARGS:
        - url <STRING>: URL of the database from which to request the models
        - outdir <STRING>: Output directory to which save the fetched models
        - organisms: <LIST(STRING)>: Organism names for which to get the appropriate models  
    """
    # Template URL for model download given from BIGG MODELS
    DOWNLOAD_URL_TEMPLATE = 'http://bigg.ucsd.edu/api/v2/models/MODEL_ID/download'
     
    outdir = Path(outdir or "./models")
    outdir.mkdir(parents=True, exist_ok=True)

    response = make_request(url)
    if not response:
        print("Failed to retrieve model list.")
        return

    models = response.json().get("results", [])
    if not models:
        print("No models found in the response.")
        return

    selected_models = [m for m in models if any(org.lower() in m["organism"].lower() for org in organisms)] if organisms else models
    
    if not selected_models:
        print("No models matched the specified organisms.")
        return

    info = []

    for model in selected_models:
        model_id = model['bigg_id']
        organism = model['organism']
        download_url = DOWNLOAD_URL_TEMPLATE.replace('MODEL_ID', model_id)
        file_path = outdir / f"{model_id}.json"

        print(f"Downloading model: {model_id} ({organism})")

        try:
            req = requests.get(download_url)
            req.raise_for_status()
            
            with open(file_path, "wb") as f: f.write(req.content)
            print(f"Saved: {file_path}")

            info.append((model_id,organism))
            break
        except requests.exceptions.RequestException as e:
            print(f"Failed to download {model_id}: {e}")

    description_file = "models_description.txt"
    with open(outdir / description_file, 'w') as f:
        for model_id, org in info:
            f.write(f'{model_id} -> {org}\n')

    print(f'Successfully written {description_file}')
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
    fetch_models_based_on_org(url=url, outdir='models', organisms=['Escherichia coli'])
