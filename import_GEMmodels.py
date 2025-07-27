import requests

def fetch_models(url:str, outdir:str, organism:[str]=None):
    response = requests.get(url)
    
    try:
        response = requests.get(url=url)
        response.raise_for_status() 
        print("Request successful!")
    except requests.exceptions.HTTPError as err:
        print(f"HTTP Error: {err}")
        
    models = response.json()["results"]
    print(models)
    
    # Filter for E. coli models
    #selected_models = [m for m in models if "coli" in m["organism"].lower()] if organism else list(m)
    #for model in selected_models:
    #    print(model["bigg_id"], "-", model["organism"])
    
if __name__=='__main__':
    url = "http://bigg.ucsd.edu/api/v2/models"
    fetch_models