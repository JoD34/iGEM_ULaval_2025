import requests

response = requests.get("http://bigg.ucsd.edu/api/v2/models")
models = response.json()["results"]

# Filter for E. coli models
ecoli_models = [m for m in models if "coli" in m["organism"].lower()]
for model in ecoli_models:
    print(model["bigg_id"], "-", model["organism"])