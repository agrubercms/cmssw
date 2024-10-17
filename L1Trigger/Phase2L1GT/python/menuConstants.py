# import yaml
# scalings_fname = "/eos/home-a/alobanov/SWAN_projects/L1T-Ph2-Menu/L1Nano/scalings_v44.yaml"
# with open(scalings_fname, "r") as f:
#     scalings = yaml.safe_load(f)
import json
scalings_fname = "./scalings_v44.json"
with open(scalings_fname, 'r') as file:
    scalings = json.load(file)

def off_to_on_thr(thr, obj, id, region, scalings = scalings):
    offset = scalings[obj][id][region]["offset"]
    slope = scalings[obj][id][region]["slope"]
    return (thr - offset)/slope

objectIDs = {
    "L1nnPuppiTau": {
        "default": 225
    }
}