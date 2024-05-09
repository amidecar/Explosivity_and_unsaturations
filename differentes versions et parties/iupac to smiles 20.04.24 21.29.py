import urllib.request
import json

def iupac_to_smiles(iupac_name):
    # URL for CIR API
    url = "https://cactus.nci.nih.gov/chemical/structure/" + urllib.parse.quote_plus(iupac_name) + "/smiles"
    
    try:
        # Send request to the CIR API
        response = urllib.request.urlopen(url)
        # Read and decode the response
        smiles = response.read().decode("utf-8").strip()
        return smiles
    except Exception as e:
        print("Error:", e)
        return "Unable to convert IUPAC name to SMILES"

iupac_name = "benzene"
smiles = iupac_to_smiles(iupac_name)
print("SMILES for", iupac_name, ":", smiles)

