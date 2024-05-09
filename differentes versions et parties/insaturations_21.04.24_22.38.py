import webbrowser
from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Descriptors 
import urllib.request
import json

def canonicalize_smiles(smiles: str) -> str:
    if not isinstance(smiles, str):
        raise TypeError("Invalid type {type(smiles)}: smiles must be a string")
    
    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        raise ValueError("Could not convert smiles to mol")

    return Chem.MolToSmiles(mol)

def highlightmol(smiles,dico):
    from rdkit.Chem import Draw
    mol= Chem.MolFromSmiles(smiles)
    highlight_atoms = set()
    for matches in dico.values():
        for match in matches:
            highlight_atoms.update(match)
    img = Draw.MolToImage(mol, highlightAtoms=list(highlight_atoms))
    img.show()
    
def iupac_to_smiles(iupac_name):
    iupac_name_spaceless=iupac_name.replace(" ","%20") # this transforms sapces into the equivalent for URLs
    # URL for CIR API
    url = "https://cactus.nci.nih.gov/chemical/structure/" + urllib.parse.quote_plus(iupac_name_spaceless) + "/smiles"
    
    try:
        # Send request to the CIR API
        response = urllib.request.urlopen(url)
        # Read and decode the response
        smiles = response.read().decode("utf-8").strip()
        return smiles
    except Exception as e:
        print("Error:", e)
        return "Unable to convert IUPAC name to SMILES"

def insat(smiles):                        
    moleculee = Chem.MolFromSmiles(smiles)              # Convert the SMILES string to a molecule object
    moleculee = Chem.AddHs(moleculee)                   # Add hydrogens
    atom_counts = {"C":0, "H":0,"N":0,"F":0,"Br":0,"Cl":0,"I":0}                  # Initialise new dict
    
    # Loop through atoms in the molecule and count different atom types
    for atom in moleculee.GetAtoms():
        atom_symbol = atom.GetSymbol()
        if atom_symbol in atom_counts:
            atom_counts[atom_symbol] += 1
        else:
            atom_counts[atom_symbol] =1
        ddi=(2*atom_counts["C"] + atom_counts["N"] - atom_counts["H"] - atom_counts["F"] 
             - atom_counts["Cl"] - atom_counts["Br"] - atom_counts["I"])/2 +1
        if ddi<0:
            ddi=0
    #calculate the degree of insaturation
    return ddi

def findinsaturation(smiles:str):
    highlight_insat={}
    molecule= Chem.MolFromSmiles(smiles)
    insaturation=["[#6]#[#6]",           #alkynes
                  "[#6]=[#6]",           #alkenes
                  "[#6]=[#7]",           #imine
                  "[#6]#[#7]",          #cyanide
                  "[#6]=[#8]",]          #C=O double bonds
    for group in insaturation:
        pattern= Chem.MolFromSmarts(group)
        matches= molecule.GetSubstructMatches(pattern)
        highlight_insat[group]= matches
    return highlight_insat



    
def main():
    input1 = input("type [n] for entering name, [s] for smiles: ")
    if(input1=="s"):
        input2 = input("Need a sketcher to find smiles ? \n type[y]/[n]: ")
        if(input2=="y"):
            webbrowser.open_new("https://pubchem.ncbi.nlm.nih.gov/edit3/index.html")
            smiles = input("Type the smiles found using the sketcher: ")
        elif(input2=="n"):
            smiles = input("Type the smiles: ")
        else:
            print("invalid answer")
            return
    elif(input1=="n"):
        molecule_name= input("Type the IUPAC name:")
        smiles = iupac_to_smiles(molecule_name)
    else:
        print("invalid answer")
        return
    if smiles=="":
        print("empty molecule")
        return
    smiles=canonicalize_smiles(smiles) #FONCTION QUI CANONISE LES SMILES ET FAIT UNE ERREUR SI LE SMILES EST MAUVAIS
    print("the degree of insaturation of the molecule is:",insat(smiles))
    dico=findinsaturation(smiles)
    highlightmol(smiles,dico)
    
    
if __name__ == '__main__':
    smiles = ""
    main()