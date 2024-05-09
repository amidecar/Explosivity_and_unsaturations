from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
IPythonConsole.ipython_useSVG=True

jaaj=str(input("What is the smiles of your molecule?")) #input of the molecule smiles
#still need to had the case of the iupac name
#function below canonicalize the smile
def canonicalize_smiles(smiles: str) -> str:
    if not isinstance(smiles, str):
        raise TypeError("Invalid type {type(smiles)}: smiles must be a string")
    
    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        raise ValueError("Could not convert smiles to mol")

    return Chem.MolToSmiles(mol)
jaaj=canonicalize_smiles(jaaj)
molecule= Chem.MolFromSmiles(jaaj)
molecule
highlight_boom={}

explosive_groups=["C#C","C=CC=C","C[Mg]","C[Li]","C[Cu]","N[Li]","N=N=N","CN=NC","C[N+]#N","NN"
                  "OO","OOO","NO","N(O)O","NI","NBr","NCl","NF","Cl=O","F=O","I=O","Br=O"]
for group in explosive_groups:
    pattern= Chem.MolFromSmiles(group)
    matches= molecule.GetSubstructMatches(pattern)
    highlight_boom[group]= matches
print(highlight_boom)
    
    
    
