import webbrowser
from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Descriptors 
import urllib.request
import json
import tkinter as tk
from PIL import Image, ImageTk
import math

def balox(smiles):                        
    moleculee = Chem.MolFromSmiles(smiles)              # Convert the SMILES string to a molecule object
    moleculee = Chem.AddHs(moleculee)                   # Add hydrogens
    atom_counts = {"C":0, "H":0,"O":0}                  # Initialise new dict
    
    # Loop through atoms in the molecule and count different atom types
    for atom in moleculee.GetAtoms():
        atom_symbol = atom.GetSymbol()
        if atom_symbol in atom_counts:
            atom_counts[atom_symbol] += 1
        else:
            atom_counts[atom_symbol] =1
    molarmass = Descriptors.MolWt(moleculee)
    
    #calculate oxygen balance from Marendaz formula
    return -1600 * (2*atom_counts["C"] + atom_counts["H"]/2 - atom_counts["O"])/molarmass

def canonicalize_smiles(smiles: str) -> str:
    if not isinstance(smiles, str):
        raise TypeError("Invalid type {type(smiles)}: smiles must be a string")
    
    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        raise ValueError("Could not convert smiles to mol")

    return Chem.MolToSmiles(mol)

def findgroups(smiles:str):
    global group_check
    group_check=0
    highlight_boom={}
    molecule= Chem.MolFromSmiles(smiles)
    explosive_groups=["[#6]#[#6]",              #carbide
                      "[#6]=[#6]-[#6]=[#6]",    #diene
                      "[#6-]",                  #carbanion
                      "[#6]-[Mg]","[#6]-[Li]","[#6]-[Cu]","[#6]-[Ni]", #grignard
                      "[#7-]-[#7+]#[#7]","[#7]=[#7+]=[#7-]", #azide
                      "[#7]=[#7]",              #azo
                      "[#7-]#[#7]",             #diazonium salt
                      "[#7]-[#7]",              #hydrazine
                      "[#8]-[#8]",              #peroxyde
                      "[#8]-[#8]-[#8]",         #ozonide
                      "[#6]-1-[#8]-[#8]-[#6]-[#8]-1",       #ozonide
                      "[#7]-[#8]","[#7]=[#8]","[#7]#[#8]",   #Nitro
                      "[#7]-[I]",               #iodamine
                      "[#7]-[Br]",              #bromamine
                      "[#7]-[Cl]",              #chloramine
                      "[#7]-[F]",               #fluoramine
                      "[#7]=[I]",               
                      "[#7]=[Br]",   
                      "[#7]=[Cl]",
                      "[#7]=[F]",
                      "[#8]-[I]",               #iodate
                      "[#8]-[F]",               #fluorate
                      "[#8]-[Br]",              #bromate
                      "[#8]-[Cl]",              #chlorate
                      "[#8]=[I]",
                      "[#8]=[Br]",
                      "[#8]=[Cl]",
                      "[#8]=[F]",
                      "[#6]#[#7+]-[#8-]"]       #fulminate
    for group in explosive_groups:
        pattern= Chem.MolFromSmarts(group)
        matches= molecule.GetSubstructMatches(pattern)
        if matches != ():
            group_check+=1
        highlight_boom[group]= matches
    return highlight_boom

def highlightmol(smiles,dico):
    from rdkit.Chem import Draw
    mol= Chem.MolFromSmiles(smiles)
    highlight_atoms = set()
    for matches in dico.values():
        for match in matches:
            highlight_atoms.update(match)
    img = Draw.MolToImage(mol, highlightAtoms=list(highlight_atoms))
    
    img_tk = ImageTk.PhotoImage(img)
    
    # Create a label to display the image
    global label
    label = tk.Label(window, image=img_tk)
    label.place(relx=0.5,rely=0.5,anchor="center")
    window.mainloop()

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

def explosivity(oxygen_balance):
    global text1
    ox = math.floor(oxygen_balance)
    
    if group_check == 0:
        text1="No explosible groups, the molecule is not explosive"
    else:
        if oxygen_balance <=-200:
            text1 = "the oxygen balance is: "+str(ox) +",the molecule is not explosive"
        elif -200 <oxygen_balance <=-160:
            text1 ="the oxygen balance is "+str(ox) +",the molecule is mildly explosive"
        elif -160<oxygen_balance<=-80 or oxygen_balance>= 24:
            text1 ="the oxygen balance is "+str(ox) +",the molecule is explosive"
        elif -80 <oxygen_balance <=12:
            text1 ="the oxygen balance is "+str(ox) +",the molecule is very explosive"
        else:
            print("idk man, good luck")
    print(text1)
    label1 = tk.Label(window, text=text1, fg="black")
    label1.place(relx=0.5,rely=0.25,anchor="center")


def exclude_self():
    if name_var.get():
        smile_var.set(False)
        entryname.place(x=220,rely=0)  # Show the entry widget if "Name" is selected
        entrysmiles.pack_forget()
        window.update()
    if smile_var.get():
        name_var.set(False)
        entrysmiles.place(x=220,rely=0.05)
        entryname.pack_forget()  # Hide the entry widget if "Smile" is selected
        window.update()


def submit():
    # Retrieve the input text when the submit button is clicked
    if name_var.get():
        print("Name is selected. Entered Name:", enteredname.get())
        #label.pack_forget()
        smiles = iupac_to_smiles(enteredname.get())

        
    elif smile_var.get():
        print("Entered smile :", enteredsmiles.get())
        #label.pack_forget()
        smiles = enteredsmiles.get()


    else:
        print("Neither Name nor Smile is selected.")
    
    smiles=canonicalize_smiles(smiles) #FONCTION QUI CANONISE LES SMILES ET FAIT UNE ERREUR SI LE SMILES EST MAUVAIS
    dico=findgroups(smiles)
    oxygen_balance = balox(smiles)
    explosivity(oxygen_balance)
    highlightmol(smiles,dico)
    window.mainloop()

def open_pubchem_sketcher():
    # Open PubChem Sketcher in a new tab
    url = "https://pubchem.ncbi.nlm.nih.gov/edit3/index.html"
    webbrowser.open_new_tab(url)

window = tk.Tk()
name_var = tk.BooleanVar()
smile_var = tk.BooleanVar()
name_button = tk.Checkbutton(window, text="Check to enter the molecule name", variable=name_var, command=exclude_self)
smile_button = tk.Checkbutton(window, text="Check to enter the molecule smiles", variable=smile_var, command=exclude_self)
enteredname = tk.StringVar()
entryname = tk.Entry(window, textvariable=enteredname)
enteredsmiles = tk.StringVar()
entrysmiles = tk.Entry(window, textvariable=enteredsmiles)
    # Create a submit button
submit_button = tk.Button(window, text="Submit", command=submit)
pubchem_button = tk.Button(window, text="Open PubChem Sketcher", command=open_pubchem_sketcher)

#label1 = tk.Label(window, text=text1)

def main():
    

    window.title("Input Window")
    window.geometry("700x800")
    window.minsize(350,640)

    # Create variables to track the state of the check button
    name_button.place(x=0,y=0)


    smile_button.place(x=0,rely=0.05)
    pubchem_button.place(x=120,rely=0.125)

    # Create an entry widget for entering the name

    submit_button.place(x=0,rely=0.125)
    

    # Run the Tkinter event loop
    window.mainloop()
if __name__ == '__main__':
    smiles = ""
    main()
    
    


