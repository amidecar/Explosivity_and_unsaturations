import webbrowser
from rdkit import Chem
from rdkit.Chem import Descriptors 
import urllib.request
import tkinter as tk
from PIL import ImageTk
import math
from tkinter import ttk
from tkinter import messagebox
from tkinter import *

def ismetal(atom):
    
    """ Finds if an atom is a metal
    Input: an rdkit atom object
    Returns: True if the atom is a metal, False otherwise 
    """
    
    n = atom.GetAtomicNum()
    if((n>=19 and n<=32) or (n>=37 and n<=51) or (n>=55 and n<=84) or n==3 or n== 4 or n==11 or n==12 or n==13 or (n>=87 and n<=118)):
        return True
    return False

def balox(smiles):         
    
    """This function calculates the oxygen balance of the molecule"""
    
    moleculee = Chem.MolFromSmiles(smiles)              # Convert the SMILES string to a molecule object
    moleculee = Chem.AddHs(moleculee)                   # Add hydrogens
    atom_counts = {"C":0, "H":0,"O":0,"Metal":0}        # Initialise new dict
    
    # Loop through atoms in the molecule and count different atom types
    for atom in moleculee.GetAtoms():
        atom_symbol = atom.GetSymbol()
        if atom_symbol in atom_counts:
            atom_counts[atom_symbol] += 1
        else:
            atom_counts[atom_symbol] =1

        if(ismetal(atom)):
            atom_counts["Metal"]+=1

    molarmass = Descriptors.MolWt(moleculee)
    
    #calculate oxygen balance from Marendaz formula
    return -1600 * (2*atom_counts["C"] + atom_counts["H"]/2 + atom_counts["Metal"] - atom_counts["O"])/molarmass



def canonicalize_smiles(smiles: str, ismain=False) -> str:
    
    """this function canonicalises the smiles given. Some molecules have more than 1 smile to describe them,
    that they are recgnized by the otehr functions"""
    
    if not isinstance(smiles, str):
        if ismain:
            messagebox.showerror('Warning!', 'Invalid type {type(smiles)}: smiles must be a string')
        raise TypeError("Invalid type {type(smiles)}: smiles must be a string")
        
    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        if ismain:
            messagebox.showerror('Warning!', "Could not convert input to mol")
        raise ValueError("Could not convert input to mol")
    return Chem.MolToSmiles(mol)


def findgroups(smiles:str):
    
    """"This function recognizes the possibly explosive groups in the molecule. This is
    then used by the highlightmol function to highlight those groups. It also create
    a global valuable that checks if there is at least one possibly explosive group"""
    
    #group_check
    group_check=0
    highlight_boom={}
    molecule= Chem.MolFromSmiles(smiles)
    explosive_groups=["[#6]#[#6]",              #carbide
                      "[#6]=[#6]-[#6]=[#6]",    #diene
                      "[#6-]",                  #carbanion
                      "[#6]-[Mg]","[#6]-[Li]","[#6]-[Cu]","[#6]-[Ni]", #grignard
                      "[#7-]-[#7+]#[#7]","[#7]=[#7+]=[#7-]","[#7]-[#7+]=[#7-]" #azide
                      "[#7]=[#7]",              #azo
                      "nn",                     #aromatic nitrogen
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
                      "[#6]#[#7+]-[#8-]",       #fulminate
                      "[Ag+].[N-3]"]            #silver nitride       
    
    for group in explosive_groups:
        pattern= Chem.MolFromSmarts(group)
        matches= molecule.GetSubstructMatches(pattern)
        if matches != ():
            group_check+=1
        highlight_boom[group]= matches
    
    return (highlight_boom,group_check)


def insat(smiles):                        

    """ this function calculates the degree of insturaton of a molecule"""

    moleculee = Chem.MolFromSmiles(smiles)              # Convert the SMILES string to a molecule object
    moleculee = Chem.AddHs(moleculee)                   # Add hydrogens
    atom_counts = {}                  # Initialise new dict
    
    # Loop through atoms in the molecule and count different atom types
    for atom in moleculee.GetAtoms():
        atom_valence = atom.GetTotalValence()
        if atom_valence in atom_counts:
            atom_counts[atom_valence] += 1
        else:
            atom_counts[atom_valence] =1
    ddi=0
    for i in atom_counts.keys():
        ddi= ddi + (i-2)*atom_counts[i]
    ddi= (ddi + 2)/2
    if ddi<0:
        ddi=0
    return ddi


def findinsaturation(smiles:str):

    """This function find the insaturations"""
    
    highlight_insat={}
    molecule= Chem.MolFromSmiles(smiles)
    insaturation=["[*]=[*]",
                  "[*]#[*]",
                  "aa",
                  "[*]@[*]"
                  ]
    for group in insaturation:
        pattern= Chem.MolFromSmarts(group)
        matches= molecule.GetSubstructMatches(pattern)
        highlight_insat[group]= matches
    return highlight_insat


def highlightmol(smiles,dico):

    """this function highlights the groups found by other functions"""
    
    from rdkit.Chem import Draw
    mol= Chem.MolFromSmiles(smiles)
    highlight_atoms = set()
    highlight_bonds = []
    for matches in dico.values():
        for match in matches:
            highlight_atoms.update(match)
            for i in range(len(match) - 1):
                bond = mol.GetBondBetweenAtoms(match[i], match[i + 1])
                if bond is not None:
                    highlight_bonds.append(bond.GetIdx())  # Store bond indices
            
    img = Draw.MolToImage(mol, highlightAtoms=list(highlight_atoms), highlightBonds=highlight_bonds,kekulize = False)
    #img=img.resize((450, 450))
    #img_tk = ImageTk.PhotoImage(img)
    return img


def iupac_to_smiles(iupac_name):
    
    """this function takes the IUPAC name of the molecule as an input and transforms it
    into a smile for the other function to use."""
    
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


def explosivity(oxygen_balance,group_check):

    """this function uses the group_check and the balox function results to find a 
    molecule likelihood to explode."""
    
    ox = math.floor(oxygen_balance)

    if group_check == 0:
        text1="No explosible groups, the molecule is not explosive."
    else:
        if oxygen_balance <=-200:
            text1 = "              The oxygen balance is "+str(ox) +",the molecule is not explosive.              "
        elif -200 <oxygen_balance <=-160:
            text1 ="              The oxygen balance is "+str(ox) +",the molecule is a mildly explosive.              "
        elif -160<oxygen_balance<=-80 or 12 <oxygen_balance:
            text1 ="              The oxygen balance is "+str(ox) +",the molecule is explosive.              "
        elif -80 <oxygen_balance <=12:
            text1 ="              The oxygen balance is "+str(ox) +",the molecule is very explosive.              "
    return text1


def submitboom():

    """this function makes the interface work"""
    
    for widget in kaboomity.winfo_children():
        if ".!notebook.!frame.!label" in str(widget):
            widget.destroy()
            
    saas = int(aaa.get())
    nameee = enteredname.get()
    smilesee = enteredsmiles.get()
    # Retrieve the input text when the submit button is clicked
    if ((saas==0) and (nameee =="")):
        print("can't be null")
        if ismain:
            messagebox.showerror('Warning!', 'Error: Write something before submitting !')
        return
    elif ((smilesee=="") and (saas==1)):
        print("can't be null")
        if ismain:
           messagebox.showerror('Warning!', 'Error: Write something before submitting !')

    if saas==0:
        print("Name is selected. Entered Name:", enteredname.get())
        #label.pack_forget()
        smiles = iupac_to_smiles(enteredname.get())

        
    elif saas==1:
        print("Entered smile :", enteredsmiles.get())
        #label.pack_forget()
        smiles = enteredsmiles.get()


    else:
        print("Neither Name nor Smile is selected.")
    #print(smiles)
    
        
    smiles=canonicalize_smiles(smiles,True) #FONCTION QUI CANONISE LES SMILES ET FAIT UNE ERREUR SI LE SMILES EST MAUVAIS
    (dico,group_check)=findgroups(smiles)
    oxygen_balance = balox(smiles)
    texte1 =explosivity(oxygen_balance,group_check)
    global label1 
    label1 = tk.Label(kaboomity, text=texte1, fg="black")
    label1.pack_forget()
    label1.place(relx=0.5,rely=0.2, anchor="center")
    aa=highlightmol(smiles,dico)
    aa=aa.resize((450, 450))
    img_tk = ImageTk.PhotoImage(aa)
    global labela
    labela = tk.Label(kaboomity, image=img_tk)
    labela.place(x=45,y=140)
    
    window.mainloop()


def submitinsat():
    # Retrieve the input text when the submit button is clicked

    print("aaaaewdsgxvbansdjiohsja")
    saas = int(aaa.get())
    nameee = enteredname.get()
    smilesee = enteredsmiles.get()
    # Retrieve the input text when the submit button is clicked
    if ((saas==0) and (nameee =="")):
        print("can't be null")
        if ismain:
            messagebox.showerror('Warning!', 'Error: Write something before submitting !')
        return
    elif ((smilesee=="") and (saas==1)):
        print("can't be null")
        if ismain:
            messagebox.showerror('Warning!', 'Error: Write something before submitting !')

    if not aaa.get():
        print("Name is selected. Entered Name:", enteredname.get())
        #label.pack_forget()
        smiles = iupac_to_smiles(enteredname.get())

        
    elif aaa.get():
        print("Entered smile :", enteredsmiles.get())
        #label.pack_forget()
        smiles = enteredsmiles.get()


    else:
        print("Neither Name nor Smile is selected.")
    smiles=canonicalize_smiles(smiles,True) #FONCTION QUI CANONISE LES SMILES ET FAIT UNE ERREUR SI LE SMILES EST MAUVAIS
    dico=findinsaturation(smiles)
    texte2 = "The degree of insaturation is " + str(insat(smiles)) + "."
    global label2 
    label2 = tk.Label(insaturation, text=texte2, fg="black")
    label2.pack_forget()
    label2.place(relx=0.5,rely=0.2, anchor="center")
    bb=highlightmol(smiles,dico)
    bb=bb.resize((450, 450))
    img_tk = ImageTk.PhotoImage(bb)
    global labelb
    labelb = tk.Label(insaturation, image=img_tk)
    labelb.place(x=45,y=140)
    
    window.mainloop()


def open_pubchem_sketcher():

    """This function makes the open puchem sketcher button work"""
    
    # Open PubChem Sketcher in a new tab
    url = "https://pubchem.ncbi.nlm.nih.gov/edit3/index.html"
    webbrowser.open_new_tab(url)




def main():

    """this is the main function,it initialises the tkinter window"""
    
    global window,tabs_container,kaboomity,insaturation,aaa ,name_button,name_button2,smile_button,smile_button2,enteredname,entryname
    global entryname2,enteredsmiles,entrysmiles,entrysmiles2,submit_button,submit_button2,pubchem_button,pubchem_button2
    
    global ismain
    ismain = True #to replace if __name__ == True even if main() is called from elsewere
    
    
    
    
    window = tk.Tk()
    tabs_container=ttk.Notebook(window)
    tabs_container.pack(fill="both", expand=True)
    kaboomity= ttk.Frame(tabs_container)
    insaturation = ttk.Frame(tabs_container)
    tabs_container.add(kaboomity, text= "explosivity")
    tabs_container.add(insaturation, text= "degree of insaturation")
    tabs_container.pack()
    aaa = IntVar()

    name_button = tk.Radiobutton(kaboomity, text="Check to enter the molecule IUPAC name", variable=aaa,value=0)
    name_button2 = tk.Radiobutton(insaturation, text="Check to enter the molecule IUPAC name", variable=aaa,value=0)
    smile_button = tk.Radiobutton(kaboomity, text="Check to enter the molecule smiles", variable=aaa,value=1)
    smile_button2 = tk.Radiobutton(insaturation, text="Check to enter the molecule smiles", variable=aaa,value=1)
    enteredname = tk.StringVar()
    entryname = tk.Entry(kaboomity, textvariable=enteredname)
    entryname2 = tk.Entry(insaturation, textvariable = enteredname)
    enteredsmiles = tk.StringVar()
    entrysmiles = tk.Entry(kaboomity, textvariable=enteredsmiles)
    entrysmiles2 = tk.Entry(insaturation, textvariable=enteredsmiles)
    submit_button = tk.Button(kaboomity, text="Submit", command=submitboom)
    submit_button2 = tk.Button(insaturation, text="Submit", command=submitinsat)
    pubchem_button = tk.Button(kaboomity, text="Open PubChem Sketcher", command=open_pubchem_sketcher)
    pubchem_button2 = tk.Button(insaturation, text="Open PubChem Sketcher", command=open_pubchem_sketcher)
    
    window.title("explositivity and degree of unsaturation")
    window.geometry("540x630")
    window.minsize(540,630)
    window.resizable(False, False)
    

    # Create variables to track the state of the check button
    name_button.place(x=0,y=20)
    name_button2.place(x=0,y=20)
    smile_button.place(x=0,y=40)
    smile_button2.place(x=0,y=40)
    
    pubchem_button.place(x=70,y=70)
    pubchem_button2.place(x=70,y=70)

    entryname.place(x=250,y=25)  # Show the entry widget if "Name" is selected
    entryname2.place(x=250,y=25)


    entrysmiles.place(x=250,y=45)
    entrysmiles2.place(x=250,y=45)

    # Create an entry widget for entering the name

    submit_button.place(x=0,y=70)
    submit_button2.place(x=0,y=70)
    
    #binds the submit buttons to the return key
    kaboomity.bind('<Return>', lambda event=None: submitboom())
    insaturation.bind('<Return>', lambda event=None: submitinsat())
    entryname.bind('<Return>', lambda event=None: submitboom())  
    entryname2.bind('<Return>', lambda event=None: submitinsat())
    entrysmiles.bind('<Return>', lambda event=None: submitboom())
    entrysmiles2.bind('<Return>', lambda event=None: submitinsat())

    # Run the Tkinter event loop
    window.mainloop()


if __name__ == '__main__':
    smiles = ""
    main()
    
    


