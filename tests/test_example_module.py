from explosivity_and_unsaturations.example_module import balox, insat, ismetal, iupac_to_smiles
from explosivity_and_unsaturations.example_module import explosivity,canonicalize_smiles, findinsaturation, findgroups, main
from rdkit import Chem


# Test the function
def test_balox():
    assert balox("O=C=O") == 0, "Test failed: the wrong value was returned"
    assert balox("C") == -398.92788131895543, "Test failed: the wrong value was returned"
    assert balox("[O-][O+]=O") == 100.00625039064941, "Test failed: the wrong value was returned"

def test_insat():
    assert insat("O=O") == 1, "Test failed: the wrong value was returned"
    assert insat("C") == 0, "Test failed: the wrong value was returned"
    assert insat("C#C") == 2, "Test failed: the wrong value was returned"
    assert insat("C1CCCCC1") == 1,"Test failed: the wrong value was returned"
    assert insat("[NH4+].[N+](=O)([O-])[O-]") == 1, "Test failed: the wrong value was returned"
    assert insat("[Al].[Al].[Al].[Al].[Al].[Al].[Al].[Al].[Al].[Al].[Al].[Al].[Al].[Al].[Al]") == 0, "Test failed: the wrong value was returned"
    
def test_ismetal():
    assert ismetal(Chem.Atom("C")) == False,"Test failed: the wrong value was returned"
    assert ismetal(Chem.Atom("Co")) == True,"Test failed: the wrong value was returned"
    assert ismetal(Chem.Atom("Ge")) == True, "Test failed: the wrong value was returned"
    assert ismetal(Chem.Atom("I")) == False, "Test failed: the wrong value was returned"
    
def test_iupac_to_smiles():
    assert iupac_to_smiles("methane") == "C","Test failed: the wrong SMILES was returned"
    assert iupac_to_smiles("trinitrotoluene") == "Cc1ccc(c(c1[N+]([O-])=O)[N+]([O-])=O)[N+]([O-])=O","Test failed: the wrong SMILES was returned"
    assert iupac_to_smiles(("water")) == "O","Test failed: the wrong SMILES was returned"

def test_explosivity():
    assert explosivity(0,0) == "No explosible groups, the molecule is not explosive.", "Test failed: the wrong value was returned"
    assert explosivity(0,1) ==  "              The oxygen balance is 0,the molecule is very explosive.              ", "Test failed: the wrong value was returned"
    assert explosivity(100,1) == "              The oxygen balance is 100,the molecule is explosive.              ","Test failed: the wrong value was returned"
    assert explosivity(-180,567890) == "              The oxygen balance is -180,the molecule is a mildly explosive.              ","Test failed: the wrong value was returned"
    assert explosivity(-500,2) == "              The oxygen balance is -500,the molecule is not explosive.              ","Test failed: the wrong value was returned"
    assert explosivity(-100.567,9) == "              The oxygen balance is -101,the molecule is explosive.              ","Test failed: the wrong value was returned"
    
def test_canonicalize_smiles():
    assert canonicalize_smiles("C") == "C", "Test failed: the wrong SMILES was returned"
    assert canonicalize_smiles("C(=O)=O") == "O=C=O", "Test failed: the wrong SMILES was returned"
    assert canonicalize_smiles("C1C=CC=CC=1") == "c1ccccc1", "Test failed: the wrong SMILES was returned"

def test_findinsaturation():
    assert findinsaturation("C") == {'[*]=[*]': (), '[*]#[*]': (), 'aa': (), '[*]@[*]': ()}, "Test failed: the wrong groups were returned"
    assert findinsaturation("c1ccccc1") == {'[*]=[*]': (), '[*]#[*]': (),
                                            'aa': ((0, 1), (0, 5), (1, 2), (2, 3), (3, 4), (4, 5)),
                                            '[*]@[*]': ((0, 1), (0, 5), (1, 2), (2, 3), (3, 4), (4, 5))}, "Test failed: the wrong groups were returned"
    assert findinsaturation("C#C") == {'[*]=[*]': (), '[*]#[*]': ((0, 1),), 'aa': (), '[*]@[*]': ()}, "Test failed: the wrong groups were returned"
    assert findinsaturation("O=O") == {'[*]=[*]': ((0, 1),), '[*]#[*]': (), 'aa': (), '[*]@[*]': ()}, "Test failed: the wrong groups were returned" 
    
def test_findgroups():
    assert findgroups("C") == ({'[#6]#[#6]': (), '[#6]=[#6]-[#6]=[#6]': (), '[#6-]': (), '[#6]-[Mg]': (),
                                '[#6]-[Li]': (), '[#6]-[Cu]': (), '[#6]-[Ni]': (), '[#7-]-[#7+]#[#7]': (),
                                '[#7]=[#7+]=[#7-]': (), '[#7]-[#7+]=[#7-][#7]=[#7]': (), 'nn': (),
                                '[#7-]#[#7]': (), '[#7]-[#7]': (), '[#8]-[#8]': (),
                                '[#8]-[#8]-[#8]': (), '[#6]-1-[#8]-[#8]-[#6]-[#8]-1': (),
                                '[#7]-[#8]': (), '[#7]=[#8]': (), '[#7]#[#8]': (), '[#7]-[I]': (),
                                '[#7]-[Br]': (), '[#7]-[Cl]': (), '[#7]-[F]': (), '[#7]=[I]': (),
                                '[#7]=[Br]': (), '[#7]=[Cl]': (), '[#7]=[F]': (), '[#8]-[I]': (),
                                '[#8]-[F]': (), '[#8]-[Br]': (), '[#8]-[Cl]': (), '[#8]=[I]': (),
                                '[#8]=[Br]': (), '[#8]=[Cl]': (), '[#8]=[F]': (), '[#6]#[#7+]-[#8-]': (),
                                '[Ag+].[N-3]': ()}, 0), "Test failed: the wrong groups were returned"
    assert findgroups("C(C(CCCC(C(C(CC[N+](=O)[O-])[Cu])C(C[Li])[Ni])[Mg])C(C#C)C[C-])=CC=C") == ({'[#6]#[#6]': ((20, 21),),
                                                                                                  '[#6]=[#6]-[#6]=[#6]': ((0, 24, 25, 26),),
                                                                                                  '[#6-]': ((23,),), '[#6]-[Mg]': ((5, 18),), '[#6]-[Li]': ((15, 16),), '[#6]-[Cu]': ((7, 13),), '[#6]-[Ni]': ((14, 17),), '[#7-]-[#7+]#[#7]': (), '[#7]=[#7+]=[#7-]': (), '[#7]-[#7+]=[#7-][#7]=[#7]': (), 'nn': (), '[#7-]#[#7]': (), '[#7]-[#7]': (), '[#8]-[#8]': (), '[#8]-[#8]-[#8]': (), '[#6]-1-[#8]-[#8]-[#6]-[#8]-1': (), '[#7]-[#8]': ((10, 12),), '[#7]=[#8]': ((10, 11),), '[#7]#[#8]': (), '[#7]-[I]': (), '[#7]-[Br]': (), '[#7]-[Cl]': (), '[#7]-[F]': (), '[#7]=[I]': (), '[#7]=[Br]': (), '[#7]=[Cl]': (), '[#7]=[F]': (), '[#8]-[I]': (), '[#8]-[F]': (), '[#8]-[Br]': (), '[#8]-[Cl]': (), '[#8]=[I]': (), '[#8]=[Br]': (), '[#8]=[Cl]': (), '[#8]=[F]': (), '[#6]#[#7+]-[#8-]': (),
                                                                                                  '[Ag+].[N-3]': ()}, 9), "Test failed: the wrong groups were returned"
    assert findgroups("C(C(C(C(C(C(C(C(C(C(C(C(C(C(C(C(C(C)N=N)OO)OF)OC)OCl)OC)OBr)OC)OI)ON)NF)NO)NCl)[N+]#N)NBr)NN)(C(C)C#[N+][O-])NI") == ({'[#6]#[#6]': (),
                                                                                                                                              '[#6]=[#6]-[#6]=[#6]': (), '[#6-]': (), '[#6]-[Mg]': (), '[#6]-[Li]': (), '[#6]-[Cu]': (), '[#6]-[Ni]': (), '[#7-]-[#7+]#[#7]': (), '[#7]=[#7+]=[#7-]': (), '[#7]-[#7+]=[#7-][#7]=[#7]': (), 'nn': (), '[#7-]#[#7]': (), '[#7]-[#7]': ((48, 49),), '[#8]-[#8]': ((20, 21),), '[#8]-[#8]-[#8]': (), '[#6]-1-[#8]-[#8]-[#6]-[#8]-1': (), '[#7]-[#8]': ((37, 36), (40, 41), (53, 54)), '[#7]=[#8]': (), '[#7]#[#8]': (), '[#7]-[I]': ((55, 56),), '[#7]-[Br]': ((46, 47),), '[#7]-[Cl]': ((42, 43),), '[#7]-[F]': ((38, 39),), '[#7]=[I]': (), '[#7]=[Br]': (), '[#7]=[Cl]': (), '[#7]=[F]': (), '[#8]-[I]': ((34, 35),), '[#8]-[F]': ((22, 23),), '[#8]-[Br]': ((30, 31),), '[#8]-[Cl]': ((26, 27),), '[#8]=[I]': (), '[#8]=[Br]': (), '[#8]=[Cl]': (), '[#8]=[F]': (), '[#6]#[#7+]-[#8-]': ((52, 53, 54),), '[Ag+].[N-3]': ()},
                                                                                                                                             12), "Test failed: the wrong groups were returned"
    assert findgroups("[Ag+].[N-3]") == ({'[#6]#[#6]': (), '[#6]=[#6]-[#6]=[#6]': (), '[#6-]': (), '[#6]-[Mg]': (), '[#6]-[Li]': (), '[#6]-[Cu]': (), '[#6]-[Ni]': (), '[#7-]-[#7+]#[#7]': (), '[#7]=[#7+]=[#7-]': (), '[#7]-[#7+]=[#7-][#7]=[#7]': (), 'nn': (), '[#7-]#[#7]': (), '[#7]-[#7]': (), '[#8]-[#8]': (), '[#8]-[#8]-[#8]': (), '[#6]-1-[#8]-[#8]-[#6]-[#8]-1': (), '[#7]-[#8]': (), '[#7]=[#8]': (), '[#7]#[#8]': (), '[#7]-[I]': (), '[#7]-[Br]': (), '[#7]-[Cl]': (), '[#7]-[F]': (), '[#7]=[I]': (),
                                          '[#7]=[Br]': (), '[#7]=[Cl]': (), '[#7]=[F]': (), '[#8]-[I]': (), '[#8]-[F]': (), '[#8]-[Br]': (), '[#8]-[Cl]': (), '[#8]=[I]': (), '[#8]=[Br]': (), '[#8]=[Cl]': (), '[#8]=[F]': (), '[#6]#[#7+]-[#8-]': (), '[Ag+].[N-3]': ((0, 1),)},
                                         1), "Test failed: the wrong groups were returned"

    
def test_main():
    assert main() == None, "Test failed, should return nothing"
    
