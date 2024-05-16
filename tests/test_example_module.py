from explosivity_and_unsaturations.example_module import balox, insat


# Test the function
def test_hello_smiles():
    assert balox("O=C=O") == 0, "Test failed: SMILES input"

def test_jaaj():
    assert insat("O=C=O") == 2, "Test failed: wrong value"