from explosivity_and_unsaturations.example_module import balox


# Test the function
def test_hello_smiles():
    assert balox("O=C=O") == 0, "Test failed: SMILES input"
