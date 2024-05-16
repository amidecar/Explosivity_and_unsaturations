from explosivity_and_unsaturations.example_module import hello_smiles


# Test the function
def test_hello_smiles():
    assert hello_smiles("C(=O)O") == "Hello, C(=O)O!", "Test failed: SMILES input"
