import sys
sys.path.append('/Users/wang-work/matagent')

from mcp_server.tools.basic.structure_utils import filter_molecules_by_substructure

try:
    smiles_list = ["CCO", "CCC", "CC(=O)O", "C1=CC=CC=C1"]

    # Test 1: Filter for alcohol group (O)
    smarts_alcohol = "[OH]"
    result_alcohol = filter_molecules_by_substructure(smiles_list, smarts_alcohol)
    print(f"Molecules with alcohol group ({smarts_alcohol}): {result_alcohol}")

    # Test 2: Filter for benzene ring (c1ccccc1)
    smarts_benzene = "c1ccccc1"
    result_benzene = filter_molecules_by_substructure(smiles_list, smarts_benzene)
    print(f"Molecules with benzene ring ({smarts_benzene}): {result_benzene}")

    # Test 3: Filter for carboxylic acid group (C(=O)O)
    smarts_carboxylic = "C(=O)O"
    result_carboxylic = filter_molecules_by_substructure(smiles_list, smarts_carboxylic)
    print(f"Molecules with carboxylic acid group ({smarts_carboxylic}): {result_carboxylic}")

    # Test 4: Invalid SMARTS pattern
    invalid_smarts = "[XYZ]"
    try:
        filter_molecules_by_substructure(smiles_list, invalid_smarts)
    except ValueError as e:
        print(f"Caught expected error for invalid SMARTS: {e}")

    # Test 5: Empty smiles_list
    empty_list_result = filter_molecules_by_substructure([], smarts_alcohol)
    print(f"Molecules from empty list: {empty_list_result}")

    # Test 6: smiles_list with invalid SMILES
    smiles_list_with_invalid = ["CCO", "invalid_smiles", "CC(=O)O"]
    result_with_invalid = filter_molecules_by_substructure(smiles_list_with_invalid, smarts_alcohol)
    print(f"Molecules with alcohol group (with invalid SMILES in list): {result_with_invalid}")

except Exception as e:
    print(f"Test failed with error: {e}")
