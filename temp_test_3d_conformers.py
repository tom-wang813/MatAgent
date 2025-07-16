import sys
sys.path.append('/Users/wang-work/matagent')

from mcp_server.tools.basic.structure_utils import generate_3d_conformers

try:
    smiles = "CCO" # Ethanol
    sdf_content = generate_3d_conformers(smiles, num_conformers=3)
    print(f"Generated 3D Conformers SDF (first 500 chars):\n{sdf_content[:500]}...")

    smiles_benzene = "C1=CC=CC=C1" # Benzene
    sdf_content_benzene = generate_3d_conformers(smiles_benzene, num_conformers=1)
    print(f"\nGenerated 3D Conformers SDF for Benzene (first 500 chars):\n{sdf_content_benzene[:500]}...")

    invalid_smiles = "invalid_smiles"
    error_conformers = generate_3d_conformers(invalid_smiles, num_conformers=1)
    print(f"Error Conformers: {error_conformers}")
except Exception as e:
    print(f"Test failed with error: {e}")
