from rdkit import Chem
from rdkit.Chem import AllChem
from pydantic import BaseModel, Field
from typing import List

from mcp_server.core.registry import register_basic_tool
from mcp_server.core.schemas import BasicToolMetadata

def _get_mol(smiles_string: str):
    """Helper function to convert SMILES to a Mol object, with error handling."""
    mol = Chem.MolFromSmiles(smiles_string)
    if mol is None:
        raise ValueError(f"Invalid SMILES string provided: {smiles_string}")
    return mol

class GenerateMorganFingerprintInput(BaseModel):
    smiles_string: str = Field(..., description="The SMILES representation of the molecule.")

@register_basic_tool(
    metadata=BasicToolMetadata(
        name="generate_morgan_fingerprint",
        description="Generates a Morgan fingerprint (ECFP2) for a given SMILES string.",
        tags=["rdkit", "chemistry", "fingerprint"],
        author="RDKit Wrapper",
        parameters=GenerateMorganFingerprintInput.model_json_schema()
    )
)
def generate_morgan_fingerprint(input: GenerateMorganFingerprintInput) -> List[int]:
    """
    Generates a Morgan fingerprint (ECFP2) for a given SMILES string.
    Returns the fingerprint as a list of integers (bit vector).
    """
    mol = _get_mol(input.smiles_string)
    # Generate Morgan fingerprint with radius 2 and 2048 bits
    fingerprint = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
    return list(fingerprint.ToBitString())