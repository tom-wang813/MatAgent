from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem
from pydantic import BaseModel, Field
from typing import List

from mcp_server.core.registry import register_basic_tool
from mcp_server.core.schemas import BasicToolMetadata


class GetElementPropertiesInput(BaseModel):
    symbol: str = Field(..., description="The chemical symbol of the element (e.g., 'C', 'Fe').")

@register_basic_tool(
    metadata=BasicToolMetadata(
        name="get_element_properties",
        description="Get properties for a chemical element.",
        parameters=GetElementPropertiesInput.model_json_schema()
    )
)
def get_element_properties(input: GetElementPropertiesInput) -> dict:
    """Returns a dictionary of properties for a given element symbol."""
    try:
        periodic_table = Chem.GetPeriodicTable()
        atomic_num = periodic_table.GetAtomicNumber(input.symbol)
        return {
            "atomic_number": atomic_num,
            "atomic_weight": periodic_table.GetAtomicWeight(input.symbol),
            "element_name": periodic_table.GetElementName(atomic_num),
        }
    except RuntimeError:
        raise ValueError(f"Invalid element symbol: '{input.symbol}'")

# This is a placeholder for a more complex tool.
# In a real scenario, this would query a database.
_VIRTUAL_MOLECULE_DB = {
    "CCO": "Ethanol",
    "C1CCCCC1": "Cyclohexane",
    "c1ccccc1": "Benzene",
    "CC(=O)O": "Acetic Acid",
}

class FindSimilarMoleculesInput(BaseModel):
    smiles_string: str = Field(..., description="The SMILES string of the molecule to search for.")
    similarity_threshold: float = Field(0.7, description="The Tanimoto similarity threshold (0.0 to 1.0).")

@register_basic_tool(
    metadata=BasicToolMetadata(
        name="find_similar_molecules",
        description="Finds molecules in a virtual database similar to the input molecule.",
        parameters=FindSimilarMoleculesInput.model_json_schema()
    )
)
def find_similar_molecules(input: FindSimilarMoleculesInput) -> list:
    """Finds similar molecules using Tanimoto similarity on Morgan fingerprints."""
    try:
        input_mol = Chem.MolFromSmiles(input.smiles_string)
        if input_mol is None:
            raise ValueError("Invalid input SMILES string.")
        input_fp = AllChem.GetMorganFingerprintAsBitVect(input_mol, 2, nBits=2048)

        similar_molecules = []
        for db_smiles, name in _VIRTUAL_MOLECULE_DB.items():
            db_mol = Chem.MolFromSmiles(db_smiles)
            if db_mol is None:
                continue
            db_fp = AllChem.GetMorganFingerprintAsBitVect(db_mol, 2, nBits=2048)
            
            similarity = Chem.DataStructs.TanimotoSimilarity(input_fp, db_fp)
            
            if similarity >= input.similarity_threshold:
                similar_molecules.append({
                    "name": name,
                    "smiles": db_smiles,
                    "similarity": round(similarity, 4),
                })
        
        # Sort by similarity, descending
        similar_molecules.sort(key=lambda x: x['similarity'], reverse=True)
        return similar_molecules

    except Exception as e:
        raise ValueError(f"An error occurred during similarity search: {e}")

class CountFunctionalGroupsInput(BaseModel):
    smiles_string: str = Field(..., description="The SMILES string of the molecule.")

@register_basic_tool(
    metadata=BasicToolMetadata(
        name="count_functional_groups",
        description="Counts common functional groups in a molecule given its SMILES string.",
        parameters=CountFunctionalGroupsInput.model_json_schema()
    )
)
def count_functional_groups(input: CountFunctionalGroupsInput) -> dict:
    """
    Counts common functional groups in a molecule.
    """
    mol = Chem.MolFromSmiles(input.smiles_string)
    if mol is None:
        raise ValueError("Invalid SMILES string.")

    # Define SMARTS patterns for common functional groups
    # This is a simplified set, can be expanded.
    functional_groups = {
        "alcohol": "[OX2H]",
        "amine": "[NX3;H2,H1;!$(NC=O)]", # Primary/secondary amine, not amide
        "aldehyde": "[CX3H1](=O)[#6]",
        "ketone": "[CX3](=O)(-[#6])-[#6]",
        "carboxylic_acid": "[CX3](=O)[OX2H]",
        "ester": "[CX3](=O)[OX2][#6]",
        "amide": "[CX3](=O)[NX3]",
        "ether": "[OX2](-[#6])-[#6]",
        "nitrile": "[NX1]#[CX2]",
        "nitro": "[NX3](=O)=O",
        "sulfonic_acid": "[SX4](=O)(=O)[OX2H]",
        "thiol": "[SX2H]"
    }

    counts = {}
    for name, smarts in functional_groups.items():
        patt = Chem.MolFromSmarts(smarts)
        if patt is not None:
            counts[name] = len(mol.GetSubstructMatches(patt))
    return counts

class CheckAromaticityInput(BaseModel):
    smiles_string: str = Field(..., description="The SMILES string of the molecule.")

@register_basic_tool(
    metadata=BasicToolMetadata(
        name="check_aromaticity",
        description="Checks if a molecule represented by a SMILES string is aromatic.",
        parameters=CheckAromaticityInput.model_json_schema()
    )
)
def check_aromaticity(input: CheckAromaticityInput) -> bool:
    """
    Checks if a molecule is aromatic.
    """
    mol = Chem.MolFromSmiles(input.smiles_string)
    if mol is None:
        raise ValueError("Invalid SMILES string.")
    
    # RDKit automatically perceives aromaticity upon molecule creation.
    # We can check if any ring is aromatic.
    for ring in mol.GetRingInfo().AtomRings():
        is_aromatic_ring = True
        for atom_idx in ring:
            if not mol.GetAtomWithIdx(atom_idx).GetIsAromatic():
                is_aromatic_ring = False
                break
        if is_aromatic_ring:
            return True
    return False

class CountChiralCentersInput(BaseModel):
    smiles_string: str = Field(..., description="The SMILES string of the molecule.")

@register_basic_tool(
    metadata=BasicToolMetadata(
        name="count_chiral_centers",
        description="Counts the number of chiral centers in a molecule given its SMILES string.",
        parameters=CountChiralCentersInput.model_json_schema()
    )
)
def count_chiral_centers(input: CountChiralCentersInput) -> int:
    """
    Counts the number of chiral centers in a molecule.
    """
    mol = Chem.MolFromSmiles(input.smiles_string)
    if mol is None:
        raise ValueError("Invalid SMILES string.")
    
    # Assign stereochemistry if not already present (e.g., from non-isomeric SMILES)
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    
    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    return len(chiral_centers)