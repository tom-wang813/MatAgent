from rdkit import Chem
from rdkit.Chem.Draw import MolToImage, MolToSVG
from rdkit.Chem import AllChem # Import AllChem for 3D conformer generation
import base64
from io import BytesIO
from pydantic import BaseModel, Field
from typing import Literal, List

from mcp_server.core.registry import register_basic_tool
from mcp_server.core.schemas import BasicToolMetadata

def _get_mol(smiles_string: str):
    """Helper function to convert SMILES to a Mol object, with error handling."""
    mol = Chem.MolFromSmiles(smiles_string)
    if mol is None:
        raise ValueError(f"Invalid SMILES string provided: {smiles_string}")
    return mol

class GenerateMoleculeImageInput(BaseModel):
    smiles_string: str = Field(..., description="The SMILES representation of the molecule.")
    image_format: Literal["png", "svg"] = Field("png", description="The desired image format (png or svg).")

# --- Tool: Generate Molecule Image ---
@register_basic_tool(
    metadata=BasicToolMetadata(
        name="generate_molecule_image",
        description="Generates a 2D image of a molecule from its SMILES string. Returns the image as a Base64 encoded string.",
        tags=["rdkit", "chemistry", "visualization"],
        author="RDKit Wrapper",
        parameters=GenerateMoleculeImageInput.model_json_schema()
    )
)
def generate_molecule_image(input: GenerateMoleculeImageInput) -> str:
    mol = _get_mol(input.smiles_string)
    if input.image_format == "png":
        img = MolToImage(mol)
        buffered = BytesIO()
        img.save(buffered, format="PNG")
        img_str = base64.b64encode(buffered.getvalue()).decode("utf-8")
        return f"data:image/png;base64,{img_str}"
    elif input.image_format == "svg":
        img_str = MolToSVG(mol)
        img_str_b64 = base64.b64encode(img_str.encode("utf-8")).decode("utf-8")
        return f"data:image/svg+xml;base64,{img_str_b64}"
    else:
        raise ValueError("Unsupported image format. Choose 'png' or 'svg'.")

class ConvertSmilesToFormatInput(BaseModel):
    smiles_string: str = Field(..., description="The SMILES representation of the molecule.")
    target_format: Literal["inchi", "inchikey", "mol", "sdf"] = Field(..., description="The target format (inchi, inchikey, mol, or sdf).")

# --- Tool: Convert SMILES to Format ---
@register_basic_tool(
    metadata=BasicToolMetadata(
        name="convert_smiles_to_format",
        description="Converts a SMILES string to various chemical file formats (InChI, InChIKey, Mol, SDF).",
        tags=["rdkit", "chemistry", "format-conversion"],
        author="RDKit Wrapper",
        parameters=ConvertSmilesToFormatInput.model_json_schema()
    )
)
def convert_smiles_to_format(input: ConvertSmilesToFormatInput) -> str:
    mol = _get_mol(input.smiles_string)
    if input.target_format == "inchi":
        return Chem.MolToInchi(mol)
    elif input.target_format == "inchikey":
        return Chem.MolToInchiKey(mol)
    elif input.target_format == "mol":
        return Chem.MolToMolBlock(mol)
    elif input.target_format == "sdf":
        # For a single molecule, MolToMolBlock is often used for SDF-like output
        # For multi-molecule SDF, a different tool would be needed.
        return Chem.MolToMolBlock(mol) # RDKit's MolToMolBlock generates a V2000 mol block, which is often used within SDF files.
    else:
        raise ValueError(f"Unsupported target format: {input.target_format}")

class GenerateSDFFileFromListInput(BaseModel):
    smiles_list: List[str] = Field(..., description="A list of SMILES representations of the molecules.")

# --- Tool: Generate SDF File from SMILES List ---
@register_basic_tool(
    metadata=BasicToolMetadata(
        name="generate_sdf_file",
        description="Generates an SDF file containing multiple molecules from a list of SMILES strings.",
        tags=["rdkit", "chemistry", "file-generation"],
        author="RDKit Wrapper",
        parameters=GenerateSDFFileFromListInput.model_json_schema()
    )
)
def generate_sdf_file(input: GenerateSDFFileFromListInput) -> str:
    sdf_content = ""
    for smiles_string in input.smiles_list:
        mol = Chem.MolFromSmiles(smiles_string)
        if mol is None:
            # Optionally, handle invalid SMILES more gracefully, e.g., skip or log
            continue
        sdf_content += Chem.MolToMolBlock(mol) + "\n$\n"
    return sdf_content

class FilterMoleculesBySubstructureInput(BaseModel):
    smiles_list: List[str] = Field(..., description="A list of SMILES representations of the molecules to filter.")
    smarts_pattern: str = Field(..., description="The SMARTS pattern representing the substructure to search for.")

# --- Tool: Filter Molecules by Substructure ---
@register_basic_tool(
    metadata=BasicToolMetadata(
        name="filter_molecules_by_substructure",
        description="Filters a list of SMILES strings, returning only those that contain a specified substructure (SMARTS pattern).",
        tags=["rdkit", "chemistry", "substructure", "filtering"],
        author="RDKit Wrapper",
        parameters=FilterMoleculesBySubstructureInput.model_json_schema()
    )
)
def filter_molecules_by_substructure(input: FilterMoleculesBySubstructureInput) -> List[str]:
    try:
        substructure_mol = Chem.MolFromSmarts(input.smarts_pattern)
        if substructure_mol is None:
            raise ValueError(f"Invalid SMARTS pattern provided: {input.smarts_pattern}")
    except Exception as e:
        raise ValueError(f"Error parsing SMARTS pattern '{input.smarts_pattern}': {e}")

    matching_smiles = []
    for smiles_string in input.smiles_list:
        try:
            mol = _get_mol(smiles_string)
            if mol.HasSubstructMatch(substructure_mol):
                matching_smiles.append(smiles_string)
        except ValueError:
            # Skip invalid SMILES strings in the list
            continue
    return matching_smiles

class Generate3DConformersInput(BaseModel):
    smiles_string: str = Field(..., description="The SMILES representation of the molecule.")
    num_conformers: int = Field(10, description="The number of conformers to generate.")
    max_attempts: int = Field(100, description="Maximum attempts to generate each conformer.")

# --- Tool: Generate 3D Conformers ---
@register_basic_tool(
    metadata=BasicToolMetadata(
        name="generate_3d_conformers",
        description="Generates a set of 3D conformers for a given molecule and returns them as an SDF file string.",
        tags=["rdkit", "chemistry", "3d", "conformer", "file-generation"],
        author="RDKit Wrapper",
        parameters=Generate3DConformersInput.model_json_schema()
    )
)
def generate_3d_conformers(input: Generate3DConformersInput) -> str:
    mol = _get_mol(input.smiles_string)
    mol = Chem.AddHs(mol) # Add hydrogens for 3D generation

    AllChem.EmbedMultipleConfs(mol, numConfs=input.num_conformers, maxAttempts=input.max_attempts)
    
    # Optimize conformers and write to SDF
    sdf_content = ""
    for i, conf in enumerate(mol.GetConformers()):
        AllChem.MMFFOptimizeMolecule(mol, confId=i)
        sdf_content += Chem.MolToMolBlock(mol, confId=i) + "\n$\n"
    return sdf_content

class ConvertSmilesToInchiInput(BaseModel):
    smiles_string: str = Field(..., description="The SMILES representation of the molecule.")

# --- Tool: Convert SMILES to InChI ---
@register_basic_tool(
    metadata=BasicToolMetadata(
        name="convert_smiles_to_inchi",
        description="Converts a standard SMILES string into its corresponding International Chemical Identifier (InChI) string.",
        tags=["rdkit", "chemistry", "format-conversion"],
        author="RDKit Wrapper",
        parameters=ConvertSmilesToInchiInput.model_json_schema()
    )
)
def convert_smiles_to_inchi(input: ConvertSmilesToInchiInput) -> str:
    mol = _get_mol(input.smiles_string)
    return Chem.MolToInchi(mol)

class ValidateSmilesInput(BaseModel):
    smiles_string: str = Field(..., description="The SMILES representation of the molecule.")

# --- Tool: Validate SMILES String ---
@register_basic_tool(
    metadata=BasicToolMetadata(
        name="validate_smiles",
        description="Checks if a SMILES string represents a valid chemical structure that can be parsed by RDKit.",
        tags=["rdkit", "chemistry", "validation"],
        author="RDKit Wrapper",
        parameters=ValidateSmilesInput.model_json_schema()
    )
)
def validate_smiles(input: ValidateSmilesInput) -> bool:
    mol = Chem.MolFromSmiles(input.smiles_string)
    return mol is not None