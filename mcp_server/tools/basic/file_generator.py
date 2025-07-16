import os
from rdkit import Chem
from rdkit.Chem import AllChem
from pydantic import BaseModel, Field
from typing import List

from mcp_server.core.registry import register_basic_tool, BasicToolMetadata

class GenerateTextInput(BaseModel):
    file_name: str = Field(..., description="The name of the file to generate.")
    content: str = Field(..., description="The content to write to the file.")

@register_basic_tool(
    metadata=BasicToolMetadata(
        name="generate_text_file",
        description="Generates a text file with the given content and returns its path.",
        parameters=GenerateTextInput.model_json_schema()
    )
)
def generate_text_file(input: GenerateTextInput) -> dict:
    """
    Generates a text file in the temporary files directory.
    """
    temp_dir = "/Users/wang-work/matagent/mcp_server/temp_files"
    file_path = os.path.join(temp_dir, input.file_name)
    
    with open(file_path, "w") as f:
        f.write(input.content)
        
    return {"file_path": file_path, "message": f"File '{input.file_name}' generated successfully."}

class GenerateSDFInput(BaseModel):
    smiles_list: List[str] = Field(..., description="A list of SMILES representations of the molecules.")

@register_basic_tool(
    metadata=BasicToolMetadata(
        name="generate_sdf_file",
        description="Generates an SDF file containing multiple molecules from a list of SMILES strings.",
        parameters=GenerateSDFInput.model_json_schema()
    )
)
def generate_sdf_file(input: GenerateSDFInput) -> dict:
    """
    Generates an SDF file containing multiple molecules from a list of SMILES strings.
    """
    temp_dir = "/Users/wang-work/matagent/mcp_server/temp_files"
    file_name = f"molecules_{os.urandom(4).hex()}.sdf"
    file_path = os.path.join(temp_dir, file_name)

    writer = Chem.SDWriter(file_path)
    for smiles in input.smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            AllChem.Compute2DCoords(mol) # Generate 2D coordinates for visualization
            writer.write(mol)
    writer.close()

    return {"file_path": file_path, "message": f"SDF file '{file_name}' generated successfully."}
