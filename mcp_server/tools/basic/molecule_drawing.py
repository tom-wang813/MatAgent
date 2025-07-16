from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
import base64
from io import BytesIO
from typing import List, Tuple

from pydantic import BaseModel, Field

from mcp_server.core.registry import register_basic_tool, BasicToolMetadata

class SmilesToSvgInput(BaseModel):
    smiles: str = Field(..., description="The SMILES string of the molecule.")
    mol_size: Tuple[int, int] = Field((300, 300), description="The size of the molecule image (width, height).")

@register_basic_tool(
    metadata=BasicToolMetadata(
        name="smiles_to_svg",
        description="Converts a SMILES string to an SVG image of the molecule.",
        parameters=SmilesToSvgInput.model_json_schema()
    )
)
def smiles_to_svg(input: SmilesToSvgInput) -> str:
    """
    Converts a SMILES string to an SVG image of the molecule.
    """
    mol = Chem.MolFromSmiles(input.smiles)
    if mol is None:
        return "<svg width='300' height='300'><text x='10' y='20'>Invalid SMILES</text></svg>"

    drawer = rdMolDraw2DSVG(*input.mol_size)
    drawer.drawOptions().clearBackground = False
    drawer.drawMolecule(mol)
    drawer.finishDrawing()
    svg = drawer.GetDrawingText()
    return svg

class SmilesToPngBase64Input(BaseModel):
    smiles: str = Field(..., description="The SMILES string of the molecule.")
    mol_size: Tuple[int, int] = Field((300, 300), description="The size of the molecule image (width, height).")

@register_basic_tool(
    metadata=BasicToolMetadata(
        name="smiles_to_png_base64",
        description="Converts a SMILES string to a Base64 encoded PNG image of the molecule.",
        parameters=SmilesToPngBase64Input.model_json_schema()
    )
)
def smiles_to_png_base64(input: SmilesToPngBase64Input) -> str:
    """
    Converts a SMILES string to a Base64 encoded PNG image of the molecule.
    """
    mol = Chem.MolFromSmiles(input.smiles)
    if mol is None:
        # Return a placeholder or error image if SMILES is invalid
        return "" 

    drawer = rdMolDraw2DCairo(*input.mol_size)
    drawer.drawOptions().clearBackground = False
    drawer.drawMolecule(mol)
    
    # Get PNG as bytes
    png_bytes = drawer.GetDrawingText()
    
    # Encode to Base64
    base64_png = base64.b64encode(png_bytes).decode('utf-8')
    return base64_png
