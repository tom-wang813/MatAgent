from Bio.SeqUtils.ProtParam import ProteinAnalysis
from pydantic import BaseModel, Field

from mcp_server.core.registry import register_basic_tool
from mcp_server.core.schemas import BasicToolMetadata

class CalculateProteinIsoelectricPointInput(BaseModel):
    amino_acid_sequence: str = Field(..., description="The amino acid sequence of the protein or peptide (e.g., 'AGCTY').")

@register_basic_tool(
    metadata=BasicToolMetadata(
        name="calculate_protein_isoelectric_point",
        description="Calculates the isoelectric point (pI) of a protein or peptide sequence.",
        parameters=CalculateProteinIsoelectricPointInput.model_json_schema()
    )
)
def calculate_protein_isoelectric_point(input: CalculateProteinIsoelectricPointInput) -> float:
    """
    Calculates the isoelectric point (pI) of a protein or peptide sequence.
    """
    try:
        analysed_seq = ProteinAnalysis(input.amino_acid_sequence)
        return analysed_seq.isoelectric_point()
    except Exception as e:
        raise ValueError(f"Error calculating isoelectric point: {e}")

class CalculateProteinExtinctionCoefficientInput(BaseModel):
    amino_acid_sequence: str = Field(..., description="The amino acid sequence of the protein or peptide (e.g., 'AGCTY').")

@register_basic_tool(
    metadata=BasicToolMetadata(
        name="calculate_protein_extinction_coefficient",
        description="Calculates the extinction coefficient of a protein or peptide sequence at 280nm.",
        parameters=CalculateProteinExtinctionCoefficientInput.model_json_schema()
    )
)
def calculate_protein_extinction_coefficient(input: CalculateProteinExtinctionCoefficientInput) -> float:
    """
    Calculates the extinction coefficient of a protein or peptide sequence at 280nm.
    """
    try:
        analysed_seq = ProteinAnalysis(input.amino_acid_sequence)
        # The extinction_coefficient method returns (extinction_coefficient_reduced, extinction_coefficient_oxidized)
        # We usually care about the reduced form for typical protein solutions.
        return analysed_seq.extinction_coefficient()[0]
    except Exception as e:
        raise ValueError(f"Error calculating extinction coefficient: {e}")
