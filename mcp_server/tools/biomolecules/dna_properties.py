from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
from pydantic import BaseModel, Field

from mcp_server.core.registry import register_basic_tool
from mcp_server.core.schemas import BasicToolMetadata

class PredictDnaMeltingTemperatureInput(BaseModel):
    dna_sequence: str = Field(..., description="The DNA sequence (e.g., 'ATGCGT').")

@register_basic_tool(
    metadata=BasicToolMetadata(
        name="predict_dna_melting_temperature",
        description="Predicts the melting temperature (Tm) of a short DNA duplex using the nearest-neighbor method.",
        parameters=PredictDnaMeltingTemperatureInput.model_json_schema()
    )
)
def predict_dna_melting_temperature(input: PredictDnaMeltingTemperatureInput) -> float:
    """
    Predicts the melting temperature (Tm) of a short DNA duplex.
    """
    try:
        # Default parameters for nearest-neighbor method are usually suitable for basic prediction.
        # Can add more parameters (e.g., Na+, K+, Mg2+, DNA_conc, dNTP_conc) if needed for more advanced use cases.
        return mt.Tm_NN(input.dna_sequence)
    except Exception as e:
        raise ValueError(f"Error predicting DNA melting temperature: {e}")

class ReverseComplementDnaInput(BaseModel):
    dna_sequence: str = Field(..., description="The DNA sequence (e.g., 'ATGCGT').")

@register_basic_tool(
    metadata=BasicToolMetadata(
        name="reverse_complement_dna",
        description="Generates the reverse complement sequence of a given DNA sequence.",
        parameters=ReverseComplementDnaInput.model_json_schema()
    )
)
def reverse_complement_dna(input: ReverseComplementDnaInput) -> str:
    """
    Generates the reverse complement sequence of a DNA sequence.
    """
    try:
        seq_obj = Seq(input.dna_sequence)
        return str(seq_obj.reverse_complement())
    except Exception as e:
        raise ValueError(f"Error generating reverse complement: {e}")

class TranslateDnaToProteinInput(BaseModel):
    dna_sequence: str = Field(..., description="The DNA sequence (e.g., 'ATGCGT').")

@register_basic_tool(
    metadata=BasicToolMetadata(
        name="translate_dna_to_protein",
        description="Translates a DNA sequence into its corresponding protein sequence.",
        parameters=TranslateDnaToProteinInput.model_json_schema()
    )
)
def translate_dna_to_protein(input: TranslateDnaToProteinInput) -> str:
    """
    Translates a DNA sequence into its corresponding protein sequence.
    """
    try:
        seq_obj = Seq(input.dna_sequence)
        return str(seq_obj.translate())
    except Exception as e:
        raise ValueError(f"Error translating DNA to protein: {e}")
