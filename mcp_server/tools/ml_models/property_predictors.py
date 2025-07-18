from typing import Optional
from pydantic import BaseModel, Field

from mcp_server.core.registry import register_ml_model
from mcp_server.core.schemas import MLModelMetadata
from shared_utils.logging_config import get_logger
from mcp_server.tools.ml_models.model_loader import _get_model, MockModel

logger = get_logger(__name__)

class PredictAqueousSolubilityInput(BaseModel):
    smiles_string: str = Field(..., description="The SMILES representation of the molecule.")

class PredictRadiusOfGyrationInput(BaseModel):
    smiles_string: str = Field(..., description="The SMILES representation of the monomer or polymer unit.")

class PredictHeatCapacityInput(BaseModel):
    smiles_string: str = Field(..., description="The SMILES representation of the molecule.")

# --- Tool: Predict Aqueous Solubility ---
@register_ml_model(
    metadata=MLModelMetadata(
        name="predict_aqueous_solubility",
        description="Predicts the aqueous solubility (LogS) of a molecule using a pre-trained GNN model.",
        tags=["machine-learning", "prediction", "solubility"],
        author="In-house ML Team",
        model_name="solubility_model_v1",
        model_version="1.0.0",
        parameters=PredictAqueousSolubilityInput.model_json_schema()
    )
)
def predict_solubility(input: PredictAqueousSolubilityInput, trace_id: str = "N/A") -> str:
    from mcp_server.tasks import predict_aqueous_solubility_task # Delayed import
    task = predict_aqueous_solubility_task.delay(input.smiles_string, trace_id)
    return task.id

# --- Tool: Predict Radius of Gyration (Rg) ---
@register_ml_model(
    metadata=MLModelMetadata(
        name="predict_radius_of_gyration",
        description="Predicts the Radius of Gyration (Rg) for a polymer chain from its monomer SMILES.",
        tags=["machine-learning", "prediction", "polymer", "conformation"],
        author="In-house ML Team",
        model_name="rg_model_v1",
        model_version="1.0.0",
        parameters=PredictRadiusOfGyrationInput.model_json_schema()
    )
)
def predict_rg(input: PredictRadiusOfGyrationInput, trace_id: str = "N/A") -> str:
    from mcp_server.tasks import predict_radius_of_gyration_task # Delayed import
    task = predict_radius_of_gyration_task.delay(input.smiles_string, trace_id)
    return task.id

# --- Tool: Predict Heat Capacity (Cp) ---
@register_ml_model(
    metadata=MLModelMetadata(
        name="predict_heat_capacity",
        description="Predicts the heat capacity (Cp) of a molecule at a standard temperature.",
        tags=["machine-learning", "prediction", "thermodynamics"],
        author="In-house ML Team",
        model_name="cp_model_v1",
        model_version="1.0.0",
        parameters=PredictHeatCapacityInput.model_json_schema()
    )
)
def predict_cp(input: PredictHeatCapacityInput, trace_id: str = "N/A") -> str:
    from mcp_server.tasks import predict_heat_capacity_task # Delayed import
    task = predict_heat_capacity_task.delay(input.smiles_string, trace_id)
    return task.id

class PredictAqueousSolubilityInput(BaseModel):
    smiles_string: str = Field(..., description="The SMILES representation of the molecule.")

# # --- Tool: Predict Aqueous Solubility ---
# @register_ml_model(
#     metadata=MLModelMetadata(
#         name="predict_aqueous_solubility",
#         description="Predicts the aqueous solubility (LogS) of a molecule using a pre-trained GNN model.",
#         tags=["machine-learning", "prediction", "solubility"],
#         author="In-house ML Team",
#         model_name="solubility_model_v1",
#         model_version="1.0.0",
#         parameters=PredictAqueousSolubilityInput.model_json_schema()
#     )
# )
# def predict_solubility(input: PredictAqueousSolubilityInput, trace_id: str = "N/A") -> str:
#     task = predict_aqueous_solubility_task.delay(input.smiles_string, trace_id)
#     return task.id

# class PredictRadiusOfGyrationInput(BaseModel):
#     smiles_string: str = Field(..., description="The SMILES representation of the monomer or polymer unit.")

# # --- Tool: Predict Radius of Gyration (Rg) ---
# @register_ml_model(
#     metadata=MLModelMetadata(
#         name="predict_radius_of_gyration",
#         description="Predicts the Radius of Gyration (Rg) for a polymer chain from its monomer SMILES.",
#         tags=["machine-learning", "prediction", "polymer", "conformation"],
#         author="In-house ML Team",
#         model_name="rg_model_v1",
#         model_version="1.0.0",
#         parameters=PredictRadiusOfGyrationInput.model_json_schema()
#     )
# )
# def predict_rg(input: PredictRadiusOfGyrationInput, trace_id: str = "N/A") -> str:
#     task = predict_radius_of_gyration_task.delay(input.smiles_string, trace_id)
#     return task.id

# class PredictHeatCapacityInput(BaseModel):
#     smiles_string: str = Field(..., description="The SMILES representation of the molecule.")

# # --- Tool: Predict Heat Capacity (Cp) ---
# @register_ml_model(
#     metadata=MLModelMetadata(
#         name="predict_heat_capacity",
#         description="Predicts the heat capacity (Cp) of a molecule at a standard temperature.",
#         tags=["machine-learning", "prediction", "thermodynamics"],
#         author="In-house ML Team",
#         model_name="cp_model_v1",
#         model_version="1.0.0",
#         parameters=PredictHeatCapacityInput.model_json_schema()
#     )
# )
# def predict_cp(input: PredictHeatCapacityInput, trace_id: str = "N/A") -> str:
#     task = predict_heat_capacity_task.delay(input.smiles_string, trace_id)
#     return task.id