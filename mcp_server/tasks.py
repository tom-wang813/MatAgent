from mcp_server.celery_app import celery_app
from mcp_server.tools.ml_models.property_predictors import _get_model
from mcp_server.tools.ml_models.property_predictors import PredictAqueousSolubilityInput, PredictRadiusOfGyrationInput, PredictHeatCapacityInput
from shared_utils.logging_config import get_logger

logger = get_logger(__name__)

@celery_app.task(bind=True)
def predict_aqueous_solubility_task(self, smiles_string: str, trace_id: str = "N/A") -> float:
    logger.info(f"Starting predict_aqueous_solubility_task for {smiles_string}", extra={'trace_id': trace_id, 'task_id': self.request.id})
    try:
        model = _get_model("solubility_model_v1", trace_id)
        result = model.predict(smiles_string)
        logger.info(f"Finished predict_aqueous_solubility_task for {smiles_string}", extra={'trace_id': trace_id, 'task_id': self.request.id, 'result': result})
        return result
    except Exception as e:
        logger.error(f"Error in predict_aqueous_solubility_task for {smiles_string}: {e}", extra={'trace_id': trace_id, 'task_id': self.request.id})
        raise

@celery_app.task(bind=True)
def predict_radius_of_gyration_task(self, smiles_string: str, trace_id: str = "N/A") -> float:
    logger.info(f"Starting predict_radius_of_gyration_task for {smiles_string}", extra={'trace_id': trace_id, 'task_id': self.request.id})
    try:
        model = _get_model("rg_model_v1", trace_id)
        result = model.predict(smiles_string)
        logger.info(f"Finished predict_radius_of_gyration_task for {smiles_string}", extra={'trace_id': trace_id, 'task_id': self.request.id, 'result': result})
        return result
    except Exception as e:
        logger.error(f"Error in predict_radius_of_gyration_task for {smiles_string}: {e}", extra={'trace_id': trace_id, 'task_id': self.request.id})
        raise

@celery_app.task(bind=True)
def predict_heat_capacity_task(self, smiles_string: str, trace_id: str = "N/A") -> float:
    logger.info(f"Starting predict_heat_capacity_task for {smiles_string}", extra={'trace_id': trace_id, 'task_id': self.request.id})
    try:
        model = _get_model("cp_model_v1", trace_id)
        result = model.predict(smiles_string)
        logger.info(f"Finished predict_heat_capacity_task for {smiles_string}", extra={'trace_id': trace_id, 'task_id': self.request.id, 'result': result})
        return result
    except Exception as e:
        logger.error(f"Error in predict_heat_capacity_task for {smiles_string}: {e}", extra={'trace_id': trace_id, 'task_id': self.request.id})
        raise
