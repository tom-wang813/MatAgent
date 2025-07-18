import os
import random
from typing import Any
from shared_utils.logging_config import get_logger
from mcp_server.core.config import get_settings

settings = get_settings()
logger = get_logger(__name__)

# A dictionary to cache all loaded models
_model_cache: dict[str, Any] = {}

class MockModel:
    """A mock model class to simulate a real ML model."""
    def __init__(self, model_path: str):
        self.model_path = model_path
        # In a real scenario, you would load weights and configs here
        logger.info(f"--- MockModel: Initializing from {model_path} ---")

    def predict(self, smiles: str) -> float:
        # Fake prediction logic for demonstration purposes
        return (len(smiles) + random.random()) * hash(self.model_path) % 100

def _get_model(model_name: str, trace_id: str = "N/A") -> MockModel:
    """Generic lazy-loading model getter."""
    if model_name not in _model_cache:
        model_path = os.path.join(
            settings.MODEL_WEIGHTS_BASE_DIR,
            'property_predictors',
            model_name,
            'model.pth'
        )
        if not os.path.exists(model_path):
            logger.error(f"Model weights not found at {model_path}", extra={'trace_id': trace_id, 'event': 'model_load_error', 'model_name': model_name})
            raise FileNotFoundError(f"Model weights not found at {model_path}")

        logger.info(f"Loading model: {model_name}", extra={'trace_id': trace_id, 'event': 'model_load_start', 'model_name': model_name})
        _model_cache[model_name] = MockModel(model_path)
        logger.info(f"Model loaded: {model_name}", extra={'trace_id': trace_id, 'event': 'model_load_end', 'model_name': model_name})
    return _model_cache[model_name]
