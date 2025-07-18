import os
from typing import List, Dict, Any
import numpy as np
from sentence_transformers import SentenceTransformer

from backend.core.config import get_settings
from shared_utils.logging_config import get_logger

settings = get_settings()
logger = get_logger(__name__)

class EmbeddingServiceError(Exception):
    """Custom exception for Embedding service errors."""
    pass

class EmbeddingService:
    def __init__(self):
        self.model = None # Initialize model as None
        self.model_name = settings.OPENROUTER_EMBEDDING_MODEL
        self._load_model() # Load model during initialization

    def _load_model(self):
        """
        Loads the SentenceTransformer model from Hugging Face Hub.
        This method is called during initialization.
        """
        logger.info(f"Loading SentenceTransformer model: {self.model_name}")
        try:
            # Load the model from Hugging Face Hub
            self.model = SentenceTransformer(self.model_name)
            logger.info(f"SentenceTransformer model {self.model_name} loaded successfully.")
        except Exception as e:
            logger.critical(f"Failed to load SentenceTransformer model {self.model_name}: {e}")
            raise EmbeddingServiceError(f"Failed to load embedding model: {e}") from e

    def get_embedding(self, text: str, trace_id: str = "N/A") -> List[float]:
        """
        Generates an embedding for the given text using the loaded SentenceTransformer model.
        """
        if not self.model:
            raise EmbeddingServiceError("Embedding model not loaded.")

        # Ensure text is a string
        if not isinstance(text, str):
            text = str(text) if text is not None else ""

        logger.info(
            "Generating embedding using SentenceTransformer",
            extra={
                'trace_id': trace_id,
                'event': 'embedding_generation_start',
                'model': self.model_name,
                'input_length': len(text)
            }
        )

        try:
            # Encode the text to get the embedding
            embedding = self.model.encode(text).tolist()
            
            logger.info(
                "Embedding generated successfully",
                extra={
                    'trace_id': trace_id,
                    'event': 'embedding_generation_end',
                    'model': self.model_name,
                    'embedding_dim': len(embedding)
                }
            )
            return embedding
        except Exception as e:
            logger.error(
                f"Error generating embedding: {e}",
                extra={
                    'trace_id': trace_id,
                    'event': 'embedding_generation_error',
                    'error_message': str(e)
                }
            )
            raise EmbeddingServiceError(f"Error generating embedding: {e}") from e