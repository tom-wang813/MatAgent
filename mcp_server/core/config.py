from pydantic_settings import BaseSettings, SettingsConfigDict
from functools import lru_cache

class Settings(BaseSettings):
    MODEL_WEIGHTS_BASE_DIR: str = "./model_weights"
    CELERY_BROKER_URL: str = "redis://localhost:6379/0"
    CELERY_RESULT_BACKEND: str = "redis://localhost:6379/0"
    CELERY_ALWAYS_EAGER: bool = False # Set to True for local testing without a worker

    model_config = SettingsConfigDict(env_file=".env", extra="ignore")

@lru_cache()
def get_settings():
    """
    Returns a cached instance of the Settings class.
    This ensures settings are loaded only once.
    """
    return Settings()