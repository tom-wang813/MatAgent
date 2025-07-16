from celery import Celery
from mcp_server.core.config import get_settings

settings = get_settings()

celery_app = Celery(
    'mcp_server',
    broker=settings.CELERY_BROKER_URL,
    backend=settings.CELERY_RESULT_BACKEND,
    include=['mcp_server.tasks'] # Tasks will be defined in mcp_server/tasks.py
)

celery_app.conf.update(
    task_track_started=True,
    task_acks_late=True,
    worker_prefetch_multiplier=1,
    task_serializer='json',
    result_serializer='json',
    accept_content=['json'],
    task_always_eager=settings.CELERY_ALWAYS_EAGER, # For local testing without a running worker
)
