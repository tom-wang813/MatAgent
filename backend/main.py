from fastapi import FastAPI, Request, Response, status
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse
from pydantic import ValidationError
from contextlib import asynccontextmanager
import uuid
import httpx
from sentence_transformers import SentenceTransformer

from shared_utils.logging_config import get_logger
from backend.core.database import Base, engine
from backend.api.routers import chat, health, files
from backend.services.tool_manager import ToolManager, ToolManagerError
from backend.services.embedding_service import EmbeddingService
from backend.core.config import get_settings
from backend.api.dependencies import ErrorResponse

logger = get_logger(__name__)
settings = get_settings()

@asynccontextmanager
async def lifespan(app: FastAPI):
    # Startup events
    logger.info("Backend application startup.")
    
    # Create database tables
    Base.metadata.create_all(bind=engine)
    logger.info("Database tables created/checked.")

    # Initialize httpx.AsyncClient for external API calls
    app.state.httpx_client = httpx.AsyncClient(base_url=settings.OPENROUTER_BASE_URL, timeout=settings.LLM_TIMEOUT)
    logger.info("httpx.AsyncClient initialized.")

    # Initialize ToolManager and discover tools once globally
    try:
        # EmbeddingService manages its own model loading internally
        app.state.embedding_service = EmbeddingService()
        app.state.tool_manager = ToolManager(
            embedding_service=app.state.embedding_service,
            client=app.state.httpx_client,
            mcp_server_url=settings.MCP_SERVER_URL,
            settings=settings
        )
        await app.state.tool_manager.discover_tools()
        logger.info("Tool discovery initiated and completed during startup.")
    except Exception as e:
        logger.critical(f"Failed to discover tools on startup: {e}", extra={'event': 'tool_discovery_startup_failure'})
        # Depending on criticality, you might want to exit here or just log and continue

    # Store frontend URL for LLMService HTTP-Referer header
    app.state.frontend_url = "http://localhost:3000" # This should be configurable via settings in a real app

    yield

    # Shutdown events
    logger.info("Backend application shutdown.")
    # Close httpx.AsyncClient
    if hasattr(app.state, 'httpx_client') and app.state.httpx_client:
        await app.state.httpx_client.aclose()
        logger.info("httpx.AsyncClient closed.")

app = FastAPI(
    title="MatAgent Backend",
    description="Core backend service for MatAgent, handling LLM orchestration and tool integration.",
    version="1.0.0",
    lifespan=lifespan
)

# CORS Middleware
app.add_middleware(
    CORSMiddleware,
    allow_origins=["http://localhost:3000"], # Allow your frontend origin
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Middleware to add trace_id to request state and log incoming requests
@app.middleware("http")
async def add_trace_id_middleware(request: Request, call_next):
    trace_id = request.headers.get("X-Trace-ID", str(uuid.uuid4()))
    request.state.trace_id = trace_id
    
    logger.info(
        f"Incoming request: {request.method} {request.url.path}",
        extra={'trace_id': trace_id, 'event': 'http_request_received', 'method': request.method, 'path': request.url.path}
    )
    
    response = await call_next(request)
    
    logger.info(
        f"Outgoing response: {response.status_code}",
        extra={'trace_id': trace_id, 'event': 'http_response_sent', 'status_code': response.status_code}
    )
    return response

@app.exception_handler(ValidationError)
async def validation_exception_handler(request: Request, exc: ValidationError):
    logger.error(f"Validation error: {exc.errors()}", extra={'trace_id': getattr(request.state, 'trace_id', 'N/A')})
    return JSONResponse(
        status_code=status.HTTP_422_UNPROCESSABLE_ENTITY,
        content=ErrorResponse(
            message="Validation error",
            code="VALIDATION_ERROR",
            details=exc.errors()
        ).model_dump()
    )

@app.exception_handler(Exception)
async def generic_exception_handler(request: Request, exc: Exception):
    logger.critical(f"Unhandled exception: {exc}", extra={'trace_id': getattr(request.state, 'trace_id', 'N/A')})
    return JSONResponse(
        status_code=status.HTTP_500_INTERNAL_SERVER_ERROR,
        content=ErrorResponse(
            message="An unexpected error occurred",
            code="INTERNAL_SERVER_ERROR",
            details=str(exc)
        ).model_dump()
    )

# Include routers
app.include_router(chat.router, prefix="/api", tags=["Chat"])
app.include_router(health.router, prefix="/api", tags=["Health"])
app.include_router(files.router, prefix="/api", tags=["Files"])
