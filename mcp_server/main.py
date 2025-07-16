from fastapi import FastAPI, Request, status
from fastapi.responses import JSONResponse
from pydantic import ValidationError

# Import the API router from the api directory
from mcp_server.api.endpoints import router as api_router
from mcp_server import tools
from mcp_server.core.responses import ErrorResponse
from shared_utils.logging_config import get_logger

# Get a logger for this module
logger = get_logger(__name__)

app = FastAPI(
    title="MatAgent MCP Server",
    description="A dedicated server for Material Chemistry Properties (MCP) calculations and ML model inference.",
    version="1.0.0",
)

# Include the API endpoints
app.include_router(api_router, prefix="/api", tags=["Tools"])

@app.middleware("http")
async def add_trace_id_middleware(request: Request, call_next):
    trace_id = request.headers.get("X-Trace-ID")
    if trace_id:
        request.state.trace_id = trace_id
        logger.info("Received request with trace_id", extra={'trace_id': trace_id, 'path': request.url.path})
    response = await call_next(request)
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

@app.get("/", tags=["Status"])
def read_root():
    """A simple health check endpoint to confirm the server is running."""
    logger.info("Health check requested.")
    return {"status": "MCP Server is running"}
