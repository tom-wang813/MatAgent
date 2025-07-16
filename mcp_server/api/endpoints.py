from fastapi import APIRouter, HTTPException, Request, status
from fastapi.responses import JSONResponse
from pydantic import BaseModel
from typing import Any, Dict, List
from celery.result import AsyncResult

from mcp_server.celery_app import celery_app
from mcp_server.core.registry import get_all_tools_metadata, run_tool
from mcp_server.core.schemas import AnyToolMetadata
from mcp_server.core.responses import SuccessResponse, ErrorResponse
from shared_utils.logging_config import get_logger

# Get a logger for this module
logger = get_logger(__name__)

# Create a new router for our API endpoints
router = APIRouter()

class ToolRunRequest(BaseModel):
    """Request model for running a tool."""
    tool_name: str
    params: Dict[str, Any] = {}

class ToolRunResponse(BaseModel):
    """Response model for a successful tool run."""
    result: Dict[str, Any]

@router.get("/tools", response_model=SuccessResponse[List[AnyToolMetadata]])
def list_tools(request: Request):
    """
    Returns a list of metadata for all available tools, including basic and ML models.
    This endpoint allows the backend to discover the capabilities of this server.
    """
    trace_id = getattr(request.state, 'trace_id', 'N/A')
    logger.info("Listing all available tools.", extra={'trace_id': trace_id, 'event': 'list_tools_request'})
    tools_metadata = [tool.model_dump() for tool in get_all_tools_metadata()]
    return SuccessResponse(data=tools_metadata)

@router.post("/run_tool", response_model=SuccessResponse[Dict[str, Any]])
def execute_tool(request: ToolRunRequest, http_request: Request):
    """
    Executes a specified tool with the given parameters.
    
    - **tool_name**: The name of the tool to execute (e.g., 'calculate_molecular_weight').
    - **params**: A dictionary of parameters to pass to the tool (e.g., {"smiles_string": "CCO"}).
    """
    trace_id = getattr(http_request.state, 'trace_id', 'N/A')
    
    logger.info(
        f"Tool execution requested: {request.tool_name}",
        extra={
            'trace_id': trace_id,
            'event': 'tool_execution_start',
            'tool_name': request.tool_name,
            'parameters': request.params
        }
    )
    
    try:
        result = run_tool(request.tool_name, trace_id=trace_id, **request.params)
        logger.info(
            f"Tool execution successful: {request.tool_name}",
            extra={
                'trace_id': trace_id,
                'event': 'tool_execution_end',
                'tool_name': request.tool_name,
                'result': result # Be careful with logging large results
            }
        )
        return SuccessResponse(data=result)
    except ValueError as e:
        logger.error(
            f"Tool execution failed (ValueError): {request.tool_name} - {e}",
            extra={
                'trace_id': trace_id,
                'event': 'tool_execution_error',
                'tool_name': request.tool_name,
                'error_type': 'ValueError',
                'error_message': str(e)
            }
        )
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail=str(e))
    except Exception as e:
        logger.critical(
            f"Tool execution failed (Critical Error): {request.tool_name} - {e}",
            extra={
                'trace_id': trace_id,
                'event': 'tool_execution_critical_error',
                'tool_name': request.tool_name,
                'error_type': type(e).__name__,
                'error_message': str(e)
            }
        )
        raise HTTPException(status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, detail=f"An unexpected error occurred: {str(e)}")

@router.get("/tasks/{task_id}/status", response_model=SuccessResponse[Dict[str, Any]])
async def get_task_status(task_id: str, request: Request):
    trace_id = getattr(request.state, 'trace_id', 'N/A')
    logger.info(f"Checking status for task_id: {task_id}", extra={'trace_id': trace_id, 'task_id': task_id})
    task = AsyncResult(task_id, app=celery_app)
    
    response_data = {
        "task_id": task.id,
        "status": task.status,
        "result": None,
        "error": None
    }

    if task.ready():
        if task.successful():
            response_data["result"] = task.result
            logger.info(f"Task {task_id} successful.", extra={'trace_id': trace_id, 'task_id': task_id, 'status': task.status})
        else:
            response_data["error"] = str(task.result)
            logger.error(f"Task {task_id} failed: {task.result}", extra={'trace_id': trace_id, 'task_id': task_id, 'status': task.status, 'error': str(task.result)})
    
    return SuccessResponse(data=response_data)

