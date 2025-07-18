from typing import Callable, Dict, Any, List
from mcp_server.core.schemas import BasicToolMetadata, MLModelMetadata, AnyToolMetadata
from shared_utils.logging_config import get_logger

logger = get_logger(__name__)

# --- In-Memory Registries ---
# As you suggested, we use separate dictionaries ("sub dictionary") for each tool type.
_basic_tool_registry: Dict[str, Callable] = {}
_basic_tool_metadata: Dict[str, BasicToolMetadata] = {}

_ml_model_registry: Dict[str, Callable] = {}
_ml_model_metadata: Dict[str, MLModelMetadata] = {}


# --- Registration Logic ---

def register_basic_tool(metadata: BasicToolMetadata):
    """A decorator to register a basic, deterministic tool."""
    def decorator(func: Callable):
        logger.info(f"Registering basic tool: {metadata.name}")
        if metadata.name in _basic_tool_registry or metadata.name in _ml_model_registry:
            logger.warning(f"Tool with name '{metadata.name}' is already registered. Skipping re-registration.")
            return func # Return the original function to avoid breaking the decorator chain
        _basic_tool_registry[metadata.name] = func
        _basic_tool_metadata[metadata.name] = metadata
        return func
    return decorator

def register_ml_model(metadata: MLModelMetadata):
    """A decorator to register a machine learning model as a tool."""
    def decorator(func: Callable):
        logger.info(f"Registering ML model tool: {metadata.name}")
        if metadata.name in _basic_tool_registry or metadata.name in _ml_model_registry:
            logger.warning(f"Tool with name '{metadata.name}' is already registered. Skipping re-registration.")
            return func # Return the original function to avoid breaking the decorator chain
        _ml_model_registry[metadata.name] = func
        _ml_model_metadata[metadata.name] = metadata
        return func
    return decorator


# --- Tool Access and Execution Logic ---

def get_all_tools_metadata() -> List[AnyToolMetadata]:
    """Returns a list of metadata for all registered tools (both basic and ML)."""
    all_metadata = list(_basic_tool_metadata.values()) + list(_ml_model_metadata.values())
    return all_metadata

def find_tool(tool_name: str) -> tuple[Callable, AnyToolMetadata]:
    """Finds a tool's function and metadata by its name across all registries."""
    if tool_name in _basic_tool_registry:
        return _basic_tool_registry[tool_name], _basic_tool_metadata[tool_name]
    if tool_name in _ml_model_registry:
        return _ml_model_registry[tool_name], _ml_model_metadata[tool_name]
    raise ValueError(f"Tool '{tool_name}' not found.")

def run_tool(tool_name: str, trace_id: str = "N/A", **kwargs: Any) -> Any:
    """
    Finds and runs a tool by its name with the given arguments.
    If the tool is an ML model, it returns a task ID.
    """
    tool_function, tool_metadata = find_tool(tool_name)
    
    # Get the function signature to understand parameter structure
    import inspect
    sig = inspect.signature(tool_function)
    parameters = list(sig.parameters.keys())
    
    # Check if the function expects a single input model parameter
    if len(parameters) == 1 or (len(parameters) == 2 and 'trace_id' in parameters):
        # Find the input parameter name (not trace_id)
        input_param = next((p for p in parameters if p != 'trace_id'), None)
        if input_param:
            # Get the parameter annotation to create the input object
            param_annotation = sig.parameters[input_param].annotation
            if hasattr(param_annotation, 'model_validate'):
                # Create the input object from kwargs
                input_obj = param_annotation.model_validate(kwargs)
                if 'trace_id' in parameters:
                    return tool_function(input_obj, trace_id=trace_id)
                else:
                    return tool_function(input_obj)
    
    # Fallback to original behavior for functions with multiple parameters
    if tool_metadata.tool_type == "ml_model":
        # The ML model functions are designed to return task IDs directly
        # They already call .delay() internally
        if 'trace_id' in parameters:
            return tool_function(trace_id=trace_id, **kwargs)
        else:
            return tool_function(**kwargs)
    else:
        # For basic tools, execute directly
        if 'trace_id' in parameters:
            return tool_function(trace_id=trace_id, **kwargs)
        else:
            return tool_function(**kwargs)
