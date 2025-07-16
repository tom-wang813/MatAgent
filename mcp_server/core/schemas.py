from pydantic import BaseModel, Field, ConfigDict
from typing import List, Literal, Union, Dict, Any

class ToolMetadata(BaseModel):
    """The base model for all tools, containing common metadata."""
    name: str = Field(..., description="The unique name of the tool.")
    description: str = Field(..., description="A brief description of what the tool does.")
    version: str = Field(default="1.0.0", description="The version of the tool.")
    author: str = Field(default="MatAgent Core", description="The author or origin of the tool.")
    tags: List[str] = Field(default_factory=list, description="Tags for categorizing the tool.")
    parameters: Dict[str, Any] = Field(
        default_factory=lambda: {"type": "object", "properties": {}, "required": []},
        description="JSON Schema for the tool's input parameters."
    )

class BasicToolMetadata(ToolMetadata):
    """Metadata specific to basic, deterministic tools."""
    tool_type: Literal["basic"] = "basic"

class MLModelMetadata(ToolMetadata):
    """Metadata specific to machine learning model-based tools."""
    model_config = ConfigDict(protected_namespaces=())
    tool_type: Literal["ml_model"] = "ml_model"
    model_name: str = Field(..., description="The name of the underlying machine learning model.")
    model_version: str = Field(..., description="The version of the machine learning model.")

# A union type to represent any possible tool metadata, which will be useful for API responses.
AnyToolMetadata = Union[BasicToolMetadata, MLModelMetadata]