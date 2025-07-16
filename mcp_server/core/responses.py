from pydantic import BaseModel, Field
from typing import TypeVar, Generic, Optional, Any

T = TypeVar('T')

class BaseResponse(BaseModel):
    status: str = Field(..., description="Status of the API response (e.g., 'success', 'error').")

class SuccessResponse(BaseResponse, Generic[T]):
    status: str = "success"
    data: T = Field(..., description="The actual data returned by the API.")

class ErrorResponse(BaseResponse):
    status: str = "error"
    message: str = Field(..., description="A human-readable error message.")
    code: Optional[str] = Field(None, description="An optional error code for programmatic handling.")
    details: Optional[Any] = Field(None, description="Optional additional details about the error.")
