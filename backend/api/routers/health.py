from fastapi import APIRouter
from backend.api.dependencies import SuccessResponse

router = APIRouter()

@router.get("/health", response_model=SuccessResponse[dict])
async def health_check():
    return SuccessResponse(data={"status": "ok"})
