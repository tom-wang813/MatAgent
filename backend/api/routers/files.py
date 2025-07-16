from fastapi import APIRouter, HTTPException, status
from fastapi.responses import FileResponse
import os

router = APIRouter()

@router.get("/download/mcp_file")
async def download_mcp_file(file_path: str):
    # Ensure the file path is within the allowed temporary directory
    # This is crucial for security to prevent directory traversal attacks
    base_dir = "/Users/wang-work/matagent/mcp_server/temp_files"
    abs_file_path = os.path.abspath(os.path.join(base_dir, file_path))

    # Check if the requested file path is indeed within the base_dir
    if not abs_file_path.startswith(base_dir):
        raise HTTPException(status_code=status.HTTP_400_BAD_REQUEST, detail="Invalid file path.")

    if not os.path.exists(abs_file_path):
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="File not found.")

    file_name = os.path.basename(abs_file_path)
    return FileResponse(path=abs_file_path, filename=file_name, media_type="application/octet-stream")
