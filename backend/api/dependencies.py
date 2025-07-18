from typing import Any, Optional, TypeVar, Generic
from backend.core.config import get_settings, Settings
from backend.core.database import get_db
from backend.services.tool_manager import ToolManager, ToolManagerError
from backend.services.embedding_service import EmbeddingService
from backend.services.llm_service import LLMService
from backend.services.conversation_service import ConversationService
from backend.agent.orchestrator import AgentOrchestrator
from backend.core.responses import SuccessResponse, ErrorResponse # Import from new location
from sqlalchemy.orm import Session
from fastapi import Depends

settings = get_settings()

T = TypeVar('T')

# Removed BaseResponse, SuccessResponse, ErrorResponse definitions

# Dependency for getting settings
def get_settings_dependency() -> Settings:
    return settings

# Dependencies for services (simplified for example, actual injection might be more complex)
def get_embedding_service() -> EmbeddingService:
    return EmbeddingService()

def get_llm_service() -> LLMService:
    import httpx
    return LLMService(client=httpx.AsyncClient(), frontend_url=settings.FRONTEND_URL)

def get_tool_manager() -> ToolManager:
    # This will be initialized once and reused
    # In a real app, you might use FastAPI's Depends or a singleton pattern
    # For simplicity, we'll create it here. Ensure it's async-aware if needed.
    import httpx
    return ToolManager(
        embedding_service=get_embedding_service(),
        client=httpx.AsyncClient(),
        mcp_server_url=settings.MCP_SERVER_URL,
        settings=settings
    )

def get_conversation_service(db: Session = Depends(get_db)) -> ConversationService:
    return ConversationService(db=db)

def get_agent_orchestrator(
    db: Session = Depends(get_db)
) -> AgentOrchestrator:
    from backend.agent.orchestrator import AgentOrchestrator
    return AgentOrchestrator(
        llm_service=get_llm_service(),
        tool_manager=get_tool_manager(),
        conversation_service=get_conversation_service(db),
        settings=get_settings()
    )