from fastapi import APIRouter, Depends, HTTPException, status, Request
from fastapi.responses import StreamingResponse
from typing_extensions import Annotated
from fastapi import Query
import json
from typing import Dict, List, Optional, AsyncGenerator

from shared_utils.logging_config import get_logger
from backend.agent.orchestrator import AgentOrchestrator
from backend.services.conversation_service import ConversationService
from backend.api.dependencies import get_agent_orchestrator, get_conversation_service, SuccessResponse, ErrorResponse

logger = get_logger(__name__)

router = APIRouter()



async def stream_agent_response(agent_generator: AsyncGenerator[Dict, None], conversation_uuid: Optional[str] = None) -> AsyncGenerator[str, None]:
    """Consumes the agent's async generator and yields JSON strings for streaming."""
    if conversation_uuid:
        yield f"data: {json.dumps({'type': 'conversation_info', 'conversation_uuid': conversation_uuid})}\n\n"
    async for step in agent_generator:
        yield f"data: {json.dumps(step)}\n\n"

@router.get("/chat")
async def chat_endpoint(
    message: Annotated[str, Query(..., description="The user's message to the agent.")],
    http_request: Request,
    conversation_uuid: Annotated[Optional[str], Query(description="The UUID of the ongoing conversation.")] = None,
    orchestrator: AgentOrchestrator = Depends(get_agent_orchestrator),
    conversation_service: ConversationService = Depends(get_conversation_service),
):
    trace_id = getattr(http_request.state, 'trace_id', 'N/A')
    logger.info(
        "Streaming chat request received.",
        extra={'trace_id': trace_id, 'event': 'chat_stream_request_received', 'conversation_uuid': conversation_uuid}
    )

    conversation = None
    new_conversation_uuid = None # To store the new UUID if created
    if conversation_uuid:
        conversation = conversation_service.get_conversation_by_uuid(conversation_uuid, trace_id=trace_id)
        if not conversation:
            raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="Conversation not found.")
    else:
        conversation = conversation_service.create_conversation(trace_id=trace_id)
        new_conversation_uuid = conversation.uuid # Store the new UUID
        logger.info("New conversation created for stream.", extra={'trace_id': trace_id, 'event': 'new_conversation_created', 'uuid': conversation.uuid})

    try:
        agent_generator = orchestrator.run_agent_loop(
            conversation_id=conversation.id,
            user_message=message,
            trace_id=trace_id
        )
        # Pass the new_conversation_uuid to the streaming function
        return StreamingResponse(stream_agent_response(agent_generator, conversation_uuid=new_conversation_uuid), media_type="text/event-stream")

    except Exception as e:
        logger.error(
            f"Error processing streaming chat request: {e}",
            extra={'trace_id': trace_id, 'event': 'chat_stream_request_error'}
        )
        # Note: Cannot raise HTTPException here as the response has already started streaming for some cases.
        # The error will be logged, and the stream will close.
        return StreamingResponse([f'data: {"type": "error", "content": "{str(e)}"}\n\n'], media_type="text/event-stream", status_code=500)
