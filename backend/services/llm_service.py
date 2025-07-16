from typing import Optional, Any
import httpx

from backend.core.config import get_settings
from shared_utils.logging_config import get_logger
from backend.services.llm_providers import BaseLLMProvider, OpenRouterLLM, LLMServiceError

settings = get_settings()
logger = get_logger(__name__)

class LLMService:
    def __init__(self, client: httpx.AsyncClient, frontend_url: str):
        self.llm_provider: BaseLLMProvider = OpenRouterLLM(client, frontend_url)

    async def call_llm(
        self,
        messages: list[dict],
        tools: Optional[list[dict]] = None,
        max_tokens: Optional[int] = None,
        temperature: float = 0.7,
        trace_id: str = "N/A",
    ) -> dict[str, Any]:
        """
        Calls the LLM API using the configured provider.
        """
        return await self.llm_provider.call_llm(
            messages=messages,
            tools=tools,
            max_tokens=max_tokens,
            temperature=temperature,
            trace_id=trace_id,
        )
