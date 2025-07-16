import os
from abc import ABC, abstractmethod
from typing import List, Dict, Optional, Any
import httpx

from backend.core.config import get_settings
from shared_utils.logging_config import get_logger

settings = get_settings()
logger = get_logger(__name__)

class LLMServiceError(Exception):
    """Custom exception for LLM service errors."""
    pass

class BaseLLMProvider(ABC):
    """Abstract base class for LLM providers."""
    def __init__(self, client: httpx.AsyncClient, frontend_url: str):
        self.client = client
        self.frontend_url = frontend_url

    @abstractmethod
    async def call_llm(
        self,
        messages: List[Dict],
        tools: Optional[List[Dict]] = None,
        max_tokens: Optional[int] = None,
        temperature: float = 0.7,
        trace_id: str = "N/A",
    ) -> Dict[str, Any]:
        """Abstract method to call the LLM API."""
        pass

class OpenRouterLLM(BaseLLMProvider):
    """Concrete implementation for OpenRouter LLM provider."""
    def __init__(self, client: httpx.AsyncClient, frontend_url: str):
        super().__init__(client, frontend_url)
        self.api_key = settings.OPENROUTER_API_KEY
        self.model = settings.OPENROUTER_MODEL
        self.base_url = settings.OPENROUTER_BASE_URL

    async def call_llm(
        self,
        messages: List[Dict],
        tools: Optional[List[Dict]] = None,
        max_tokens: Optional[int] = None,
        temperature: float = 0.7,
        trace_id: str = "N/A",
    ) -> Dict[str, Any]:
        """
        Calls the OpenRouter LLM API with the given messages and tools.
        Returns the raw LLM response.
        """
        headers = {
            "Authorization": f"Bearer {self.api_key}",
            "HTTP-Referer": self.frontend_url, # Use the configured frontend URL
            "X-Title": "MatAgent",
        }

        payload = {
            "model": self.model,
            "messages": messages,
            "temperature": temperature,
        }
        if tools:
            payload["tools"] = tools
            payload["tool_choice"] = "auto"
        if max_tokens:
            payload["max_tokens"] = max_tokens

        logger.info(
            "Sending request to LLM",
            extra={
                'trace_id': trace_id,
                'event': 'llm_request_sent',
                'model': self.model,
                'messages_count': len(messages),
                'tools_count': len(tools) if tools else 0,
                'payload': payload # Log the full payload
            }
        )

        try:
            response = await self.client.post(f"{self.base_url}/chat/completions", json=payload, headers=headers)
            response.raise_for_status()  # Raise an exception for 4xx or 5xx responses
            
            llm_response = response.json()
            
            logger.info(
                "Received response from LLM",
                extra={
                    'trace_id': trace_id,
                    'event': 'llm_response_received',
                    'model': self.model,
                    'usage': llm_response.get("usage"),
                    'finish_reason': llm_response.get("choices", [{}])[0].get("finish_reason"),
                    'response_id': llm_response.get("id"),
                    'raw_response': llm_response # Be careful with logging large responses
                }
            )
            return llm_response
        except httpx.HTTPStatusError as e:
            logger.error(
                f"LLM API HTTP error: {e.response.status_code} - {e.response.text}",
                extra={
                    'trace_id': trace_id,
                    'event': 'llm_api_http_error',
                    'status_code': e.response.status_code,
                    'response_text': e.response.text
                }
            )
            raise LLMServiceError(f"LLM API HTTP error: {e.response.status_code} - {e.response.text}") from e
        except httpx.RequestError as e:
            logger.error(
                f"LLM API request error: {e}",
                extra={
                    'trace_id': trace_id,
                    'event': 'llm_api_request_error',
                    'error_message': str(e) if str(e) else repr(e),
                    'request_url': str(e.request.url) if e.request else "N/A",
                    'request_method': e.request.method if e.request else "N/A"
                }
            )
            raise LLMServiceError(f"LLM API request error: {e}") from e
        except Exception as e:
            logger.critical(
                f"An unexpected error occurred in LLMService: {e}",
                extra={
                    'trace_id': trace_id,
                    'event': 'llm_service_critical_error',
                    'error_message': str(e)
                }
            )
            raise LLMServiceError(f"An unexpected error occurred in LLMService: {e}") from e
