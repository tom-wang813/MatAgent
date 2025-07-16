import httpx
from typing import List, Dict, Any, Optional
import asyncio
import numpy as np
from tenacity import retry, stop_after_attempt, wait_fixed, retry_if_exception_type

from shared_utils.logging_config import get_logger
from backend.services.embedding_service import EmbeddingService
from backend.core.config import Settings, get_settings
from backend.core.responses import SuccessResponse, ErrorResponse # Import from new location

settings = get_settings()

logger = get_logger(__name__)

class ToolManagerError(Exception):
    """Custom exception for ToolManager errors."""
    pass

class ToolManager:
    def __init__(self, embedding_service: EmbeddingService, client: httpx.AsyncClient, mcp_server_url: str, settings: Settings):
        self.mcp_server_url = mcp_server_url
        self.settings = settings
        self.client = client
        self.embedding_service = embedding_service
        self._available_tools: List[Dict] = []
        self._tool_embeddings: Optional[np.ndarray] = None
        self._tool_names: List[str] = []

    @retry(stop=stop_after_attempt(settings.MCP_TOOL_DISCOVERY_RETRIES), wait=wait_fixed(settings.MCP_TOOL_DISCOVERY_RETRY_DELAY), retry=retry_if_exception_type(httpx.RequestError), reraise=True)
    async def discover_tools(self, trace_id: str = "N/A"):
        logger.info(
            "Attempting to discover tools from MCP server.",
            extra={'trace_id': trace_id, 'event': 'tool_discovery_start'}
        )
        try:
            response = await self.client.get(f"{self.mcp_server_url}/api/tools")
            response.raise_for_status()
            # Expecting SuccessResponse from MCP server
            parsed_response = SuccessResponse[List[Dict]].model_validate(response.json())
            raw_tools = parsed_response.data

            self._available_tools = []
            self._tool_names = []
            tool_descriptions = []

            for tool_data in raw_tools:
                # Convert MCP server tool metadata to OpenAI function calling format
                openai_tool_format = {
                    "type": "function",
                    "function": {
                        "name": tool_data["name"],
                        "description": tool_data["description"],
                        "parameters": tool_data.get("parameters", {"type": "object", "properties": {}, "required": []})
                    }
                }
                # Add tool_type to the OpenAI format for internal use
                openai_tool_format["tool_type"] = tool_data.get("tool_type", "basic")
                self._available_tools.append(openai_tool_format)
                self._tool_names.append(tool_data["name"])
                tool_descriptions.append(tool_data.get("description", tool_data["name"]))
            
            # Generate embeddings for tool descriptions
            if tool_descriptions:
                # For simplicity, we'll generate embeddings one by one. 
                # For many tools, consider batching or pre-calculating.
                embeddings = []
                for desc in tool_descriptions:
                    embeddings.append(self.embedding_service.get_embedding(desc, trace_id=trace_id))
                self._tool_embeddings = np.array(embeddings)
                logger.info(f"Tool embeddings shape: {self._tool_embeddings.shape}", extra={'trace_id': trace_id, 'event': 'tool_embeddings_generated'})

            logger.info(
                f"Successfully discovered {len(self._available_tools)} tools from MCP server.",
                extra={'trace_id': trace_id, 'event': 'tool_discovery_end', 'num_tools': len(self._available_tools)}
            )

        except httpx.HTTPStatusError as e:
            logger.error(
                f"Failed to discover tools from MCP server: HTTP error {e.response.status_code} - {e.response.text}",
                extra={'trace_id': trace_id, 'event': 'tool_discovery_error', 'status_code': e.response.status_code}
            )
            raise ToolManagerError(f"Failed to discover tools from MCP server: {e}") from e
        except Exception as e:
            logger.critical(
                f"An unexpected error occurred during tool discovery: {e}",
                extra={'trace_id': trace_id, 'event': 'tool_discovery_critical_error', 'error_message': str(e)}
            )
            raise ToolManagerError(f"An unexpected error occurred during tool discovery: {e}") from e

    def get_all_tool_definitions(self) -> List[Dict]:
        """
        Returns all discovered tools in OpenAI function calling format.
        """
        return self._available_tools

    async def get_relevant_tool_definitions(
        self,
        query_text: str,
        trace_id: str = "N/A",
        top_k: int = None # Change default to None, will use self.settings.TOP_K_SIMILAR_TOOLS if None
    ) -> List[Dict]:
        """
        Retrieves the top_k most relevant tool definitions based on a query text.
        """
        if top_k is None:
            top_k = self.settings.TOP_K_SIMILAR_TOOLS
        if not self._available_tools or self._tool_embeddings is None:
            logger.warning(
                "No tools discovered or tool embeddings not generated. Returning all tools.",
                extra={'trace_id': trace_id, 'event': 'tool_retrieval_warning'}
            )
            return self._available_tools

        query_embedding = self.embedding_service.get_embedding(query_text, trace_id=trace_id)
        query_embedding = np.array(query_embedding).reshape(1, -1)

        # Calculate cosine similarity
        # (A . B) / (||A|| * ||B||)
        similarities = np.dot(query_embedding, self._tool_embeddings.T) / \
                       (np.linalg.norm(query_embedding) * np.linalg.norm(self._tool_embeddings, axis=1))
        
        # Get top_k indices
        top_k_indices = np.argsort(similarities[0])[::-1][:top_k]

        relevant_tools = [self._available_tools[i] for i in top_k_indices]

        logger.info(
            f"Retrieved {len(relevant_tools)} relevant tools for query.",
            extra={'trace_id': trace_id, 'event': 'tool_retrieval_success', 'query': query_text, 'relevant_tools': [t["function"]["name"] for t in relevant_tools]}
        )
        return relevant_tools

    async def execute_tool(
        self,
        tool_name: str,
        tool_args: Dict[str, Any],
        trace_id: str = "N/A"
    ) -> Any:
        """
        Executes a tool on the MCP server.
        If the tool is an ML model, it polls for the result.
        """
        logger.info(
            f"Executing tool '{tool_name}' on MCP server.",
            extra={'trace_id': trace_id, 'event': 'tool_execution_start', 'tool_name': tool_name, 'tool_args': tool_args}
        )
        try:
            # Find the tool's metadata to check its type
            tool_metadata = next((t for t in self._available_tools if t["function"]["name"] == tool_name), None)
            if not tool_metadata:
                raise ToolManagerError(f"Tool '{tool_name}' not found in discovered tools.")

            payload = {"tool_name": tool_name, "params": tool_args}
            headers = {"X-Trace-ID": trace_id} # Pass trace_id to MCP server
            response = await self.client.post(f"{self.mcp_server_url}/api/run_tool", json=payload, headers=headers)
            response.raise_for_status()
            
            # Expecting SuccessResponse from MCP server
            parsed_response = SuccessResponse[Any].model_validate(response.json())
            result_data = parsed_response.data

            if tool_metadata.get("tool_type") == "ml_model":
                task_id = result_data  # For ML models, run_tool returns the task_id
                logger.info(f"ML model tool '{tool_name}' submitted as task {task_id}. Polling for result.", extra={'trace_id': trace_id, 'event': 'ml_tool_submitted', 'tool_name': tool_name, 'task_id': task_id})
                
                # Poll for task status
                while True:
                    await asyncio.sleep(self.settings.CELERY_POLLING_INTERVAL) # Poll every X seconds
                    task_status_response = await self.client.get(f"{self.mcp_server_url}/api/tasks/{task_id}/status", headers=headers)
                    task_status_response.raise_for_status()
                    
                    # Expecting SuccessResponse from MCP server
                    parsed_task_status = SuccessResponse[Dict[str, Any]].model_validate(task_status_response.json())
                    status_data = parsed_task_status.data

                    status = status_data["status"]
                    logger.debug(f"Task {task_id} status: {status}", extra={'trace_id': trace_id, 'event': 'task_polling', 'task_id': task_id, 'status': status})

                    if status == "SUCCESS":
                        final_result = status_data["result"]
                        logger.info(f"ML model tool '{tool_name}' (task {task_id}) completed successfully.", extra={'trace_id': trace_id, 'event': 'ml_tool_completed', 'tool_name': tool_name, 'task_id': task_id})
                        return final_result
                    elif status == "FAILURE":
                        error_message = status_data.get("error", "Unknown error")
                        logger.error(f"ML model tool '{tool_name}' (task {task_id}) failed: {error_message}", extra={'trace_id': trace_id, 'event': 'ml_tool_failed', 'tool_name': tool_name, 'task_id': task_id, 'error': error_message})
                        raise ToolManagerError(f"ML model tool '{tool_name}' (task {task_id}) failed: {error_message}")
                    # Continue polling if status is PENDING, STARTED, etc.
            else:
                # For basic tools, the result is returned directly
                result = result_data
                logger.info(
                    f"Tool '{tool_name}' executed successfully.",
                    extra={'trace_id': trace_id, 'event': 'tool_execution_end', 'tool_name': tool_name, 'result': result}
                )
                return result

        except httpx.HTTPStatusError as e:
            logger.error(
                f"Tool execution failed on MCP server: HTTP error {e.response.status_code} - {e.response.text}",
                extra={'trace_id': trace_id, 'event': 'tool_execution_error', 'tool_name': tool_name, 'status_code': e.response.status_code, 'response_text': e.response.text}
            )
            raise ToolManagerError(f"Tool '{tool_name}' execution failed: {e.response.status_code} - {e.response.text}") from e
        except httpx.RequestError as e:
            logger.error(
                f"Failed to connect to MCP server for tool execution: {e}",
                extra={'trace_id': trace_id, 'event': 'tool_execution_error', 'tool_name': tool_name, 'error_message': str(e)}
            )
            raise ToolManagerError(f"Failed to connect to MCP server for tool execution: {e}") from e
        except Exception as e:
            logger.critical(
                f"An unexpected error occurred during tool execution for '{tool_name}': {e}",
                extra={'trace_id': trace_id, 'event': 'tool_execution_critical_error', 'tool_name': tool_name, 'error_message': str(e)}
            )
            raise ToolManagerError(f"An unexpected error occurred during tool execution for '{tool_name}': {e}") from e

async def test_ml_tool_execution():
    """A simple test function to demonstrate ML tool execution and polling."""
    from backend.services.embedding_service import EmbeddingService
    from backend.core.config import get_settings
    import httpx

    settings = get_settings()
    embedding_service = EmbeddingService()
    client = httpx.AsyncClient()
    tool_manager = ToolManager(
        embedding_service=embedding_service,
        client=client,
        mcp_server_url=settings.MCP_SERVER_URL,
        settings=settings
    )

    try:
        await tool_manager.discover_tools()
        logger.info("Tools discovered for testing.")

        # Test an ML model tool
        smiles = "CCO"
        logger.info(f"Attempting to predict solubility for {smiles}")
        result = await tool_manager.execute_tool("predict_aqueous_solubility", {"smiles_string": smiles})
        logger.info(f"Prediction result for {smiles}: {result}")

    except ToolManagerError as e:
        logger.error(f"Test failed: {e}")
    except Exception as e:
        logger.error(f"An unexpected error occurred during test: {e}")
    finally:
        await client.aclose()

if __name__ == "__main__":
    # This block will only run if the script is executed directly
    # For testing within the Docker environment, you'd typically run this via pytest or a custom script
    asyncio.run(test_ml_tool_execution())
