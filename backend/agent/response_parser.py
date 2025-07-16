import json
from typing import List, Dict, Any, Tuple, Optional

from shared_utils.logging_config import get_logger

logger = get_logger(__name__)

class LLMResponseParser:
    def __init__(self):
        pass

    def parse_llm_response(self, llm_response: Dict[str, Any], trace_id: str = "N/A") -> Tuple[str, Optional[List[Dict]], Optional[str]]:
        """
        Parses the raw LLM response to determine the agent's action.
        Returns a tuple: (action_type, tool_calls, content)
        action_type can be 'tool_calls', 'final_answer', 'thought', 'direct_response', 'error'.
        """
        # TODO: [BUG-001] LLM Tool Calling Issue: The LLM is not consistently generating 'tool_calls' finish_reason.
        # Attempts to parse tool calls from text were unsuccessful after 5 tries. This bug is skipped.
        # The agent will currently not execute tools if the LLM does not use native tool_calls mechanism.
        choice = llm_response["choices"][0]
        finish_reason = choice["finish_reason"]
        message = choice["message"]

        if finish_reason == "tool_calls":
            tool_calls = message.get("tool_calls")
            logger.info(
                "LLM requested tool calls.",
                extra={'trace_id': trace_id, 'event': 'llm_tool_call_request', 'tool_calls_count': len(tool_calls)}
            )
            return "tool_calls", tool_calls, None

        elif finish_reason == "stop":
            content = message.get("content")
            if content and content.strip().startswith("Final Answer:"):
                final_answer = content.replace("Final Answer:", "").strip()
                logger.info("LLM provided final answer.", extra={'trace_id': trace_id, 'event': 'llm_final_answer', 'answer': final_answer})
                return "final_answer", None, final_answer
            elif content and content.strip().startswith("Thought:"):
                thought = content.strip()
                logger.info("LLM provided a thought.", extra={'trace_id': trace_id, 'event': 'llm_thought', 'thought': thought})
                return "thought", None, thought
            else:
                # Direct response, treat as final answer
                direct_response = content.strip() if content else "No response from AI."
                logger.info("LLM provided direct response.", extra={'trace_id': trace_id, 'event': 'llm_direct_response', 'response': direct_response})
                return "direct_response", None, direct_response

        else:
            logger.warning(
                f"LLM stopped with unexpected finish reason: {finish_reason}",
                extra={'trace_id': trace_id, 'event': 'llm_unexpected_stop', 'reason': finish_reason}
            )
            return "error", None, f"The AI stopped due to: {finish_reason}. Please try again."