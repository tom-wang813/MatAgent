import tiktoken
from typing import List, Dict, Optional

from shared_utils.logging_config import get_logger
from backend.core.config import Settings

logger = get_logger(__name__)

class TokenCostMonitor:
    def __init__(self, settings: Settings):
        self.settings = settings
        # Explicitly use cl100k_base for gpt-4o and similar models
        # If self.settings.OPENROUTER_MODEL can be other models, a more robust solution might be needed.
        self.tokenizer = tiktoken.get_encoding("cl100k_base")

    def estimate_tokens_from_messages(self, messages: List[Dict]) -> int:
        """
        Estimates the number of tokens in a list of messages using tiktoken.
        """
        num_tokens = 0
        for message in messages:
            num_tokens += 4  # every message follows <im_start>{role/name}\n{content}<im_end>\n
            for key, value in message.items():
                if key == "content" and value is not None:
                    # Ensure value is a string
                    content_str = str(value) if not isinstance(value, str) else value
                    num_tokens += len(self.tokenizer.encode(content_str))
                elif key == "tool_calls" and value is not None:
                    for tool_call in value:
                        name_str = str(tool_call["function"]["name"])
                        args_str = str(tool_call["function"]["arguments"])
                        num_tokens += len(self.tokenizer.encode(name_str))
                        num_tokens += len(self.tokenizer.encode(args_str))
                elif key == "name" and value is not None:
                    name_str = str(value)
                    num_tokens += len(self.tokenizer.encode(name_str))
                elif key in ["role", "tool_call_id"] and value is not None:
                    # Convert any other string-like values to string before encoding
                    value_str = str(value)
                    num_tokens += len(self.tokenizer.encode(value_str))
                # Add other keys if they contribute to token count (e.g., tool_call_id for tool messages)

        num_tokens += 2  # every reply is primed with <im_start>assistant
        return num_tokens

    def calculate_cost(self, prompt_tokens: int, completion_tokens: int) -> float:
        """
        Calculates the estimated cost based on token usage.
        """
        cost = (
            prompt_tokens * self.settings.OPENROUTER_COST_PER_PROMPT_TOKEN_USD +
            completion_tokens * self.settings.OPENROUTER_COST_PER_COMPLETION_TOKEN_USD
        )
        return cost

    def check_limits(
        self,
        current_cost: float,
        current_tokens: int,
        turn_count: int,
        trace_id: str = "N/A"
    ) -> Optional[str]:
        """
        Checks if any conversation limits have been reached.
        Returns an error message if a limit is reached, otherwise None.
        """
        if current_cost >= self.settings.MAX_CONVERSATION_COST_USD:
            logger.warning("Conversation reached max cost limit.", extra={'trace_id': trace_id, 'event': 'cost_limit_reached', 'cost': current_cost})
            return "I'm sorry, this conversation has reached its maximum cost limit. Please try a simpler query or start a new conversation."
        
        if turn_count >= self.settings.MAX_AGENT_TURNS:
            logger.warning("Agent reached max turns limit.", extra={'trace_id': trace_id, 'event': 'turn_limit_reached', 'turns': turn_count})
            return "I'm sorry, the AI reached its maximum thinking steps. Please try to simplify your request."
        
        # Note: MAX_CONVERSATION_TOKENS is checked before LLM call in orchestrator
        # This check is for the overall conversation, not just the prompt.
        # If we want to check total tokens here, we need to pass estimated_prompt_tokens + completion_tokens
        # For now, it's handled in orchestrator before LLM call.

        return None
