from typing import List, Dict

from shared_utils.logging_config import get_logger
from backend.core.config import Settings

logger = get_logger(__name__)

class PromptBuilder:
    def __init__(self, settings: Settings):
        self.settings = settings

    def build_system_prompt(self, relevant_tool_definitions: List[Dict]) -> str:
        """
        Builds the system prompt for the LLM, combining base prompt and tool definitions.
        """
        system_prompt = self.settings.BASE_SYSTEM_PROMPT

        if relevant_tool_definitions:
            tool_instructions = "\n\nAvailable tools:\n"
            for tool_def in relevant_tool_definitions:
                tool_instructions += f'- {tool_def["function"]["name"]}: {tool_def["function"]["description"]}\n'
            system_prompt += tool_instructions

        return system_prompt

    # You can add more methods here for building other parts of the prompt,
    # e.g., for summarizing history if needed.

