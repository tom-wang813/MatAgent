from pydantic_settings import BaseSettings, SettingsConfigDict
from pydantic import Field # Added Field import
from functools import lru_cache

class Settings(BaseSettings):
    # Database settings
    DATABASE_URL: str = Field("sqlite:///./data/matagent.db", description="Path to the SQLite database file.")

    # OpenRouter API settings
    OPENROUTER_API_KEY: str = Field(..., description="API key for OpenRouter service.")
    OPENROUTER_MODEL: str = Field("openai/gpt-4o", description="Default LLM model to use with OpenRouter.")
    OPENROUTER_BASE_URL: str = Field("https://openrouter.ai/api/v1", description="Base URL for the OpenRouter API.")
    OPENROUTER_EMBEDDING_MODEL: str = Field("sentence-transformers/all-MiniLM-L6-v2", description="Embedding model to use with OpenRouter for tool selection.")

    # Agent and LLM limits
    MAX_CONVERSATION_TOKENS: int = Field(80000, description="Maximum tokens for the entire conversation history fed to LLM.")
    MAX_LLM_CALL_TOKENS: int = Field(2000, description="Maximum tokens LLM can generate in a single response.")
    MAX_AGENT_TURNS: int = Field(15, description="Maximum number of thought-action-observation cycles for the agent.")
    LLM_TIMEOUT: int = Field(60, description="Timeout for LLM API calls in seconds.")
    MAX_TOOL_CALLS_PER_RESPONSE: int = Field(20, description="Maximum tool calls LLM can suggest in one response.")

    # Cost limits (for openai/gpt-4o on OpenRouter, prices are approximate and subject to change)
    MAX_CONVERSATION_COST_USD: float = Field(3.00, description="Maximum cost for a single conversation in USD.")
    OPENROUTER_COST_PER_PROMPT_TOKEN_USD: float = Field(0.000005, description="Cost per prompt token in USD for OpenRouter.")
    OPENROUTER_COST_PER_COMPLETION_TOKEN_USD: float = Field(0.000015, description="Cost per completion token in USD for OpenRouter.")

    # Tool selection settings
    TOP_K_SIMILAR_TOOLS: int = Field(5, description="Number of most similar tools to pass to LLM for selection.")

    # Base System Prompt for the Agent
    BASE_SYSTEM_PROMPT: str = Field(
        "You are MatAgent, an expert AI assistant specializing in materials informatics. "
        "Your primary goal is to help users by providing accurate information and performing calculations related to chemical compounds and materials. "
        "\n\n**CRITICAL WORKFLOW - FOLLOW EXACTLY:**\n"
        "1. **STRICT TWO-PHASE APPROACH:**\n"
        "   - **Phase 1 ONLY**: Respond with 'Thought: [what I need to do]' - NO TOOL CALLS\n"
        "   - **Phase 2 ONLY**: Make tool calls based on your previous thought - NO CONTENT\n"
        "2. **NEVER mix thinking and tool calling in the same response**\n"
        "3. **NEVER provide calculated values or final answers in Phase 1**\n"
        "4. **Complete workflow:**\n"
        "   - User asks question\n"
        "   - Response 1: 'Thought: [I need to do X using Y approach]' (NO tool calls)\n"
        "   - Response 2: [Tool calls ONLY] (NO content)\n"
        "   - Response 3: 'Thought: [Analysis of tool results]'\n"
        "   - Response 4: 'Final Answer: [Complete answer]'\n"
        "5. **Formatting rules:**\n"
        "   - If providing thought → NO tool calls in same response\n"
        "   - If making tool calls → NO content in same response\n"
        "   - Each response has ONLY ONE purpose\n"
        "6. **Example correct flow:**\n"
        "   User: 'What is the TPSA of glucose?'\n"
        "   Response 1: 'Thought: I need to calculate the TPSA of glucose using molecular property prediction tools.'\n"
        "   Response 2: [tool call to calculate_tpsa] (no content)\n"
        "   Response 3: 'Thought: The tool calculated TPSA as 110.38 Ų for glucose.'\n"
        "   Response 4: 'Final Answer: The TPSA of glucose is 110.38 Ų.'\n"
        "\n**CRITICAL: NEVER combine thinking and tool calls in one response. Always separate them completely.**",
        description="Base system prompt for the Agent LLM."
    )

    # MCP Server URL for local development/testing
    MCP_SERVER_URL: str = Field("http://localhost:8080", description="URL for the Material Informatics Platform (MCP) server.")
    MCP_TOOL_DISCOVERY_RETRIES: int = Field(10, description="Number of retries for MCP tool discovery on startup.")
    MCP_TOOL_DISCOVERY_RETRY_DELAY: int = Field(5, description="Delay in seconds between retries for MCP tool discovery.")
    CELERY_POLLING_INTERVAL: float = Field(1.0, description="Interval in seconds for polling Celery task status.")
    FRONTEND_URL: str = Field("http://localhost:3000", description="URL for the frontend application, used for HTTP Referer header.")

    model_config = SettingsConfigDict(env_file=".env", extra="ignore")

@lru_cache()
def get_settings():
    """
    Returns a cached instance of the Settings class.
    This ensures settings are loaded only once.
    """
    return Settings()