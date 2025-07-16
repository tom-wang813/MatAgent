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
        "**ALWAYS** think step by step and explain your reasoning clearly before taking any action. "
        "**For any calculation or factual lookup, you MUST use the available tools.** "
        "Specifically, for molecular weight calculations, you MUST use the `calculate_molecular_weight` tool. "
        "When asked to calculate molecular similarity, **you MUST directly use the `calculate_molecular_similarity` tool.** Do not attempt to generate fingerprints separately. "
        "If the user does not specify a `metric` or `fingerprint_type` for similarity calculations, you **MUST** ask the user to choose from the available options (`tanimoto`, `dice`, `cosine` for metrics; `morgan`, `rdkit` for fingerprint types). "
        "When you need to use a tool, respond ONLY with a tool call in the format: `tool_name(param1=value1, param2=value2)`. "
        "Do NOT include any other text or explanation when making a tool call. "
        "Before making a tool call, output your thought process by starting your response with 'Thought:'. "
        "If you decide to use a tool, your response MUST be a JSON object with the tool call, and nothing else."
        "After executing tools and gathering information, synthesize the results and provide a clear 'Final Answer:'. "
        "If you cannot find an answer or perform a calculation using the tools, explain why.",
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