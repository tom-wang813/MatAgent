import json
from typing import List, Dict, Any, Optional
import httpx # Added httpx import

from shared_utils.logging_config import get_logger
from backend.core.config import Settings
from backend.services.llm_service import LLMService, LLMServiceError
from backend.services.tool_manager import ToolManager, ToolManagerError
from backend.services.conversation_service import ConversationService
from backend.agent.prompt_builder import PromptBuilder
from backend.agent.response_parser import LLMResponseParser
from backend.agent.token_cost_monitor import TokenCostMonitor

logger = get_logger(__name__)

class AgentOrchestrator:
    def __init__(
        self,
        llm_service: LLMService,
        tool_manager: ToolManager,
        conversation_service: ConversationService,
        settings: Settings
    ):
        self.llm_service = llm_service
        self.tool_manager = tool_manager
        self.conversation_service = conversation_service
        self.settings = settings
        self.prompt_builder = PromptBuilder(settings)
        self.response_parser = LLMResponseParser()
        self.token_cost_monitor = TokenCostMonitor(settings)

    async def run_agent_loop(
        self,
        conversation_id: int,
        user_message: str,
        trace_id: str = "N/A"
    ):
        """
        Runs the main Agent loop and yields each step as it occurs.
        """
        current_cost = 0.0
        current_tokens = 0
        turn_count = 0
        
        db_messages = self.conversation_service.get_messages(conversation_id, trace_id=trace_id)
        messages: List[Dict] = []
        
        # This list will temporarily hold messages for the current turn,
        # allowing us to group assistant tool_calls with their tool_outputs.
        current_turn_messages: List[Dict] = []
        
        for msg in db_messages:
            if msg.role == "assistant":
                # If we encounter a new assistant message, and there are messages in current_turn_messages,
                # it means the previous turn is complete. Append them to the main messages list.
                if current_turn_messages:
                    messages.extend(current_turn_messages)
                    current_turn_messages = []

                try:
                    parsed_content = json.loads(msg.content)
                    assistant_message_dict = {"role": msg.role}
                    if "tool_calls" in parsed_content:
                        assistant_message_dict["tool_calls"] = parsed_content["tool_calls"]
                        current_turn_messages.append(assistant_message_dict)
                    else:
                        if "content" in parsed_content:
                            assistant_message_dict["content"] = parsed_content["content"]
                        elif "thought" in parsed_content:
                            assistant_message_dict["content"] = parsed_content["thought"]
                        
                        if not assistant_message_dict.get("tool_calls") and not assistant_message_dict.get("content"):
                            assistant_message_dict["content"] = msg.content
                        
                        messages.append(assistant_message_dict)

                except json.JSONDecodeError:
                    messages.append({"role": msg.role, "content": msg.content})
            elif msg.role == "tool":
                # Tool messages always follow an assistant message with tool_calls.
                # Append them to current_turn_messages.
                current_turn_messages.append({"role": msg.role, "tool_call_id": msg.tool_call_id, "content": msg.content})
            else:
                # If we encounter a user message or any other role, and there are messages in current_turn_messages,
                # it means the previous turn (assistant + tools) is complete.
                if current_turn_messages:
                    messages.extend(current_turn_messages)
                    current_turn_messages = []
                messages.append({"role": msg.role, "content": msg.content})

        # After the loop, append any remaining messages in current_turn_messages (for the last turn).
        if current_turn_messages:
            messages.extend(current_turn_messages)

        self.conversation_service.add_message(conversation_id, "user", user_message, trace_id=trace_id)
        messages.append({"role": "user", "content": user_message})

        while True:
            turn_count += 1
            logger.info(f"Agent loop turn {turn_count} started.", extra={'trace_id': trace_id, 'event': 'agent_turn_start', 'turn': turn_count})

            limit_exceeded_message = self.token_cost_monitor.check_limits(
                current_cost=current_cost,
                current_tokens=current_tokens,
                turn_count=turn_count,
                trace_id=trace_id
            )
            if limit_exceeded_message:
                yield {"type": "final_answer", "content": limit_exceeded_message}
                break

            query_for_tools = messages[-1]["content"]
            relevant_tool_definitions = await self.tool_manager.get_relevant_tool_definitions(
                query_text=query_for_tools,
                trace_id=trace_id,
                top_k=self.settings.TOP_K_SIMILAR_TOOLS
            )
            logger.info(f"Selected {len(relevant_tool_definitions)} relevant tools.", extra={'trace_id': trace_id, 'event': 'tool_selection_complete', 'tools': [t["function"]["name"] for t in relevant_tool_definitions]})

            system_prompt = self.prompt_builder.build_system_prompt(relevant_tool_definitions)
            llm_messages = [{"role": "system", "content": system_prompt}] + messages

            estimated_prompt_tokens = self.token_cost_monitor.estimate_tokens_from_messages(llm_messages)
            if estimated_prompt_tokens > self.settings.MAX_CONVERSATION_TOKENS:
                yield {"type": "final_answer", "content": "Conversation history too long. Please start a new conversation."}
                logger.warning("Conversation history exceeded max tokens.", extra={'trace_id': trace_id, 'event': 'token_limit_exceeded', 'estimated_tokens': estimated_prompt_tokens})
                break

            try:
                llm_response = await self.llm_service.call_llm(
                    messages=llm_messages,
                    tools=relevant_tool_definitions,
                    max_tokens=self.settings.MAX_LLM_CALL_TOKENS,
                    trace_id=trace_id
                )
                
                usage = llm_response.get("usage", {})
                prompt_tokens = usage.get("prompt_tokens", 0)
                completion_tokens = usage.get("completion_tokens", 0)
                
                if "total_cost_usd" in usage:
                    current_cost += usage["total_cost_usd"]
                else:
                    current_cost += self.token_cost_monitor.calculate_cost(prompt_tokens, completion_tokens)
                current_tokens += usage.get("total_tokens", 0)

                action_type, tool_calls, content = self.response_parser.parse_llm_response(llm_response, trace_id)
                
                # TODO: [BUG-001] LLM Tool Calling Issue: The LLM is not consistently generating 'tool_calls' finish_reason.
                # Attempts to parse tool calls from text were unsuccessful after 5 tries. This bug is skipped.
                # The agent will currently not execute tools if the LLM does not use native tool_calls mechanism.
                if action_type == "tool_calls":
                    yield {"type": "tool_calls", "content": tool_calls}
                    if len(tool_calls) > self.settings.MAX_TOOL_CALLS_PER_RESPONSE:
                        yield {"type": "final_answer", "content": "The AI tried to call too many tools at once."}
                        logger.warning("LLM attempted too many tool calls.", extra={'trace_id': trace_id, 'event': 'too_many_tool_calls', 'count': len(tool_calls)})
                        break

                    messages.append(llm_response["choices"][0]["message"])
                    self.conversation_service.add_message(conversation_id, "assistant", json.dumps(llm_response["choices"][0]["message"]), trace_id=trace_id)

                    tool_outputs = []
                    for tool_call in tool_calls:
                        tool_name = tool_call["function"]["name"]
                        tool_args = json.loads(tool_call["function"]["arguments"])
                        
                        logger.info(f"Executing tool: {tool_name}", extra={'trace_id': trace_id, 'event': 'tool_execution_request', 'tool_name': tool_name, 'tool_args': tool_args})
                        try:
                            tool_result = await self.tool_manager.execute_tool(
                                tool_name=tool_name,
                                tool_args=tool_args,
                                trace_id=trace_id
                            )
                            tool_outputs.append({"tool_call_id": tool_call["id"], "output": tool_result})
                            logger.info(f"Tool {tool_name} executed successfully.", extra={'trace_id': trace_id, 'event': 'tool_execution_success', 'tool_name': tool_name})
                        except (ToolManagerError, Exception) as e:
                            tool_outputs.append({"tool_call_id": tool_call["id"], "output": f"Error: {str(e)}"})
                            logger.error(f"Tool {tool_name} execution failed: {e}", extra={'trace_id': trace_id, 'event': 'tool_execution_failure', 'tool_name': tool_name, 'error': str(e)})
                    
                    for output in tool_outputs:
                        yield {"type": "tool_output", "content": json.dumps(output["output"])}
                        messages.append({"role": "tool", "tool_call_id": output["tool_call_id"], "content": output["output"]})
                        self.conversation_service.add_message(conversation_id, "tool", json.dumps(output["output"]), tool_call_id=output["tool_call_id"], trace_id=trace_id)

                    continue

                elif action_type == "final_answer" or action_type == "direct_response":
                    yield {"type": "final_answer", "content": content}
                    messages.append(llm_response["choices"][0]["message"])
                    self.conversation_service.add_message(conversation_id, "assistant", json.dumps(llm_response["choices"][0]["message"]), trace_id=trace_id)
                    logger.info("Agent provided final answer.", extra={'trace_id': trace_id, 'event': 'agent_final_answer'})
                    break

                elif action_type == "thought":
                    yield {"type": "thought", "content": content}
                    messages.append(llm_response["choices"][0]["message"])
                    self.conversation_service.add_message(conversation_id, "assistant", json.dumps(llm_response["choices"][0]["message"]), trace_id=trace_id)
                    logger.info("Agent provided a thought.", extra={'trace_id': trace_id, 'event': 'agent_thought', 'thought': content})
                    continue

                elif action_type == "error":
                    yield {"type": "final_answer", "content": content}
                    break

            except LLMServiceError as e:
                yield {"type": "final_answer", "content": f"An error occurred with the AI service: {e}"}
                logger.error(f"LLM service error: {e}", extra={'trace_id': trace_id, 'event': 'llm_service_error'})
                break
            except Exception as e:
                yield {"type": "final_answer", "content": f"An unexpected error occurred: {e}"}
                logger.critical(f"Critical error in agent loop: {e}", extra={'trace_id': trace_id, 'event': 'agent_critical_error'})
                break
        
        logger.info("Agent loop finished.", extra={'trace_id': trace_id, 'event': 'agent_loop_end', 'total_cost_usd': current_cost, 'turns_taken': turn_count})
