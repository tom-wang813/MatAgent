import React, { useState, useRef, useEffect } from 'react';
import { Box, TextField, Button, Typography, CircularProgress, Paper, Avatar, Chip, Accordion, AccordionSummary, AccordionDetails } from '@mui/material';
import SendIcon from '@mui/icons-material/Send';
import SmartToyIcon from '@mui/icons-material/SmartToy';
import PersonIcon from '@mui/icons-material/Person';
import ErrorOutlineIcon from '@mui/icons-material/ErrorOutline';
import InfoOutlinedIcon from '@mui/icons-material/InfoOutlined';
import ExpandMoreIcon from '@mui/icons-material/ExpandMore';
import ExpandLessIcon from '@mui/icons-material/ExpandLess';
import SettingsIcon from '@mui/icons-material/Settings';
import EnhancedMarkdown from './EnhancedMarkdown';
import ToolCallDisplay from './ToolCallDisplay';
import { sendMessage } from '../utils/api'; // Import sendMessage from api.js

const MAX_TOOL_OUTPUT_LENGTH = 200; // Define max length for tool output display

// Helper function to group agent steps into logical turns
const groupAgentSteps = (steps) => {
  if (!steps || !Array.isArray(steps)) {
    return [];
  }
  
  const grouped = [];
  let currentGroup = null;
  let toolOutputs = []; // Track multiple tool outputs for a single tool_calls

  steps.forEach(step => {
    if (step.type === 'thought') {
      // If we have a previous group, push it before starting a new one
      if (currentGroup) {
        if (toolOutputs.length > 0) {
          currentGroup.tool_outputs = toolOutputs;
          toolOutputs = [];
        }
        grouped.push(currentGroup);
      }
      currentGroup = { thought: step.content, tool_calls: [], tool_outputs: [] };
    } else if (step.type === 'tool_calls') {
      if (!currentGroup) {
        currentGroup = { thought: null, tool_calls: [], tool_outputs: [] };
      }
      // Concatenate all tool calls for this group
      currentGroup.tool_calls = currentGroup.tool_calls.concat(step.content);
    } else if (step.type === 'tool_output') {
      // Collect tool outputs - there might be multiple for one tool_calls
      try {
        const outputData = typeof step.content === 'string' ? JSON.parse(step.content) : step.content;
        toolOutputs.push(outputData);
      } catch (e) {
        toolOutputs.push(step.content);
      }
    } else if (step.type === 'final_answer') {
      // Final answer ends any current group
      if (currentGroup) {
        if (toolOutputs.length > 0) {
          currentGroup.tool_outputs = toolOutputs;
          toolOutputs = [];
        }
        grouped.push(currentGroup);
        currentGroup = null;
      }
    }
  });

  // Push any remaining group
  if (currentGroup) {
    if (toolOutputs.length > 0) {
      currentGroup.tool_outputs = toolOutputs;
    }
    grouped.push(currentGroup);
  }

  return grouped;
};

const ToolOutputContent = ({ output }) => {
  const [showFull, setShowFull] = useState(false);
  const isLong = output && output.length > MAX_TOOL_OUTPUT_LENGTH;
  const displayedOutput = isLong && !showFull ? output.substring(0, MAX_TOOL_OUTPUT_LENGTH) + '...' : output;

  return (
    <>
      <Typography component="div" sx={{ mt: 1 }}>
        <EnhancedMarkdown>{displayedOutput}</EnhancedMarkdown>
      </Typography>
      {isLong && (
        <Button 
          onClick={() => setShowFull(!showFull)}
          size="small"
          endIcon={showFull ? <ExpandLessIcon /> : <ExpandMoreIcon />}
          sx={{ mt: 0.5 }} // Adjusted margin top for consistency
        >
          {showFull ? 'Show Less' : 'Show More'}
        </Button>
      )}
    </>
  );
};

// Component for collapsible agent steps
const AgentStepsDisplay = ({ agentSteps, defaultExpanded = false }) => {
  const [expanded, setExpanded] = useState(defaultExpanded);
  
  if (!agentSteps || agentSteps.length === 0) {
    return null;
  }
  
  const groupedSteps = groupAgentSteps(agentSteps);
  
  if (!groupedSteps || groupedSteps.length === 0) {
    return null;
  }

  const getSummary = () => {
    const thoughtCount = groupedSteps.filter(step => step.thought).length;
    const toolCallCount = groupedSteps.reduce((acc, step) => acc + (step.tool_calls?.length || 0), 0);
    const outputCount = groupedSteps.reduce((acc, step) => acc + (step.tool_outputs?.length || 0), 0);
    
    return `${groupedSteps.length} 步骤 • ${thoughtCount} 思考 • ${toolCallCount} 工具调用 • ${outputCount} 输出`;
  };

  return (
    <Accordion 
      expanded={expanded} 
      onChange={(event, isExpanded) => setExpanded(isExpanded)}
      sx={{ mb: 1, '&:before': { display: 'none' } }}
    >
      <AccordionSummary
        expandIcon={<ExpandMoreIcon />}
        sx={{ 
          bgcolor: 'grey.50', 
          minHeight: 40,
          '& .MuiAccordionSummary-content': { margin: '8px 0' }
        }}
      >
        <Box sx={{ display: 'flex', alignItems: 'center', gap: 1 }}>
          <SettingsIcon fontSize="small" color="action" />
          <Typography variant="body2" color="text.secondary">
            {getSummary()}
          </Typography>
        </Box>
      </AccordionSummary>
      <AccordionDetails sx={{ pt: 1 }}>
        {groupedSteps.map((stepGroup, groupIndex) => (
          <Box key={groupIndex} sx={{ mb: 1, p: 1, bgcolor: 'grey.50', borderRadius: 1 }}>
            <Typography variant="subtitle2">步骤 {groupIndex + 1}</Typography>
            {stepGroup.thought && (
              <Typography variant="body2" sx={{ mt: 0.5 }}>
                <strong>Agent Thought:</strong> <EnhancedMarkdown>{stepGroup.thought}</EnhancedMarkdown>
              </Typography>
            )}
            {stepGroup.tool_calls && stepGroup.tool_calls.length > 0 && (
              <Box sx={{ mt: 0.5 }}>
                <Typography variant="body2"><strong>Tool Calls:</strong></Typography>
                <ToolCallDisplay toolCalls={stepGroup.tool_calls} />
              </Box>
            )}
            {stepGroup.tool_outputs && stepGroup.tool_outputs.length > 0 && (
              <Box sx={{ mt: 0.5 }}>
                <Typography variant="body2"><strong>Tool Outputs:</strong></Typography>
                {stepGroup.tool_outputs.map((output, outputIndex) => (
                  <Box key={outputIndex} sx={{ mt: 0.5, pl: 1, borderLeft: '2px solid', borderColor: 'grey.300' }}>
                    <Typography variant="caption" color="text.secondary">Output {outputIndex + 1}:</Typography>
                    <ToolOutputContent output={JSON.stringify(output)} />
                  </Box>
                ))}
              </Box>
            )}
          </Box>
        ))}
      </AccordionDetails>
    </Accordion>
  );
};

const ChatWindow = ({
  conversation,
  onSendMessage,
  onUpdateTitle,
  isLoading,
  setIsLoading,
  onUpdateConversationUuid // New prop to update conversation_uuid
}) => {
  const [inputValue, setInputValue] = useState('');
  const [currentAssistantMessage, setCurrentAssistantMessage] = useState(null);
  const eventSourceRef = useRef(null); // To store the EventSource instance
  const messagesEndRef = useRef(null);
  const textareaRef = useRef(null);
  const [currentConversationUuid, setCurrentConversationUuid] = useState(conversation?.uuid || null); // New state for conversation_uuid

  useEffect(() => {
    if (conversation?.uuid) {
      setCurrentConversationUuid(conversation.uuid);
    }
  }, [conversation?.uuid]);

  const scrollToBottom = () => {
    messagesEndRef.current?.scrollIntoView({ behavior: 'smooth' });
  };

  useEffect(() => {
    scrollToBottom();
  }, [conversation?.messages, currentAssistantMessage]);

  const handleKeyPress = (e) => {
    if (e.key === 'Enter' && !e.shiftKey) {
      e.preventDefault();
      handleSubmit(e);
    }
  };

  const handleSubmit = async (e) => {
    e.preventDefault();
    if (!inputValue.trim() || isLoading) return;

    // 如果没有conversation，自动创建一个新的
    if (!conversation) {
      console.warn('No conversation found, this should not happen in normal flow');
      return;
    }

    const userMessage = {
      id: Date.now().toString(),
      type: 'user',
      content: inputValue.trim(),
      timestamp: new Date().toISOString()
    };

    onSendMessage(conversation.id, userMessage);

    if (conversation.messages.length === 0) {
      const title = inputValue.trim().length > 30
        ? inputValue.trim().substring(0, 30) + '...'
        : inputValue.trim();
      onUpdateTitle(conversation.id, title);
    }

    const currentInput = inputValue;
    setInputValue('');
    setIsLoading(true);

    // Initialize a new assistant message for streaming
    const newAssistantMessage = {
      id: (Date.now() + 1).toString(),
      type: 'assistant',
      content: '',
      timestamp: new Date().toISOString(),
      agent_steps: [],
      status: 'pending',
    };
    setCurrentAssistantMessage(newAssistantMessage);

    // Close any existing EventSource connection before starting a new one
    if (eventSourceRef.current) {
      eventSourceRef.current.close();
    }

    try {
      eventSourceRef.current = sendMessage(
        currentInput,
        currentConversationUuid,
        (step) => { // onNewStep callback
          setCurrentAssistantMessage(prev => {
            // Ensure prev is not null before spreading
            const updatedMessage = prev ? { ...prev } : {
              id: (Date.now() + 1).toString(),
              type: 'assistant',
              content: '',
              timestamp: new Date().toISOString(),
              agent_steps: [],
              status: 'pending',
            };

            if (step.type === 'conversation_info') {
              if (step.conversation_uuid) {
                setCurrentConversationUuid(step.conversation_uuid);
                onUpdateConversationUuid(conversation.id, step.conversation_uuid);
              }
            } else if (step.type === 'text_chunk') {
              updatedMessage.content += step.content;
            } else if (step.type === 'final_answer') {
              updatedMessage.content = step.content;
              updatedMessage.status = 'success';
            } else if (step.type === 'error') {
              updatedMessage.content = step.content;
              updatedMessage.status = 'error';
            } else if (step.type === 'end') {
              // This marks the end of the stream
              updatedMessage.status = updatedMessage.status === 'error' ? 'error' : 'success'; // Preserve error status if any
              setIsLoading(false);
              if (eventSourceRef.current) {
                eventSourceRef.current.close();
              }
            } else {
              // For thought, tool_calls, tool_output, append to agent_steps
              updatedMessage.agent_steps = [...updatedMessage.agent_steps, step];
            }
            return updatedMessage;
          });
        },
        (error) => { // onError callback
          console.error('Error during streaming message:', error);
          let errorMessage = '連接錯誤';
          
          if (error.message) {
            if (error.message.includes('SSE連接錯誤')) {
              errorMessage = '與服務器的連接中斷，請檢查網絡連接並重試';
            } else if (error.message.includes('fetch')) {
              errorMessage = '無法連接到服務器，請確保服務正在運行';
            } else if (error.message.includes('Invalid URL')) {
              errorMessage = 'API 配置錯誤，請聯繫管理員';
            } else {
              errorMessage = error.message;
            }
          }
          
          setCurrentAssistantMessage(prev => ({
            ...prev,
            content: `錯誤: ${errorMessage}`,
            status: 'error',
          }));
          setIsLoading(false); // Stop loading on error
        },
        () => { // onOpen callback
          // Optional: handle connection open
          console.log('SSE connection opened successfully.');
        }
      );
    } catch (error) {
      console.error('Error initiating streaming message:', error);
      let errorMessage = '無法啟動對話';
      
      if (error.message) {
        if (error.message.includes('Invalid URL')) {
          errorMessage = 'API 配置錯誤，請聯繫管理員';
        } else if (error.message.includes('fetch')) {
          errorMessage = '無法連接到服務器，請確保服務正在運行';
        } else {
          errorMessage = error.message;
        }
      }
      
      setCurrentAssistantMessage(prev => ({
        ...prev,
        content: `抱歉，无法启动对话：${errorMessage}。请稍后再试。`,
        status: 'error',
      }));
      setIsLoading(false);
    } finally {
      // The EventSource will be closed by the onError or when the stream naturally ends.
      // We set isLoading to false when the 'end' signal is received or on error.
      // The 'end' signal handling needs to be inside onNewStep for type 'end'.
    }
  };

  // Effect to close EventSource when component unmounts
  useEffect(() => {
    return () => {
      if (eventSourceRef.current) {
        eventSourceRef.current.close();
      }
    };
  }, []);

  // Effect to add the completed assistant message to conversation messages
  useEffect(() => {
    // Only add if currentAssistantMessage exists, is not loading, and has a final status
    if (currentAssistantMessage && !isLoading && (currentAssistantMessage.status === 'success' || currentAssistantMessage.status === 'error')) {
      onSendMessage(conversation.id, currentAssistantMessage);
      setCurrentAssistantMessage(null); // Clear current message
    }
  }, [isLoading, currentAssistantMessage, conversation, onSendMessage]);

  const adjustTextareaHeight = () => {
    const textarea = textareaRef.current;
    if (textarea) {
      textarea.style.height = 'auto';
      textarea.style.height = Math.min(textarea.scrollHeight, 120) + 'px';
    }
  };

  useEffect(() => {
    adjustTextareaHeight();
  }, [inputValue]);

  const formatTime = (timestamp) => {
    return new Date(timestamp).toLocaleTimeString('zh-CN', {
      hour: '2-digit',
      minute: '2-digit'
    });
  };

  if (!conversation) {
    console.log('ChatWindow: No conversation provided, showing welcome screen');
    return (
      <Box sx={{
        flexGrow: 1,
        display: 'flex',
        flexDirection: 'column',
        justifyContent: 'center',
        alignItems: 'center',
        p: 3,
        textAlign: 'center',
        color: 'text.secondary',
      }}>
        <SmartToyIcon sx={{ fontSize: 80, mb: 2 }} />
        <Typography variant="h5" component="h3" gutterBottom>
          欢迎使用 BioaGen Agent
        </Typography>
        <Typography variant="body1" sx={{ maxWidth: 500 }}>
          这是一个基于AI的生物分子分析助手。您可以询问有关分子结构、性质预测、
          化学反应等相关问题。点击左侧"新建对话"开始您的分析之旅。
        </Typography>
      </Box>
    );
  }

  return (
    <Box sx={{
      flexGrow: 1,
      display: 'flex',
      flexDirection: 'column',
      height: '100%',
      overflow: 'hidden',
    }}>
      <Box sx={{ p: 2, borderBottom: '1px solid', borderColor: 'divider' }}>
        <Typography variant="h6">{conversation.title}</Typography>
      </Box>

      {conversation.available_tools && conversation.available_tools.length > 0 && (
        <Box sx={{ p: 2, bgcolor: 'action.hover', borderBottom: '1px solid', borderColor: 'divider' }}>
          <Typography variant="subtitle2" gutterBottom>Available Tools:</Typography>
          <Box component="ul" sx={{ listStyle: 'none', p: 0, m: 0 }}>
            {conversation.available_tools.map((tool, index) => (
              <li key={index}>
                <Chip label={`${tool.function.name}: ${tool.function.description}`} size="small" sx={{ mb: 0.5 }} />
              </li>
            ))}
          </Box>
        </Box>
      )}

      <Box sx={{ flexGrow: 1, overflowY: 'auto', p: 2 }}>
        {conversation.messages.map((message) => (
          <Box key={message.id} sx={{
            display: 'flex',
            mb: 2,
            justifyContent: message.type === 'user' ? 'flex-end' : 'flex-start',
          }}>
            {message.type === 'assistant' && (
              <Avatar sx={{ bgcolor: 'primary.main', mr: 1 }}><SmartToyIcon fontSize="small" /></Avatar>
            )}
            <Paper elevation={1} sx={{
              p: 1,
              borderRadius: 1,
              maxWidth: '85%', // 增加宽度从70%到85%
              bgcolor: message.type === 'user' ? 'primary.light' : 'background.paper',
              color: message.type === 'user' ? 'primary.contrastText' : 'text.primary',
            }}>
              <AgentStepsDisplay agentSteps={message.agent_steps} defaultExpanded={true} />
              <Typography component="div">
                <EnhancedMarkdown>{message.content}</EnhancedMarkdown>
              </Typography>

              {/* Display status indicator */}
              {message.status && message.status !== 'success' && (
                <Chip
                  icon={message.status === 'error' ? <ErrorOutlineIcon /> : <InfoOutlinedIcon />}
                  label={message.status === 'error' ? '执行失败' : '部分成功'}
                  size="small"
                  color={message.status === 'error' ? 'error' : 'warning'}
                  sx={{ mt: 1 }}
                />
              )}

              <Typography variant="caption" color="text.secondary" sx={{ display: 'block', mt: 0.5 }}>
                {formatTime(message.timestamp)}
              </Typography>
            </Paper>
            {message.type === 'user' && (
              <Avatar sx={{ bgcolor: 'secondary.main', ml: 1 }}><PersonIcon fontSize="small" /></Avatar>
            )}
          </Box>
        ))}

        {/* Display the current streaming message */}
        {currentAssistantMessage && (
          <Box sx={{
            display: 'flex',
            mb: 2,
            justifyContent: 'flex-start',
          }}>
            <Avatar sx={{ bgcolor: 'primary.main', mr: 1 }}><SmartToyIcon fontSize="small" /></Avatar>
            <Paper elevation={1} sx={{
              p: 1,
              borderRadius: 1,
              maxWidth: '85%', // 增加宽度从70%到85%
              bgcolor: 'background.paper',
              color: 'text.primary',
            }}>
              {currentAssistantMessage.agent_steps && currentAssistantMessage.agent_steps.length > 0 && (
                <Box sx={{ mb: 1 }}>
                  {groupAgentSteps(currentAssistantMessage.agent_steps).map((stepGroup, groupIndex) => (
                    <Box key={groupIndex} sx={{ mb: 1, p: 1, bgcolor: 'grey.50', borderRadius: 1 }}>
                      <Typography variant="subtitle2">步骤 {groupIndex + 1}</Typography>
                      {stepGroup.thought && (
                        <Typography variant="body2" sx={{ mt: 0.5 }}>
                          <strong>Agent Thought:</strong> <EnhancedMarkdown>{stepGroup.thought}</EnhancedMarkdown>
                        </Typography>
                      )}
                      {stepGroup.tool_calls && stepGroup.tool_calls.length > 0 && (
                        <Box sx={{ mt: 0.5 }}>
                          <Typography variant="body2"><strong>Tool Calls:</strong></Typography>
                          <ToolCallDisplay toolCalls={stepGroup.tool_calls} />
                        </Box>
                      )}
                      {stepGroup.tool_outputs && stepGroup.tool_outputs.length > 0 && (
                        <Box sx={{ mt: 0.5 }}>
                          <Typography variant="body2"><strong>Tool Outputs:</strong></Typography>
                          {stepGroup.tool_outputs.map((output, outputIndex) => (
                            <Box key={outputIndex} sx={{ mt: 0.5, pl: 1, borderLeft: '2px solid', borderColor: 'grey.300' }}>
                              <Typography variant="caption" color="text.secondary">Output {outputIndex + 1}:</Typography>
                              <ToolOutputContent output={JSON.stringify(output)} />
                            </Box>
                          ))}
                        </Box>
                      )}
                    </Box>
                  ))}
                </Box>
              )}
              <Typography component="div">
                <EnhancedMarkdown>{currentAssistantMessage.content}</EnhancedMarkdown>
              </Typography>
              {isLoading && (
                <CircularProgress size={20} sx={{ mt: 1 }} />
              )}
            </Paper>
          </Box>
        )}

        {isLoading && !currentAssistantMessage && (
          <Box sx={{
            display: 'flex',
            mb: 2,
            justifyContent: 'flex-start',
          }}>
            <Avatar sx={{ bgcolor: 'primary.main', mr: 1 }}><SmartToyIcon fontSize="small" /></Avatar>
            <Paper elevation={1} sx={{
              p: 1,
              borderRadius: 1,
              maxWidth: '85%', // 增加宽度从70%到85%
              bgcolor: 'background.paper',
              color: 'text.primary',
            }}>
              <CircularProgress size={20} />
              <Typography variant="body2" sx={{ mt: 1 }}>正在思考中</Typography>
            </Paper>
          </Box>
        )}
        <div ref={messagesEndRef} />
      </Box>

      <Box sx={{ p: 2, borderTop: '1px solid', borderColor: 'divider', display: 'flex', alignItems: 'flex-end', gap: 1 }}>
        <TextField
          fullWidth
          multiline
          maxRows={4}
          value={inputValue}
          onChange={(e) => setInputValue(e.target.value)}
          onKeyPress={handleKeyPress}
          placeholder="输入您的问题...例如：分析乙醇分子的性质"
          variant="outlined"
          size="small"
          disabled={isLoading}
          inputRef={textareaRef}
          sx={{ flexGrow: 1 }}
        />
        <Button
          variant="contained"
          endIcon={<SendIcon />}
          onClick={handleSubmit}
          disabled={!inputValue.trim() || isLoading}
          sx={{ height: '40px' }} // Adjust height to match TextField
        >
          发送
        </Button>
      </Box>
    </Box>
  );
};

export default ChatWindow;
