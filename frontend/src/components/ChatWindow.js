import React, { useState, useRef, useEffect } from 'react';
import { Box, TextField, Button, Typography, CircularProgress, Paper, Avatar, Chip } from '@mui/material';
import SendIcon from '@mui/icons-material/Send';
import SmartToyIcon from '@mui/icons-material/SmartToy';
import PersonIcon from '@mui/icons-material/Person';
import ErrorOutlineIcon from '@mui/icons-material/ErrorOutline';
import InfoOutlinedIcon from '@mui/icons-material/InfoOutlined';
import ExpandMoreIcon from '@mui/icons-material/ExpandMore';
import ExpandLessIcon from '@mui/icons-material/ExpandLess';
import ReactMarkdown from 'react-markdown';
import remarkMath from 'remark-math';
import rehypeKatex from 'rehype-katex';
import 'katex/dist/katex.min.css'; // Import KaTeX CSS
import ToolCallDisplay from './ToolCallDisplay';
import { sendMessage } from '../utils/api'; // Import sendMessage from api.js

const MAX_TOOL_OUTPUT_LENGTH = 200; // Define max length for tool output display

// Helper function to group agent steps into logical turns
const groupAgentSteps = (steps) => {
  const grouped = [];
  let currentGroup = null;

  steps.forEach(step => {
    if (step.type === 'thought') {
      if (currentGroup && (currentGroup.tool_calls.length > 0 || currentGroup.tool_output)) {
        grouped.push(currentGroup);
      }
      currentGroup = { thought: step.content, tool_calls: [], tool_output: null };
    } else if (step.type === 'tool_calls') {
      if (!currentGroup) {
        currentGroup = { thought: null, tool_calls: [], tool_output: null };
      }
      currentGroup.tool_calls = currentGroup.tool_calls.concat(step.content); // step.content is now the tool_calls array
    } else if (step.type === 'tool_output') {
      if (!currentGroup) {
        currentGroup = { thought: null, tool_calls: [], tool_output: null };
      }
      currentGroup.tool_output = step.content.output; // step.content is now the tool_output object
      grouped.push(currentGroup); // Tool output usually marks the end of a step
      currentGroup = null; // Reset for next group
    } else if (step.type === 'final_answer') {
      if (currentGroup) {
        grouped.push(currentGroup);
      }
      // Final answer is handled separately, not grouped as a step
    }
  });

  if (currentGroup) {
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
        <ReactMarkdown>{displayedOutput}</ReactMarkdown>
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

  // Effect to add the completed assistant message to conversation messages
  useEffect(() => {
    if (currentAssistantMessage && !isLoading && conversation) {
      onSendMessage(conversation.id, currentAssistantMessage);
      setCurrentAssistantMessage(null); // Clear current message
    }
  }, [isLoading, currentAssistantMessage, conversation, onSendMessage]);

  const handleKeyPress = (e) => {
    if (e.key === 'Enter' && !e.shiftKey) {
      e.preventDefault();
      handleSubmit(e);
    }
  };

  const handleSubmit = async (e) => {
    e.preventDefault();
    if (!inputValue.trim() || isLoading || !conversation) return;

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
          setCurrentAssistantMessage(prev => ({
            ...prev,
            content: `抱歉，发生了错误：${error.message || '未知错误'}。请稍后再试。`,
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
      setCurrentAssistantMessage(prev => ({
        ...prev,
        content: `抱歉，无法启动对话：${error.message || '未知错误'}。请稍后再试。`,
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
              maxWidth: '70%',
              bgcolor: message.type === 'user' ? 'primary.light' : 'background.paper',
              color: message.type === 'user' ? 'primary.contrastText' : 'text.primary',
            }}>
              {message.agent_steps && message.agent_steps.length > 0 && (
                <Box sx={{ mb: 1 }}>
                  {groupAgentSteps(message.agent_steps).map((stepGroup, groupIndex) => (
                    <Box key={groupIndex} sx={{ mb: 1, p: 1, bgcolor: 'grey.50', borderRadius: 1 }}>
                      <Typography variant="subtitle2">步骤 {groupIndex + 1}</Typography>
                      {stepGroup.thought && (
                        <Typography variant="body2" sx={{ mt: 0.5 }}>
                          <strong>Agent Thought:</strong> <ReactMarkdown>{stepGroup.thought}</ReactMarkdown>
                        </Typography>
                      )}
                      {stepGroup.tool_calls && stepGroup.tool_calls.length > 0 && (
                        <Box sx={{ mt: 0.5 }}>
                          <Typography variant="body2"><strong>Tool Calls:</strong></Typography>
                          <ToolCallDisplay toolCalls={stepGroup.tool_calls} />
                        </Box>
                      )}
                      {stepGroup.tool_output && (
                        <Box sx={{ mt: 0.5 }}>
                          <Typography variant="body2"><strong>Tool Output:</strong></Typography>
                          <ToolOutputContent output={JSON.stringify(stepGroup.tool_output)} />
                        </Box>
                      )}
                    </Box>
                  ))}
                </Box>
              )}
              <Typography component="div">
                <ReactMarkdown
                  remarkPlugins={[remarkMath]}
                  rehypePlugins={[rehypeKatex]}
                >{message.content}</ReactMarkdown>
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
              maxWidth: '70%',
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
                          <strong>Agent Thought:</strong> <ReactMarkdown>{stepGroup.thought}</ReactMarkdown>
                        </Typography>
                      )}
                      {stepGroup.tool_calls && stepGroup.tool_calls.length > 0 && (
                        <Box sx={{ mt: 0.5 }}>
                          <Typography variant="body2"><strong>Tool Calls:</strong></Typography>
                          <ToolCallDisplay toolCalls={stepGroup.tool_calls} />
                        </Box>
                      )}
                      {stepGroup.tool_output && (
                        <Box sx={{ mt: 0.5 }}>
                          <Typography variant="body2"><strong>Tool Output:</strong></Typography>
                          <ToolOutputContent output={JSON.stringify(stepGroup.tool_output)} />
                        </Box>
                      )}
                    </Box>
                  ))}
                </Box>
              )}
              <Typography component="div">
                <ReactMarkdown
                  remarkPlugins={[remarkMath]}
                  rehypePlugins={[rehypeKatex]}
                >{currentAssistantMessage.content}</ReactMarkdown>
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
              maxWidth: '70%',
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