import React, { useState, useEffect } from 'react';
import ChatSidebar from './components/ChatSidebar';
import ChatWindow from './components/ChatWindow';
import { Box } from '@mui/material';
// import './App.css'; // Removed as MUI will handle styling

function App() {
  const [conversations, setConversations] = useState([]);
  const [activeConversationId, setActiveConversationId] = useState(null);
  const [isLoading, setIsLoading] = useState(false);

  // Load conversations from localStorage on mount
  useEffect(() => {
    let parsedConversations = [];
    try {
      const savedConversations = localStorage.getItem('bioagen_conversations');
      if (savedConversations) {
        parsedConversations = JSON.parse(savedConversations);
      }
    } catch (error) {
      console.error("Error parsing conversations from localStorage, clearing it.", error);
      localStorage.removeItem('bioagen_conversations');
      parsedConversations = [];
    }

    if (parsedConversations.length > 0) {
      setConversations(parsedConversations);
      setActiveConversationId(parsedConversations[0].id);
    } else {
      // If no conversations found or parsing failed, create a new one automatically
      createNewConversation();
    }
  }, []);

  // Save conversations to localStorage whenever they change
  useEffect(() => {
    localStorage.setItem('bioagen_conversations', JSON.stringify(conversations));
  }, [conversations]);

  const createNewConversation = () => {
    const newConversation = {
      id: Date.now().toString(),
      title: '新对话',
      messages: [],
      createdAt: new Date().toISOString(),
      updatedAt: new Date().toISOString()
    };
    
    setConversations(prev => [newConversation, ...prev]);
    setActiveConversationId(newConversation.id);
  };

  const updateConversationTitle = (conversationId, title) => {
    setConversations(prev => prev.map(conv => 
      conv.id === conversationId 
        ? { ...conv, title, updatedAt: new Date().toISOString() }
        : conv
    ));
  };

  const addMessage = (conversationId, message) => {
    setConversations(prev => prev.map(conv => 
      conv.id === conversationId 
        ? { 
            ...conv, 
            messages: [...conv.messages],
            updatedAt: new Date().toISOString()
          }
        : conv
    ));
  };

  const deleteConversation = (conversationId) => {
    setConversations(prev => {
      const filtered = prev.filter(conv => conv.id !== conversationId);
      if (activeConversationId === conversationId) {
        setActiveConversationId(filtered.length > 0 ? filtered[0].id : null);
      }
      return filtered;
    });
  };

  const updateConversationUuid = (conversationId, uuid) => {
    setConversations(prev => prev.map(conv => 
      conv.id === conversationId 
        ? { ...conv, uuid, updatedAt: new Date().toISOString() }
        : conv
    ));
  };

  const activeConversation = conversations.find(conv => conv.id === activeConversationId);

  return (
    <Box sx={{
      display: 'flex',
      height: '100vh',
      width: '100vw',
      overflow: 'hidden',
    }}>
      <ChatSidebar
        conversations={conversations}
        activeConversationId={activeConversationId}
        onSelectConversation={setActiveConversationId}
        onNewConversation={createNewConversation}
        onDeleteConversation={deleteConversation}
      />
      <ChatWindow
        conversation={activeConversation}
        onSendMessage={addMessage}
        onUpdateTitle={updateConversationTitle}
        isLoading={isLoading}
        setIsLoading={setIsLoading}
        onUpdateConversationUuid={updateConversationUuid}
      />
    </Box>
  );
}

export default App;
