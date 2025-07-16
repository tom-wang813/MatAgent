import React from 'react';
import { Box, Typography, Button, List, ListItem, ListItemText, IconButton, Paper } from '@mui/material';
import AddIcon from '@mui/icons-material/Add';
import DeleteIcon from '@mui/icons-material/Delete';

const ChatSidebar = ({
  conversations,
  activeConversationId,
  onSelectConversation,
  onNewConversation,
  onDeleteConversation
}) => {
  const formatDate = (dateString) => {
    const date = new Date(dateString);
    const now = new Date();
    const diff = now - date;
    const days = Math.floor(diff / (1000 * 60 * 60 * 24));
    
    if (days === 0) {
      return '今天';
    } else if (days === 1) {
      return '昨天';
    } else if (days < 7) {
      return `${days}天前`;
    } else {
      return date.toLocaleDateString('zh-CN');
    }
  };

  return (
    <Box sx={{
      width: 280,
      bgcolor: 'background.paper',
      borderRight: '1px solid', 
      borderColor: 'divider',
      display: 'flex',
      flexDirection: 'column',
      height: '100%',
    }}>
      <Box sx={{ p: 2, borderBottom: '1px solid', borderColor: 'divider', display: 'flex', alignItems: 'center', justifyContent: 'space-between' }}>
        <Typography variant="h6" sx={{ fontWeight: 'bold' }}>🧬 BioaGen Agent</Typography>
        <Button
          variant="contained"
          startIcon={<AddIcon />}
          onClick={onNewConversation}
          size="small"
        >
          新建对话
        </Button>
      </Box>
      
      <List sx={{ flexGrow: 1, overflowY: 'auto' }}>
        {conversations.length === 0 ? (
          <Typography variant="body2" color="text.secondary" sx={{ p: 2, textAlign: 'center' }}>
            暂无对话历史
          </Typography>
        ) : (
          conversations.map((conversation) => (
            <ListItem
              key={conversation.id}
              button
              selected={conversation.id === activeConversationId}
              onClick={() => onSelectConversation(conversation.id)}
              sx={{
                borderBottom: '1px solid', 
                borderColor: 'divider',
                '&.Mui-selected': {
                  bgcolor: 'action.selected',
                },
                '&:hover': {
                  bgcolor: 'action.hover',
                },
              }}
            >
              <ListItemText
                primary={<Typography variant="subtitle1" noWrap>{conversation.title}</Typography>}
                secondary={<Typography variant="body2" color="text.secondary">{formatDate(conversation.updatedAt)}</Typography>}
              />
              <IconButton
                edge="end"
                aria-label="delete"
                onClick={(e) => {
                  e.stopPropagation();
                  onDeleteConversation(conversation.id);
                }}
                size="small"
              >
                <DeleteIcon fontSize="small" />
              </IconButton>
            </ListItem>
          ))
        )}
      </List>
    </Box>
  );
};

export default ChatSidebar;
