import React from 'react';
import { render, screen, fireEvent } from '@testing-library/react';
import '@testing-library/jest-dom';
import ChatSidebar from '../components/ChatSidebar';

describe('ChatSidebar', () => {
  const mockConversations = [
    { id: '1', title: 'Conversation 1', updatedAt: new Date().toISOString() },
    { id: '2', title: 'Conversation 2', updatedAt: new Date().toISOString() },
  ];

  const mockOnSelectConversation = jest.fn();
  const mockOnNewConversation = jest.fn();
  const mockOnDeleteConversation = jest.fn();

  beforeEach(() => {
    jest.clearAllMocks();
  });

  test('renders sidebar title', () => {
    render(
      <ChatSidebar
        conversations={[]}
        activeConversationId={null}
        onSelectConversation={mockOnSelectConversation}
        onNewConversation={mockOnNewConversation}
        onDeleteConversation={mockOnDeleteConversation}
      />
    );
    expect(screen.getByText('🧬 BioaGen Agent')).toBeInTheDocument();
  });

  test('renders "新建对话" button', () => {
    render(
      <ChatSidebar
        conversations={[]}
        activeConversationId={null}
        onSelectConversation={mockOnSelectConversation}
        onNewConversation={mockOnNewConversation}
        onDeleteConversation={mockOnDeleteConversation}
      />
    );
    expect(screen.getByRole('button', { name: /新建对话/i })).toBeInTheDocument();
  });

  test('calls onNewConversation when "新建对话" button is clicked', () => {
    render(
      <ChatSidebar
        conversations={[]}
        activeConversationId={null}
        onSelectConversation={mockOnSelectConversation}
        onNewConversation={mockOnNewConversation}
        onDeleteConversation={mockOnDeleteConversation}
      />
    );
    fireEvent.click(screen.getByRole('button', { name: /新建对话/i }));
    expect(mockOnNewConversation).toHaveBeenCalledTimes(1);
  });

  test('renders conversation list when conversations are provided', () => {
    render(
      <ChatSidebar
        conversations={mockConversations}
        activeConversationId={null}
        onSelectConversation={mockOnSelectConversation}
        onNewConversation={mockOnNewConversation}
        onDeleteConversation={mockOnDeleteConversation}
      />
    );
    expect(screen.getByText('Conversation 1')).toBeInTheDocument();
    expect(screen.getByText('Conversation 2')).toBeInTheDocument();
  });

  test('displays "暂无对话历史" when no conversations', () => {
    render(
      <ChatSidebar
        conversations={[]}
        activeConversationId={null}
        onSelectConversation={mockOnSelectConversation}
        onNewConversation={mockOnNewConversation}
        onDeleteConversation={mockOnDeleteConversation}
      />
    );
    expect(screen.getByText('暂无对话历史')).toBeInTheDocument();
  });

  test('calls onSelectConversation when a conversation item is clicked', () => {
    render(
      <ChatSidebar
        conversations={mockConversations}
        activeConversationId={null}
        onSelectConversation={mockOnSelectConversation}
        onNewConversation={mockOnNewConversation}
        onDeleteConversation={mockOnDeleteConversation}
      />
    );
    fireEvent.click(screen.getByText('Conversation 1'));
    expect(mockOnSelectConversation).toHaveBeenCalledTimes(1);
    expect(mockOnSelectConversation).toHaveBeenCalledWith('1');
  });

  test('applies "active" class to the active conversation', () => {
    render(
      <ChatSidebar
        conversations={mockConversations}
        activeConversationId="1"
        onSelectConversation={mockOnSelectConversation}
        onNewConversation={mockOnNewConversation}
        onDeleteConversation={mockOnDeleteConversation}
      />
    );
    expect(screen.getByText('Conversation 1').closest('.conversation-item')).toHaveClass('active');
    expect(screen.getByText('Conversation 2').closest('.conversation-item')).not.toHaveClass('active');
  });

  test('calls onDeleteConversation when delete button is clicked', () => {
    render(
      <ChatSidebar
        conversations={mockConversations}
        activeConversationId={null}
        onSelectConversation={mockOnSelectConversation}
        onNewConversation={mockOnNewConversation}
        onDeleteConversation={mockOnDeleteConversation}
      />
    );
    fireEvent.click(screen.getAllByTitle('删除对话')[0]); // Click the first delete button
    expect(mockOnDeleteConversation).toHaveBeenCalledTimes(1);
    expect(mockOnDeleteConversation).toHaveBeenCalledWith('1');
  });
});
