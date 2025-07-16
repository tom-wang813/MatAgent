from sqlalchemy.orm import Session
from typing import List, Optional, Dict
from datetime import datetime
import uuid

from backend.core.models import Conversation, Message
from shared_utils.logging_config import get_logger

logger = get_logger(__name__)

class ConversationService:
    def __init__(self, db: Session):
        self.db = db

    def create_conversation(self, title: Optional[str] = None, trace_id: str = "N/A", created_at: Optional[datetime] = None) -> Conversation:
        """
        Creates a new conversation in the database.
        """
        new_uuid = str(uuid.uuid4())
        conversation = Conversation(uuid=new_uuid, title=title if title else "New Chat", created_at=created_at)
        self.db.add(conversation)
        self.db.commit()
        self.db.refresh(conversation)
        logger.info(f"Created new conversation: {conversation.uuid}", extra={'trace_id': trace_id, 'event': 'conversation_created', 'conversation_id': conversation.id, 'uuid': conversation.uuid})
        return conversation

    def get_conversation(self, conversation_id: int, trace_id: str = "N/A") -> Optional[Conversation]:
        """
        Retrieves a conversation by its ID.
        """
        conversation = self.db.query(Conversation).filter(Conversation.id == conversation_id).first()
        if conversation:
            logger.debug(f"Retrieved conversation: {conversation.uuid}", extra={'trace_id': trace_id, 'event': 'conversation_retrieved', 'conversation_id': conversation.id})
        return conversation

    def get_conversation_by_uuid(self, conversation_uuid: str, trace_id: str = "N/A") -> Optional[Conversation]:
        """
        Retrieves a conversation by its UUID.
        """
        conversation = self.db.query(Conversation).filter(Conversation.uuid == conversation_uuid).first()
        if conversation:
            logger.debug(f"Retrieved conversation by UUID: {conversation.uuid}", extra={'trace_id': trace_id, 'event': 'conversation_retrieved_by_uuid', 'conversation_id': conversation.id, 'uuid': conversation.uuid})
        return conversation

    def add_message(self, conversation_id: int, role: str, content: str, tool_call_id: Optional[str] = None, trace_id: str = "N/A") -> Message:
        """
        Adds a new message to an existing conversation.
        """
        message = Message(conversation_id=conversation_id, role=role, content=content, tool_call_id=tool_call_id)
        self.db.add(message)
        self.db.commit()
        self.db.refresh(message)
        logger.info(f"Added message to conversation {conversation_id} (role: {role})", extra={'trace_id': trace_id, 'event': 'message_added', 'conversation_id': conversation_id, 'role': role})
        return message

    def get_messages(self, conversation_id: int, trace_id: str = "N/A") -> List[Message]:
        """
        Retrieves all messages for a given conversation, ordered by creation time.
        """
        messages = self.db.query(Message).filter(Message.conversation_id == conversation_id).order_by(Message.created_at).all()
        logger.debug(f"Retrieved {len(messages)} messages for conversation {conversation_id}", extra={'trace_id': trace_id, 'event': 'messages_retrieved', 'conversation_id': conversation_id, 'count': len(messages)})
        return messages

    def update_conversation_title(self, conversation_id: int, new_title: str, trace_id: str = "N/A") -> Optional[Conversation]:
        """
        Updates the title of a conversation.
        """
        conversation = self.get_conversation(conversation_id, trace_id=trace_id)
        if conversation:
            conversation.title = new_title
            self.db.commit()
            self.db.refresh(conversation)
            logger.info(f"Updated title for conversation {conversation_id} to '{new_title}'", extra={'trace_id': trace_id, 'event': 'conversation_title_updated', 'conversation_id': conversation.id, 'new_title': new_title})
        return conversation

    def delete_conversation(self, conversation_id: int, trace_id: str = "N/A"):
        """
        Deletes a conversation and all its associated messages.
        """
        conversation = self.get_conversation(conversation_id, trace_id=trace_id)
        if conversation:
            self.db.delete(conversation)
            self.db.commit()
            logger.info(f"Deleted conversation {conversation_id}", extra={'trace_id': trace_id, 'event': 'conversation_deleted', 'conversation_id': conversation.id})

    def get_all_conversations(self) -> List[Conversation]:
        """
        Retrieves all conversations.
        """
        conversations = self.db.query(Conversation).order_by(Conversation.created_at.desc()).all()
        logger.debug(f"Retrieved {len(conversations)} conversations.", extra={'event': 'all_conversations_retrieved', 'count': len(conversations)})
        return conversations
