import httpx
import json
import time
import pytest
import asyncio
import pytest_asyncio

# Fixture for httpx.AsyncClient
@pytest_asyncio.fixture(scope="module")
async def async_client():
    async with httpx.AsyncClient() as client:
        yield client

# Helper function to wait for backend health check
async def wait_for_backend_health(client: httpx.AsyncClient, max_retries=30, delay=2):
    health_url = "http://localhost:8000/api/health"
    for i in range(max_retries):
        try:
            print(f"Attempt {i+1}/{max_retries}: Waiting for backend health check...")
            response = await client.get(health_url, timeout=5)
            if response.status_code == 200 and response.json().get("status") == "ok":
                print("Backend is healthy!")
                return True
        except httpx.RequestError as e:
            print(f"Backend not ready yet: {e}")
        await asyncio.sleep(delay)
    print("Backend did not become healthy within the expected time.")
    return False

# Helper function to send a chat message and extract conversation_uuid if present
async def send_chat_message(client: httpx.AsyncClient, message, conversation_uuid=None):
    url = "http://localhost:8000/api/chat"
    payload = {"message": message}
    if conversation_uuid:
        payload["conversation_uuid"] = conversation_uuid
    
    print(f"Sending message: {message} with conversation_uuid: {conversation_uuid}")
    response = await client.post(url, json=payload, timeout=120) # Increased timeout
    response.raise_for_status()
    
    response_text = response.text
    print(f"HTTP Response Status: {response.status_code}")
    print(f"HTTP Response Content (first 500 chars):\n{response_text[:500]}...")

    extracted_uuid = None
    for line in response_text.splitlines():
        try:
            if line.strip(): # Ensure line is not empty
                data = json.loads(line)
                if data.get("type") == "conversation_info" and "conversation_uuid" in data:
                    extracted_uuid = data["conversation_uuid"]
                    print(f"Extracted conversation_uuid from stream: {extracted_uuid}")
                    break
        except json.JSONDecodeError:
            continue # Not a valid JSON line, skip

    await asyncio.sleep(2) # Add a small delay after sending message
    return extracted_uuid, response_text

@pytest.mark.asyncio
async def test_chat_flow(async_client):
    # Step 0: Wait for backend to be healthy
    assert await wait_for_backend_health(async_client), "Backend did not become healthy."

    # Step 1: Send the first message and get conversation_uuid from response
    first_message_text = "分析一下乙醇的性质"
    conversation_uuid, _ = await send_chat_message(async_client, first_message_text)
    
    assert conversation_uuid is not None, "Failed to retrieve conversation_uuid from HTTP response."
    print(f"Retrieved conversation_uuid: {conversation_uuid}")

    # Step 2: Send the second message with the conversation_uuid
    second_message_text = "然后和葡萄糖进行详细的对比"
    _, response_text = await send_chat_message(async_client, second_message_text, conversation_uuid)
    assert "Final Answer" in response_text or "tool_calls" in response_text, "Second message did not result in expected response."
    print(f"Second message response: {response_text}")
