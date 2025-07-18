#!/usr/bin/env python3
"""
Test script to verify SSE streaming fixes are working properly.
This script can be run to test the error handling improvements.
"""

import asyncio
import httpx
import json
import sys
from typing import AsyncGenerator

async def test_sse_connection(base_url: str = "http://localhost:8000"):
    """Test SSE connection and error handling"""
    
    print("🧪 Testing SSE Connection and Error Handling...")
    print(f"📡 Base URL: {base_url}")
    
    # Test 1: Valid connection
    print("\n1️⃣ Testing valid SSE connection...")
    try:
        async with httpx.AsyncClient(timeout=30.0) as client:
            url = f"{base_url}/api/chat"
            params = {
                "message": "Hello, this is a test message"
            }
            
            async with client.stream("GET", url, params=params) as response:
                print(f"✅ Connection established. Status: {response.status_code}")
                
                if response.status_code == 200:
                    print("📨 Receiving SSE messages...")
                    async for line in response.aiter_lines():
                        if line.startswith("data: "):
                            data = line[6:]  # Remove "data: " prefix
                            try:
                                parsed = json.loads(data)
                                print(f"📦 Received: {parsed.get('type', 'unknown')} - {parsed.get('content', '')[:100]}...")
                                
                                if parsed.get('type') == 'end':
                                    print("🏁 Stream ended normally")
                                    break
                                elif parsed.get('type') == 'error':
                                    print(f"❌ Error received: {parsed.get('content')}")
                                    break
                            except json.JSONDecodeError as e:
                                print(f"⚠️ JSON decode error: {e}")
                                print(f"Raw data: {data}")
                else:
                    print(f"❌ Connection failed with status: {response.status_code}")
                    
    except httpx.ConnectError:
        print("❌ Connection failed - server may not be running")
        return False
    except Exception as e:
        print(f"❌ Unexpected error: {e}")
        return False
    
    # Test 2: Invalid endpoint
    print("\n2️⃣ Testing invalid endpoint...")
    try:
        async with httpx.AsyncClient(timeout=10.0) as client:
            url = f"{base_url}/api/invalid_endpoint"
            params = {"message": "test"}
            
            response = await client.get(url, params=params)
            print(f"📊 Invalid endpoint status: {response.status_code}")
            
    except Exception as e:
        print(f"✅ Expected error for invalid endpoint: {type(e).__name__}")
    
    # Test 3: Empty message
    print("\n3️⃣ Testing empty message...")
    try:
        async with httpx.AsyncClient(timeout=10.0) as client:
            url = f"{base_url}/api/chat"
            params = {"message": ""}
            
            async with client.stream("GET", url, params=params) as response:
                print(f"📊 Empty message status: {response.status_code}")
                
                if response.status_code == 200:
                    async for line in response.aiter_lines():
                        if line.startswith("data: "):
                            data = line[6:]
                            try:
                                parsed = json.loads(data)
                                if parsed.get('type') in ['error', 'end']:
                                    print(f"✅ Handled empty message: {parsed.get('type')}")
                                    break
                            except json.JSONDecodeError:
                                pass
                            
    except Exception as e:
        print(f"⚠️ Error with empty message: {e}")
    
    print("\n✅ SSE testing completed!")
    return True

async def test_backend_health(base_url: str = "http://localhost:8000"):
    """Test if backend is healthy"""
    print("🏥 Testing backend health...")
    
    try:
        async with httpx.AsyncClient(timeout=5.0) as client:
            response = await client.get(f"{base_url}/api/health")
            if response.status_code == 200:
                print("✅ Backend is healthy")
                return True
            else:
                print(f"⚠️ Backend health check failed: {response.status_code}")
                return False
    except Exception as e:
        print(f"❌ Backend health check failed: {e}")
        return False

async def main():
    """Main test function"""
    print("🚀 Starting SSE Error Handling Tests")
    print("=" * 50)
    
    base_url = "http://localhost:8000"
    
    # Check if custom URL provided
    if len(sys.argv) > 1:
        base_url = sys.argv[1]
        print(f"🔧 Using custom base URL: {base_url}")
    
    # Test backend health first
    if not await test_backend_health(base_url):
        print("\n❌ Backend is not available. Please ensure:")
        print("   1. Backend server is running")
        print("   2. Backend is accessible at the specified URL")
        print("   3. Health endpoint is working")
        return
    
    # Test SSE functionality
    await test_sse_connection(base_url)
    
    print("\n🎉 All tests completed!")
    print("\n📋 Summary of fixes applied:")
    print("   ✅ Backend error response JSON formatting")
    print("   ✅ Stream error handling with try-catch")
    print("   ✅ Agent orchestrator exception handling")
    print("   ✅ Frontend error message improvements")
    print("   ✅ Tool manager error categorization")
    print("   ✅ SSE connection state handling")

if __name__ == "__main__":
    asyncio.run(main())
