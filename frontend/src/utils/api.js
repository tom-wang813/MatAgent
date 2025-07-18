export const sendMessage = (message, conversationId, onNewStep, onError, onOpen) => {
  const backendUrl = process.env.REACT_APP_BACKEND_URL || 'http://localhost:8000/api';

  // Construct URL with query parameters for GET request with EventSource
  // Handle relative URLs properly
  let fullUrl;
  if (backendUrl.startsWith('http')) {
    fullUrl = `${backendUrl}/chat`;
  } else {
    // For relative URLs like '/api', construct full URL based on current location
    fullUrl = `${window.location.protocol}//${window.location.host}${backendUrl}/chat`;
  }
  
  let url;
  try {
    url = new URL(fullUrl);
  } catch (error) {
    console.error('Failed to construct URL from:', fullUrl, 'Error:', error);
    if (onError) onError(new Error(`Invalid URL: ${fullUrl}`));
    return;
  }
  url.searchParams.append('message', message);
  if (conversationId) {
    url.searchParams.append('conversation_uuid', conversationId);
  }
  // EventSource does not support custom headers directly for GET requests.
  // If X-Trace-ID is critical, it needs to be passed as a query parameter or handled differently.
  // For now, we'll omit it as it's primarily for backend logging.

  const eventSource = new EventSource(url.toString());

  eventSource.onopen = () => {
    console.log('SSE connection opened.');
    if (onOpen) onOpen();
  };

  eventSource.onmessage = (event) => {
    console.log('SSE message received:', event.data); // Debug log
    if (event.data) {
      try {
        const step = JSON.parse(event.data);
        console.log('Parsed SSE step:', step); // Debug log
        onNewStep(step);
      } catch (e) {
        console.error('Error parsing SSE message:', e, 'Data:', event.data);
        // Don't trigger error for parsing issues, just log them
        // if (onError) onError(e);
      }
    }
  };

  eventSource.onerror = (err) => {
    console.error('SSE error:', err);
    eventSource.close();
    
    // Create a more descriptive error based on the readyState
    let errorMessage = 'SSE連接錯誤，請檢查網絡連接';
    
    if (eventSource.readyState === EventSource.CONNECTING) {
      errorMessage = '正在嘗試連接服務器...';
    } else if (eventSource.readyState === EventSource.CLOSED) {
      errorMessage = '與服務器的連接已關閉，請重試';
    } else {
      // Check if it's a network error or server error
      if (err.target && err.target.readyState === EventSource.CLOSED) {
        errorMessage = '服務器連接中斷，請檢查服務器狀態';
      } else {
        errorMessage = '網絡連接錯誤，請檢查網絡設置';
      }
    }
    
    const error = new Error(errorMessage);
    if (onError) onError(error);
  };

  // Return the EventSource instance so it can be closed by the caller if needed
  return eventSource;
};
