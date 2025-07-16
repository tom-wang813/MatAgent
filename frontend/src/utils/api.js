export const sendMessage = (message, conversationId, onNewStep, onError, onOpen) => {
  const backendUrl = process.env.REACT_APP_BACKEND_URL || 'http://localhost:8000/api';
  const traceId = `trace-${Date.now()}-${Math.random().toString(36).substring(2, 15)}`;

  // Construct URL with query parameters for GET request with EventSource
  const url = new URL(`${backendUrl}/chat`);
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
    if (event.data) {
      try {
        const step = JSON.parse(event.data);
        onNewStep(step);
      } catch (e) {
        console.error('Error parsing SSE message:', e, 'Data:', event.data);
        if (onError) onError(e);
      }
    }
  };

  eventSource.onerror = (err) => {
    console.error('SSE error:', err);
    eventSource.close();
    if (onError) onError(err);
  };

  // Return the EventSource instance so it can be closed by the caller if needed
  return eventSource;
};
