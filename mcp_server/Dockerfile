# Use the official Python 3.11 slim image as the base
FROM python:3.11-slim

# Set working directory inside the container
WORKDIR /app

# Install system-level build dependencies (cached unless you change this line)
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
         build-essential \
         libxrender1 \
         libxext6 \
    && rm -rf /var/lib/apt/lists/*

# Copy Python dependency list into the image
COPY ./requirements.txt .

# Upgrade pip and install all Python packages (including Celery, which should be in requirements.txt)
RUN pip install --no-cache-dir --upgrade pip     && pip install --no-cache-dir -r requirements.txt

# Copy only the necessary application source code
COPY ./mcp_server /app/mcp_server
COPY ./shared_utils /app/shared_utils

# Ensure Python can find your modules
ENV PYTHONPATH=/app

# Expose port 8080 for incoming HTTP
EXPOSE 8080

# Default command: start Uvicorn serving your FastAPI app
CMD ["uvicorn", "mcp_server.main:app", "--host", "0.0.0.0", "--port", "8080"]
