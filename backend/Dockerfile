# Use official Python 3.11 slim image as base
FROM python:3.11-slim

# Set working directory
WORKDIR /app

# Install system build dependencies and curl for health checks
RUN apt-get update \
    && apt-get install -y --no-install-recommends build-essential curl \
    && rm -rf /var/lib/apt/lists/*

# Copy and install Python dependencies
# Make sure your requirements.txt already includes celery, tenacity, python-json-logger, etc.
COPY ./requirements.txt .
RUN pip install --no-cache-dir --upgrade pip \
    && pip install --no-cache-dir -r requirements.txt

# Copy only necessary application source code
COPY ./backend /app/backend
COPY ./shared_utils /app/shared_utils

# Ensure Python can locate your app modules
ENV PYTHONPATH=/app

# Expose HTTP port
EXPOSE 8000

# Default command to run the FastAPI app via Uvicorn
CMD ["uvicorn", "backend.main:app", "--host", "0.0.0.0", "--port", "8000"]
