# Use an official Python runtime as a parent image
FROM python:3.11-slim

# Install build dependencies
RUN apt-get update && apt-get install -y --no-install-recommends build-essential libxrender1 libxext6 && rm -rf /var/lib/apt/lists/*

# Set the working directory in the container
WORKDIR /app

# Copy the unified requirements file from the project root
COPY ../requirements.txt .

# Install any needed packages specified in requirements.txt
RUN pip install --no-cache-dir --upgrade pip &&     pip install --no-cache-dir -r requirements.txt &&     pip install celery

# Copy the entire project into the container
COPY . .

# Set PYTHONPATH to include the current working directory (/app) so imports work correctly
ENV PYTHONPATH=/app

# Make port 8080 available to the world outside this container
EXPOSE 8080

# Run main.py when the container launches
CMD ["uvicorn", "mcp_server.main:app", "--host", "0.0.0.0", "--port", "8080"]
