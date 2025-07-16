# Use an official Python runtime as a parent image
FROM registry.hub.docker.com/library/python:3.11-slim

# Install build dependencies
RUN apt-get update && apt-get install -y --no-install-recommends build-essential && rm -rf /var/lib/apt/lists/*

# Set the working directory in the container
WORKDIR /app

# Copy the unified requirements file from the project root
COPY ../requirements.txt .

# Install any needed packages specified in requirements.txt
RUN pip install --no-cache-dir --upgrade pip &&     pip install --no-cache-dir -r requirements.txt &&     pip show python-json-logger &&     ls -l /usr/local/lib/python3.11/site-packages/ &&     pip install tenacity &&     pip install tenacity

# Copy the entire project into the container
COPY . .

# Set PYTHONPATH to include the current working directory (/app) so imports work correctly
ENV PYTHONPATH=/app

# Make port 8000 available to the world outside this container
EXPOSE 8000

# Run main.py when the container launches
CMD ["uvicorn", "backend.main:app", "--host", "0.0.0.0", "--port", "8000"]