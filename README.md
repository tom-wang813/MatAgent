# MatAgent: AI-Powered Materials Science Assistant

## Project Introduction

MatAgent is an innovative AI-powered platform designed to revolutionize materials science and chemistry research. It acts as an intelligent assistant, leveraging the power of large language models (LLMs) to interact with specialized machine learning models and computational tools. The primary goal of MatAgent is to empower researchers by providing on-demand predictions for various material properties based on molecular structures (e.g., SMILES strings), streamlining the research and development process.

### Key Features:
- **Intelligent Agent Orchestration:** An advanced AI agent understands user queries, plans actions, and orchestrates the use of appropriate tools and ML models.
- **Material Property Prediction:** Integrates with specialized machine learning models to predict critical material properties such as aqueous solubility, radius of gyration (Rg), and heat capacity (Cp).
- **User-Friendly Interface:** A responsive web frontend provides an intuitive chat-based interface for seamless interaction with the AI assistant.
- **Modular and Scalable Architecture:** Built with a microservices approach using FastAPI for the backend and MCP server, ensuring maintainability and scalability.
- **Containerized Deployment:** Utilizes Docker and Docker Compose for easy, consistent, and reproducible deployment across different environments.

### Architecture Overview:

The MatAgent platform is composed of three main services, orchestrated using Docker Compose:

1.  **Frontend (React.js):**
    *   A modern web application built with React.js.
    *   Provides the graphical user interface (GUI) for users to input queries and view responses from the AI agent.
    *   Communicates with the Backend API to send user messages and receive AI-generated responses, including tool outputs.

2.  **Backend (FastAPI):**
    *   The central intelligence hub of MatAgent, developed using FastAPI (Python).
    *   Hosts the core AI agent logic, including:
        *   **LLM Orchestration:** Manages interactions with external Large Language Models (e.g., via OpenRouter).
        *   **Conversation Management:** Stores and retrieves conversation history.
        *   **Tool Management:** Discovers, selects, and executes specialized tools.
        *   **Embedding Service:** Utilizes `sentence-transformers` for generating embeddings, crucial for tool selection and context understanding.
    *   Connects to the MCP Server to utilize its specialized ML models as tools.

3.  **MCP Server (FastAPI):**
    *   A dedicated Microservice for Material Chemistry Properties (MCP) calculations, also built with FastAPI (Python).
    *   Houses various pre-trained machine learning models (e.g., for solubility, Rg, Cp prediction).
    *   Exposes API endpoints that the Backend's AI agent can call as "tools" to perform specific material property predictions.
    *   Designed to be extensible, allowing for the easy addition of new computational tools and ML models.

## How It Works:

1.  A user enters a query (e.g., "What is the aqueous solubility of ethanol (CCO)?") into the Frontend.
2.  The Frontend sends the query to the Backend API.
3.  The Backend's AI agent processes the query using the LLM.
4.  The LLM, guided by the agent's logic, identifies the need for a specific tool (e.g., `predict_aqueous_solubility`).
5.  The agent calls the corresponding API endpoint on the MCP Server, passing the necessary parameters (e.g., SMILES string "CCO").
6.  The MCP Server executes its internal ML model to predict the property and returns the result to the Backend.
7.  The Backend's AI agent integrates this result into a coherent response and sends it back to the Frontend.
8.  The Frontend displays the AI's response, including the predicted property value, to the user.

## How to Deploy

This project uses Docker and Docker Compose for easy setup and deployment.

### Prerequisites

Before you begin, ensure you have the following installed on your system:

*   **Docker Desktop:** Includes Docker Engine and Docker Compose.
    *   [Download Docker Desktop](https://www.docker.com/products/docker-desktop/)

### Environment Variables

You need to set up an API key for the LLM service (e.g., OpenRouter). Create a `.env` file in the `backend/` directory with the following content:

```
OPENROUTER_API_KEY=your_openrouter_api_key_here
```

Replace `your_openrouter_api_key_here` with your actual API key.

### Steps to Deploy

1.  **Clone the Repository:**
    If you haven't already, clone the MatAgent repository to your local machine:
    ```bash
    git clone https://github.com/your-repo/matagent.git
    cd matagent
    ```
    (Note: Replace `https://github.com/your-repo/matagent.git` with the actual repository URL if different.)

2.  **Build and Run with Docker Compose:**
    Navigate to the root directory of the cloned repository (where `docker-compose.yml` is located) and run the following command:

    ```bash
    docker-compose up --build
    ```
    *   `--build`: This flag ensures that Docker rebuilds the images for all services (frontend, backend, mcp_server) if there are any changes in the Dockerfiles or dependencies. It's good practice to include this the first time you run it or after pulling updates.

    This command will:
    *   Build the Docker images for the `frontend`, `backend`, and `mcp_server` services.
    *   Start all three services in detached mode (you'll see logs in your terminal).
    *   The `backend` service will wait for the `mcp_server` to be ready.
    *   The `frontend` service will wait for the `backend` to be ready.

3.  **Access the Application:**
    Once all services are up and running, you can access the MatAgent frontend in your web browser:

    *   **Frontend:** `http://localhost`

    The frontend is configured to run on port 80, and Docker Compose maps this to your host's port 80.

4.  **Verify Services (Optional):**
    You can check the status of individual services:
    *   **Backend Health Check:** `http://localhost:8000/api/health` (You should see a success message)
    *   **MCP Server Health Check:** `http://localhost:8080/` (You should see "MCP Server is running")

### Stopping the Application

To stop all running services and remove the containers, networks, and volumes created by `docker-compose up`:

```bash
docker-compose down
```

If you want to stop and remove volumes (e.g., for a clean start of the database), use:

```bash
docker-compose down -v
```

This will remove the `matagent.db` file created by the backend, allowing for a fresh database on the next `docker-compose up`.