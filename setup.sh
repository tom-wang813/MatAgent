#!/bin/bash

# MatAgent One-Click Setup Script
# Automated environment configuration and project startup

set -e  # Exit on any error

echo "🚀 Welcome to MatAgent One-Click Setup Script!"
echo "=============================================="

# Check if Docker is installed
check_docker() {
    echo "🔍 Checking Docker environment..."
    if ! command -v docker &> /dev/null; then
        echo "❌ Error: Docker is not installed"
        echo "Please install Docker Desktop first: https://www.docker.com/products/docker-desktop/"
        exit 1
    fi
    
    if ! docker info &> /dev/null; then
        echo "❌ Error: Docker is not running"
        echo "Please start Docker Desktop and try again"
        exit 1
    fi
    
    if ! command -v docker-compose &> /dev/null; then
        echo "❌ Error: Docker Compose is not installed"
        echo "Please install Docker Compose or use a newer version of Docker Desktop"
        exit 1
    fi
    
    echo "✅ Docker environment check passed"
}

# Configure API Key
setup_api_key() {
    echo ""
    echo "🔑 Configure OpenRouter API Key"
    echo "------------------------------"
    
    if [ -f ".env" ]; then
        echo "Found existing .env file"
        if grep -q "OPENROUTER_API_KEY" .env; then
            existing_key=$(grep "OPENROUTER_API_KEY" .env | cut -d'=' -f2)
            if [ -n "$existing_key" ] && [ "$existing_key" != "your_openrouter_api_key_here" ]; then
                echo "✅ API Key is already configured"
                read -p "Do you want to update the existing API Key? (y/N): " update_key
                if [ "$update_key" != "y" ] && [ "$update_key" != "Y" ]; then
                    return 0
                fi
            fi
        fi
    fi
    
    echo ""
    echo "📝 Please enter your OpenRouter API Key:"
    echo "   If you don't have one, visit: https://openrouter.ai/keys"
    echo "   Press Enter to skip this step (you can configure manually later)"
    echo ""
    read -p "API Key: " api_key
    
    if [ -z "$api_key" ]; then
        echo "⚠️  Skipping API Key configuration"
        echo "Please manually create .env file later and add: OPENROUTER_API_KEY=your_key_here"
        # Create example .env file
        cat > .env << EOF
# MatAgent Environment Configuration
# Please replace with your actual OpenRouter API Key
OPENROUTER_API_KEY=your_openrouter_api_key_here

# Optional: Customize other settings
# OPENROUTER_MODEL=openai/gpt-4o
# OPENROUTER_BASE_URL=https://openrouter.ai/api/v1
EOF
        echo "✅ Created example .env file"
    else
        # Create .env file
        cat > .env << EOF
# MatAgent Environment Configuration
OPENROUTER_API_KEY=$api_key

# Optional: Customize other settings
# OPENROUTER_MODEL=openai/gpt-4o
# OPENROUTER_BASE_URL=https://openrouter.ai/api/v1
EOF
        echo "✅ API Key configuration completed"
    fi
}

# Create necessary directories
setup_directories() {
    echo ""
    echo "📁 Creating necessary directories..."
    mkdir -p data
    echo "✅ Directory creation completed"
}

# Build and start services
start_services() {
    echo ""
    echo "🏗️  Building and starting MatAgent services..."
    echo "This may take a few minutes..."
    echo ""
    
    # Stop existing services
    echo "🛑 Stopping existing services..."
    docker-compose down --remove-orphans 2>/dev/null || true
    
    # Build and start
    echo "🚀 Starting services..."
    if docker-compose up --build -d; then
        echo ""
        echo "✅ Services started successfully!"
        
        # Wait for health checks
        echo "⏳ Waiting for service initialization..."
        sleep 10
        
        # Check service status
        echo ""
        echo "📊 Service status:"
        docker-compose ps
        
        echo ""
        echo "🎉 MatAgent deployment completed!"
        echo ""
        echo "🌐 Access URLs:"
        echo "   Frontend UI: http://localhost"
        echo "   Backend API: http://localhost:8000"
        echo "   Health check: http://localhost:8000/api/health"
        echo ""
        echo "📖 Usage instructions:"
        echo "   - Open http://localhost in your browser to start using"
        echo "   - Ask materials-related questions, such as molecular property prediction"
        echo "   - Supports SMILES molecular formula input"
        echo ""
        echo "🛠️  Management commands:"
        echo "   View logs: docker-compose logs -f"
        echo "   Stop services: docker-compose down"
        echo "   Restart services: docker-compose restart"
        
    else
        echo "❌ Service startup failed"
        echo "Please check the error messages and try again"
        exit 1
    fi
}

# Show help information
show_help() {
    echo ""
    echo "🆘 Troubleshooting:"
    echo ""
    echo "1. If you encounter port conflicts:"
    echo "   sudo lsof -i :80 :8000 :8080"
    echo "   kill -9 <PID>"
    echo ""
    echo "2. If Docker build fails:"
    echo "   docker system prune -f"
    echo "   docker-compose build --no-cache"
    echo ""
    echo "3. If you need to reconfigure API Key:"
    echo "   Edit .env file"
    echo "   docker-compose restart"
    echo ""
    echo "4. View detailed logs:"
    echo "   docker-compose logs backend"
    echo "   docker-compose logs mcp_server"
    echo ""
    echo "📞 Need help? Please check project documentation or submit an Issue"
}

# Main execution flow
main() {
    check_docker
    setup_api_key
    setup_directories
    start_services
    show_help
}

# Handle command line arguments
case "${1:-}" in
    --help|-h)
        echo "MatAgent One-Click Setup Script"
        echo ""
        echo "Usage: $0 [options]"
        echo ""
        echo "Options:"
        echo "  --help, -h     Show this help message"
        echo "  --stop         Stop all services"
        echo "  --restart      Restart all services"
        echo "  --logs         Show service logs"
        echo "  --status       Show service status"
        echo ""
        exit 0
        ;;
    --stop)
        echo "🛑 Stopping MatAgent services..."
        docker-compose down
        echo "✅ Services stopped"
        exit 0
        ;;
    --restart)
        echo "🔄 Restarting MatAgent services..."
        docker-compose restart
        echo "✅ Services restarted"
        docker-compose ps
        exit 0
        ;;
    --logs)
        echo "📋 Showing service logs..."
        docker-compose logs -f
        exit 0
        ;;
    --status)
        echo "📊 Service status:"
        docker-compose ps
        exit 0
        ;;
    "")
        main
        ;;
    *)
        echo "❌ Unknown option: $1"
        echo "Use $0 --help for help"
        exit 1
        ;;
esac
