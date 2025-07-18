# Stage 1: build
FROM node:20-alpine AS builder
WORKDIR /app

# 1) Copy only dependency manifests, install, then copy the rest
COPY package.json package-lock.json ./
RUN npm ci

# 2) Copy source and build
COPY . .
ARG REACT_APP_BACKEND_URL
ENV REACT_APP_BACKEND_URL=$REACT_APP_BACKEND_URL
RUN npm run build

# Stage 2: serve
FROM nginx:alpine
COPY --from=builder /app/build /usr/share/nginx/html
COPY nginx.conf /etc/nginx/conf.d/default.conf
EXPOSE 80
CMD ["nginx", "-g", "daemon off;"]
