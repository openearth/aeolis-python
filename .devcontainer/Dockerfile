FROM python:3.11-slim

# Set the working directory
WORKDIR /usr/src/app

# Copy the rest of the application code
COPY . .

RUN pip install --no-cache-dir hatchling && \
    pip install --no-cache-dir .
