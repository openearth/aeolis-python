FROM python:3.11-slim

# Set the working directory
WORKDIR /usr/src/app

# Copy the pyproject.toml and install dependencies
COPY pyproject.toml ./
RUN pip install --no-cache-dir hatchling && \
    pip install --no-cache-dir .

# Copy the rest of the application code
COPY . .

# Command to run when the container starts
CMD ["python", "-m", "aeolis.console"]