# Build stage
FROM rust:1.83-slim-bookworm AS builder

WORKDIR /usr/src/scnado
COPY . .

RUN apt-get update && apt-get install -y pkg-config libssl-dev python3-dev python3-pip python3-venv
RUN pip3 install --break-system-packages maturin

# Build the Rust binary and Python package
RUN maturin build --release

# Runtime stage
FROM python:3.11-slim-bookworm

WORKDIR /app

# Install runtime dependencies
RUN apt-get update && apt-get install -y \
    libssl3 \
    && rm -rf /var/lib/apt/lists/*

# Copy the built wheel from the builder stage and install it
COPY --from=builder /usr/src/scnado/target/wheels/*.whl /tmp/
RUN pip install /tmp/*.whl && rm /tmp/*.whl

# Set the entrypoint to the scnado binary (if it was installed as a script)
# Or just keep it as a base image for Snakemake
CMD ["scnado", "--help"]
