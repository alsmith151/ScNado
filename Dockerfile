# Build stage
FROM python:3.12-bookworm AS builder

WORKDIR /usr/src/scnado
COPY . .

# Install Rust nightly and build dependencies
RUN apt-get update && apt-get install -y curl pkg-config libssl-dev build-essential cmake
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y --default-toolchain nightly
ENV PATH="/root/.cargo/bin:${PATH}"

RUN pip install maturin

# Build the Rust binary and Python package
RUN maturin build --release --interpreter python3.12

# Create a directory for all wheels (including dependencies)
RUN mkdir /wheels && \
    pip wheel /usr/src/scnado/target/wheels/*.whl -w /wheels

# Runtime stage
FROM python:3.12-slim-bookworm

WORKDIR /app

# Install minimal runtime libraries
RUN apt-get update && apt-get install -y \
    libssl3 \
    && rm -rf /var/lib/apt/lists/*

# Copy all pre-built wheels from the builder stage
COPY --from=builder /wheels /wheels
RUN pip install --no-cache-dir /wheels/*.whl && rm -rf /wheels

# Set the entrypoint to the scnado binary (if it was installed as a script)
# Or just keep it as a base image for Snakemake
CMD ["scnado", "--help"]
