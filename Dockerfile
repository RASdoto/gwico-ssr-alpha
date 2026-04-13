FROM python:3.12-slim

LABEL maintainer="GWICO-SSR Team"
LABEL version="0.1.0a1"
LABEL description="GWICO-SSR: In-silico SSR identification and characterization platform"

WORKDIR /app

# Install system dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends gcc && \
    rm -rf /var/lib/apt/lists/*

# Copy project files
COPY pyproject.toml README.md ./
COPY src/ src/
COPY config/ config/
COPY data/ data/
COPY examples/ examples/

# Install the package
RUN pip install --no-cache-dir .

# Create output and data directories
RUN mkdir -p /app/outputs /app/data

# Default config
ENV GWICO_SSR_DB_URL=sqlite:////app/data/gwico_ssr.db
ENV GWICO_SSR_LOG_FORMAT=json

ENTRYPOINT ["python", "-m", "gwico_ssr"]
CMD ["--help"]
