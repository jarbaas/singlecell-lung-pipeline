FROM python:3.12-slim AS builder
RUN python -m venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"
RUN pip install --no-cache-dir --only-binary=:all: \
    anndata==0.10.8 \
    scanpy==1.10.2 \
    pandas==2.2.3 \
    celltypist==1.6.3 \
    harmonypy==0.1.0 \
    scikit-image==0.22.0

FROM python:3.12-slim
RUN apt-get update && apt-get install -y --no-install-recommends \
    procps \
    && rm -rf /var/lib/apt/lists/*
RUN useradd -m -s /bin/bash bio_user
RUN mkdir -p /workspace && chown -R bio_user:bio_user /workspace
WORKDIR /workspace
COPY --from=builder /opt/venv /opt/venv
ENV PATH="/opt/venv/bin:$PATH"
USER bio_user
CMD ["bash"]