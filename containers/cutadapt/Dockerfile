FROM python:3.11-slim

LABEL software="cutadapt"
LABEL version="1.11"

# Install cutadapt with fixed version for reproducibility
RUN pip install --no-cache-dir cutadapt==1.11

ENTRYPOINT ["cutadapt"]
