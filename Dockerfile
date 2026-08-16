FROM python:3.13-slim

WORKDIR /app

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY src/ ./src/
COPY data/reference/Genomes/ ./data/reference/Genomes/

ENTRYPOINT ["python", "src/pipeline.py"]
