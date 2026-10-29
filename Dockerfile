# Shared/classroom deployment of Protein Analyzer (hand-crafted model; ESM-2 left out to keep the image small).
#   docker build -t protein-analyzer .
#   docker run -p 8000:8000 protein-analyzer
FROM python:3.12-slim

WORKDIR /app
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY app.py serve.py ./
COPY protein_analyzer ./protein_analyzer
COPY templates ./templates
COPY static ./static
# Ship the trained hand-crafted model, reports and figures (not the ESM-2 model or caches).
COPY results/dataset.csv results/dataset_summary.json results/report*.md results/metrics_handcrafted.json \
     results/model_handcrafted.joblib results/*.png ./results/

ENV PA_PUBLIC=1 PA_HOST=0.0.0.0 PA_PORT=8000 PA_OPEN_BROWSER=0 PA_RATE_PER_MIN=30
EXPOSE 8000
RUN useradd --create-home appuser && chown -R appuser /app
USER appuser
CMD ["python", "serve.py"]
