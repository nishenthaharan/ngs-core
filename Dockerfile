FROM python:3.12-slim AS runtime

ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

RUN groupadd --system ngscore && useradd --system --gid ngscore --create-home ngscore

WORKDIR /app
COPY pyproject.toml README.md LICENSE ./
COPY src ./src
RUN python -m pip install --no-cache-dir .

USER ngscore
WORKDIR /data
ENTRYPOINT ["ngs-core"]
CMD ["--help"]
