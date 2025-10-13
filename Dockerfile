FROM python:3.12-slim

# Install system dependencies and Chromium
RUN apt-get update && apt-get install -y \
    wget \
    fonts-liberation \
    libasound2 \
    libatk-bridge2.0-0 \
    libatk1.0-0 \
    libc6 \
    libcairo2 \
    libcups2 \
    libdbus-1-3 \
    libdrm2 \
    libgbm1 \
    libgtk-3-0 \
    libnspr4 \
    libnss3 \
    libx11-6 \
    libxcomposite1 \
    libxdamage1 \
    libxext6 \
    libxfixes3 \
    libxrandr2 \
    libxrender1 \
    libxss1 \
    libxtst6 \
    xdg-utils \
    chromium \
    --no-install-recommends \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY . /app
RUN pip install -e . --no-cache-dir

EXPOSE 80
ENV NAME=TelomereHunter2
ENV PYTHONWARNINGS="ignore::FutureWarning"
CMD ["/bin/bash"]