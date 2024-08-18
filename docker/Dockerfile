# Please install docker according to the guideline provided here: https://docs.docker.com/engine/install/ubuntu/

# Use Debian as the base image
FROM debian:latest

# Set non-interactive frontend for APT
ARG DEBIAN_FRONTEND=noninteractive

# Install basic packages and SageMath
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    build-essential \
    git \
    cmake \
    wget \
    curl \
    python3 \
    python3-full \
    sagemath && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Install MiniZinc
WORKDIR /home 
RUN wget https://github.com/MiniZinc/MiniZincIDE/releases/download/2.8.5/MiniZincIDE-2.8.5-bundle-linux-x86_64.tgz && \
    tar -xvzf MiniZincIDE-2.8.5-bundle-linux-x86_64.tgz && \
    mv MiniZincIDE-2.8.5-bundle-linux-x86_64 minizinc && \
    rm MiniZincIDE-2.8.5-bundle-linux-x86_64.tgz && \
    ln -s /home/minizinc/bin/minizinc /usr/local/bin/minizinc

# Clone DL repository
RUN git clone https://github.com/hadipourh/DL

# Set working directory
WORKDIR /home/DL

# Create a virtual environment
RUN python3 -m venv myenv

# Install required Python packages
RUN myenv/bin/python3 -m pip install --no-cache-dir pyyaml minizinc gurobipy numpy

# Clean up APT cache
RUN apt-get clean && rm -rf /var/lib/apt/lists/*
