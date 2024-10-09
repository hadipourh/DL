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
RUN LATEST_MINIZINC_VERSION=$(curl -s https://api.github.com/repos/MiniZinc/MiniZincIDE/releases/latest | grep -oP '"tag_name": "\K(.*)(?=")') \
    && wget "https://github.com/MiniZinc/MiniZincIDE/releases/download/$LATEST_MINIZINC_VERSION/MiniZincIDE-$LATEST_MINIZINC_VERSION-bundle-linux-x86_64.tgz" \
    && mkdir MiniZinc \
    && tar -xvzf "MiniZincIDE-$LATEST_MINIZINC_VERSION-bundle-linux-x86_64.tgz" -C MiniZinc --strip-components=1 \
    && rm "MiniZincIDE-$LATEST_MINIZINC_VERSION-bundle-linux-x86_64.tgz"
RUN ln -s /home/MiniZinc/bin/minizinc /usr/local/bin/minizinc

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
