#!/bin/bash

set -e

# 1. Check if `uv` is installed
if ! command -v uv &> /dev/null; then
    echo "uv is not installed. Installing uv..."
    curl -LsSf https://astral.sh/uv/install.sh | sh    
    source "$HOME/.local/bin/env"
else
    echo "uv is already installed."
fi

# 2. Check if Python 3.12 is installed, otherwise use `uv` to install it
if ! python3.12 --version &> /dev/null; then
    echo "Python 3.12 is not installed. Installing Python 3.12 with uv..."
    uv python install 3.12
else
    echo "Python 3.12 is already installed."
fi

# 3. Check if `.venv` exists, otherwise create one
rm -rf .venv
echo "Creating virtual environment..."
uv venv --python 3.12

source .venv/bin/activate

CHECKPOINTS_DIR="checkpoints"
if [ ! -d "$CHECKPOINTS_DIR" ]; then
    echo "Downloading checkpoints from BOX"
    curl -L https://brandeis.box.com/shared/static/r16dalx2vodecg93da1gb9micrtmcayt --output checkpoints.zip
    unzip checkpoints.zip -d checkpoints
else
    echo "Checkpoints folder already exists. Skipping download."
fi

if [ -f "requirements.in" ]; then
    echo "Installing dependencies from requirements.in..."
    uv pip install -r requirements.in
else
    echo "No requirements.in found. Skipping pip install."
fi
