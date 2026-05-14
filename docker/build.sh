#!/bin/bash
set -ex

# Run from the repo root regardless of where the script is invoked from
cd "$(dirname "$0")/.."

# Extract version
version=$(cat ./VERSION | sed -nre 's/^[^0-9]*(([0-9]+\.)*[0-9]+).*/\1/p')

# Build the image (Forcing linux/amd64 for compatibility)
docker build --platform linux/amd64 -t danielnguyener/riboflow:latest -f docker/Dockerfile .

# Generate lists (optional, but useful)
docker run --rm danielnguyener/riboflow:latest apt list --installed 2>/dev/null > ./docker/apt.list
docker run --rm danielnguyener/riboflow:latest bash -c "source activate ribo_genome && conda list" > ./docker/conda.list

echo "Build complete. Version: $version"
