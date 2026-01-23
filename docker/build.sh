#!/bin/bash
set -ex

# Copy necessary files to the docker folder so the build context is small
cp ../VERSION ./VERSION
cp ../environment.yaml ./environment.yaml

# Extract version
version=$(cat ./VERSION | sed -nre 's/^[^0-9]*(([0-9]+\.)*[0-9]+).*/\1/p')

function cleanup {
    rm ./VERSION
    rm ./environment.yaml
}
trap cleanup EXIT

# Build the image (Forcing linux/amd64 for compatibility)
docker build --platform linux/amd64 -t danielnguyener/riboflow:latest .

# Generate lists (optional, but useful)
docker run -it danielnguyener/riboflow:latest --rm apt list --installed 2>/dev/null > ./apt.list
docker run -it danielnguyener/riboflow:latest --rm danielnguyener/riboflow:latest bash -c "source activate ribo_genome && conda list" > ./conda.list

echo "Build complete. Version: $version"
