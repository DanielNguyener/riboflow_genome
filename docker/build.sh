#!/bin/bash

set -e

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

IMAGE_NAME="riboflow"
TAG="latest"
BUILD_CONTEXT="$PROJECT_ROOT"
DOCKERFILE="$SCRIPT_DIR/Dockerfile"

while [[ $# -gt 0 ]]; do
    case $1 in
        -n|--name)
            IMAGE_NAME="$2"
            shift 2
            ;;
        -t|--tag)
            TAG="$2"
            shift 2
            ;;
        -h|--help)
            echo "Usage: $0 [OPTIONS]"
            echo "Options:"
            echo "  -n, --name NAME    Image name (default: riboflow)"
            echo "  -t, --tag TAG       Image tag (default: latest)"
            echo "  -h, --help         Show this help message"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use -h or --help for usage information"
            exit 1
            ;;
    esac
done

VERSION=$(grep "version =" "$PROJECT_ROOT/VERSION" | cut -d"'" -f2)

if [ -z "$VERSION" ]; then
    echo "Warning: Could not extract version from $PROJECT_ROOT/VERSION"
    VERSION="unknown"
fi

FULL_IMAGE_NAME="${IMAGE_NAME}:${TAG}"
VERSION_IMAGE_NAME="${IMAGE_NAME}:${VERSION}"

echo "Building Docker image: $FULL_IMAGE_NAME"
echo "Detected Version: $VERSION"
echo "Build context: $BUILD_CONTEXT"
echo "Dockerfile: $DOCKERFILE"
echo ""

cd "$BUILD_CONTEXT"

docker build \
    -f "$DOCKERFILE" \
    -t "$FULL_IMAGE_NAME" \
    .

if [ "$VERSION" != "unknown" ]; then
    echo "Tagging image with version: $VERSION"
    docker tag "$FULL_IMAGE_NAME" "$VERSION_IMAGE_NAME"
fi

echo ""
echo "Generating package lists in $SCRIPT_DIR..."

docker run --rm "$FULL_IMAGE_NAME" apt list --installed 2>/dev/null | sed 's/\x1b\[[0-9;]*m//g' > "$SCRIPT_DIR/apt.list" || echo "Warning: Could not generate apt.list"

docker run --rm "$FULL_IMAGE_NAME" bash -c "source /opt/conda/etc/profile.d/conda.sh && conda list -n ribo_genome" > "$SCRIPT_DIR/conda.list" || echo "Warning: Could not generate conda.list"

echo "Lists generated: apt.list, conda.list"
echo ""
echo "Build complete!"
echo "Images available:"
echo "  - $FULL_IMAGE_NAME"
if [ "$VERSION" != "unknown" ]; then
    echo "  - $VERSION_IMAGE_NAME"
fi
echo ""
echo "To run the pipeline:"
echo "  docker run --rm -v \$(pwd):/riboflow $FULL_IMAGE_NAME nextflow run RiboFlow.groovy -params-file project.yaml"
