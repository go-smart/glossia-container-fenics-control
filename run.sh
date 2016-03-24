#!/bin/sh

IMAGE_ID="gosmart/glossia-goosefoot"

if [ -f input/glossia-simimage.txt ];
then
    if [ "$1" = "--simimage" ];
    then
        IMAGE_ID=$(cat input/glossia-simimage.txt)
        echo "Using image ${IMAGE_ID}"
    elif [ "$1" = "--no-simimage" ];
    else
        echo <<EOM
Diagnostic input indicates a specific image ID (i.e. version of the simulation software) "
Use --simimage to pull it in, or --no-simimage to ignore. While this image may be on your
computer already, it is recommended that you only allow pulling of remote images if the
simulation comes from a trusted source
EOM
    fi
fi

export COMPOSE_API_VERSION=$(docker version | grep 'Server API' | awk '{ print $NF }')
IMAGE_ID="${IMAGE_ID}" docker-compose up
