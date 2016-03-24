#!/bin/sh

IMAGE_ID="gosmart/glossia-goosefoot"

if [ -f input/glossia-simimage.txt ];
then
    IMAGE_ID=$(cat input/glossia-simimage.txt)
fi

export COMPOSE_API_VERSION=$(docker version | grep 'Server API' | awk '{ print $NF }')
IMAGE_ID="${IMAGE_ID}" docker-compose up
