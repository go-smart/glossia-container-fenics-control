#!/bin/sh

export COMPOSE_API_VERSION=$(docker version | grep 'Server API' | awk '{ print $NF }')
docker-compose up
