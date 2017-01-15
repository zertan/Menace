#!/bin/bash

if [[ "$TRAVIS_OS_NAME" == "osx" ]]; then
	brew update
    export CONDA_OS=MacOSX
elif [[ "$TRAVIS_OS_NAME" == "linux" ]]; then
	sudo apt-get update
	export CONDA_OS=Linux
fi