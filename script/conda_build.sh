#!/bin/bash

if [[ "$TRAVIS_OS_NAME" == "osx" ]]; then
	brew update
    export CONDA_OS=MacOSX
elif [[ "$TRAVIS_OS_NAME" == "linux" ]]; then
	sudo apt-get update
	export CONDA_OS=Linux
fi

if [[ "$TRAVIS_PYTHON_VERSION" == "2.7" ]]; then
	wget https://repo.continuum.io/miniconda/Miniconda2-latest-$CONDA_OS-x86_64.sh -O miniconda.sh;
else
	wget https://repo.continuum.io/miniconda/Miniconda3-latest-$CONDA_OS-x86_64.sh -O miniconda.sh;
fi

bash miniconda.sh -b -p $HOME/miniconda

export PATH="$HOME/miniconda/bin:$PATH"
hash -r
conda config --set always_yes yes --set changeps1 no --set anaconda_upload no
conda update -q conda
conda info -a # Useful for debugging any issues with conda
conda install conda-build anaconda-client

conda config --add channels conda-forge
conda config --add channels bioconda
conda config --add channels zertan

conda create -q -n test-environment python=$TRAVIS_PYTHON_VERSION
source activate test-environment
conda build recipe/menace
conda install -f --use-local -y menace