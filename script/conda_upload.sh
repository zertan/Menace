# Only need to change these two variables
PKG_NAME=menace
USER=zertan

OS=$TRAVIS_OS_NAME-64
export CONDA_BLD_PATH=~/conda-bld
export VERSION=`python -c "exec(\"execfile('menace/version.py')\nprint __version__\")"`
anaconda -t $CONDA_UPLOAD_TOKEN upload -u $USER $CONDA_BLD_PATH/$OS/$PKG_NAME-$VERSION-py27_0.tar.bz2 --force