#! /bin/bash

# This script will install extra packages needed for 2D limit plot:
# --> draw2DScan.py <--
#
# The script will install additional python packages needed.
#
# NOTE: after the installation you may still need to set paths manually:
#
#   export PATH=$PWD/python_packages/bin:$PATH
#   export PYTHONPATH=$PWD/python_packages/lib/python2.7/site-packages:$PYTHONPATH
#

function setup_env {
    export PATH=$PWD/python_packages/bin:$PATH
    export PYTHONPATH=$PWD/python_packages/lib/python2.7/site-packages:$PYTHONPATH
}

if [ ! -d python_packages ]; then
    echo "Installing needed python packages..."

    mkdir python_packages

    prefix=$PWD/python_packages

    curl -O https://bootstrap.pypa.io/get-pip.py
    python get-pip.py --prefix=$prefix
    rm get-pip.py

    setup_env

    # Install numpy
    pip install numpy --prefix=$prefix
    pip install python-dateutil --prefix=$prefix
    pip install matplotlib --prefix=$prefix --global-option=build_ext --global-option="-L$(python2.7-config --prefix)/lib"
    pip install six --prefix=$prefix --global-option=build_ext --global-option="-L$(python2.7-config --prefix)/lib"
else
    setup_env
fi