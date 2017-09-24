# Setup

The files in this directory are used to run SpyBug on Riss during development.
The relative path in this directory are set such that SpyBug is executed from
the directory it has been cloned into (tools/spybug/example).

## Download and Installation

The following steps assume none of the packages are already installed. SpyBug
will be installed into the tools directory.

pushd tools
git clone https://bitbucket.org/mlindauer/spybug.git
pip install numpy --user
pip install git+https://github.com/automl/ParameterConfigSpace.git --user
pip install git+https://github.com/automl/pysmac --user
popd

## Preparing the PCS File

As Riss parameters change over time, make sure to generate an up-to-date PCS
file before running SpyBug. In case the release binary should be used, change
the link accordingly!

pushd tools/spybug-files
ln -s ../../build/bin/riss-core riss || true
./riss -config= -pcs-file=riss.pcs
popd

## Running

SpyBug will be executed from the directory it has been installed into.

pushd tools/spybug/example/
python ../SpyBug.py --scenario ../../spybug-files/spybug_scenario.txt --n 100000 --time_limit 100000
popd
