#!/usr/bin/env bash
#
# simple script to test SMAC functionality
# note: will install additional python packages via pip --user

# make sure we fail as soon as a command fails
set -e

ROOTDIR="$(git rev-parse --show-toplevel)"
if [ -z "$ROOTDIR" ]
then
	echo "error: not called from within a git repository"
	exit 1
fi
cd "$ROOTDIR"

if [ ! -x "$ROOTDIR"/build/bin/riss-core ]
then
	echo "error: cannot find a riss executable 'build/bin/riss-core'"
	exit 1
fi

# make sure python3 is available and has all the packages we want to have
echo
echo "Install SMAC3 dependencies ..."
sudo apt-get install python3-pip swig

pushd tools

# get SMAC3
if [ -d SMAC3 ]
then
	cd SMAC3
	git pull origin development 
	cd ..
else
	git clone https://github.com/automl/SMAC3.git
	git checkout development 
fi

# install required modules via pip3 for the current user
cat requirements.txt | xargs -n 1 -L 1 pip3 install --user
popd

# copy from spear example
echo
echo "Setup Riss example ..."
pushd tools/SMAC3/examples/

# "create" a directory for riss, based on spear, but use our scenario file
ln -sf spear_qcp riss
cp "$ROOTDIR"/tools/smac-files/scenario_riss.txt riss

# get all the other files right for riss
cd riss/target_algorithm
mkdir -p riss-python
cd riss-python

# get link to riss binary, and get configuration file here
echo
echo "Create Riss configuration ..."
cp "$ROOTDIR"/tools/smac-files/* .
ln -sf "$ROOTDIR"/build/bin/riss-core riss
./riss -config= -pcs-file=riss.pcs
popd

# give smac a try
echo
echo "Start running SMAC ..."
pushd tools/SMAC3/examples/riss
python3 ../../scripts/smac --scenario scenario_riss.txt --verbose DEBUG
popd
