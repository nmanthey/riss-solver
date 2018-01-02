#!/usr/bin/env bash
# test-SMAC3.sh, Norbert Manthey, 2018, LGPL v2, see LICENSE
#
# simple python script to run SMAC3 tuning example
#

# fail early
set -e

# check for system dependencies
for tool in python3 pip3 swig
do
	if ! command -v "$tool" &> /dev/null
	then
		echo "error: did not find $tool, please install and run again!"
		exit 1
	fi
done

# check for python dependencies
if ! command -v virtualenv &> /dev/null
then
	echo "error: did not find virtualenv, trying to install via pip3 --user"
	pip3 install virtualenv --user
fi

# go to target directory
cd tools

# get SMAC3 source if not present already
PRESENT=
[ ! -d SMAC3 ] || PRESENT=1
if [ -z "$PRESENT" ]
then
	git clone https://github.com/automl/SMAC3.git
fi

# enter SMAC3 environment
pushd SMAC3

# setup python environment for testing, and install dependencies
if [ -z "$PRESENT" ]
then
	virtualenv SMAC_python_env
	source SMAC_python_env/bin/activate
	cat requirements.txt | xargs -n 1 -L 1 pip3 install
	python3 setup.py install
fi

# run spear example
pushd examples/spear_qcp/
bash -x run.sh 
popd

# tear down
deactivate
popd
