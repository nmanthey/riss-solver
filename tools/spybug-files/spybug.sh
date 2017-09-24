#!/bin/bash
# spybug.sh, Norbert Manthey, 2017, LGPL v2, see LICENSE
#
# This script runs Riss in combination with SpyBug to check for problems with
# parameter combinations. The first parameter to the script is the timeout to
# run the analysis. It defaults to 15 minutes.

# fail early
set -e -u

# by default, run 900 seconds
TIMEOUT=900
[ "$#" -lt 1 ] || TIMEOUT=$1

ROOTDIR="$(git rev-parse --show-toplevel)"
if [ -z "$ROOTDIR" ]
then
	echo "error: not called from within a git repository"
	exit 1
fi

cd "$ROOTDIR"

if [ ! -x build/bin/riss-core ]
then
	echo "error: cannot find a riss executable 'build/bin/riss-core'"
	exit 1
fi


# check whether all tools are installed
echo
echo "Install tool dependencies ..."
pushd tools
	if [ -d spybug ]
	then
		cd spybug
		git pull origin master
		cd ..
	else
		git clone https://bitbucket.org/mlindauer/spybug.git
	fi
	pip install numpy --user
	pip install git+https://github.com/automl/ParameterConfigSpace.git --user
	pip install git+https://github.com/automl/pysmac --user
popd

echo
echo "Prepare the PCS File ..."
pushd tools/spybug-files
	ln -sf ../../build/bin/riss-core riss || true
	./riss -config= -pcs-file=riss.pcs
popd

# Running (from the example directory of SpyBug, because there all files are
# present already.
echo
echo "Run SpyBug for $TIMEOUT seconds ..."
pushd tools/spybug/example/
	python ../SpyBug.py --scenario ../../spybug-files/spybug_scenario.txt --n 100000 --time_limit $TIMEOUT
popd
