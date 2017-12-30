#!/usr/bin/env bash
#
# simple script to test ipasir functionality

# make sure we fail as soon as a command fails
set -e

echo "test ipasir"

# test whether we are in the correct path
[ -f scripts/ipasir-makefile ] || exit 1
[ -x scripts/make-ipasir.sh ] || exit 1
[ -f regression/cnfs/sat.cnf ] || exit 1
[ -f regression/cnfs/unsat.cnf ] || exit 1

# build the ipasir package
./scripts/make-ipasir.sh riss_7

# setup the ipasir environment
cd tools/

# get the ipasir environment from github
IPASIRNAME=ipasirgithub
if [ ! -d "$IPASIRNAME" ]; then
    git clone https://github.com/biotomas/ipasir.git "$IPASIRNAME"
else
    pushd "$IPASIRNAME"
    git pull origin master
    popd
fi

cd "$IPASIRNAME"

# setup everything for riss
mkdir -p sat/riss_7
rm -rf sat/riss_7/*
mv ../../riss_7.tar.gz sat/riss_7/
cp ../../scripts/ipasir-makefile sat/riss_7/makefile

# build the dynamic library
make -C sat/riss_7 libipasirriss_7.so
ls sat/riss_7/libriss-coprocessor.so
rm sat/riss_7/libriss-coprocessor.so

# build one example tool
./scripts/mkone.sh genipafolio riss_7

# test two simple cases
STATUS=0
bin/genipafolio-riss_7 ../../regression/cnfs/sat.cnf || STATUS=$?
[ $STATUS -eq 10 ] || exit 1
bin/genipafolio-riss_7 ../../regression/cnfs/unsat.cnf || STATUS=$?
[ $STATUS -eq 20 ] || exit 1

# report
echo "success"
