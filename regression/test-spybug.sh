#!/usr/bin/env bash
#
# simple script to test spybug functionality
# note: will install additional python packages via pip --user

# make sure we fail as soon as a command fails
set -e -x

echo "test spybug"
script_dir=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

# test whether script is there
[ -x "$script_dir"/../tools/spybug-files/spybug.sh ] || exit 1

# run the spybug script with 180 seconds search
bash -x "$script_dir"/../tools/spybug-files/spybug.sh 180
