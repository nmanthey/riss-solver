#!/bin/bash
# LGPL v2, see LICENSE
# 
# Run astyle for all C++ headers and files in riss-toolbox
#
# Usage: astyle.sh [directory]

if [[ "$1x" == "x" ]]; then
    root=$PWD
else
    root="$1"
fi

find $root \( -path ./tools -o -path ./build -o -path ./release \) -prune -o \( -name '*.cc' -or -name '*.c' -or -name '*.h' \) -exec astyle --options=.astylerc {} \;
# clean up backup files
find $root \( -path ./tools -o -path ./build -o -path ./release \) -prune -o \( -name '*.cc.orig' -or -name '*.c.orig' -or -name '*.h.orig' \) -exec rm -v {} \;
