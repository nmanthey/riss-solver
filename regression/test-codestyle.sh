#!/usr/bin/env bash
#
# script to check code style of all source files

# abort in case a command fails
set -e

# get all files
FILES=$(find \( -path ./tools -o -path ./build -o -path ./release -o -path ./coverity-dir \) -prune -o \( -name '*.cc' -or -name '*.c' -or -name '*.h' \))

# run astyle with dry-run on these files
RET=0
CHECKED=0
for file in $FILES
do
    [ -f $file ] || continue

    OUTPUT=$(astyle --options=.astylerc --dry-run --formatted $file)
    if [ -n "$OUTPUT" ]; then
        echo "style failure for file $file"
        RET=1
    fi
    CHECKED=$(($CHECKED + 1))
done

# report
[ "$RET" -eq 0 ] && echo "success -- checked $CHECKED files"
exit $RET
