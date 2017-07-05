# a few simple commands that should work once the solvers are built

# make sure each command succeeds
set -e

echo "test whether we execute from repo directory"
[ -x regression/useful-commands.sh ] || exit 1

LOG=$(mktemp)

echo "count models with portfolio solver"
# count models with pfolio solver, should result in 7 models
./build/bin/pfolio regression/cnfs/sat.cnf -models=7 2>&1 | tee $LOG | awk '/c found models: / { if ($4 >= 7) {exit 0} else {exit 1}}'
