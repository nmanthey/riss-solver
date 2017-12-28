# Riss Solver Tutorial

This file gives a simple example on how to retrieve the solver, build it on a
Linux system, and run the solver with a basic example. All steps listed in this
tutorial can be executed consecutively.

To have a look how the solver could be used, check the file "scripts/ci.sh",
which calls many different use cases of the solver and implemented tools. If
you only want to use Riss to solver SAT problems from input files, this
tutorial is the right starting point.

## Retrieving the Source

The source code of the solver is available on github. Just get a copy to build
the tool:

```bash
git clone https://github.com/nmanthey/riss-solver.git
cd riss-solver
```

## Building the Solver

The build system uses cmake, where some variables can be used to build the
solver with support for logging proofs or other features. Here, we just use
the default setup:

```bash
mkdir -p build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make riss-core -j $(nproc)
```

After the build completes, there should be a binary riss-core in the directory
bin.

```bash
cd bin
ls -l riss-core
```

## Using the solver

### Accepted Input Format

Riss, as well as many other SAT solvers, read CNF formulas. These formulas are
usually encoded in the DIMACS format, which is a simple plain text file format
which is based on numbers. The file has a header "p cnf #vars #clauses", and
comments can be added on extra lines that start with a 'c'. All other lines
contain clauses, where each clause is terminated with a '0', where a clause
might be spread over multiple lines.

Riss also accepts files that do not have the header. Furthermore, Riss can
parse gzipped DIMACS files.

### An Actual Example

In the following, we create a simple satisfiable formula that encodes that
clause "1 or (not 2) or 3", and store it in the file "input.cnf". Next,
we run Riss and store the output in the file "output.txt". The exit code of
the solver indicates whether the formula was satisfiable (10) or
unsatisfiable (20). If the exit code is 0, then the formula was not solved and
Riss cannot determine the status of the formula.

```bash
echo "1 -2 3 0" > input.cnf
./riss-core input.cnf output.txt
STATUS=$?
echo $STATUS
```

The solver is pretty verbose, and prints lots of information on stderr. This
output can be reduced by adding "-verb=0".

### Parsing the Output

The output of the solver has two important lines. The line starting with 's'
indicates whether the formula has been satisfiable ("s SATISFIABLE") or
unsatisfiable ("s UNSATISFIABLE"). In the former case, a model for the formula
is printed on line(s) starting with a 'v'. For the above input, the following
output is generated.

```bash
cat output.txt
grep "^s" output.txt
grep "^v" output.txt
```

## Other Features

Riss and the implemented tools can be used to solver multiple other problems
(run "riss-core --help -helpLevel=0 --help-verb" for details):

 * Riss can enumerate multiple/all models of a formula
 * Riss can print DRAT proofs to prove unsatisfiability
 * Riss can be used as a library, as it implements the IPASIR interface
 * Riss can print its parameter specification to automatically tune it for a
   benchmark using tools like ParamILS or SMAC
 * Coprocessor can be used to simplify CNF formulas
