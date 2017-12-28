# Riss tool collection. 2010-2017, Norbert Manthey, LGPL v2, see LICENSE

You can receive the latest copy of the Riss and Coprocessor tool collection from
http://tools.computational-logic.org

Please send bug reports to norbert.manthey@tu-dresden.de

All tools in this package are highly configurable. For a simple comparison, most 
techniques are disabled by default. Hence, you should specify options to use. For
both Riss and Coprocessor, using -config=Riss427 is a fair starting point.

## Getting Started

For a quick example on how to retrieve, build and use the solver on a Linux
machine, have a look at doc/TUTORIAL.md!

## Components

This software package contains the SAT solver Riss, and might contain related 
tools:

| Tool        | Description                                 |
| ----------- | ------------------------------------------- |
| Coprocessor |  CNF simplifier                             |
| Priss       |  simple portfolio SAT solver                |


## Building

The tool box uses [CMake](http://cmake.org/) as build tool chain. In-source
builds are not recommended. It is better to build the different build types
(release, debug) in a separate directory.

### Dependencies

To compile the project, libz is required, which can be installed for example
with

```bash
sudo apt-get install zlib1g-dev
```

To test the ipasir package, curl is required, which can be installed for
example with

```bash
sudo apt-get install curl
```

### Building an example configuration


```bash
# Create a directory for the build
mkdir debug
cd debug/

# Create a debug build configuration
cmake -D CMAKE_BUILD_TYPE=Debug ..

# Create release build Intel compiler
cmake -D CMAKE_BUILD_TYPE=Release -D CMAKE_CXX_COMPILER=icpc ..

# If you want DRAT proof support, configure cmake the following
cmake -D DRATPROOF=ON ..

# You get a list of all targets with
make help

# Two different configurations of the riss solver
# (with and without simplification enabled)
make riss-simp
make riss-core

# Copy all scripts to the bin/ directory
make scripts
```

### Options

To configure your build, pass the described options to cmake like this

```bash
cmake -D OPTION_NAME=value ..
```

| Option             | Description                                            | Default |
| ------------------ | ------------------------------------------------------ | ------- |
| CMAKE_BUILD_TYPE   | "Debug" or "Release"                                   | Release |
| CMAKE_CXX_COMPILER | C++ compiler that should be used                       |     g++ |
| STATIC_BINARIES    | Build fully statically linked binaries                 |      ON |
| WARNINGS           | Set verbose warning flags                              |     OFF |
| QUIET              | Disable Wcpp                                           |     OFF |


## Incremental Solving

Riss supports two different C interfaces, where one is the IPASIR interface, which has been
set up for incremental track of the SAT Race in 2015. The actual interface of Riss supports
a few more routines. Furthermore, Coprocessor's simplification can be used via a C interface.

### Build the solver for the IPASIR interface

After running cmake (see above), build the following library:

```bash
make riss-coprocessor-lib-static
```
Then, include the header file "riss/ipasir.h" into your project, and link against the library.

## Common Usage

The available parameters can be listed for each tool by calling:

    bin/<tool> --help

Due to the large number of parameters, a more helpful alternative is:
    
    bin/<tool> --help -helpLevel=0 --help-verb

The configuration specification can be written to a pcs file automatically. 
Before using this file with an automated configuration framework, please check
that only necessary parameters appear in the file. The procedure will include 
all parameters, also the cpu or memory limit parameters, and is currently considered
experimental.

  bin/<tool> -pcs-file=<pcs-filename>
  

Using Riss to solve a formula <input.cnf> use

    bin/riss <input.cnf> [solution] [solverOptions]

Using Riss to solve a formula <input.cnf> and generating a DRUP proof, the 
following call can be used. Depending on the proof verification tool the option
`-no-proofFormat` should be specified. Note, not all formula simplification
techniques support generating a DRAT proof.

    bin/riss <input.cnf> -proof=input.proof -no-proofFormat [solution] [solverOptions]

The script `cp.sh` provides a simple setup for using Coprocessor as a formula
simplification tool and afterwards running a SAT solver. The used SAT solver can
be exchanged by changing the variable `satsolver` in the head of the script.

    bin/cp.sh <input.cnf> [coprocessorOptions]

The parallel portfolio solver priss uses incarnations of riss and executes them in
parallel. To obtain a version that executes exact copies of the solver, issue the 
following call, and add the CNF formula as well. Furthermore, you might want to specify
the number of used threads by adding "-threads=X"

    bin/priss -ppconfig= -no-addSetup -no-pr -no-ps -psetup=PLAIN -pAllSetup=-independent
