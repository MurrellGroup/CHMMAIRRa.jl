# CHMMAIRRa.jl

[![Documentation](https://img.shields.io/badge/docs-latest-blue.svg)](https://murrellgroup.github.io/CHMMAIRRa.jl/)
[![MIT license](https://img.shields.io/badge/license-MIT-green.svg)](https://opensource.org/license/MIT)
[![Coverage](https://codecov.io/gh/MurrellGroup/CHMMAIRRa.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/MurrellGroup/CHMMAIRRa.jl)


CHMMAIRRa.jl is a method for reference-based detection of chimeric adaptive immune receptor repertoire (AIRR) sequences. 

We use a Hidden Markov Model (HMM) to represent sequences as being generated from a single reference with mutations or from multiple references with mutations. The core HMM functionality is implemented in the [CHMMera.jl](https://github.com/MurrellGroup/CHMMera.jl) package, while CHMMAIRRa.jl provides an AIRR-specific wrapper for the HMM.


# Installation

The tool is primarily distributed as a an executable for Linux and MacOS from the [releases](https://github.com/MurrellGroup/CHMMAIRRa.jl/releases) page. Once you have downloaded the executable, you will need to unpack it and add it to your PATH to allow it to be run from any directory in the terminal.

```
tar -xvf CHMMAIRRa.tar.gz # unpack the executable
export PATH=$PATH:$(pwd)/CHMMAIRRa/bin/ # add the executable to your PATH
```

The tool is also available from the [Julia package manager](https://pkg.julialang.org/).

# Usage

For chimera detection with a minimal output, run the following command:

```
CHMMAIRRa --receptor IG --V-fasta V.fasta --assignments AIRR.tsv --out CHMMAIRRa_out.tsv
```