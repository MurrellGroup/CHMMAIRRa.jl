# CHMMAIRRa.jl

[![Documentation](https://img.shields.io/badge/docs-latest-blue.svg)](https://murrellgroup.github.io/CHMMAIRRa.jl/)
[![MIT license](https://img.shields.io/badge/license-MIT-green.svg)](https://opensource.org/license/MIT)
[![Coverage](https://codecov.io/gh/MurrellGroup/CHMMAIRRa.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/MurrellGroup/CHMMAIRRa.jl)


CHMMAIRRa.jl is a method for reference-based detection of chimeric adaptive immune receptor repertoire (AIRR) sequences.

We use a Hidden Markov Model (HMM) to represent sequences as being generated from a single reference with mutations or from multiple references with mutations. The core HMM functionality is implemented in the [CHMMera.jl](https://github.com/MurrellGroup/CHMMera.jl) package, while CHMMAIRRa.jl provides an AIRR-specific wrapper for the HMM.

## Quick start

### Using the executable

The tool is available as an [executable](https://github.com/MurrellGroup/CHMMAIRRa.jl/releases). You will need to unpack it and add it to your PATH to allow it to be run from any directory in the terminal.
```bash
unzip CHMMAIRRa_ubuntu-linux_x86-64_v0.0.1.zip # unpack the executable
export PATH=$PATH:$(pwd)/CHMMAIRRa_ubuntu-linux_x86-64_v0.0.1/bin/ # add the executable to your PATH. Place this line in your ~/.bashrc or ~/.zshrc file to make it permanent.
CHMMAIRRa --receptor IG --V-fasta V.fasta --assignments MiAIRR.tsv --out CHMMAIRRa_out.tsv
```

### Using the Julia package

```julia
using Pkg; Pkg.add("CHMMAIRRa")
using CHMMAIRRa
detect_chimeras_from_files("V.fasta", "MiAIRR.tsv", "CHMMAIRRa_out.tsv", receptor = "IG")
```

If you are analyzing TCRs, you can use the ```--receptor TCR``` argument.

### Filter Chimeras from AIRR Dataset
To get your input AIRR file back with chimeras removed (preserving all original columns):
```bash
CHMMAIRRa --receptor IG --V-fasta V.fasta --assignments MiAIRR.tsv \
          --out CHMMAIRRa_out.tsv --non-chimeric-MiAIRR filtered_MiAIRR.tsv
```

