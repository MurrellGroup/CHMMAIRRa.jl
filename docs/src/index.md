# CHMMAIRRa.jl

[![MIT license](https://img.shields.io/badge/license-MIT-green.svg)](https://opensource.org/license/MIT)
[![Coverage](https://codecov.io/gh/MurrellGroup/CHMMAIRRa.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/MurrellGroup/CHMMAIRRa.jl)


CHMMAIRRa.jl is a method for reference-based detection of chimeric adaptive immune receptor repertoires (AIRR) sequences.

Chimera detection is performed using a Hidden Markov Model (HMM) that models sequences as being generated from a single reference with mutations or from multiple references with mutations. The core HMM functionality is implemented in the [CHMMera.jl](https://github.com/MurrellGroup/CHMMera.jl) package, while CHMMAIRRa.jl provides an AIRR-specific wrapper for the HMM.


# Installation

The tool is primarily distributed as a an executable for Linux and MacOS from the [releases](https://github.com/MurrellGroup/CHMMAIRRa.jl/releases) page. Once you have downloaded the executable, you will need to unpack it and add it to your PATH to allow it to be run from any directory in the terminal.

```
tar -xvf CHMMAIRRa.tar.gz # unpack the executable
export PATH=$PATH:$(pwd)/CHMMAIRRa/bin/ # add the executable to your PATH
```

The tool is also available from the [Julia package manager](https://pkg.julialang.org/), and can be installed with the following command:

```
using Pkg; Pkg.add("CHMMAIRRa")
```

See the [Julia usage](#Julia-usage) section for more details.

# Usage

For chimera detection with a minimal output, run the following command:

```
CHMMAIRRa --receptor IG --V-fasta V.fasta --assignments AIRR.tsv --out CHMMAIRRa_out.tsv
```

For chimera detection with a detailed output including recombination path information, run the following command:

```
CHMMAIRRa --receptor IG --V-fasta V.fasta --assignments AIRR.tsv --out CHMMAIRRa_out.tsv --detailed
```

If you want the output to include all columns from the input file, use the ```--non-chimeric-MiAIRR``` and ```--chimeric-MiAIRR``` arguments:

```
CHMMAIRRa --receptor IG --V-fasta V.fasta --assignments AIRR.tsv --out CHMMAIRRa_out.tsv --non-chimeric-MiAIRR nonchimeric.tsv --chimeric-MiAIRR chimeric.tsv
```


# Arguments
## Mandatory arguments
1. ```--receptor``` : The receptor type, either ```IG``` or ```TCR```.
2. ```--V-fasta``` : The database V sequences in FASTA format.
3. ```--assignments``` : The VDJ assignments in MiAIRR format. Must have been assigned to the ```--V-fasta``` database.
4. ```--out``` : The output file path. See [CHMMAIRRa output columns](#chmmairra-output-columns) for more details.

## Optional arguments
### Output options
1. ```--detailed``` : Whether to output detailed information about the recombination path. See [CHMMAIRRa output columns](#chmmairra-output-columns) for more details.
2. ```--non-chimeric-MiAIRR``` : The output file path for non-chimeric sequences, including all columns from the input file.
3. ```--chimeric-MiAIRR``` : The output file path for chimeric sequences, including all columns from the input file.

### Extra outputs
1. ```--recombfreqplot``` : A plot of chimerism per recombination event, normalized by allele frequency. Can be useful for identifying missing database sequences.
2. ```--chimeric-alignments``` : Fasta file containing chimeric sequences and the germline alignments they matched to. This also provides a hamming_distances column in the chimeric assignments file.

### Technical parameters
1. ```--p-threshold``` : Posterior threshold of switching templates (default: 0.95).
2. ```--HMM-parameters``` : If set, use HMM parameters in this file. Overrides parameter presets from --receptor. See src/IG_parameters.tsv or src/TCR_parameters.tsv for formatting.
3. ```--mafft``` : The path to the MAFFT executable. Automatically found if not provided.
4. ```--align-database``` : Whether to align the database V sequences. Will be done automatically if database V sequences are not the same length.
5. ```--subsample``` : Randomly subsample this number of unique sequences to run chimera detection on.
6. ```--deduplicate``` : Deduplicate input by v_sequence_alignment column. Note that if subsample is given, subsampling is performed after deduplication.
7. ```--chunk-size``` : Read in chunk-size lines of input at a time, write output for that chunk, and move to the next chunk. This decreases memory consumption in exchange for increasing runtime.
8. ```--count-chimeric-segments``` : Whether to count the number of times chimeric segments appear in nonchimeric sequences. Adds ```chimera_segment_counts``` column.
9. ```--trim``` : How many nucleotides to trim off the inner ends of the chimeric segments when searching for matches in nonchimeric sequences.
10. ```--seed``` : Seed to use for random subsampling.
11. ```--quiet``` : Whether to suppress non-error messages.

# CHMMAIRRa output columns
## Mandatory columns
- ```threaded``` : The query V sequence threaded onto the database alignment.
- ```chimera_probability``` : The probability of the V sequence being a chimera.
- ```chimeric``` : Whether the V sequence is chimeric.

## Optional columns generated when using the ```--detailed``` arg
- ```startingpoint``` : The nonchimeric starting state of the viterbi path (i.e. the HMM's v_call).
- ```recombinations``` : A list of recombination points according to the viterbi path of the query. Format: (v_call 1, v_call 2, recombination position in threaded query sequence).
- ```recombinations_degapped``` : Same as recombinations, but with the non-gapped position of the recombination.
- ```pathevaluation``` : The probability of the second likeliest reference in the HMM. This gives us a measure of confidence for the recombination path itself.
- ```nrecombinations``` : The number of recombination points in the viterbi path.

## Optional outputs generated when using the ```--count-chimeric-segments``` arg
- ```segment_1_count``` : The number of times the left chimeric segment appears in nonchimeric sequences.
- ```segment_2_count``` : The number of times the right chimeric segment appears in nonchimeric sequences.


# Julia usage

CHMMAIRRa.jl's Julia interface provides two functions: ```detect_chimeras``` and ```detect_chimeras_from_files```; the first takes paths to the input files and writes to a file, while the second takes a DataFrame and returns a DataFrame.

```
using CHMMAIRRa

# run from DataFrame
refnames, refseqs = read_fasta("V.fasta")
assignments = CSV.read("MiAIRR.tsv", DataFrame, delim = '\t')
chmmairra_output = detect_chimeras(refnames, refseqs, assignments, receptor = "IG")

# run from file paths
detect_chimeras_from_files("V.fasta", "MiAIRR.tsv", "CHMMAIRRa_out.tsv", receptor = "IG")
```

# Cite

If you use CHMMAIRRa.jl in your work, please cite the following paper:

```
```
