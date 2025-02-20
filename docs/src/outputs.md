### Output columns
#### Mandatory columns
- ```threaded``` : The query V sequence threaded onto the database alignment.
- ```chimera_probability``` : The probability of the V sequence being a chimera.
- ```chimeric``` : Whether ```chimera_probability``` is greater than the threshold.

#### Optional columns generated when using the ```--detailed``` arg
- ```startingpoint``` : The nonchimeric starting state of the viterbi path (i.e. the HMM's v_call).
- ```recombinations``` : A list of recombination points according to the viterbi path of the query. Format: (v_call 1, v_call 2, recombination position in threaded query sequence).
- ```recombinations_degapped``` : Same as recombinations, but with the non-gapped position of the recombination.
- ```pathevaluation``` : The probability of the second likeliest reference in the HMM. This gives us a measure of confidence for the recombination path itself.
- ```nrecombinations``` : The number of recombination points in the viterbi path.

#### Optional outputs generated when using the ```--count-chimeric-segments``` arg
- ```segment_1_count``` : The number of times the left chimeric segment appears in nonchimeric sequences.
- ```segment_2_count``` : The number of times the right chimeric segment appears in nonchimeric sequences.

### Chimeric alignments

The tool can also output a fasta file containing alignments between chimeric sequences and their correspond germline genes. Also adds hamming_distances columns to the output file, indicating the number of mismatches between the chimeric and germline sequences.

```bash
CHMMAIRRa --receptor IG --V-fasta V.fasta --assignments MiAIRR.tsv \
          --out CHMMAIRRa_out.tsv --chimeric-alignments chimeric_alignments.fasta
```

### Recombination occurrences plot

Sometimes missing database sequences can cause overrepresented recombinations. To generate a barplot of the top 10 recombinations normalized by their expected frequency, use the ```--recombfreqplot``` argument.

```bash
CHMMAIRRa --receptor IG --V-fasta V.fasta --assignments MiAIRR.tsv \
          --out CHMMAIRRa_out.tsv --recombfreqplot recombfreqplot.txt
```

By default, the plot is produced using unicode. To produce an svg/pdf using Plots.jl, you'll need to load it with ```using Plots``` and run CHMMAIRRa as a Julia package.
