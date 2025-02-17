var documenterSearchIndex = {"docs":
[{"location":"api/#API-Reference","page":"API","title":"API Reference","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"Modules = [CHMMAIRRa]","category":"page"},{"location":"api/#CHMMAIRRa.detect_chimeras-Tuple{Vector{String}, Vector{String}, DataFrames.DataFrame}","page":"API","title":"CHMMAIRRa.detect_chimeras","text":"detect_chimeras(refnames::Vector{String}, refseqs::Vector{String}, assignments::DataFrame;\n                HMM_parameters::Union{Dict, Nothing} = nothing,\n                p_threshold::Float64 = 0.95,\n                receptor::String = \"IG\",\n\n                chimeric_alignments::Bool = false,\n                recombfreqplot::Bool = false,\n                detailed::Bool = false,\n                subsample::Union{Int, Nothing} = nothing,\n                deduplicate::Bool = false,\n\n                count_chimeric_segments::Bool = false,\n                trim::Int = 5,\n                seed::Int = 123,\n\n                mafft::Union{String, Nothing} = nothing,\n                align_database::Bool = false,\n\n                quiet::Bool = false,\n                disable_internal_dedup::Bool = false\n                )::NamedTuple{(:out, :chimeric_alignments, :recombfreqplot),\n                            Tuple{DataFrame,\n                            Union{Nothing, Tuple{Vector{String}, Vector{String}}},\n                            Union{Nothing, Plots.Plot}\n                            }\n                }\n\nRun CHMMAIRRa on parsed arguments, returning output as a Tuple containing relevant outputs. Output tuple fields:\n\nout: the input assignments dataframe with chimerism probabilities.\nchimericalignments: names and sequences for alignments showing which reference sequences recombined at which positions to form detected chimeras. Set to nothing if chimericalignments is false.\nrecombfreqplot: barplot of most frequent recombinations, normalized by their expected frequencies. Set to nothing if recombfreqplot is false.\n\nArguments:\n\nrefnames::Vector{String}: Names of the reference sequences.\nrefseqs::Vector{String}: Sequences of the reference sequences.\nassignments::DataFrame: VDJ assignments, in MiAIRR format.\n\nOptional arguments:\n\nHMM_parameters::Union{Dict, Nothing}: If set, use HMM parameters in this dictionary. Overrides parameter presets from receptor.\np_threshold::Float64: Posterior threshold of switching templates.\nreceptor::String: Which adaptive immune receptor to run on. Options: IG, TCR.\nchimeric_alignments::Bool: Whether to generate and return chimeric sequences and the germline alignments they matched to.\nrecombfreqplot::Bool: Whether to generate and return a plot of chimerism per recombination event.\ndetailed::Bool: Whether to include details about the chimeras. Calculate viterbi path for chimeric sequences to add startingpoint, recombinations, recombinationsdegapped, and maxlog_prob columns.\nsubsample::Union{Int, Nothing}: Randomly subsample this number of unique sequences to run chimera detection on.\ndeduplicate::Bool: Whether to deduplicate input by vsequencealignment column. Note that if subsample is given, subsampling is performed after deduplication.\ncount_chimeric_segments::Bool: If set, count the number of times chimeric segments appear in nonchimeric sequences. Adds chimerasegmentcounts column.\ntrim::Int: How many nucleotides to trim off the inner ends of the chimeric segments when searching for matches in nonchimeric sequences.\nseed::Int: Seed to use for random subsampling.\nmafft::Union{String, Nothing}: Path to MAFFT binary. If not provided, will try to find mafft in system PATH.\nalign_database::Bool: Whether to align database V sequences.\nquiet::Bool: If set, hides non-error messages.\ndisable_internal_dedup::Bool: If set, does not deduplicate sequences internally. This was added for timing benchmark purposes and has no effect on the output.\n\n\n\n\n\n","category":"method"},{"location":"api/#CHMMAIRRa.detect_chimeras_from_files-Tuple{String, String, String}","page":"API","title":"CHMMAIRRa.detect_chimeras_from_files","text":"detect_chimeras_from_files(V_fasta::String, assignments::String, out::String;\n\n        p_threshold::Float64 = 0.95,\n        receptor::String = \"IG\",\n\n        non_chimeric_MiAIRR::Union{String, Nothing} = nothing,\n        chimeric_MiAIRR::Union{String, Nothing} = nothing,\n        chimeric_alignments::Union{String, Nothing} = nothing,\n        recombfreqplot::Union{String, Nothing} = nothing,\n        detailed::Bool = false,\n        subsample::Union{Int, Nothing} = nothing,\n        deduplicate::Bool = false,\n        chunk_size::Union{Int, Nothing} = nothing,\n\n        HMM_parameters::Union{String, Nothing} = nothing,\n        mafft::Union{String, Nothing} = nothing,\n        align_database::Bool = false,\n        count_chimeric_segments::Bool = false,\n        trim::Int = 5,\n        seed::Int = 123,\n\n        quiet::Bool = false,\n        disable_internal_dedup::Bool = false\n        )\n\nRun CHMMAIRRa from input files to writing of output files, taking care of I/O.\n\nArguments:\n\nV_fasta::String: V database in fasta format.\nassignments::String: VDJ assignments, in MiAIRR format.\nout::String: Write CHMMAIRRa output.\n\nOptional arguments:\n\np_threshold::Float64: Posterior threshold of switching templates.\nreceptor::String: Which adaptive immune receptor to run on. Options: IG, TCR.\nnon_chimeric_MiAIRR::Union{String, Nothing}: Write assignments of sequences determined to be non-chimeric.\nchimeric_MiAIRR::Union{String, Nothing}: Write assignments of sequences determined to be chimeric.\nchimeric_alignments::Union{String, Nothing}: Write fasta file containing chimeric sequences and the germline alignments they matched to.\nrecombfreqplot::Union{String, Nothing}: Generate a plot of chimerism per recombination event.\ndetailed::Bool: Include details about the chimeras. Calculate viterbi path for chimeric sequences to add startingpoint, recombinations, recombinationsdegapped, and maxlog_prob columns.\nsubsample::Union{Int, Nothing}: Randomly subsample this number of unique sequences to run chimera detection on.\ndeduplicate::Bool: Deduplicate input by vsequencealignment column. Note that if subsample is given, subsampling is performed after deduplication.\nchunk_size::Union{Int, Nothing}: Read in chunk-size lines of input at a time, write output for that chunk, and move to the next chunk.\nHMM_parameters::Union{String, Nothing}: If set, use HMM parameters in this file. Overrides parameter presets from –receptor. See src/IGparameters.tsv or src/TCRparameters.tsv for formatting.\nmafft::Union{String, Nothing}: Path to MAFFT binary. If not provided, will try to find mafft in system PATH.\nalign_database::Bool: Whether to align database V sequences.\ncount_chimeric_segments::Bool: If set, count the number of times chimeric segments appear in nonchimeric sequences. Adds chimerasegmentcounts column.\ntrim::Int: How many nucleotides to trim off the inner ends of the chimeric segments when searching for matches in nonchimeric sequences.\nseed::Int: Seed to use for random subsampling.\nquiet::Bool: If set, hides non-error messages.\ndisable_internal_dedup::Bool: If set, does not deduplicate sequences internally. This was added for timing benchmark purposes and has no effect on the output.\n\n\n\n\n\n","category":"method"},{"location":"api/#CHMMAIRRa.julia_main-Tuple{}","page":"API","title":"CHMMAIRRa.julia_main","text":"julia_main()::Cint\n\nEntry point C call for main function for Julia\n\n\n\n\n\n","category":"method"},{"location":"#CHMMAIRRa.jl","page":"Overview","title":"CHMMAIRRa.jl","text":"","category":"section"},{"location":"","page":"Overview","title":"Overview","text":"(Image: MIT license) (Image: Coverage)","category":"page"},{"location":"","page":"Overview","title":"Overview","text":"CHMMAIRRa.jl is a method for reference-based detection of chimeric adaptive immune receptor repertoires (AIRR) sequences.","category":"page"},{"location":"","page":"Overview","title":"Overview","text":"Chimera detection is performed using a Hidden Markov Model (HMM) that models sequences as being generated from a single reference with mutations or from multiple references with mutations. The core HMM functionality is implemented in the CHMMera.jl package, while CHMMAIRRa.jl provides an AIRR-specific wrapper for the HMM.","category":"page"},{"location":"#Installation","page":"Overview","title":"Installation","text":"","category":"section"},{"location":"","page":"Overview","title":"Overview","text":"The tool is primarily distributed as a an executable for Linux and MacOS from the releases page. Once you have downloaded the executable, you will need to unpack it and add it to your PATH to allow it to be run from any directory in the terminal.","category":"page"},{"location":"","page":"Overview","title":"Overview","text":"tar -xvf CHMMAIRRa.tar.gz # unpack the executable\nexport PATH=$PATH:$(pwd)/CHMMAIRRa/bin/ # add the executable to your PATH","category":"page"},{"location":"","page":"Overview","title":"Overview","text":"The tool is also available from the Julia package manager, and can be installed with the following command:","category":"page"},{"location":"","page":"Overview","title":"Overview","text":"using Pkg; Pkg.add(\"CHMMAIRRa\")","category":"page"},{"location":"","page":"Overview","title":"Overview","text":"See the Julia usage section for more details.","category":"page"},{"location":"#Usage","page":"Overview","title":"Usage","text":"","category":"section"},{"location":"","page":"Overview","title":"Overview","text":"For chimera detection with a minimal output, run the following command:","category":"page"},{"location":"","page":"Overview","title":"Overview","text":"CHMMAIRRa --receptor IG --V-fasta V.fasta --assignments AIRR.tsv --out CHMMAIRRa_out.tsv","category":"page"},{"location":"","page":"Overview","title":"Overview","text":"For chimera detection with a detailed output including recombination path information, run the following command:","category":"page"},{"location":"","page":"Overview","title":"Overview","text":"CHMMAIRRa --receptor IG --V-fasta V.fasta --assignments AIRR.tsv --out CHMMAIRRa_out.tsv --detailed","category":"page"},{"location":"","page":"Overview","title":"Overview","text":"If you want the output to include all columns from the input file, use the --non-chimeric-MiAIRR and --chimeric-MiAIRR arguments:","category":"page"},{"location":"","page":"Overview","title":"Overview","text":"CHMMAIRRa --receptor IG --V-fasta V.fasta --assignments AIRR.tsv --out CHMMAIRRa_out.tsv --non-chimeric-MiAIRR nonchimeric.tsv --chimeric-MiAIRR chimeric.tsv","category":"page"},{"location":"#Mandatory-arguments","page":"Overview","title":"Mandatory arguments","text":"","category":"section"},{"location":"#Inputs","page":"Overview","title":"Inputs","text":"","category":"section"},{"location":"","page":"Overview","title":"Overview","text":"--receptor : The receptor type, either IG or TCR.\n--V-fasta : The database V sequences in FASTA format.\n--assignments : The VDJ assignments in MiAIRR format. Must have been assigned to the --V-fasta database.","category":"page"},{"location":"#Outputs","page":"Overview","title":"Outputs","text":"","category":"section"},{"location":"","page":"Overview","title":"Overview","text":"--out : The output file path. See CHMMAIRRa output columns for more details.","category":"page"},{"location":"#Optional-arguments","page":"Overview","title":"Optional arguments","text":"","category":"section"},{"location":"#Output-options","page":"Overview","title":"Output options","text":"","category":"section"},{"location":"","page":"Overview","title":"Overview","text":"--detailed : Whether to output detailed information about the recombination path. See CHMMAIRRa output columns for more details.\n--non-chimeric-MiAIRR : The output file path for non-chimeric sequences, including all columns from the input file.\n--chimeric-MiAIRR : The output file path for chimeric sequences, including all columns from the input file.","category":"page"},{"location":"#Extra-outputs","page":"Overview","title":"Extra outputs","text":"","category":"section"},{"location":"","page":"Overview","title":"Overview","text":"--recombfreqplot : A plot of chimerism per recombination event, normalized by allele frequency. Can be useful for identifying missing database sequences.\n--chimeric-alignments : Fasta file containing chimeric sequences and the germline alignments they matched to. This also provides a hamming_distances column in the chimeric assignments file.","category":"page"},{"location":"#Technical-parameters","page":"Overview","title":"Technical parameters","text":"","category":"section"},{"location":"","page":"Overview","title":"Overview","text":"--p-threshold : Posterior threshold of switching templates (default: 0.95).\n--HMM-parameters : If set, use HMM parameters in this file. Overrides parameter presets from –receptor. See src/IGparameters.tsv or src/TCRparameters.tsv for formatting.\n--mafft : The path to the MAFFT executable. Automatically found if not provided.\n--align-database : Whether to align the database V sequences. Will be done automatically if database V sequences are not the same length.\n--subsample : Randomly subsample this number of unique sequences to run chimera detection on.\n--deduplicate : Deduplicate input by vsequencealignment column. Note that if subsample is given, subsampling is performed after deduplication.\n--chunk-size : Read in chunk-size lines of input at a time, write output for that chunk, and move to the next chunk. This decreases memory consumption in exchange for increasing runtime.\n--count-chimeric-segments : Whether to count the number of times chimeric segments appear in nonchimeric sequences. Adds chimera_segment_counts column.\n--trim : How many nucleotides to trim off the inner ends of the chimeric segments when searching for matches in nonchimeric sequences.\n--seed : Seed to use for random subsampling.\n--quiet : Whether to suppress non-error messages.","category":"page"},{"location":"#CHMMAIRRa-output-columns","page":"Overview","title":"CHMMAIRRa output columns","text":"","category":"section"},{"location":"#Mandatory-columns","page":"Overview","title":"Mandatory columns","text":"","category":"section"},{"location":"","page":"Overview","title":"Overview","text":"threaded : The query V sequence threaded onto the database alignment.\nchimera_probability : The probability of the V sequence being a chimera.\nchimeric : Whether the V sequence is chimeric.","category":"page"},{"location":"#Optional-columns-generated-when-using-the-detailed-arg","page":"Overview","title":"Optional columns generated when using the --detailed arg","text":"","category":"section"},{"location":"","page":"Overview","title":"Overview","text":"startingpoint : The nonchimeric starting state of the viterbi path (i.e. the HMM's v_call).\nrecombinations : A list of recombination points according to the viterbi path of the query. Format: (vcall 1, vcall 2, recombination position in threaded query sequence).\nrecombinations_degapped : Same as recombinations, but with the non-gapped position of the recombination.\npathevaluation : The probability of the second likeliest reference in the HMM. This gives us a measure of confidence for the recombination path itself.\nnrecombinations : The number of recombination points in the viterbi path.","category":"page"},{"location":"#Optional-outputs-generated-when-using-the-count-chimeric-segments-arg","page":"Overview","title":"Optional outputs generated when using the --count-chimeric-segments arg","text":"","category":"section"},{"location":"","page":"Overview","title":"Overview","text":"segment_1_count : The number of times the left chimeric segment appears in nonchimeric sequences.\nsegment_2_count : The number of times the right chimeric segment appears in nonchimeric sequences.","category":"page"},{"location":"#Julia-usage","page":"Overview","title":"Julia usage","text":"","category":"section"},{"location":"","page":"Overview","title":"Overview","text":"CHMMAIRRa.jl's Julia interface provides two functions: detect_chimeras and detect_chimeras_from_files; the first takes paths to the input files and writes to a file, while the second takes a DataFrame and returns a DataFrame.","category":"page"},{"location":"","page":"Overview","title":"Overview","text":"using CHMMAIRRa\n\n# run from DataFrame\nrefnames, refseqs = read_fasta(\"V.fasta\")\nassignments = CSV.read(\"MiAIRR.tsv\", DataFrame, delim = '\\t')\nchmmairra_output = detect_chimeras(refnames, refseqs, assignments, receptor = \"IG\")\n\n# run from file paths\ndetect_chimeras_from_files(\"V.fasta\", \"MiAIRR.tsv\", \"CHMMAIRRa_out.tsv\", receptor = \"IG\")","category":"page"},{"location":"#Cite","page":"Overview","title":"Cite","text":"","category":"section"},{"location":"","page":"Overview","title":"Overview","text":"If you use CHMMAIRRa.jl in your work, please cite the following paper:","category":"page"},{"location":"","page":"Overview","title":"Overview","text":"","category":"page"}]
}
