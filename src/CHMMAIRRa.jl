module CHMMAIRRa

using CHMMera, ArgParse, CSV, DataFrames, FASTX, Logging, Random, Requires

# declare so we can override it in the extension
include("utils.jl")
include("plot.jl")



function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--p-threshold"
            help = "Posterior threshold of switching templates."
            arg_type = Float64
            default = 0.95

        # output options
        "--non-chimeric-MiAIRR"
            help = "Write assignments of sequences determined to be non-chimeric."
            arg_type = String
            default = nothing
        "--chimeric-MiAIRR"
            help = "Write assignments of sequences determined to be chimeric."
            arg_type = String
            default = nothing
        "--chimeric-alignments"
            help = "Write fasta file containing chimeric sequences and the germline alignments they matched to.
                    This also provides a hamming_distances columns in the chimeric assignments file."
            arg_type = String
            default = nothing
        "--recombfreqplot"
            help = "Generate a plot of chimerism per recombination event."
            arg_type = String
            default = nothing
        "--detailed"
            help = "Include details about the chimeras. Calculate viterbi path for chimeric sequences to add startingpoint, recombinations, recombinations_degapped, and max_log_prob columns."
            action = :store_true
        "--subsample"
            help = "Randomly subsample this number of unique sequences to run chimera detection on."
            arg_type = Int
            default = nothing
        "--deduplicate"
            help = "Deduplicate input by v_sequence_alignment column. Note that if subsample is given, subsampling is performed after deduplication."
            action = :store_true
        "--chunk-size"
            help = "Read in chunk-size lines of input at a time, write output for that chunk, and move to the next chunk.
                    This decreases memory consumption in exchange for increasing runtime."
            arg_type = Int
            default = nothing

        # technical parameters
        "--HMM-parameters"
            help = "If set, use HMM parameters in this file. Overrides parameter presets from --receptor. See src/IG_parameters.tsv or src/TCR_parameters.tsv for formatting."
            arg_type = String
            default = nothing
        "--mafft"
            help = "Path to MAFFT binary. If not provided, will try to find mafft in system PATH."
            arg_type = String
            default = nothing
        "--align-database"
            help = "Whether to align database V sequences."
            action = :store_true
        "--count-chimeric-segments"
            help = "If set, count the number of times chimeric segments appear in nonchimeric sequences. Adds chimera_segment_counts column.
                    Finding the segments of a chimera in high quantities among nonchimeric sequences further validates the chimeric event."
            action = :store_true
        "--trim"
            help = "How many nucleotides to trim off the inner ends of the chimeric segments when searching for matches in nonchimeric sequences."
            arg_type = Int
            default = 5
        "--seed"
            help = "Seed to use for random subsampling."
            arg_type = Int
            default = 123

        "--quiet"
            help = "If set, hides non-error messages."
            action = :store_true

        "--disable-internal-dedup"
            help = "If set, does not deduplicate sequences internally. This was added for timing benchmark purposes and has no effect on the output."
            action = :store_true

        # required arguments
        "--receptor"
            help = "Which adaptive immune receptor to run on. Options: IG, TCR."
            arg_type = String
            required = true
        "--V-fasta"
            help = "V database in fasta format."
            arg_type = String
            required = true
        "--assignments"
            help = "VDJ assignments, in MiAIRR format."
            arg_type = String
            required = true
        "--out"
            help = "Write CHMMAIRRa output."
            arg_type = String
            required = true

        end
    return parse_args(s)
end

# Assignments with missing values in these AIRR columns will get thrown out
REQUIRED_COLUMNS = ["sequence_id", "v_call", "v_sequence_alignment", "v_germline_alignment"]
REQUIRED_COLUMNS_TYPES = Dict("sequence_id" => String, "v_call" => String, "v_sequence_alignment" => String, "v_germline_alignment" => String)

"""
    detect_chimeras_from_files(V_fasta::String, assignments::String, out::String;

            p_threshold::Float64 = 0.95,
            receptor::String = "IG",

            non_chimeric_MiAIRR::Union{String, Nothing} = nothing,
            chimeric_MiAIRR::Union{String, Nothing} = nothing,
            chimeric_alignments::Union{String, Nothing} = nothing,
            recombfreqplot::Union{String, Nothing} = nothing,
            detailed::Bool = false,
            subsample::Union{Int, Nothing} = nothing,
            deduplicate::Bool = false,
            chunk_size::Union{Int, Nothing} = nothing,

            HMM_parameters::Union{String, Nothing} = nothing,
            mafft::Union{String, Nothing} = nothing,
            align_database::Bool = false,
            count_chimeric_segments::Bool = false,
            trim::Int = 5,
            seed::Int = 123,

            quiet::Bool = false,
            disable_internal_dedup::Bool = false
            )
Run CHMMAIRRa from input files to writing of output files, taking care of I/O.

# Arguments:
- `V_fasta::String`: V database in fasta format.
- `assignments::String`: VDJ assignments, in MiAIRR format.
- `out::String`: Write CHMMAIRRa output.
# Optional arguments:
- `p_threshold::Float64`: Posterior threshold of switching templates.
- `receptor::String`: Which adaptive immune receptor to run on. Options: IG, TCR.
- `non_chimeric_MiAIRR::Union{String, Nothing}`: Write assignments of sequences determined to be non-chimeric.
- `chimeric_MiAIRR::Union{String, Nothing}`: Write assignments of sequences determined to be chimeric.
- `chimeric_alignments::Union{String, Nothing}`: Write fasta file containing chimeric sequences and the germline alignments they matched to.
- `recombfreqplot::Union{String, Nothing}`: Generate a plot of chimerism per recombination event.
- `detailed::Bool`: Include details about the chimeras. Calculate viterbi path for chimeric sequences to add startingpoint, recombinations, recombinations_degapped, and max_log_prob columns.
- `subsample::Union{Int, Nothing}`: Randomly subsample this number of unique sequences to run chimera detection on.
- `deduplicate::Bool`: Deduplicate input by v_sequence_alignment column. Note that if subsample is given, subsampling is performed after deduplication.
- `chunk_size::Union{Int, Nothing}`: Read in chunk-size lines of input at a time, write output for that chunk, and move to the next chunk.
- `HMM_parameters::Union{String, Nothing}`: If set, use HMM parameters in this file. Overrides parameter presets from --receptor. See src/IG_parameters.tsv or src/TCR_parameters.tsv for formatting.
- `mafft::Union{String, Nothing}`: Path to MAFFT binary. If not provided, will try to find mafft in system PATH.
- `align_database::Bool`: Whether to align database V sequences.
- `count_chimeric_segments::Bool`: If set, count the number of times chimeric segments appear in nonchimeric sequences. Adds chimera_segment_counts column.
- `trim::Int`: How many nucleotides to trim off the inner ends of the chimeric segments when searching for matches in nonchimeric sequences.
- `seed::Int`: Seed to use for random subsampling.
- `quiet::Bool`: If set, hides non-error messages.
- `disable_internal_dedup::Bool`: If set, does not deduplicate sequences internally. This was added for timing benchmark purposes and has no effect on the output.
"""
function detect_chimeras_from_files(V_fasta::String, assignments::String, out::String;

            p_threshold::Float64 = 0.95,
            receptor::String = "IG",

            non_chimeric_MiAIRR::Union{String, Nothing} = nothing,
            chimeric_MiAIRR::Union{String, Nothing} = nothing,
            chimeric_alignments::Union{String, Nothing} = nothing,
            recombfreqplot::Union{String, Nothing} = nothing,
            detailed::Bool = false,
            subsample::Union{Int, Nothing} = nothing,
            deduplicate::Bool = false,
            chunk_size::Union{Int, Nothing} = nothing,

            HMM_parameters::Union{String, Nothing} = nothing,
            mafft::Union{String, Nothing} = nothing,
            align_database::Bool = false,
            count_chimeric_segments::Bool = false,
            trim::Int = 5,
            seed::Int = 123,

            quiet::Bool = false,
            disable_internal_dedup::Bool = false
            )

    if quiet
        disable_logging(LogLevel(10))
    end

    if (! isnothing(chunk_size)) .&& (! isnothing(subsample))
        error("The chunk-size and subsample arguments are incompatible. Please use one at a time.")
    end

    # you need to run viterbi to get the segments these options depend on
    if (! isnothing(chimeric_alignments)) | (count_chimeric_segments) | (! isnothing(recombfreqplot))
        detailed = true
    end

    # get HMM parameters from a file if provided or from one of the preset options
    if ! isnothing(HMM_parameters)
        HMM_parameters = parse_parameters(HMM_parameters)
    end

    # read in V database
    refnames, refseqs = read_fasta(V_fasta)
    @assert length(refnames) == length(refseqs)

    # if we are planning to write an AIRR file, we need to read in all of the columns. Otherwise, read only what's required
    read_cols = (isnothing(chimeric_MiAIRR)) && isnothing(non_chimeric_MiAIRR) ? REQUIRED_COLUMNS : get_columns(assignments)
    # no chunking
    if isnothing(chunk_size)
        assignments = CSV.read(assignments, DataFrame, select=read_cols, types=REQUIRED_COLUMNS_TYPES, delim='\t')
        chmmairra_output = detect_chimeras(refnames, refseqs, assignments;
                                                            HMM_parameters = HMM_parameters,
                                                            p_threshold = p_threshold,
                                                            receptor = receptor,

                                                            chimeric_alignments = ! isnothing(chimeric_alignments),
                                                            recombfreqplot = ! isnothing(recombfreqplot),
                                                            detailed = detailed,
                                                            subsample = subsample,
                                                            deduplicate = deduplicate,

                                                            count_chimeric_segments = count_chimeric_segments,
                                                            trim = trim,
                                                            seed = seed,

                                                            mafft = mafft,
                                                            align_database = align_database,

                                                            quiet = quiet,
                                                            disable_internal_dedup = disable_internal_dedup
                                                            )
        if ! isnothing(chimeric_alignments)
            write_fasta(chimeric_alignments, chmmairra_output.chimeric_alignments[2], seq_names = chmmairra_output.chimeric_alignments[1])
        end
        write_cols = vcat("sequence_id", setdiff(names(chmmairra_output.out), read_cols))
        CSV.write(out, chmmairra_output.out[!,write_cols], delim = '\t', compress = is_gz_path(out))

        if ! isnothing(chimeric_MiAIRR)
            chimeric_out = chmmairra_output.out[chmmairra_output.out.chimeric,:]
            CSV.write(chimeric_MiAIRR, chimeric_out, delim = '\t', compress = is_gz_path(chimeric_MiAIRR))
        end

        if ! isnothing(non_chimeric_MiAIRR)
            non_chimeric_out = chmmairra_output.out[.! chmmairra_output.out.chimeric,:]
            CSV.write(non_chimeric_MiAIRR, non_chimeric_out, delim = '\t', compress = is_gz_path(non_chimeric_MiAIRR))
        end



        if typeof(chmmairra_output.recombfreqplot) == Vector{Vector{Char}}
            write_recombfreqplot(chmmairra_output.recombfreqplot, recombfreqplot)
        elseif ! isnothing(chmmairra_output.recombfreqplot)
            Base.get_extension(CHMMAIRRa, :PlotExt).my_savefig(chmmairra_output.recombfreqplot, recombfreqplot)
        end

    # read input and write output one chunk at a time
    else
        # delete files where we are trying to write to and open IO streams as needed
        out_s = clear_and_open_for_append(out)
        chimeric_s = isnothing(chimeric_MiAIRR) ? nothing : clear_and_open_for_append(chimeric_MiAIRR)
        non_chimeric_s = isnothing(non_chimeric_MiAIRR) ? nothing : clear_and_open_for_append(non_chimeric_MiAIRR)
        chimeric_alignments_s = isnothing(chimeric_alignments) ? nothing : clear_and_open_for_append(chimeric_alignments)

        assignments_iterator, col_inds = setup_eachline_iterator(assignments, read_cols)

        chunk_n = 1
        while true
            @info "------------------------------ Processing chunk $(chunk_n) ------------------------------"

            assignments = get_next_chunk(assignments_iterator, col_inds, read_cols, chunk_size, REQUIRED_COLUMNS_TYPES)

            chmmairra_output = detect_chimeras(refnames, refseqs, assignments;
                                                                HMM_parameters = HMM_parameters,
                                                                p_threshold = p_threshold,
                                                                receptor = receptor,

                                                                chimeric_alignments = ! isnothing(chimeric_alignments),
                                                                recombfreqplot = false, # recombfreqplot not handled for chunked mode
                                                                detailed = detailed,
                                                                subsample = nothing,
                                                                deduplicate = deduplicate,

                                                                count_chimeric_segments = count_chimeric_segments,
                                                                trim = trim,
                                                                seed = seed,

                                                                mafft = mafft,
                                                                align_database = align_database,

                                                                quiet = quiet,
                                                                disable_internal_dedup = disable_internal_dedup
                                                                )

            if ! isnothing(chimeric_alignments_s)
                append_to_fasta(chimeric_alignments_s, chmmairra_output.chimeric_alignments[1], chmmairra_output.chimeric_alignments[2])
            end

            if chunk_n == 1
                write_cols = vcat("sequence_id", setdiff(names(chmmairra_output.out), read_cols))
                # write column names to output file
                CSV.write(out_s, chmmairra_output.out[[],write_cols], delim = '\t', append = false)
            end

            if size(assignments)[1] == 0 # covers scenario when the input file length is divisible by chunk-size
                break
            end

            CSV.write(out_s, chmmairra_output.out[!,write_cols], delim = '\t', append = true)

            if ! isnothing(chimeric_MiAIRR)
                chimeric_out = chmmairra_output.out[chmmairra_output.out.chimeric,:]
                CSV.write(chimeric_s, chimeric_out, delim = '\t', append = true)
            end

            if ! isnothing(non_chimeric_MiAIRR)
                non_chimeric_out = chmmairra_output.out[.! chmmairra_output.out.chimeric,:]
                CSV.write(non_chimeric_s, non_chimeric_out, delim = '\t', append = true)
            end

            # if we didn't read in chunk-size, it means we reached the end of the file
            if size(assignments)[1] < chunk_size
                break
            end
            chunk_n += 1

        end
        close(out_s)
        isnothing(chimeric_s) ? nothing : close(chimeric_s)
        isnothing(chimeric_alignments_s) ? nothing : close(chimeric_alignments_s)
    end
    @info "Done."
    return chmmairra_output
end

"""
    detect_chimeras(refnames::Vector{String}, refseqs::Vector{String}, assignments::DataFrame;
                    HMM_parameters::Union{Dict, Nothing} = nothing,
                    p_threshold::Float64 = 0.95,
                    receptor::String = "IG",

                    chimeric_alignments::Bool = false,
                    recombfreqplot::Bool = false,
                    detailed::Bool = false,
                    subsample::Union{Int, Nothing} = nothing,
                    deduplicate::Bool = false,

                    count_chimeric_segments::Bool = false,
                    trim::Int = 5,
                    seed::Int = 123,

                    mafft::Union{String, Nothing} = nothing,
                    align_database::Bool = false,

                    quiet::Bool = false,
                    disable_internal_dedup::Bool = false
                    )::NamedTuple{(:out, :chimeric_alignments, :recombfreqplot),
                                Tuple{DataFrame,
                                Union{Nothing, Tuple{Vector{String}, Vector{String}}},
                                Any}
                                }
                    }



Run CHMMAIRRa on parsed arguments, returning output as a Tuple containing relevant outputs.
Output tuple fields:
1. out: the input assignments dataframe with chimerism probabilities.
2. chimeric_alignments: names and sequences for alignments showing which reference sequences recombined at which positions to form detected chimeras. Set to nothing if chimeric_alignments is false.
3. recombfreqplot: barplot of most frequent recombinations, normalized by their expected frequencies. Set to nothing if recombfreqplot is false.

# Arguments:
- `refnames::Vector{String}`: Names of the reference sequences.
- `refseqs::Vector{String}`: Sequences of the reference sequences.
- `assignments::DataFrame`: VDJ assignments, in MiAIRR format.
# Optional arguments:
- `HMM_parameters::Union{Dict, Nothing}`: If set, use HMM parameters in this dictionary. Overrides parameter presets from receptor.
- `p_threshold::Float64`: Posterior threshold of switching templates.
- `receptor::String`: Which adaptive immune receptor to run on. Options: IG, TCR.
- `chimeric_alignments::Bool`: Whether to generate and return chimeric sequences and the germline alignments they matched to.
- `recombfreqplot::Bool`: Whether to generate and return a plot of chimerism per recombination event.
- `detailed::Bool`: Whether to include details about the chimeras. Calculate viterbi path for chimeric sequences to add startingpoint, recombinations, recombinations_degapped, and max_log_prob columns.
- `subsample::Union{Int, Nothing}`: Randomly subsample this number of unique sequences to run chimera detection on.
- `deduplicate::Bool`: Whether to deduplicate input by v_sequence_alignment column. Note that if subsample is given, subsampling is performed after deduplication.
- `count_chimeric_segments::Bool`: If set, count the number of times chimeric segments appear in nonchimeric sequences. Adds chimera_segment_counts column.
- `trim::Int`: How many nucleotides to trim off the inner ends of the chimeric segments when searching for matches in nonchimeric sequences.
- `seed::Int`: Seed to use for random subsampling.
- `mafft::Union{String, Nothing}`: Path to MAFFT binary. If not provided, will try to find mafft in system PATH.
- `align_database::Bool`: Whether to align database V sequences.
- `quiet::Bool`: If set, hides non-error messages.
- `disable_internal_dedup::Bool`: If set, does not deduplicate sequences internally. This was added for timing benchmark purposes and has no effect on the output.
"""
function detect_chimeras(refnames::Vector{String}, refseqs::Vector{String}, assignments::DataFrame;
                    HMM_parameters::Union{Dict, Nothing} = nothing,
                    p_threshold::Float64 = 0.95,
                    receptor::String = "IG",

                    chimeric_alignments::Bool = false,
                    recombfreqplot::Bool = false,
                    detailed::Bool = false,
                    subsample::Union{Int, Nothing} = nothing,
                    deduplicate::Bool = false,

                    count_chimeric_segments::Bool = false,
                    trim::Int = 5,
                    seed::Int = 123,

                    mafft::Union{String, Nothing} = nothing,
                    align_database::Bool = false,

                    quiet::Bool = false,
                    disable_internal_dedup::Bool = false
                    )::NamedTuple{(:out, :chimeric_alignments, :recombfreqplot),
                                Tuple{DataFrame,
                                Union{Nothing, Tuple{Vector{String}, Vector{String}}},
                                Any
                                }
                    }

    if quiet
        disable_logging(LogLevel(10))
    end

    # using a named tuple to handle multiple optional outputs
    output = (
        out = assignments,
        chimeric_alignments = nothing,
        recombfreqplot = nothing
    )

    # align the database if the length of the sequences is not the same
    if align_database  | (length(unique(length.(refseqs))) != 1)
        @info "Aligning V database sequences using mafft"
        refnames, refseqs = mafft_wrapper(refseqs, refnames, mafft = mafft, threads = Base.Threads.nthreads())
    end

    refseqs = uppercase.(refseqs)

    # no HMM_parameters provided, use default parameters for receptor
    if isnothing(HMM_parameters) & ((receptor == "IG") | (receptor == "TCR"))
        HMM_parameters = parse_parameters(joinpath(@__DIR__, "..", "data", receptor * "_parameters.tsv"))
    elseif isnothing(HMM_parameters)
        error(string("Invalid receptor ", receptor, ". Receptor must be either IG or TCR."))
    end

    # check for missing values in required columns
    to_remove = ismissing.(assignments.v_call) .| ismissing.(assignments.v_sequence_alignment) .| ismissing.(assignments.v_germline_alignment)
    req_str = join(REQUIRED_COLUMNS, ", ")
    @info "Removing $(sum(to_remove)) rows with missing values in required columns $(req_str)"
    assignments = assignments[.!to_remove,:]

    if deduplicate
        assignments = assignments[.!nonunique(assignments[!,["v_sequence_alignment"]]),:]
        @info "Number of unique sequences after deduplication: $(size(assignments)[1])"
    end

    if (! isnothing(subsample)) && (size(assignments)[1] > subsample)
        assignments = assignments[shuffle(MersenneTwister(seed), 1:size(assignments)[1])[1:subsample],:]
    end

    queries = uppercase.(assignments[!,"v_sequence_alignment"])
    q2refs = uppercase.(assignments[!,"v_germline_alignment"])
    q2ref_names = String.(map(x->split(x, ",")[1], assignments[!,"v_call"]))
    degapped_refs = degap.(refseqs);
    refname2ind = Dict(zip(refnames,collect(eachindex(refnames))))
    ali_length = maximum(length.(refseqs))
    @assert length.(q2refs) == length.(queries)
    @info "Threading alignment"
    threaded = thread_all(queries, q2refs, q2ref_names, refseqs, degapped_refs, refname2ind, ali_length);

    # for timing benchmarks, we don't want to deduplicate internally
    if disable_internal_dedup
        unique_threaded = threaded
    else
        unique_threaded = unique(threaded)
    end

    # we only need to run once per value in threaded vector
    # the strategy is to run everything on unique threaded values, then leftjoin onto assignment table before writing to file
    @info "Running forward algorithm to get chimerism probabilities"
    unique_chimera_probs = get_chimera_probabilities(unique_threaded, refseqs, bw = HMM_parameters["method"] == "BW", mutation_probabilities = HMM_parameters["mutation_probabilities"], base_mutation_probability = HMM_parameters["base_mutation_probability"], prior_probability = HMM_parameters["prior_probability"]);
    chimera_info = DataFrame(threaded = unique_threaded, chimera_probability = unique_chimera_probs)
    # collect additional information about the chimeras using viterbi
    if detailed | chimeric_alignments
        @info "Running viterbi to get chimeric paths"
        recombination_information = get_recombination_events(unique_threaded, refseqs, bw = HMM_parameters["method"] == "BW", mutation_probabilities = HMM_parameters["mutation_probabilities"], base_mutation_probability = HMM_parameters["base_mutation_probability"], prior_probability = HMM_parameters["prior_probability"], detailed = true);
        # extract viterbi related information
        # columns to indicate the positions the recombinations occurred in
        chimeric_recombination_vecs = map(el->el.recombinations, recombination_information)
        chimeric_recombinations = [length(recombination_events) > 0 ? replace(join([(refnames[recomb.left], refnames[recomb.right], recomb.position) for recomb in recombination_events], ";"), "\""=>"") : "" for recombination_events in chimeric_recombination_vecs]
        chimeric_recombinations_degapped = [length(recombination_events) > 0 ? replace(join([(refnames[recomb.left], refnames[recomb.right], recomb.position - ngaps(refseqs[recomb.left], recomb.position), recomb.position - ngaps(refseqs[recomb.right], recomb.position)) for recomb in recombination_events], ";"), "\""=>"") : "" for recombination_events in chimeric_recombination_vecs ]
        # starting point in the viterbi path - also serves as a CHMMera allele "assignment"
        unique_starting_points = map(el->refnames[el.startingpoint], recombination_information)
        unique_pathevaluations = map(el->el.pathevaluation, recombination_information)
        chimera_info = hcat(chimera_info, DataFrame(startingpoint = unique_starting_points, recombinations = chimeric_recombinations, recombinations_degapped = chimeric_recombinations_degapped, pathevaluation = unique_pathevaluations))
    end

    # map chimera information back onto original sequences
    assignments[!,"threaded"] = threaded

    # for timing benchmarks, we don't want to deduplicate internally
    if disable_internal_dedup
        assignments = hcat(assignments, chimera_info[:,setdiff(names(chimera_info), names(assignments))])
    else
        assignments = leftjoin(assignments, chimera_info, on = :threaded)
    end

    @info "Number of chimeric sequences: $(sum(assignments.chimera_probability .> p_threshold)) ($(round(sum(assignments.chimera_probability .> p_threshold)/nrow(assignments) * 100, digits = 3))%)"
    assignments[!,"chimeric"] = assignments.chimera_probability .> p_threshold
    # find the number of times chimeric segments are seen in nonchimeric sequences
    # only works for single recombination chimeras
    if count_chimeric_segments
        @info "Counting chimeric segments"
        assignments = add_segment_count_columns(assignments, p_threshold, trim)
    end

    if chimeric_alignments
        @info "Finding chimeric alignments"
        alignments = get_chimeric_alignments(assignments[assignments.chimera_probability .> p_threshold,:], refnames, refseqs)
        if length(alignments) > 0
            chimeric_alignments_obj = collect(Iterators.flatten(map(x -> x[1], alignments))), collect(Iterators.flatten(map(x -> x[2], alignments)))
        else
            chimeric_alignments_obj = (Vector{String}(undef, 0), Vector{String}(undef, 0))
        end
        assignments[!,"hamming_distances"] .= ""
        assignments[assignments.chimera_probability .> p_threshold,"hamming_distances"] .= [el[3] for el in alignments]
        output = merge(output, (chimeric_alignments = chimeric_alignments_obj,))
    end

    if recombfreqplot
        chimeras_per_recombination = get_chimerism_per_recombination(assignments)
        p = plot_chimerism_per_recombination(chimeras_per_recombination)
        output = merge(output, (recombfreqplot = p,))
    end

    output = merge(output, (out = assignments,))
    return output
end

"""
    julia_main()::Cint

Entry point C call for main function for Julia
"""
function julia_main()::Cint
    parsed_args = parse_commandline()

    if parsed_args["quiet"]
        disable_logging(LogLevel(10))
    end

    @info "Parsed args:"
    for (arg,val) in parsed_args
        @info "  $arg  =>  $val"
    end

    detect_chimeras_from_files(parsed_args["V-fasta"], parsed_args["assignments"], parsed_args["out"];

            p_threshold = parsed_args["p-threshold"],
            receptor = parsed_args["receptor"],

            non_chimeric_MiAIRR = parsed_args["non-chimeric-MiAIRR"],
            chimeric_MiAIRR = parsed_args["chimeric-MiAIRR"],
            chimeric_alignments = parsed_args["chimeric-alignments"],
            recombfreqplot = parsed_args["recombfreqplot"],
            detailed = parsed_args["detailed"],
            subsample = parsed_args["subsample"],
            deduplicate = parsed_args["deduplicate"],
            chunk_size = parsed_args["chunk-size"],

            HMM_parameters = parsed_args["HMM-parameters"],
            mafft = parsed_args["mafft"],
            align_database = parsed_args["align-database"],
            count_chimeric_segments = parsed_args["count-chimeric-segments"],
            trim = parsed_args["trim"],
            seed = parsed_args["seed"],

            quiet = parsed_args["quiet"],
            disable_internal_dedup = parsed_args["disable-internal-dedup"]
            )
    return 0
end

export detect_chimeras_from_files, detect_chimeras

end