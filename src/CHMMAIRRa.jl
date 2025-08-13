module CHMMAIRRa

using CHMMera, ArgParse, CSV, DataFrames, Random

include("utils.jl")
include("plot.jl")

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin

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
        "--p-threshold"
            help = "Posterior threshold of switching templates."
            arg_type = Float64
            default = 0.95
        "--J-fasta"
            help = "J database in fasta format. If provided, will run on both V and J and pick the maximum probability of chimerism."
            arg_type = String
            default = nothing
        "--min-DFR"
            help = "Minimum differences from reference (DFR) threshold to consider a recombination event."
            arg_type = Int
            default = 1
        "--HMM-parameters"
            help = "If set, use HMM parameters in this file. Overrides parameter presets from --receptor. See src/IG_parameters.tsv or src/TCR_parameters.tsv for formatting."
            arg_type = String
            default = nothing
        "--mafft"
            help = "Path to MAFFT binary. If not provided, will try to find mafft in system PATH."
            arg_type = String
            default = nothing
        "--align-database"
            help = "Whether to align database sequences."
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
REQUIRED_COLUMNS_TYPES = Dict("sequence_id" => String, "v_call" => String, "v_sequence_alignment" => String, "v_germline_alignment" => String,
                            "j_call" => String, "j_sequence_alignment" => String, "j_germline_alignment" => String)

"""
    detect_chimeras_from_files(V_fasta::String, assignments::String, out::String;

            receptor::String = "IG",
            J_fasta::Union{String, Nothing} = nothing,

            non_chimeric_MiAIRR::Union{String, Nothing} = nothing,
            chimeric_MiAIRR::Union{String, Nothing} = nothing,
            chimeric_alignments::Union{String, Nothing} = nothing,
            recombfreqplot::Union{String, Nothing} = nothing,
            detailed::Bool = false,
            subsample::Union{Int, Nothing} = nothing,
            deduplicate::Bool = false,
            chunk_size::Union{Int, Nothing} = nothing,
           
            p_threshold::Float64 = 0.95,
            min_DFR::Int = 1,
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
- `receptor::String`: Which adaptive immune receptor to run on. Options: IG, TCR.
- `J_fasta::Union{String, Nothing}`: J database in fasta format. If provided, will run on both V and J and the sequence will be chimeric if either V or J is chimeric.
- `non_chimeric_MiAIRR::Union{String, Nothing}`: Write assignments of sequences determined to be non-chimeric.
- `chimeric_MiAIRR::Union{String, Nothing}`: Write assignments of sequences determined to be chimeric.
- `chimeric_alignments::Union{String, Nothing}`: Write fasta file containing chimeric sequences and the germline alignments they matched to.
- `recombfreqplot::Union{String, Nothing}`: Generate a plot of chimerism per recombination event.
- `detailed::Bool`: Include details about the chimeras. Calculate viterbi path for chimeric sequences to add startingpoint, recombinations, recombinations_degapped, and max_log_prob columns.
- `subsample::Union{Int, Nothing}`: Randomly subsample this number of unique sequences to run chimera detection on.
- `deduplicate::Bool`: Deduplicate input by v_sequence_alignment column. Note that if subsample is given, subsampling is performed after deduplication.
- `chunk_size::Union{Int, Nothing}`: Read in chunk-size lines of input at a time, write output for that chunk, and move to the next chunk.
- `p_threshold::Float64`: Posterior threshold of switching templates.
- `min_DFR::Int`: Minimum differences from reference (DFR) threshold to consider a recombination event.
- `HMM_parameters::Union{String, Nothing}`: If set, use HMM parameters in this file. Overrides parameter presets from --receptor. See src/IG_parameters.tsv or src/TCR_parameters.tsv for formatting.
- `mafft::Union{String, Nothing}`: Path to MAFFT binary. If not provided, will try to find mafft in system PATH.
- `align_database::Bool`: Whether to align database sequences.
- `count_chimeric_segments::Bool`: If set, count the number of times chimeric segments appear in nonchimeric sequences. Adds chimera_segment_counts column.
- `trim::Int`: How many nucleotides to trim off the inner ends of the chimeric segments when searching for matches in nonchimeric sequences.
- `seed::Int`: Seed to use for random subsampling.
- `quiet::Bool`: If set, hides non-error messages.
- `disable_internal_dedup::Bool`: If set, does not deduplicate sequences internally. This was added for timing benchmark purposes and has no effect on the output.
"""
function detect_chimeras_from_files(V_fasta::String, assignments::String, out::String;

            receptor::String = "IG",
            J_fasta::Union{String, Nothing} = nothing,

            non_chimeric_MiAIRR::Union{String, Nothing} = nothing,
            chimeric_MiAIRR::Union{String, Nothing} = nothing,
            chimeric_alignments::Union{String, Nothing} = nothing,
            recombfreqplot::Union{String, Nothing} = nothing,
            detailed::Bool = false,
            subsample::Union{Int, Nothing} = nothing,
            deduplicate::Bool = false,
            chunk_size::Union{Int, Nothing} = nothing,
            
            p_threshold::Float64 = 0.95,
            min_DFR::Int = 1,
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
    V_refnames, V_refseqs = read_fasta(V_fasta)

    if ! isnothing(J_fasta)
        J_db = read_fasta(J_fasta)
        global REQUIRED_COLUMNS
        REQUIRED_COLUMNS = vcat(REQUIRED_COLUMNS, ["j_call", "j_sequence_alignment", "j_germline_alignment"])
    else
        J_db = nothing
    end

    # if we are planning to write an AIRR file, we need to read in all of the columns. Otherwise, read only what's required
    read_cols = (isnothing(chimeric_MiAIRR)) && isnothing(non_chimeric_MiAIRR) ? REQUIRED_COLUMNS : get_columns(assignments)
    # no chunking
    if isnothing(chunk_size)
        assignments = CSV.read(assignments, DataFrame, select=read_cols, types=REQUIRED_COLUMNS_TYPES, delim='\t')
        chmmairra_output = detect_chimeras((V_refnames, V_refseqs), assignments;
                                                            J_db = J_db,
                                                            p_threshold = p_threshold,
                                                            min_DFR = min_DFR,
                                                            HMM_parameters = HMM_parameters,
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

        if chmmairra_output.recombfreqplot isa Vector{Vector{Char}}
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
            log_info("------------------------------ Processing chunk $(chunk_n) ------------------------------", quiet)

            assignments = get_next_chunk(assignments_iterator, col_inds, read_cols, chunk_size, REQUIRED_COLUMNS_TYPES)

            if size(assignments)[1] == 0 # covers scenario when the input file length is divisible by chunk-size
                break
            end

            chmmairra_output = detect_chimeras((V_refnames, V_refseqs), assignments;
                                                                J_db = J_db,
                                                                p_threshold = p_threshold,
                                                                min_DFR = min_DFR,
                                                                HMM_parameters = HMM_parameters,

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
    log_info("Done.", quiet)
    return chmmairra_output
end

"""
    detect_chimeras(V_db::Tuple{Vector{String}, Vector{String}}, assignments::DataFrame;
                    J_db::Union{Tuple{Vector{String}, Vector{String}}, Nothing} = nothing,
                    p_threshold::Float64 = 0.95,
                    min_DFR::Int = 1,
                    HMM_parameters::Union{Dict, Nothing} = nothing,
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
- `V_db::Tuple{Vector{String}, Vector{String}}`: Tuple containing (reference names, reference sequences) for the V database.
- `assignments::DataFrame`: VDJ assignments, in MiAIRR format.
# Optional arguments:
- `J_db::Union{Tuple{Vector{String}, Vector{String}}, Nothing}`: Optional tuple containing (reference names, reference sequences) for the J database. If provided, will run on both V and J and the sequence will be chimeric if either V or J is chimeric.
- `p_threshold::Float64`: Posterior threshold of switching templates.
- `min_DFR::Int`: Minimum differences from reference (DFR) threshold to consider a recombination event.
- `HMM_parameters::Union{Dict, Nothing}`: If set, use HMM parameters in this dictionary. Overrides parameter presets from receptor.
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
- `align_database::Bool`: Whether to align database sequences.
- `quiet::Bool`: If set, hides non-error messages.
- `disable_internal_dedup::Bool`: If set, does not deduplicate sequences internally. This was added for timing benchmark purposes and has no effect on the output.
"""
function detect_chimeras(V_db::Tuple{Vector{String}, Vector{String}}, assignments::DataFrame;

                    J_db::Union{Tuple{Vector{String}, Vector{String}}, Nothing} = nothing,
                    
                    p_threshold::Float64 = 0.95,
                    min_DFR::Int = 1,
                    HMM_parameters::Union{Dict, Nothing} = nothing,
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

    V_output = detect_chimeras(V_db[1], V_db[2], assignments, 'v';
                            p_threshold = p_threshold,
                            min_DFR = min_DFR,
                            HMM_parameters = HMM_parameters,

                            receptor = receptor,

                            chimeric_alignments = chimeric_alignments,
                            recombfreqplot = recombfreqplot,
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

    if isnothing(J_db)
        V_output.out[!,"chimeric"] .= V_output.out[!,"v_chimeric"]
        return V_output
    else
        VJ_output = (
            out = nothing,
            chimeric_alignments = nothing,
            recombfreqplot = nothing
        )
        J_output = detect_chimeras(J_db[1], J_db[2], assignments[!,["sequence_id", "j_call", "j_sequence_alignment", "j_germline_alignment"]], 'j';
                                p_threshold = p_threshold,
                                min_DFR = min_DFR,
                                HMM_parameters = HMM_parameters,

                                receptor = receptor,

                                chimeric_alignments = chimeric_alignments,
                                recombfreqplot = recombfreqplot,
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

        output_df = hcat(V_output.out, J_output.out[:,setdiff(names(J_output.out), ["sequence_id", "j_call", "j_sequence_alignment", "j_germline_alignment"])])
        output_df[!,"chimeric"] .= (V_output.out[!,"v_chimeric"] .| J_output.out[!,"j_chimeric"])
        VJ_output = merge(VJ_output, (out = output_df,))
        if chimeric_alignments
            chimeric_alignments_tup = (vcat(V_output.chimeric_alignments[1], J_output.chimeric_alignments[1]), 
                                       vcat(V_output.chimeric_alignments[2], J_output.chimeric_alignments[2]))
            VJ_output = merge(VJ_output, (chimeric_alignments = chimeric_alignments_tup,))
        end


        if recombfreqplot    
            if isdefined(Base, :get_extension) && (Base.get_extension(CHMMAIRRa, :PlotExt) !== nothing)
                recombfreqplot_plot = Plots.plot(V_output.recombfreqplot, J_output.recombfreqplot, layout=(2,1))
            else
                recombfreqplot_plot = vcat(V_output.recombfreqplot, J_output.recombfreqplot)
            end
            
            VJ_output = merge(VJ_output, (recombfreqplot = recombfreqplot_plot,))
        end
        return VJ_output
    end
end

"""
    detect_chimeras(refnames::Vector{String}, refseqs::Vector{String}, assignments::DataFrame, gene::Char;
                    p_threshold::Float64 = 0.95,
                    min_DFR::Int = 1,
                    HMM_parameters::Union{Dict, Nothing} = nothing,
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
- `gene::Char`: Which gene to run on. Options: 'v', 'j'.
# Optional arguments:
- `p_threshold::Float64`: Posterior threshold of switching templates.
- `min_DFR::Int`: Minimum differences from reference (DFR) threshold to consider a recombination event.
- `HMM_parameters::Union{Dict, Nothing}`: If set, use HMM parameters in this dictionary. Overrides parameter presets from receptor.
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
- `align_database::Bool`: Whether to align database sequences.
- `quiet::Bool`: If set, hides non-error messages.
- `disable_internal_dedup::Bool`: If set, does not deduplicate sequences internally. This was added for timing benchmark purposes and has no effect on the output.
"""
function detect_chimeras(refnames::Vector{String}, refseqs::Vector{String}, assignments::DataFrame, gene::Char;
                    p_threshold::Float64 = 0.95,
                    min_DFR::Int = 1,
                    HMM_parameters::Union{Dict, Nothing} = nothing,
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

    # using a named tuple to handle multiple optional outputs
    output = (
        out = assignments,
        chimeric_alignments = nothing,
        recombfreqplot = nothing
    )

    # align the database if the length of the sequences is not the same
    if align_database  | (length(unique(length.(refseqs))) != 1)
        log_info("Aligning $(gene) database sequences using mafft", quiet)
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
    to_remove = ismissing.(assignments[!,"$(gene)_call"]) .| ismissing.(assignments[!,"$(gene)_sequence_alignment"]) .| ismissing.(assignments[!,"$(gene)_sequence_alignment"])
    req_str = join(REQUIRED_COLUMNS, ", ")
    log_info("Removing $(sum(to_remove)) rows with missing values in required columns $(req_str)", quiet)
    assignments = assignments[.!to_remove,:]
    
    if deduplicate
        assignments = assignments[.!nonunique(assignments[!,["$(gene)_sequence_alignment"]]),:]
        log_info("Number of unique sequences after deduplication: $(size(assignments)[1])", quiet)
    end

    if (! isnothing(subsample)) && (size(assignments)[1] > subsample)
        assignments = assignments[shuffle(MersenneTwister(seed), 1:size(assignments)[1])[1:subsample],:]
    end

    # for timing benchmarks, we don't want to deduplicate internally
    if disable_internal_dedup
        unique_assignments = copy(assignments[!,["$(gene)_sequence_alignment", "$(gene)_germline_alignment", "$(gene)_call"]])
    else
        unique_assignments = unique(assignments, ["$(gene)_sequence_alignment", "$(gene)_germline_alignment", "$(gene)_call"])[!,["$(gene)_sequence_alignment", "$(gene)_germline_alignment", "$(gene)_call"]]
        log_info("Number sequences before and after internal deduplication: $(nrow(assignments)) -> $(nrow(unique_assignments))", quiet)
    end

    add_DFR_column!(unique_assignments, gene)
    above_DFR_threshold_mask = unique_assignments[!,"$(gene)_DFR"] .>= min_DFR
    log_info("$(sum(.! above_DFR_threshold_mask)) sequences below DFR threshold $(min_DFR) will be assigned chimera probability 0.0", quiet)
    
    log_info("Threading alignments", quiet)
    threaded = get_threaded(unique_assignments, refnames, refseqs, gene = gene, quiet = quiet)

    threaded_above_DFR_threshold = threaded[above_DFR_threshold_mask]
    # we only need to run once per value in threaded vector
    # the strategy is to run everything on unique threaded values, then leftjoin onto assignment table before writing to file
    log_info("Running forward algorithm to get chimerism probabilities", quiet)
    chimera_probs = get_chimera_probabilities(threaded_above_DFR_threshold, refseqs, bw = HMM_parameters["method"] == "BW", mutation_probabilities = HMM_parameters["mutation_probabilities"], base_mutation_probability = HMM_parameters["base_mutation_probability"], prior_probability = HMM_parameters["prior_probability"]);
    # Add threaded sequences and assign chimera probabilities based on DFR threshold
    unique_assignments[!,"$(gene)_threaded"] = threaded
    unique_assignments[!,"$(gene)_chimera_probability"] = zeros(Float64, nrow(unique_assignments))
    unique_assignments[above_DFR_threshold_mask,"$(gene)_chimera_probability"] = chimera_probs
    # collect additional information about the chimeras using viterbi
    if detailed | chimeric_alignments
        log_info("Running viterbi to get chimeric paths", quiet)
        recombination_information = get_recombination_events(unique_assignments[!,"$(gene)_threaded"], refseqs, bw = HMM_parameters["method"] == "BW", mutation_probabilities = HMM_parameters["mutation_probabilities"], base_mutation_probability = HMM_parameters["base_mutation_probability"], prior_probability = HMM_parameters["prior_probability"], detailed = true);
        # extract viterbi related information
        # columns to indicate the positions the recombinations occurred in
        chimeric_recombination_vecs = map(el->el.recombinations, recombination_information)
        chimeric_recombinations = [length(recombination_events) > 0 ? replace(join([(refnames[recomb.left], refnames[recomb.right], recomb.position) for recomb in recombination_events], ";"), "\""=>"") : "" for recombination_events in chimeric_recombination_vecs]
        chimeric_recombinations_degapped = [length(recombination_events) > 0 ? replace(join([(refnames[recomb.left], refnames[recomb.right], recomb.position - ngaps(refseqs[recomb.left], recomb.position), recomb.position - ngaps(refseqs[recomb.right], recomb.position)) for recomb in recombination_events], ";"), "\""=>"") : "" for recombination_events in chimeric_recombination_vecs ]
        # starting point in the viterbi path - also serves as a CHMMera allele "assignment"
        unique_starting_points = map(el->refnames[el.startingpoint], recombination_information)
        unique_pathevaluations = map(el->el.pathevaluation, recombination_information)
        unique_assignments = hcat(unique_assignments, DataFrame(Symbol("$(gene)_startingpoint") => unique_starting_points, 
                                                                Symbol("$(gene)_recombinations") => chimeric_recombinations, 
                                                                Symbol("$(gene)_recombinations_degapped") => chimeric_recombinations_degapped, 
                                                                Symbol("$(gene)_pathevaluation") => unique_pathevaluations))
    end

    # for timing benchmarks, we don't want to deduplicate internally
    if disable_internal_dedup
        assignments = hcat(assignments, unique_assignments[:,setdiff(names(unique_assignments), names(assignments))])
    else
        assignments = leftjoin(assignments, unique_assignments, on = ["$(gene)_sequence_alignment", "$(gene)_germline_alignment", "$(gene)_call"])
    end

    assignments[!,"$(gene)_chimeric"] = assignments[!,"$(gene)_chimera_probability"] .> p_threshold
    nchimeric = sum(assignments[!,"$(gene)_chimeric"])
    percent_chimeric = round(sum(assignments[!,"$(gene)_chimeric"] .> p_threshold)/nrow(assignments) * 100, digits = 3)
    log_info("Number of chimeric sequences: $(nchimeric) ($(percent_chimeric)%)", quiet)
    log_info("Total rows $(nrow(assignments))", quiet)
    # find the number of times chimeric segments are seen in nonchimeric sequences
    # only works for single recombination chimeras
    if count_chimeric_segments
        log_info("Counting chimeric segments", quiet)
        assignments = add_segment_count_columns(assignments, p_threshold, trim, gene)
    end

    if chimeric_alignments
        log_info("Finding chimeric alignments", quiet)
        alignments = get_chimeric_alignments(assignments[assignments[!,"$(gene)_chimera_probability"] .> p_threshold,:], refnames, refseqs, gene)
        if length(alignments) > 0
            chimeric_alignments_obj = collect(Iterators.flatten(map(x -> x[1], alignments))), collect(Iterators.flatten(map(x -> x[2], alignments)))
        else
            chimeric_alignments_obj = (Vector{String}(undef, 0), Vector{String}(undef, 0))
        end
        assignments[!,"$(gene)_hamming_distances"] .= ""
        assignments[assignments[!,"$(gene)_chimera_probability"] .> p_threshold,"$(gene)_hamming_distances"] .= [el[3] for el in alignments]
        output = merge(output, (chimeric_alignments = chimeric_alignments_obj,))
    end

    if recombfreqplot
        chimeras_per_recombination = get_chimerism_per_recombination(assignments, gene)
        if nrow(chimeras_per_recombination) > 0
            p = plot_chimerism_per_recombination(chimeras_per_recombination, gene)
        else
            p = nothing
        end
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

    log_info("Parsed args:", quiet)
    for (arg,val) in parsed_args
        log_info("  $arg  =>  $val", quiet)
    end

    detect_chimeras_from_files(parsed_args["V-fasta"], parsed_args["assignments"], parsed_args["out"];
            receptor = parsed_args["receptor"],

            non_chimeric_MiAIRR = parsed_args["non-chimeric-MiAIRR"],
            chimeric_MiAIRR = parsed_args["chimeric-MiAIRR"],
            chimeric_alignments = parsed_args["chimeric-alignments"],
            recombfreqplot = parsed_args["recombfreqplot"],
            detailed = parsed_args["detailed"],
            subsample = parsed_args["subsample"],
            deduplicate = parsed_args["deduplicate"],
            chunk_size = parsed_args["chunk-size"],

            p_threshold = parsed_args["p-threshold"],
            J_fasta = parsed_args["J-fasta"],
            min_DFR = parsed_args["min-DFR"],
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