using CSV, DataFrames, CodecZlib

#These functions take an IgBlast alignment, and thread it onto the reference alignment, transposing gaps
#This currently discards insertions relative to the sequence that IgBlast says its closest to.
#If other references sequences had similar insertions, then this is discarded signal
#I don't think this is going to make much practical difference though
function thread_seq(v_sequence_alignment::AbstractString,
    v_germline_alignment::AbstractString,
    v_call::AbstractString,
    refseqs::Vector{String},
    degapped_refs::Vector{String},
    refname2ind::Dict{String, Int64},
    ali_length::Int)::String
    dg = degap(v_germline_alignment)
    gapped_full_ref = refseqs[refname2ind[v_call]]
    # Because igblast gives a local alignment, sometimes missing ends:
    matchrange = findfirst(dg, degapped_refs[refname2ind[v_call]])
    # discard insertions relative to the sequence that IgBlast says its closest to
    keep_pos = collect(v_germline_alignment) .!= '-';
    no_inserts = collect(v_sequence_alignment)[keep_pos];
    # add Ns to the start of the sequence to match the reference
    no_inserts = vcat(['N' for i in 1:matchrange[1] - 1], no_inserts)
    threaded = fill('-', ali_length);
    pos = 1
    # thread the sequence onto the reference
    for i in 1:ali_length
        # if we've run out of sequence, add Ns
        if pos <= length(no_inserts)
            # if the reference is not a gap, add the sequence
            if gapped_full_ref[i] != '-'
                threaded[i] = no_inserts[pos]
                pos += 1
            end
        end
    end
    return join(threaded)
end

function thread_all(queries::Vector{String},q2refs::Vector{String},q2ref_names::Vector{String},refseqs::Vector{String},degapped_refs::Vector{String},refname2ind::Dict{String, Int64},ali_length::Int64)::Vector{String}
    threaded = Vector{String}(undef, length(queries))
    Threads.@threads for i in eachindex(queries)
        threaded[i] = thread_seq(queries[i], q2refs[i], q2ref_names[i], refseqs, degapped_refs, refname2ind, ali_length)
    end
    return threaded
end

function add_v_DFR_column!(assignments::DataFrame)::DataFrame
    assignments[!,"v_DFR"] .= length.(assignments.v_germline_alignment) .- Int.(round.(assignments.v_identity ./ 100 .* length.(assignments.v_germline_alignment)))
    return assignments
end

function count_segment_1(segment_1::String, threaded::Vector{String}, n::Vector{Int64})::Int64 return sum(n[startswith.(threaded, segment_1)]) end
function count_segment_2(segment_2::String, threaded::Vector{String}, n::Vector{Int64})::Int64 return sum(n[endswith.(threaded, segment_2)]) end

function add_segment_count_columns(assignments::DataFrame, p_threshold::Float64, trim::Int64)::DataFrame
    # for each chimeric assignment, determine how many rows in the nonchimeric assignments contains chimeric segments 1 (left) and 2 (right)
    # add the counts to segment_1_count and segment_2_count columns
    # only works with single recombination chimeras in order to maintain simple substring search and avoid expensive operations.

    function parse_recombinations(s)
        t = strip.(split(s[2:end-1], ","))
        return (t[1], t[2], parse(Int64, t[3]))
    end

    # only use single recombination chimeras
    assignments[!,"nrecombinations"] = map(el-> (ismissing(el)) | (el == "") ? 0 : count(x->x==';', el) + 1, assignments.recombinations)
    single_recomb_selector = assignments.nrecombinations .== 1
    parsed_recombinations = parse_recombinations.(assignments[single_recomb_selector,"recombinations"])

    # get chimeric and nonchimeric threaded arrays to look for segments
    threaded_ct = combine(groupby(assignments[assignments.chimera_probability .<= p_threshold,:], :threaded), nrow => :n)
    nonchimeric_threaded = [el[1:end - trim] for el in threaded_ct.threaded]
    nonchimeric_n = threaded_ct.n
    chimeric_threaded = assignments[single_recomb_selector,:].threaded

    # extract left and right segments
    segment_1s = [threaded[1:recombination[3] - trim] for (threaded, recombination) in zip(chimeric_threaded, parsed_recombinations)]
    segment_2s = [threaded[recombination[3] + 1 + trim:end - trim] for (threaded, recombination) in zip(chimeric_threaded, parsed_recombinations)]

    # count occurrences of chimeric segments in nonchimeric data
    segment_1_cts = zeros(Int64, sum(single_recomb_selector))
    Threads.@threads for (i, segment_1) in collect(enumerate(segment_1s))
        segment_1_cts[i] = count_segment_1(segment_1, copy(nonchimeric_threaded), copy(nonchimeric_n))
    end

    segment_2_cts = zeros(Int64, sum(single_recomb_selector))
    Threads.@threads for (i, segment_2) in collect(enumerate(segment_2s))
        segment_2_cts[i] = count_segment_2(segment_2, copy(nonchimeric_threaded), copy(nonchimeric_n))
    end

    # add to chimeric_assignments dataframe
    assignments[!,"segment_1_count"] .= fill(-1, nrow(assignments))
    assignments[!,"segment_2_count"] .= fill(-1, nrow(assignments))
    assignments[single_recomb_selector,"segment_1_count"] .= segment_1_cts
    assignments[single_recomb_selector,"segment_2_count"] .= segment_2_cts
    return assignments
end

function mafft_wrapper(seqs::Vector{String}, seqnames::Vector{String}; mafft::Union{String, Nothing} = nothing, threads::Int = Base.Threads.nthreads())
    # Find mafft executable
    mafft_path = if isnothing(mafft)
        Sys.which("mafft")
    else
        mafft
    end

    if isnothing(mafft_path)
        error("Could not find mafft executable. Please ensure mafft is installed and in your PATH, or provide the path using --mafft")
    end

    mktempdir() do mydir
        input_fasta = joinpath(mydir, "sequences.fasta")
        write_fasta(input_fasta, seqs, seq_names = seqnames)
        cmd = `$(mafft_path) --thread $(threads) $(input_fasta)`
        io = IOBuffer()
        run(pipeline(cmd, stdout=io, stderr=devnull))
        stdo = String(take!(io))
        results = split.(split(stdo, ">")[2:end], "\n", limit = 2)
        names = map(x->string(x[1]), results)
        seqs = [uppercase(replace(result[2], "\n" => "")) for result in results ]
        return names, seqs
    end
end

function degap(s::String)::String
    return replace(s, Pair("-", ""))
end

function ngaps(s::String, i::Int64)::Int
    return count(j->(j=='-'), s[1:i])
end

function parse_parameters(parameter_file::String)::Dict
    parameter_df = DataFrame(CSV.File(parameter_file, delim = "\t", header = false))
    parameters = Dict()
    parameters["method"] = parameter_df[parameter_df[:,1] .== "method",2][1]
    parameters["mutation_probabilities"] = parse.(Float64, split(parameter_df[parameter_df[:,1] .== "mutation_probabilities",2][1], ","))
    parameters["base_mutation_probability"] = parse(Float64, parameter_df[parameter_df[:,1] .== "base_mutation_probability",2][1])
    parameters["prior_probability"] = parse(Float64, parameter_df[parameter_df[:,1] .== "prior_probability",2][1])
    return parameters
end

#################### Functions related to chimera alignments ####################
function parse_recombs_str(recombs_str::String)
    if recombs_str == ""
        return []
    end
    recombs_arr = [split(recomb_str[2:end-1], ", ") for recomb_str in  split(recombs_str, ";")]
    recombs_arr = [(recomb_arr[1], recomb_arr[2], parse(Int64, recomb_arr[3])) for recomb_arr in recombs_arr]
    recombs_arr = sort(recombs_arr, by = x -> x[3])
    return recombs_arr
end

function get_chimeric_alignments(chimeric_assignments::DataFrame, refnames::Vector{String}, refseqs::Vector{String})
    query_names = chimeric_assignments.sequence_id
    chimeric_threaded = chimeric_assignments.threaded
    chimeric_recombs = parse_recombs_str.(chimeric_assignments.recombinations)
    chimeric_alignments = Vector(undef, nrow(chimeric_assignments))
    Threads.@threads for i in 1:nrow(chimeric_assignments)
        chimeric_alignments[i] = get_chimeric_alignment(query_names[i], chimeric_threaded[i], refnames, refseqs, chimeric_recombs[i])
    end
    return chimeric_alignments
end

function get_chimeric_alignment(query_name::AbstractString, chimeric_threaded::String, refnames::Vector{String}, refseqs::Vector{String}, recombs::Vector)
    if length(recombs) == 0
        return Vector{String}(undef, 0),Vector{String}(undef, 0), ""
    end
    ali_length = length(chimeric_threaded)
    sequences, sequence_ids, prev_split  = Vector{String}(undef, 0), Vector{String}(undef, 0), 0
    for recomb in recombs
        indices = (prev_split + 1, recomb[3])
        sequences = push!(sequences, "-"^(prev_split)*chimeric_threaded[indices[1]:indices[2]]*"-"^(ali_length-recomb[3]))
        sequence_ids = push!(sequence_ids, string(query_name, "_", indices[1], "to", indices[2]))
        prev_split = recomb[3]
    end
    sequences = push!(sequences, "-"^(prev_split)*chimeric_threaded[prev_split + 1:end])
    sequence_ids = push!(sequence_ids, string(query_name, "_", prev_split + 1, "toEnd"))
    sequence_ids = push!(sequence_ids, [recombs[1][1], recombs[1][2]]..., [recombs[i][2] for i in 2:lastindex(recombs)]...)
    sequences = push!(sequences, [refseqs[findfirst(x->x == recombs[1][1], refnames)], refseqs[findfirst(x->x == recombs[1][2], refnames)]]...,
                                [refseqs[findfirst(x->x == recombs[i][2], refnames)] for i in 2:lastindex(recombs)]...)

    # hamming distance between each sequence section and its corresponding reference
    prev_split = 0
    hamming_distances = Vector{Int64}(undef, 0)
    for i in 1:lastindex(recombs)
        indices = (prev_split + 1, recombs[i][3])
        d = hamming_no_N(sequences[i][indices[1]:indices[2]], sequences[i + length(recombs) + 1][indices[1]:indices[2]])
        push!(hamming_distances, d)
        sequence_ids[i] = "$(sequence_ids[i])_ham$(d)"
        prev_split = recombs[i][3]
    end
    # don't count distance from N padding at the end
    last_sequence = sequences[length(recombs) + 1][prev_split + 1:end]
    reversed_last_sequence = reverse(last_sequence)
    Npadding = 0
    while reversed_last_sequence[Npadding + 1] == 'N'
        Npadding += 1
    end

    d = hamming_no_N(last_sequence[1:length(last_sequence) - Npadding], sequences[end][prev_split + 1:length(sequences[end]) - Npadding])
    push!(hamming_distances, d)
    sequence_ids[length(recombs) + 1] = "$(sequence_ids[length(recombs) + 1])_ham$(d)"
    return sequence_ids, sequences, join(hamming_distances, ",")
end

# Normal hamming distance, but Ns have a distance of 0 to anything to avoid counting the padding
function hamming_no_N(s1::String, s2::String)::Int64
    @assert length(s1) == length(s2)
    return sum([c1 != c2 for (c1, c2) in zip(s1, s2) if (! ((c1 == 'N') & (c2 != '-'))) & (! ((c2 == 'N') & (c1 != '-')))])
end

#################### fasta file handling ####################
function read_fasta(file_path::String)
    sequence_ids, sequences = Vector{String}(undef, 0), Vector{String}(undef, 0)
    current_sequence = Vector{String}(undef, 0)

    # Open the FASTA file
    io = open(file_path, "r")

    for line in eachline(io)
        if startswith(line, ">")
            # If a sequence is already being read, store it
            if !isempty(current_sequence)
                push!(sequences, join(current_sequence))
                current_sequence = Vector{String}(undef, 0)
            end
            # Extract the header (remove the ">" character)
            push!(sequence_ids, strip(line[2:end]))
        else
            # Append the sequence line to the current sequence
            push!(current_sequence, strip(line))
        end
    end

    # Store the last sequence
    if !isempty(current_sequence)
        push!(sequences, join(current_sequence))
    end

    close(io)
    return sequence_ids, sequences
end

function write_fasta(filepath::String, sequences::Vector{String}; seq_names = nothing)
    if seq_names === nothing
        seq_names = ["$(i)" for i = 1:length(sequences)]
    end
    io = open(filepath, "w")
    append_to_fasta(io, seq_names, sequences)
    close(io)
end

function append_to_fasta(io, sequence_ids::Vector{String}, sequences::Vector{String})
    if length(sequence_ids) != length(sequences)
        throw(ArgumentError("Number of headers must be equal to the number of sequences"))
    end

    for (sequence_id, sequence) in zip(sequence_ids, sequences)
        println(io, ">$sequence_id")
        write(io, sequence)
        println(io)  # Add a newline after each sequence
    end
end


#################### file handling ####################
function clear_and_open_for_append(filepath::String)
    io = open(filepath, "w")
    close(io)
    if is_gz_path(filepath)
        return GzipCompressorStream(open(filepath, "a"))
    else
        return open(filepath, "a")
    end
end

function is_gz_path(filepath::String)::Bool
    return filepath[end-1:end] == "gz" ? true : false
end


#################### chunked iterator for reading in large tsv files ####################

function get_columns(filepath::String)::Vector{String}
    # list columns for uncompressed or gzip files
    if filepath[end - 2:end] == ".gz"
        f = GzipDecompressorStream(open(filepath))
        return String.(split(readline(f), "\t"))
    else
        f = open(filepath)
        return String.(split(String(readline(f)), "\t"))
    end
end

function get_column_indices(filepath::String, required_columns::Vector{String})::Vector{Int64}
    file_columns = get_columns(filepath)
    return [findfirst(x->x == required_column, file_columns) for required_column in required_columns]
end

function setup_eachline_iterator(filepath::String, columns::Vector{String})
    if is_gz_path(filepath)
        it = eachline(GzipDecompressorStream(open(filepath)))
    else
        it = eachline(open(filepath))
    end
    col_inds = get_column_indices(filepath, columns)
    colnames = Iterators.iterate(it)
    return it, col_inds
end

function get_next_chunk(it::Base.EachLine, col_inds::Vector{Int}, col_names::Vector{String}, chunk_size::Int, type_map::Dict{String, DataType})::DataFrame
    # using an eachline iterator, get the next dataframe chunk of size chunk_size
    assignments = collect(Iterators.take(it, chunk_size))
    split_assignments = split.(assignments, '\t', limit = maximum(col_inds) + 1)
    indexed_assignments = [getindex.(split_assignments, col_ind) for col_ind in col_inds]
    df = DataFrame()
    for (i, col_name) in enumerate(col_names)
        df[!,col_name] = indexed_assignments[i]
    end
    convert_column_types!(df, type_map)
    return df
end

function convert_column_types!(df::DataFrame, type_map::Dict{String, DataType})
    for (col, new_type) in type_map
        if hasproperty(df, col)
            old_type = eltype(df[!,col])
            if old_type != new_type
                try
                    # First try direct conversion (handles numeric type conversions like Int -> Float64)
                    df[!, col] = convert.(Vector{new_type}, df[!, col])
                catch e1
                    try
                        # If direct conversion fails, try parsing (handles String -> numeric conversions)
                        if new_type <: AbstractString
                            df[!, col] = String.(df[!, col])
                        elseif (new_type <: Number) && (old_type <: AbstractString || old_type == String)
                            df[!, col] = parse.(new_type, df[!, col])
                        else
                            # For other cases, try element-wise conversion
                            df[!, col] = [convert(new_type, x) for x in df[!, col]]
                        end
                    catch e2
                        @warn "Failed to convert column $col from $old_type to $new_type. Error 1: $e1, Error 2: $e2"
                        # Keep the original column unchanged
                    end
                end
            end
        else
            @warn "Column $col not found in the DataFrame"
        end
    end
    return df
end

function get_chimerism_per_recombination(CHMMAIRRa_out::DataFrame; top_n::Int = 10)
    function recombined_alleles_freqprod(recombs, allele2freq, default_freq)
        parsed_recombs = parse_recombinations_string(recombs)
        p = get(allele2freq, parsed_recombs[1].left_allele, default_freq)
        for el in parsed_recombs
            p *= get(allele2freq, parsed_recombs[1].right_allele, default_freq)
        end
        return p
    end

    function get_frequencies(df::DataFrame, col::String)
        cts = combine(groupby(df, col), nrow => :n)
        cts[!,"frequency"] .= cts[!,"n"] ./ sum(cts[!,"n"])
        return cts
    end

    allele_freqs = get_frequencies(CHMMAIRRa_out, "v_call")
    median_freq = sort(allele_freqs[!, "frequency"])[maximum([floor(Int, nrow(allele_freqs) / 2), 1])]
    allele2freq = Dict(allele_freqs[!, "v_call"] => allele_freqs[!, "frequency"])
    chimeras_by_allele = get_frequencies(CHMMAIRRa_out[(.! ismissing.(CHMMAIRRa_out.recombinations_degapped)) .& (CHMMAIRRa_out.recombinations_degapped .!= ""),:], "recombinations_degapped")

    sort!(chimeras_by_allele, :n, rev = true)
    chimeras_by_allele = chimeras_by_allele[1:minimum([top_n, nrow(chimeras_by_allele)]),:]



    chimeras_by_allele[!,"recombined_alleles_freqprod"] = [recombined_alleles_freqprod(recombs, allele2freq, median_freq) for recombs in chimeras_by_allele[!, "recombinations_degapped"]]
    chimeras_by_allele[!,"chimeras_normalized"] .= chimeras_by_allele.n ./ chimeras_by_allele[!,"recombined_alleles_freqprod"]
    sort!(chimeras_by_allele, :chimeras_normalized, rev = true)
    return chimeras_by_allele
end

function parse_recombinations_string(recomb_str)
    return [parse_recombination_string(recomb) for recomb in split(recomb_str, ";")]
end

function parse_recombination_string(recomb)
    recomb = recomb[2:end - 1]
    recomb_split = split(recomb, ", ")
    return (left_allele = recomb_split[1], right_allele = recomb_split[2], left_pos_degapped = parse(Int, recomb_split[3]), right_pos_degapped = parse(Int, recomb_split[4]))
end


#################### Simple logging functions ####################
function log_info(msg::String, quiet::Bool = false)
    if ! quiet
        println(msg)
    end
end


