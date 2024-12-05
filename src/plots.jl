using DataFrames, Plots

function get_frequencies(df::DataFrame, col::String)
    cts = combine(groupby(df, col), nrow => :n)
    cts[!,"frequency"] .= cts[!,"n"] ./ sum(cts[!,"n"])
    return cts
end

function wrap_text(str, width=50)
    arr = []
    start_ind = 1
    while start_ind < length(str)
        push!(arr, str[start_ind:minimum([start_ind + width - 1, length(str)])])
        start_ind += width
    end
    return join(arr, "\n")
end

function parse_recombinations_string(recomb_str)
    return [parse_recombination_string(recomb) for recomb in split(recomb_str, ";")]
end

function parse_recombination_string(recomb)
    recomb = recomb[2:end - 1]
    recomb_split = split(recomb, ", ")
    return (left_allele = recomb_split[1], right_allele = recomb_split[2], left_pos_degapped = parse(Int, recomb_split[3]), right_pos_degapped = parse(Int, recomb_split[4]))
end

function get_chimerism_per_recombination(CHMMAIRRa_out::DataFrame; top_n::Int = 10)
    allele_freqs = get_frequencies(CHMMAIRRa_out, "v_call")
    median_freq = sort(allele_freqs[!, "frequency"])[maximum([floor(Int, nrow(allele_freqs) / 2), 1])]
    allele2freq = Dict(allele_freqs[!, "v_call"] => allele_freqs[!, "frequency"])
    chimeras_by_allele = get_frequencies(CHMMAIRRa_out[(.! ismissing.(CHMMAIRRa_out.recombinations_degapped)) .& (CHMMAIRRa_out.recombinations_degapped .!= ""),:], "recombinations_degapped")

    sort!(chimeras_by_allele, :n, rev = true)
    chimeras_by_allele = chimeras_by_allele[1:minimum([top_n, nrow(chimeras_by_allele)]),:]

    function recombined_alleles_freqprod(recombs, allele2freq, default_freq)
        parsed_recombs = parse_recombinations_string(recombs)
        p = get(allele2freq, parsed_recombs[1].left_allele, default_freq)
        for el in parsed_recombs
            p *= get(allele2freq, parsed_recombs[1].right_allele, default_freq)
        end
        return p
    end

    chimeras_by_allele[!,"recombined_alleles_freqprod"] = [recombined_alleles_freqprod(recombs, allele2freq, median_freq) for recombs in chimeras_by_allele[!, "recombinations_degapped"]]
    chimeras_by_allele[!,"chimeras_normalized"] .= chimeras_by_allele.n ./ chimeras_by_allele[!,"recombined_alleles_freqprod"]
    sort!(chimeras_by_allele, :chimeras_normalized, rev = true)
    return chimeras_by_allele
end

function plot_chimerism_per_recombination(chimera_by_recombination::DataFrame; textwrap_width::Int=50, fontsize::Int=8)
    wrapped_labels = [wrap_text(label, textwrap_width) for label in chimera_by_recombination[!, "recombinations_degapped"]]
    num_newlines = sum([length(findall(x->x=='\n', label)) for label in wrapped_labels])
    m = maximum(chimera_by_recombination.chimeras_normalized)

    # Set up the plot size and layout
    plot_size = (300 + (textwrap_width * 8), (40 * nrow(chimera_by_recombination)) + (20 * num_newlines) + 100)

    # Create the bar plot
    p = Plots.bar(
        1:nrow(chimera_by_recombination),
        chimera_by_recombination.chimeras_normalized,
        label = "",
        xlabel = "Expected frequency normalized count",
        ylabel = "Recombination",
        yticks = (1:nrow(chimera_by_recombination), wrapped_labels),
        xlims = (0, m * 1.3),
        ylims = (0, nrow(chimera_by_recombination) + 1),
        size = plot_size,
        color = :skyblue,
        orientation = :horizontal
    )

    # Add count labels above bars
    for i in 1:nrow(chimera_by_recombination)
        annotate!(p, chimera_by_recombination.chimeras_normalized[i] + (m / 50), i, text("n=$(chimera_by_recombination.n[i])", :black, :left, 8))
    end

    return p
end