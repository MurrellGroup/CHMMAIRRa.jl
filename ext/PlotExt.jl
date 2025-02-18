module PlotExt

using CHMMAIRRa
using DataFrames, Plots

function wrap_text(str, width=50)
    arr = []
    start_ind = 1
    while start_ind < length(str)
        push!(arr, str[start_ind:minimum([start_ind + width - 1, length(str)])])
        start_ind += width
    end
    return join(arr, "\n")
end

function plotsjl_plot_chimerism_per_recombination(chimeras_per_recombination::DataFrame; textwrap_width::Int=50, fontsize::Int=8)
    @info "Plotting chimerism per recombination"
    wrapped_labels = [wrap_text(label, textwrap_width) for label in chimeras_per_recombination[!, "recombinations_degapped"]]
    num_newlines = sum([length(findall(x->x=='\n', label)) for label in wrapped_labels])
    m = maximum(chimeras_per_recombination.chimeras_normalized)

    # Set up the plot size and layout
    plot_size = (300 + (textwrap_width * 8), (40 * nrow(chimeras_per_recombination)) + (20 * num_newlines) + 100)

    # Create the bar plot
    p = Plots.bar(
        1:nrow(chimeras_per_recombination),
        chimeras_per_recombination.chimeras_normalized,
        label = "",
        xlabel = "Expected frequency normalized count",
        ylabel = "Recombination",
        yticks = (1:nrow(chimeras_per_recombination), wrapped_labels),
        xlims = (0, m * 1.3),
        ylims = (0, nrow(chimeras_per_recombination) + 1),
        size = plot_size,
        color = :skyblue,
        orientation = :horizontal
    )

    # Add count labels above bars
    for i in 1:nrow(chimeras_per_recombination)
        annotate!(p, chimeras_per_recombination.chimeras_normalized[i] + (m / 50), i, text("n=$(chimeras_per_recombination.n[i])", :black, :left, 8))
    end

    return p
end

function my_savefig(p, filename)
    savefig(p, filename)
end

end