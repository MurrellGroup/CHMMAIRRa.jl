using DataFrames

function xval2num_blocks(x_vals::Vector{Float64}; max_num_blocks::Int = 50)
    max_x_value = maximum(x_vals)
    x_vals_normalized = x_vals ./ max_x_value
    num_blocks = floor.(Int, x_vals_normalized .* max_num_blocks)
    return num_blocks
end

function print_unicode_plot(plot::Vector{Vector{Char}})
    for row in plot
        println(join(row))
    end
end

function count2nchars(count::Int)
    return collect(string(" n=", count))
end

function write_recombfreqplot(plot::Vector{Vector{Char}}, out::String)
    open(out, "w") do f
        for row in plot
            println(f, join(row))
        end
    end
end

function plot_chimerism_per_recombination(chimeras_per_recombination::DataFrame, gene::Char)
    if isdefined(Base, :get_extension) && (Base.get_extension(CHMMAIRRa, :PlotExt) !== nothing)
        return Base.get_extension(CHMMAIRRa, :PlotExt).plotsjl_plot_chimerism_per_recombination(chimeras_per_recombination)
    else
        return unicode_plot_chimerism_per_recombination(chimeras_per_recombination, gene)
    end
end

function unicode_plot_chimerism_per_recombination(chimeras_per_recombination::DataFrame, gene::Char)
    xvalues = Vector{Float64}(chimeras_per_recombination[!, "chimeras_normalized"])
    ylabels = Vector{String}(chimeras_per_recombination[!, "$(gene)_recombinations_degapped"])
    counts = Vector{Int}(chimeras_per_recombination[!, "n"])
    return unicode_barplot(xvalues, ylabels, counts)
end

# create a very simple horizontally-oriented unicode barplot
function unicode_barplot(xvalues::Vector{Float64}, ylabels::Vector{String}, counts::Vector{Int}; max_num_blocks::Int = 50, x_axis_label::String = "Chimeras normalized", BLOCK = '█')::Vector{Vector{Char}}
    ## setup numbers we need
    num_blocks = xval2num_blocks(xvalues)
    max_label_length = maximum(length.(ylabels))
    n_labels = count2nchars.(counts)
    # determine plot width. We need to fit labels, y axis, and bars
    plot_width = max_label_length + 4 + max_num_blocks + maximum(length.(n_labels))
    # determine plot height. We need to fit labels, spaces between labels, x axis, and x axis labels
    plot_height = (length(ylabels) * 2) + 3
    # determine x axis label start and end
    x_axis_label_start = max_label_length + 4 + floor(Int, ((max_num_blocks - length(x_axis_label)) / 2))
    x_axis_label_end = x_axis_label_start + length(x_axis_label) - 1
    # set up x axis tickmark
    max_x_value_chars =  collect(string(floor(Int, maximum(xvalues))))
    x_value_tickmark_start = plot_width - length(max_x_value_chars) + 1


    ## setup plot

    plot = [fill(' ', plot_width) for _ in 1:plot_height]

    for i in 1:length(ylabels)
        # Y axis
        plot[(i * 2) - 1][max_label_length + 3] = '|'
        plot[i * 2][max_label_length + 3] = '|'
        # Y axis labels
        plot[i * 2][max_label_length - length(ylabels[i]) + 1:max_label_length] .= collect(ylabels[i])
        # Bars
        plot[i * 2][max_label_length + 4:max_label_length + 4 + num_blocks[i] - 1] .= BLOCK
        # counts
        plot[i * 2][max_label_length + 4 + num_blocks[i]:max_label_length + 4 + num_blocks[i] + length(n_labels[i]) - 1] .= n_labels[i]
    end
    plot[end - 2] .= '_'
    plot[end][x_axis_label_start:x_axis_label_end] .= collect(x_axis_label)

    # print plot
    print_unicode_plot(plot)

    return plot
end


