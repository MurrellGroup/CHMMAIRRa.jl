using CHMMAIRRa, Test, DataFrames, CSV

function compare_tsvs(file1::String, file2::String; columns::String = "all")
    df1 = DataFrame(CSV.File(file1, delim = "\t", header = false))
    df2 = DataFrame(CSV.File(file2, delim = "\t", header = false))
    if columns == "intersection"
        columns = intersect(names(df1), names(df2))
    elseif columns == "all"
        @test names(df1) == names(df2)
        columns = names(df1)
    end
    for col in columns
        if typeof(df1[!, col]) == Vector{<:AbstractFloat}
            @test all(isapprox.(df1[!, col], df2[!, col]))
        else
            @test all(isequal.(df1[!, col], df2[!, col]))
        end
    end
end

@testset "CHMMAIRRa.jl" begin
    @testset "utils.jl" begin

        v_sequence_alignments = ["ACGT"]
        v_germline_alignments = ["A-GT"]
        v_calls = ["ref1"]
        refseqs = ["A-GT-"]
        degapped_refseqs = ["AGT"]
        refname2ind = Dict("ref1" => 1)
        ali_length = 5

        result = CHMMAIRRa.thread_all(v_sequence_alignments, v_germline_alignments, v_calls, refseqs, degapped_refseqs, refname2ind, ali_length)
        @test result == ["A-GT-"]

        params = CHMMAIRRa.parse_parameters("parameters.tsv")
        @test params["method"] == "DB"
        @test params["mutation_probabilities"] == [0.1, 0.2]
        @test params["base_mutation_probability"] == 0.05
        @test params["prior_probability"] == 0.05

        @test CHMMAIRRa.parse_recombs_str("(A, B, 1);(C, D, 2)") == [("A", "B", 1), ("C", "D", 2)]
        @test CHMMAIRRa.hamming_no_N("EFGHNNNN", "QQQQQQQ-") == 5

        @test CHMMAIRRa.count_segment_1("ABC", ["ABC", "ABCD", "ABCE"], [1, 2, 3]) == 6
        @test CHMMAIRRa.count_segment_2("ABC", ["ABC", "QABCD", "ACE"], [1, 2, 3]) == 1

        n, s = CHMMAIRRa.read_fasta("V.fasta")
        @test n == ["ref1", "ref2"]
        @test s == ["ACGTACGTACGT", "ACCACCACCAAT"]

        CHMMAIRRa.write_fasta("test.fasta", ["AGT"], seq_names = ["ref1"])
        n, s = CHMMAIRRa.read_fasta("test.fasta")
        @test n == ["ref1"]
        @test s == ["AGT"]

        io = CHMMAIRRa.clear_and_open_for_append("test.fasta")
        CHMMAIRRa.append_to_fasta(io, ["ref1"], ["AGT"])
        close(io)
        n, s = CHMMAIRRa.read_fasta("test.fasta")
        @test n == ["ref1"]
        @test s == ["AGT"]

        df = DataFrame(col1 = [1, 2])
        type_map = Dict("col1" => Float64)
        @test typeof(CHMMAIRRa.convert_column_types!(df, type_map).col1) == Vector{Float64}

        it, col_inds = CHMMAIRRa.setup_eachline_iterator("assignments.tsv", CHMMAIRRa.REQUIRED_COLUMNS)
        df = CHMMAIRRa.get_next_chunk(it, col_inds, CHMMAIRRa.REQUIRED_COLUMNS, 10, CHMMAIRRa.REQUIRED_COLUMNS_TYPES)
        @test df == DataFrame(sequence_id = ["seq1", "seq2"], v_call = ["ref2", "ref2"], v_sequence_alignment = ["ACGTACACCAGGAT", "ACCACCACCAGT"], v_germline_alignment = ["ACCACCACCA--AT", "ACCACCACCAAT"])

        aligned_names, aligned_seqs = CHMMAIRRa.mafft_wrapper(["ACGT", "CCCACGTGG"], ["ref1", "ref2"])
        @test aligned_names == ["ref1", "ref2"]
        @test aligned_seqs == ["---ACGT--", "CCCACGTGG"]

        @test CHMMAIRRa.degap("AB--DDE") == "ABDDE"
        @test CHMMAIRRa.ngaps("AB--DDE", 3) == 1
        @test CHMMAIRRa.is_gz_path("file.gz")
    end

    @testset "CHMMAIRRa.jl" begin
        CHMMAIRRa_out = CHMMAIRRa.detect_chimeras_from_files("V.fasta", "assignments.tsv", "CHMMAIRRa_out_temp.tsv",
                                    non_chimeric_MiAIRR = "non_chimeric_airr_temp.tsv",
                                    chimeric_MiAIRR = "chimeric_airr_temp.tsv",
                                    chimeric_alignments = "chimeric_alignments_temp.tsv",
                                    recombfreqplot = "recombfreqplot.pdf",
                                    detailed = true,
                                    count_chimeric_segments = true)
        compare_tsvs("CHMMAIRRA_out_temp.tsv", "CHMMAIRRA_out.tsv", columns = "all")

        CHMMAIRRa.detect_chimeras_from_files("V.fasta", "assignments.tsv", "CHMMAIRRa_out_temp.tsv",
                                    non_chimeric_MiAIRR = "non_chimeric_airr_temp.tsv",
                                    chimeric_MiAIRR = "chimeric_airr_temp.tsv",
                                    chimeric_alignments = "chimeric_alignments_temp.tsv",
                                    recombfreqplot = "recombfreqplot.pdf",
                                    detailed = true,
                                    chunk_size = 1)
        # count_chimeric_segments doesn't work with chunk_size, so we compare the columns in common between the files
        compare_tsvs("CHMMAIRRA_out_temp.tsv", "CHMMAIRRA_out.tsv", columns = "intersection")
    end
end
