using Test

@testset "analysis" begin

    scripts = sort(filter(endswith(".jl"), readdir(repo("Analysis"))))

    @testset "every script parses" begin
        for f in scripts
            @test parse_error(repo("Analysis", f)) === nothing
        end
    end

    @testset "there is a script for every figure one is claimed for" begin
        for n in SCRIPTED_FIGURES
            @test isfile(repo("Analysis", "make_fig$(n).jl"))
        end
    end

    @testset "figures without a script are the documented two" begin
        # Arrange: 3 is the 1D notebook, 13 has no identified generator
        missing_scripts = [n for n in PAPER_FIGURES
                           if !isfile(repo("Analysis", "make_fig$(n).jl"))]

        # Act / Assert
        @test sort(missing_scripts) == [3, 13]
    end

    @testset "no Windows drive paths survive on the paper path" begin
        # talk_figures.jl and not_in_paper.jl are parked rather than maintained,
        # and keep the paths they were written with.
        parked = ["talk_figures.jl", "not_in_paper.jl"]

        for f in scripts
            f in parked && continue
            live = join([split(l, "#")[1] for l in eachline(repo("Analysis", f))], "\n")
            @test !occursin(r"\"[A-Z]:/", live)
        end
    end

    @testset "no figure script shadows the configured run set" begin
        # A `dir_df = "..."` inside a function silently ignores DATA_DIR.
        for n in SCRIPTED_FIGURES
            live = join([split(l, "#")[1] for l in eachline(repo("Analysis", "make_fig$(n).jl"))], "\n")
            @test !occursin(r"dir_df\s*=\s*\"", live)
            @test !occursin(r"dir_fig\s*=\s*\"", live)
        end
    end

    @testset "output goes through FIG_DIR and input through DATA_DIR" begin
        # Arrange
        cfg = read(repo("Analysis", "config.jl"), String)

        # Act / Assert
        @test occursin("DATA_DIR", cfg)
        @test occursin("FIG_DIR", cfg)
        # and the run set is resolved before the heavy package load, so a
        # missing DATA_DIR fails in a second rather than after an install
        @test findfirst("DATA_DIR", cfg)[1] < findfirst("using_pkg", cfg)[1]
    end

    @testset "config refuses to run without a run set" begin
        # Arrange / Act
        out = IOBuffer()
        ok = try
            withenv("DATA_DIR" => nothing) do
                run(pipeline(`$(Base.julia_cmd()) --startup-file=no -e
                              "include(\"$(repo("Analysis","config.jl"))\")"`,
                             stdout=devnull, stderr=out))
            end
            true
        catch
            false
        end

        # Assert: it must fail, and say what to set
        @test !ok
        @test occursin("DATA_DIR", String(take!(out)))
    end

    @testset "the stability analysis regenerates Fig5's data byte for byte" begin
        # Arrange: this is the one driver that needs no simulation output
        mktempdir() do out
            # Act
            run(pipeline(`$(Base.julia_cmd()) --startup-file=no
                          $(repo("Analysis","make_fig5.jl")) $out`,
                         stdout=devnull, stderr=devnull))

            # Assert
            for f in ("tau1.csv", "tau5.csv")
                @test isfile(joinpath(out, f))
                @test read(joinpath(out, f), String) ==
                      read(repo("figures", "Fig5", f), String)
            end
        end
    end

    @testset "the phase-diagram tables ship with the repository" begin
        for f in ("DF.csv", "DF2.csv", "DF_ph2.csv", "DF2_ph2.csv",
                  "DF_analyse.csv", "DF_tikz_norm_adjusted.csv")
            @test isfile(repo("Analysis", "PDpaper", f))
        end
    end

    @testset "PhaseDiagram finds those tables without setup" begin
        src = read(repo("Analysis", "PhaseDiagram.jl"), String)
        live = join([split(l, "#")[1] for l in split(src, "\n")], "\n")
        @test occursin("PDpaper", live)
        @test !occursin("D:/PDpaper", live)
    end

    @testset "talk and dropped-figure code is kept out of the figure scripts" begin
        @test isfile(repo("Analysis", "talk_figures.jl"))
        @test isfile(repo("Analysis", "not_in_paper.jl"))
        for n in SCRIPTED_FIGURES
            src = read(repo("Analysis", "make_fig$(n).jl"), String)
            @test !occursin("_sout", src)
        end
    end
end
