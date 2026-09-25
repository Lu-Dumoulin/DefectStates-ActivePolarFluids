using Test

@testset "linear stability analysis" begin

    @testset "both sources parse" begin
        for f in ("LSA.jl", "linear_stability.jl")
            @test parse_error(repo("LSA", f)) === nothing
        end
    end

    @testset "it needs nothing but base Julia" begin
        # No `using`, no include: it must run with a bare julia.
        src = read(repo("LSA", "linear_stability.jl"), String)
        live = join([split(l, "#")[1] for l in split(src, "\n")], "\n")
        @test !occursin(r"^\s*using ", live)
        @test !occursin("include(", live)
        @test !occursin("using_pkg", live)
    end

    @testset "it regenerates figure 5's data byte for byte" begin
        mktempdir() do out
            # Act: run it exactly as a reader would, with a bare julia
            run(pipeline(`$(Base.julia_cmd()) --startup-file=no
                          $(repo("LSA", "linear_stability.jl")) $out`,
                         stdout=devnull, stderr=devnull))

            # Assert
            for f in ("tau1.csv", "tau5.csv")
                @test isfile(joinpath(out, f))
                @test read(joinpath(out, f), String) ==
                      read(repo("figures", "Fig5", f), String)
            end
        end
    end

    @testset "its default output is figure 5's directory" begin
        src = read(repo("LSA", "linear_stability.jl"), String)
        @test occursin("\"figures\", \"Fig5\"", src)
    end
end
