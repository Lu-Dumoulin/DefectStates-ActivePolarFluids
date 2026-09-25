using Test

# Nothing published here should depend on, or disclose, the machine the runs
# were done on: its directory layout means nothing to a reader, and the scratch
# paths carry an account name. Patterns are assembled from fragments so this
# file does not match itself.

@testset "hygiene" begin

    tracked = filter(!isempty, split(read(`git -C $(repo()) ls-files`, String), "\n"))
    # skip the test suite itself, and binary assets
    sources = filter(tracked) do f
        !startswith(f, "test/") &&
        !any(endswith(f, e) for e in (".png", ".pdf", ".jld", ".csv", ".gz"))
    end

    @testset "there are files to check" begin
        @test length(sources) > 20
    end

    @testset "no account name or scratch path" begin
        # Arrange
        account = "dumou" * "lil"
        scratch = "/srv/" * "beegfs"

        # Act / Assert
        for f in sources
            txt = read(repo(f), String)
            @test !occursin(account, txt)
            @test !occursin(scratch, txt)
        end
    end

    @testset "no cluster directory layout" begin
        # Arrange: the run sets lived in directories named like this
        runset = "FFT_" * "2D_P_"

        # Act / Assert
        for f in sources
            @test !occursin(runset, read(repo(f), String))
        end
    end

    @testset "no host or site names" begin
        # Arrange
        host = "baob" * "ab"
        site = "uni" * "ge"

        # Act / Assert
        for f in sources
            txt = lowercase(read(repo(f), String))
            @test !occursin(host, txt)
            @test !occursin(site, txt)
        end
    end

    @testset "no Windows drive paths outside the parked and vendored code" begin
        # talk_figures.jl and not_in_paper.jl are parked rather than
        # maintained. Utilities/ is vendored byte-for-byte from a lab-internal
        # directory, and carries a "F:/Images/" default in a file-copy helper
        # nothing here calls; editing it would break that byte-identity for no
        # gain, and it discloses nothing about the machine that ran anything.
        exempt = ("Analysis/talk_figures.jl", "Analysis/not_in_paper.jl")
        for f in sources
            (f in exempt || startswith(f, "Utilities/")) && continue
            endswith(f, ".jl") || continue
            live = join([split(l, "#")[1] for l in eachline(repo(f))], "\n")
            @test !occursin(r"\"[A-Z]:/", live)
        end
    end
end
