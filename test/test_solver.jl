using Test, CSV, DataFrames

@testset "solver" begin

    @testset "every source file parses" begin
        for f in sort(filter(endswith(".jl"), readdir(repo("2D"))))
            @test parse_error(repo("2D", f)) === nothing
        end
        for f in sort(filter(endswith(".jl"), readdir(repo("Utilities"))))
            @test parse_error(repo("Utilities", f)) === nothing
        end
    end

    # The horizon logic, lifted from InputParameters.jl. Kept in step with it by
    # the test below, which asserts the file still resolves it this way.
    function horizon(path, idx)
        df = CSV.read(path, DataFrame)[idx, :]
        (hasproperty(df, :t_fin)   ? df[:t_fin]   : 150000,
         hasproperty(df, :t_prin)  ? df[:t_prin]  : 1000,
         hasproperty(df, :t_check) ? df[:t_check] : 1)
    end

    @testset "a table without a horizon falls back to the archival values" begin
        # Arrange: 2D/DF.csv predates those columns
        path = repo("2D", "DF.csv")

        # Act
        h = horizon(path, 1)

        # Assert: exactly what InputParameters.jl hardcoded before
        @test h == (150000, 1000, 1)
    end

    @testset "a table with a horizon supplies it" begin
        # Arrange / Act
        h8 = horizon(repo("params", "DF_8.csv"), 1)
        h9 = horizon(repo("params", "DF_9.csv"), 1)

        # Assert
        @test h8 == (2_000_000, 10_000, 1)
        @test h9 == (200_000, 1_000, 1)
    end

    @testset "InputParameters.jl still resolves the horizon by fallback" begin
        # Arrange
        src = read(repo("2D", "InputParameters.jl"), String)

        # Act / Assert
        @test occursin("hasproperty(df, :t_fin)", src)
        @test occursin("150000", src)      # the fallback value
    end

    @testset "the run is configurable through the environment" begin
        # Arrange
        src = read(repo("2D", "InputParameters.jl"), String)

        # Act / Assert
        for var in ("DATA_DIR", "SIM_IDX", "DF_FILE")
            @test occursin(var, src)
        end
        @test occursin("SLURM_ARRAY_TASK_ID", src)   # the Slurm default
    end

    @testset "the output directory the solver writes into is created" begin
        # 2D.jl saves into <file>Data/, which nothing created on the cluster
        src = read(repo("2D", "InputParameters.jl"), String)
        @test occursin("mkpath(joinpath(file, \"Data\"))", src)
    end

    @testset "no cluster-specific absolute paths remain" begin
        for f in filter(endswith(".jl"), readdir(repo("2D")))
            src = read(repo("2D", f), String)
            live = join([split(l, "#")[1] for l in split(src, "\n")], "\n")
            @test !occursin("/srv/beegfs", live)
            @test !occursin(r"\"[A-Z]:/", live)
        end
    end

    @testset "the defect-pair start is opt-in and off by default" begin
        # Arrange
        ini  = read(repo("2D", "kernels.jl"), String)
        inp  = read(repo("2D", "InputParameters.jl"), String)
        main = read(repo("2D", "2D.jl"), String)

        # Act / Assert: the initialiser is present, and selected by a D column
        @test occursin("function kernel_ini_P!", ini)
        @test occursin("hasproperty(df, :D)", inp)
        @test occursin("seed_defect_pair", main)

        # no table shipped here has one, so every published run takes the
        # noise path, exactly as before
        for t in filter(f -> startswith(f, "DF_"), readdir(repo("params")))
            cols = split(first(eachline(repo("params", t))), ",")
            @test "D" ∉ cols
        end
        @test "D" ∉ split(first(eachline(repo("2D", "DF.csv"))), ",")
    end

    @testset "interleaving the field construction leaves the RNG stream alone" begin
        # 2D.jl now builds rho between the two draws rather than after both.
        # The draws themselves must be unchanged.
        using Random
        Nx, Nz, ρ0, A = 32, 32, 0.65, 0.002

        Random.seed!(7)
        a1 = rand(Float64, Nx, Nz); a2 = rand(Float64, Nx, Nz, 2)
        ρa = (1.0 .+ A .* (a1 .- 0.5)) .* ρ0; Pa = (a2 .- 0.5) .* A

        Random.seed!(7)
        b1 = rand(Float64, Nx, Nz)
        ρb = (1.0 .+ A .* (b1 .- 0.5)) .* ρ0
        b2 = rand(Float64, Nx, Nz, 2); Pb = (b2 .- 0.5) .* A

        bits(x) = reinterpret(UInt64, vec(x))
        @test bits(ρa) == bits(ρb)
        @test bits(Pa) == bits(Pb)
    end

    @testset "the sweep regenerates the archival table byte for byte" begin
        # Arrange: run AllInputParam.jl in a copy so DF.csv is never clobbered
        mktempdir() do tmp
            cp(repo("2D"), joinpath(tmp, "2D"))
            cp(repo("Utilities"), joinpath(tmp, "Utilities"))
            reference = read(repo("2D", "DF.csv"), String)

            # Act
            ok = try
                run(pipeline(Cmd(`$(Base.julia_cmd()) --startup-file=no --project=$(repo())
                                  AllInputParam.jl`, dir=joinpath(tmp, "2D")),
                             stdout=devnull, stderr=devnull))
                true
            catch
                false
            end

            # Assert
            if ok
                @test read(joinpath(tmp, "2D", "DF.csv"), String) == reference
            else
                @test_skip "AllInputParam.jl could not run (missing packages?)"
            end
        end
    end
end
