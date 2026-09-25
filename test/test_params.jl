using Test, CSV, DataFrames

@testset "parameter tables" begin

    tables = sort(filter(f -> startswith(f, "DF_") && endswith(f, ".csv"),
                         readdir(repo("params"))))

    @testset "the expected tables are present" begin
        # Arrange
        expected = ["DF_1.csv", "DF_4.csv", "DF_6.csv", "DF_8.csv", "DF_9.csv", "DF_11.csv", "DF_12.csv"]

        # Act / Assert
        @test sort(expected) == tables
    end

    @testset "every table has the columns the solver reads" begin
        # Arrange
        required = [:fn, :L, :rho0, :kd, :rhocr, :ap, :nu1, :gamma, :kp, :zetap,
                    :zetap2, :zetarho, :ar, :xi, :M, :seed, :dx, :dz, :dtmin,
                    :t_fin, :t_prin, :t_check]

        # Act / Assert
        for t in tables
            cols = propertynames(CSV.read(repo("params", t), DataFrame))
            @test isempty(setdiff(required, cols))
        end
    end

    @testset "a is tied to the activity as Table I states" begin
        # a' = 4 zeta_rho'/3, and 4/3 when zeta_rho = 0
        for t in tables
            df = CSV.read(repo("params", t), DataFrame)
            for r in eachrow(df)
                expected = r.zetarho == 0 ? 4/3 : abs(r.zetarho)*4/3
                @test r.ar ≈ expected
            end
        end
    end

    @testset "the ar column holds the value the run uses" begin
        # The sweep, the per-figure tables and the solver all associate it the
        # same way now, so the column is not merely informational: the solver
        # reads it. See params/README.md.
        for t in tables
            df = CSV.read(repo("params", t), DataFrame)
            for r in eachrow(df)
                @test r.ar == (r.zetarho == 0 ? 4/3 : abs(r.zetarho)*4/3)
            end
        end

        src = read(repo("2D", "InputParameters.jl"), String)
        @test occursin("const ar::Float64 = Tf(df[:ar])", src)
        @test !occursin("abs(ζρ)*4/3", src)     # no longer recomputed
    end

    @testset "Table I constants hold in every row" begin
        # Arrange
        fixed = Dict(:ap => 0.1, :nu1 => 0.0, :gamma => 1.0, :kp => 1e-4,
                     :zetap => 0.0, :zetap2 => 0.0, :xi => 1.0, :M => 1e-4,
                     :rhocr => 0.0, :dx => 0.01, :dz => 0.01, :dtmin => 0.01)

        # Act / Assert
        for t in tables
            df = CSV.read(repo("params", t), DataFrame)
            for (col, val) in fixed
                @test all(df[!, col] .≈ val)
            end
        end
    end

    @testset "the renewal rate is positive and t_check is 1" begin
        for t in tables
            df = CSV.read(repo("params", t), DataFrame)
            @test all(df.kd .> 0)
            @test all(df.t_check .== 1)
            @test all(df.t_prin .> 0)
            @test all(df.t_fin .>= df.t_prin)
        end
    end

    @testset "fn numbers the rows from 1" begin
        for t in tables
            df = CSV.read(repo("params", t), DataFrame)
            @test df.fn == collect(1:nrow(df))
        end
    end

    # --- the two tables that can be checked against ones actually used --------

    # ar is derived from zetarho, and Analysis/PDpaper/DF.csv predates the fix
    # that made the column match what the solver uses, so compare the inputs.
    physical(df) = sort([(r.rho0, r.kd, r.zetarho) for r in eachrow(df)])

    @testset "DF_9 is the sweep the runs used" begin
        # The table those runs used is no longer in the repository, so the grid
        # is pinned here instead. It is a full factorial over these three.
        df = CSV.read(repo("params", "DF_9.csv"), DataFrame)

        @test sort(unique(df.rho0)) == collect(0.4:0.1:1.5)
        @test sort(unique(df.zetarho)) == [0.0, 1, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24]
        @test sort(unique(round.(df.kd, digits=10))) ==
              [0.01, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 2.0, 5.0, 10.0]
        @test nrow(df) == 12*14*10 == 1680
        @test length(unique(zip(df.rho0, df.zetarho, df.kd))) == nrow(df)
    end

    @testset "DF_11 reproduces the archival L=50 table" begin
        # Arrange
        ours  = CSV.read(repo("params", "DF_11.csv"), DataFrame)
        theirs = CSV.read(repo("2D", "DF.csv"), DataFrame)

        # Act / Assert
        @test nrow(ours) == nrow(theirs)
        @test physical(ours) == physical(theirs)
    end

    @testset "DF_4 has one row per point of Fig4's heatmap" begin
        # Arrange
        df = CSV.read(repo("params", "DF_4.csv"), DataFrame)
        ndef = CSV.read(repo("figures", "Fig4", "Ndef.csv"), DataFrame)

        # Act / Assert
        @test nrow(df) == nrow(ndef) == 1008
        @test sort(unique(df.rho0)) == collect(0.4:0.1:1.5)
        @test sort(unique(df.zetarho)) == [0.0, 1, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24]
        @test sort(unique(round.(df.kd, digits=10))) == [0.01, 0.1, 0.2, 1.0, 5.0, 10.0]
        @test length(unique(zip(df.rho0, df.zetarho, df.kd))) == nrow(df)
    end

    @testset "DF_9 covers the grid Fig9 analysed" begin
        # Arrange
        df = CSV.read(repo("params", "DF_9.csv"), DataFrame)
        densdef = CSV.read(repo("figures", "Fig9", "g34_densdef.csv"), DataFrame)

        # Act / Assert
        @test nrow(df) == nrow(densdef) == 1680
        @test length(unique(df.kd)) == 10
    end

    @testset "Fig1's densities are the ones its caption names" begin
        # Arrange
        df = CSV.read(repo("params", "DF_1.csv"), DataFrame)

        # Act / Assert
        @test sort(df.rho0) == [0.45, 0.6, 0.75]
        @test all(df.L .== 50)
        @test all(df.kd .≈ 1/5)        # tau = 5
        @test all(df.zetarho .== 4)
    end

    @testset "Fig6 is the density list the caption implies" begin
        # Arrange
        df = CSV.read(repo("params", "DF_6.csv"), DataFrame)

        # Act / Assert
        @test sort(df.rho0) == [0.4, 0.5, 0.6, 0.65, 0.7, 0.75, 0.8, 1.2]
        @test all(df.L .== 10)
        @test all(df.kd .≈ 1/5)          # tau = 5
        @test all(df.zetarho .== 4)

        # Fig11's caption says rho0 = 0.45 and 0.55 appear there but not here
        @test 0.45 ∉ df.rho0
        @test 0.55 ∉ df.rho0
        f11 = CSV.read(repo("params", "DF_11.csv"), DataFrame)
        @test 0.45 ∈ f11.rho0
        @test 0.55 ∈ f11.rho0
    end

    @testset "Fig8 and Fig12 describe the same two solutions" begin
        # Arrange
        a = CSV.read(repo("params", "DF_8.csv"), DataFrame)
        b = CSV.read(repo("params", "DF_12.csv"), DataFrame)

        # Act / Assert
        @test physical(a) == physical(b)
        @test nrow(a) == 2
        @test all(a.t_fin .== 2_000_000)
    end
end
