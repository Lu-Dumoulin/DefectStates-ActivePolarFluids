using Test, CSV, DataFrames

@testset "parameter tables" begin

    tables = sort(filter(f -> startswith(f, "DF_") && endswith(f, ".csv"),
                         readdir(repo("params"))))

    @testset "the expected tables are present" begin
        # Arrange
        expected = ["DF_1.csv", "DF_4.csv", "DF_8.csv", "DF_9.csv", "DF_11.csv", "DF_12.csv"]

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

    physical(df) = sort([(r.rho0, r.kd, r.zetarho, r.ar) for r in eachrow(df)])

    @testset "DF_9 reproduces the sweep the phase-diagram runs used" begin
        # Arrange
        ours  = CSV.read(repo("params", "DF_9.csv"), DataFrame)
        theirs = CSV.read(repo("Analysis", "PDpaper", "DF.csv"), DataFrame)

        # Act / Assert
        @test nrow(ours) == 1680
        @test physical(ours) == physical(theirs)
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
        @test length(unique(df.rho0)) == 12
        @test length(unique(df.zetarho)) == 14
        @test length(unique(df.kd)) == 6
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
