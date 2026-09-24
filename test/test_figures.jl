using Test

@testset "figures" begin

    figdirs = ["Fig$(n)" for n in PAPER_FIGURES if isdir(repo("figures", "Fig$(n)"))]

    @testset "the collected figures each have a wrapper" begin
        for d in figdirs
            @test isfile(repo("figures", d, "main.tex"))
        end
    end

    @testset "every figure is laid out the same way" begin
        # main.tex plus a fig_*.tex holding the picture
        for d in figdirs
            body = filter(f -> startswith(f, "fig_") && endswith(f, ".tex"),
                          readdir(repo("figures", d)))
            @test length(body) == 1
        end
    end

    @testset "every referenced image and table is present" begin
        for d in figdirs, tex in filter(endswith(".tex"), readdir(repo("figures", d)))
            for ref in referenced_files(repo("figures", d, tex))
                # sources sit beside the figure, flattened
                @test isfile(repo("figures", d, basename(ref))) ||
                      isfile(repo("figures", d, ref))
            end
        end
    end

    @testset "no figure still points into the article's tree" begin
        for d in figdirs, tex in filter(endswith(".tex"), readdir(repo("figures", d)))
            live = join([split(l, "%")[1] for l in eachline(repo("figures", d, tex))], "\n")
            @test !occursin("Fig/Fig_", live)
            @test !occursin("Fig/Tikzpdf", live)
        end
    end

    @testset "wrapped figures reproduce the article's text block" begin
        # Figures built inside the article inherit its \linewidth and font, so
        # their wrappers must set both. 246pt for a one-column figure, 510pt
        # for the two-column figure*.
        onecol = ["Fig2", "Fig4", "Fig8", "Fig9"]
        twocol = ["Fig6"]

        for d in vcat(onecol, twocol)
            src = read(repo("figures", d, "main.tex"), String)
            width = d in twocol ? "510.0pt" : "246.0pt"
            @test occursin("\\setlength{\\linewidth}{$width}", src)
            @test occursin("\\fontsize{9}{10.5}", src)   # revtex inside a figure
            @test occursin("\\input{abrv}", src)          # \zr and friends
        end
    end

    @testset "self-contained figures are left alone" begin
        # These were built by their own projects and never saw the article's
        # text block; imposing one changes the output.
        for d in ["Fig1", "Fig3", "Fig5", "Fig7", "Fig10"]
            src = read(repo("figures", d, "main.tex"), String)
            @test !occursin("\\fontsize{9}{10.5}", src)
        end
    end

    @testset "Fig3 plots exactly the data the 1D notebook wrote" begin
        for f in readdir(repo("1D", "figdata"))
            a = repo("1D", "figdata", f)
            b = repo("figures", "Fig3", "figdata", f)
            @test isfile(b)
            @test read(a, String) == read(b, String)
        end
    end

    @testset "the build script covers every figure" begin
        src = read(repo("figures", "build_all.sh"), String)
        for d in figdirs
            @test occursin(d, src)
        end
    end

    # --- optional: needs a TeX installation ----------------------------------

    @testset "each figure compiles" begin
        if !have("pdflatex")
            @test_skip "pdflatex not installed"
        else
            mktempdir() do tmp
                for d in figdirs
                    work = joinpath(tmp, d)
                    cp(repo("figures", d), work)
                    ok = try
                        run(pipeline(Cmd(`pdflatex -interaction=nonstopmode -halt-on-error main.tex`,
                                         dir=work), stdout=devnull, stderr=devnull))
                        true
                    catch
                        false
                    end
                    @test ok
                    ok && @test isfile(joinpath(work, "main.pdf"))
                    # a missing macro silently drops a label rather than failing
                    if isfile(joinpath(work, "main.log"))
                        log = read(joinpath(work, "main.log"), String)
                        @test !occursin("Undefined control sequence", log)
                    end
                end
            end
        end
    end
end
