using Test

# Optional. The submitted figure PDFs are not part of the repository; when a
# copy of the resubmission directory is available, every figure this repository
# builds must reproduce its PDF exactly. Point SUBMITTED_DIR at it:
#
#     SUBMITTED_DIR=~/Documents/submittedFiles julia --project=. test/runtests.jl

@testset "against the submitted article" begin

    dir = get(ENV, "SUBMITTED_DIR", "")
    if isempty(dir) || !isdir(expanduser(dir))
        @test_skip "SUBMITTED_DIR not set; skipping comparison with the article"
    elseif !have("pdflatex") || !have("gs")
        @test_skip "needs pdflatex and ghostscript"
    else
        dir = expanduser(dir)

        "Ink bounding box, in points, as (width, height)."
        function bbox(pdf)
            out = read(pipeline(`gs -q -dNOPAUSE -dBATCH -sDEVICE=bbox $pdf`,
                                stderr=`cat`), String)
            m = match(r"%%HiResBoundingBox:\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)", out)
            m === nothing && return nothing
            v = parse.(Float64, m.captures)
            (round(v[3]-v[1], digits=2), round(v[4]-v[2], digits=2))
        end

        mktempdir() do tmp
            for n in PAPER_FIGURES
                d = "Fig$(n)"
                isdir(repo("figures", d)) || continue
                submitted = joinpath(dir, "Dumoulin_etal_Fig$(n).pdf")
                isfile(submitted) || continue

                # Arrange
                work = joinpath(tmp, d)
                cp(repo("figures", d), work)

                # Act
                try
                    run(pipeline(Cmd(`pdflatex -interaction=nonstopmode -halt-on-error main.tex`,
                                     dir=work), stdout=devnull, stderr=devnull))
                catch
                end

                # Assert
                built = joinpath(work, "main.pdf")
                if isfile(built)
                    @test bbox(built) == bbox(submitted)
                else
                    @test false  # failed to build
                end
            end
        end
    end
end
