# Shared helpers. No tests here.

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))

"Path relative to the repository root."
repo(parts...) = joinpath(ROOT, parts...)

"True if `cmd` is on PATH, so optional tests can skip instead of fail."
function have(cmd)
    try
        success(pipeline(`which $cmd`, devnull))
    catch
        false
    end
end

"Parse a Julia source file, returning nothing on success or the exception."
function parse_error(path)
    try
        Meta.parseall(read(path, String))
        nothing
    catch e
        e
    end
end

"Every `\\includegraphics{...}` and `table {...}` target on an uncommented line."
function referenced_files(texpath)
    out = String[]
    for line in eachline(texpath)
        stripped = split(line, "%")[1]          # drop comments
        isempty(strip(stripped)) && continue
        for m in eachmatch(r"\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}", stripped)
            push!(out, m.captures[1])
        end
        for m in eachmatch(r"\{([A-Za-z0-9_./\-]+\.(?:csv|txt|dat))\}", stripped)
            push!(out, m.captures[1])
        end
    end
    unique(out)
end

"The figures the article contains, and where each one's picture lives."
const PAPER_FIGURES = 1:13

"Figures drawn by an Analysis script (3 is the 1D notebook; 13 is unidentified)."
const SCRIPTED_FIGURES = [1, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12]
