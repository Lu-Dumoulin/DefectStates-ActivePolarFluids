# Bitwise comparison of two snapshot directories.
#
#     julia --project=. tools/compare_snapshots.jl <dirA> <dirB>
#
# Compares the raw bits of every dataset, not the values, so -0.0 vs 0.0 and
# any last-ulp difference are both reported. Exits non-zero on any mismatch.

using JLD

bits(x::Array{Float64}) = reinterpret(UInt64, vec(x))

function compare(dirA, dirB)
    a = sort(filter(endswith(".jld"), readdir(dirA)))
    b = sort(filter(endswith(".jld"), readdir(dirB)))
    if isempty(a)
        println("no .jld files in $dirA"); return false
    end
    if a != b
        println("snapshot file lists differ:")
        println("  only in A: ", setdiff(a, b))
        println("  only in B: ", setdiff(b, a))
        return false
    end
    ok = true
    for f in a
        A = load(joinpath(dirA, f))
        B = load(joinpath(dirB, f))
        for k in ("rho", "v", "P")
            if bits(A[k]) == bits(B[k])
                continue
            end
            d = A[k] .- B[k]
            println("MISMATCH  $f  [$k]   max|Δ| = ", maximum(abs, d),
                    "   differing elements = ", count(!iszero, d), " / ", length(d))
            ok = false
        end
    end
    println(ok ? "IDENTICAL: $(length(a)) snapshots, all datasets bit-for-bit equal"
               : "DIFFERENCES FOUND")
    return ok
end

exit(compare(ARGS[1], ARGS[2]) ? 0 : 1)
