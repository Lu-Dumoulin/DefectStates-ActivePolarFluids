const ar = 4.0*4/3
const χ = 0.1
const kd = 0.2
const κ = 1e-4
const Γ = 1.0
const γ = 1e-4

A(ρ0) = 3*ar*ρ0^2-5/2*χ
a(q, ρ0, ζ) = -( kd + q^2*( ρ0^2*( (A(ρ0)-3*ρ0*ζ)/(1+2*q^2) + κ/Γ ) + γ*A(ρ0)) )
b(q, ρ0) = q^2*χ*ρ0*( ρ0^2/(1+2*q^2) + γ )
c(ρ0) = χ*ρ0/Γ
d(q, ρ0) = -ρ0^2/Γ * ( 2*χ+κ*q^2 )
ζc_simp(ρ0) = A(ρ0)/(3*ρ0) + ( sqrt(A(ρ0)*γ/(3*ρ0^2)) - sqrt(2*kd/(3*ρ0^3)) )^2 

tρ0 = Array(0.3:0.01:1.2);
tζ = Array(0.01:0.01:40);
res = zeros(length(tρ0),3);

for iρ0 in eachindex(tρ0)
    @show ρ0 = tρ0[iρ0]
    res[iρ0,1] = ρ0
    res[iρ0,2] = ζc_simp(ρ0)

    ζc = 0.0
    for ζ in tζ

        for q=0:100
            dρ = 1e-3
            dp = 1e-3
            steps = 0
            while dρ < 5e-3 && dp < 5e-3 && steps < 1000
                temp = dρ
                dρ += 0.01*( a(q,ρ0,ζ)*temp + b(q,ρ0)*dp )
                dp += 0.01*( c(ρ0)*temp + d(q,ρ0)*dp )
                steps+=1
            end
            if dρ > 1e-3 || dp > 1e-3
                ζc = ζ
                break
            end
        end

        if ζ==ζc
            res[iρ0,3] = ζc
            @show ρ0, ζc
            break
        end

    end
end

f, ax, p = lines(res[:,1], res[:,2])
lines!(ax, res[:,1], res[:,3])
f