
abstract type EField end

# """
# Incident wave:
# \$\$
# E_{inc,j} = E_0 \exp(-i k r_j - i \omega t)
# $$
# """

# https://discourse.julialang.org/t/struct-of-arrays-soa-vs-array-of-structs-aos/30015/16


# polarized:
# left and right hand circular polarized VS the orthonormal Cartesian vectors ex and ey






struct PlaneWave <: EField
    E0
    e
    k
end

field(r, E::PlaneWave) = E.E0 * exp(-im * dot(E.k, r))
# field(r::Matrix, E::PlaneWave) = [field(ri, E) for ri in eachcol(r)]

function field(r::AbstractArray, e::PlaneWave)
    out = Array{ComplexF64}(undef, size(r))
    for i in eachindex(r)
        out[i] = field(r[i], e)
    end
end

# E_inc(r, E::PlaneWave) = E.E0 * exp(-im * E.k)
# E_inc_td(t, r, E::PlaneWave, omega) = field(E) * exp(im * omega * t)
#
# k = [0 0 1]
# # u = [0,0,1] polarisibility
# E_inc(r, E_0, k) = E_0 * exp(-im * dot(k,r))
#




function E_inc(E0, kvec, r)
    Ei = zeros(ComplexF64, 3, length(r))

    for (i, ri) in enumerate(r)
        Ei[:,i] = E0 .* exp.(im * dot(kvec, ri))
    end

    return reshape(Ei,:)
end




