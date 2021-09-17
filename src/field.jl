
abstract type Field end

# """
# Incident wave:
# \$\$
# E_{inc,j} = E_0 \exp(-i k r_j - i \omega t)
# $$
# """

# https://discourse.julialang.org/t/struct-of-arrays-soa-vs-array-of-structs-aos/30015/16


# polarized:
# left and right hand circular polarized VS the orthonormal Cartesian vectors ex and ey




struct PlaneWave <: Field
    E₀::SVector{3,Float64}
    k::SVector{3,Float64}
end

function field(r, k, E₀)
    return E₀ * exp(-im * dot(k, r))
end

function field(r, E::PlaneWave)
    return field(r, E.E₀, E.k)
end

function field(r::Vector{SVector{3, Float64}}, e::PlaneWave)
    out = Array{ComplexF64}(undef, (3, length(r)))
    for i in eachindex(r)
        out[:,i] = field(r[i], e)
    end
    return out
end

# E_inc(r, E::PlaneWave) = E.E0 * exp(-im * E.k)
# E_inc_td(t, r, E::PlaneWave, omega) = field(E) * exp(im * omega * t)
#
# k = [0 0 1]
# # u = [0,0,1] polarisibility
# E_inc(r, E_0, k) = E_0 * exp(-im * dot(k,r))
#


# function E_inc(E0, kvec, r)
#     Ei = zeros(ComplexF64, 3, length(r))
#
#     for (i, ri) in enumerate(r)
#         Ei[:,i] = E0 .* exp.(im * dot(kvec, ri))
#     end
#
#     return reshape(Ei, :)
# end




