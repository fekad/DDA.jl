"""
# Incident electric field

Incident planewave:
$$ E_{inc}(r) = E_0 e^{i (kr - \omega t)} $$

Arbitrary incident electric field:
- [Jones vector](https://en.m.wikipedia.org/wiki/Jones_calculus) in the lab frame (x-y polarisation plane, z propagation),
- 3D rotation of the coordinate system

Visualisation:
- Polarisation state can be represented as a point on the Bloch (Poincare) sphere.
"""

abstract type Field end



"""
Suppose that a monochromatic plane wave of light is travelling in the positive $z$-direction, with angular frequency $\omega$ and wave vector $\mathbf{k} = (0,0,k)$, where the wavenumber k = ω/c. Then the electric and magnetic fields E and H are orthogonal to k at each point; they both lie in the plane "transverse" to the direction of motion. Furthermore, H is determined from E by 90-degree rotation and a fixed multiplier depending on the wave impedance of the medium. So the polarization of the light can be determined by studying E. The complex amplitude of E is written

$$\begin{pmatrix} E_{x}(t) \\ E_{y}(t) \\ 0 \end{pmatrix} =
\begin{pmatrix} E_{0x} e^{i(kz - \omega t + \phi_{x})} \\ E_{0y} e^{i(kz - \omega t + \phi_{y})} \\ 0 \end{pmatrix} =
\begin{pmatrix} E_{0x} e^{i\phi_{x}} \\ E_{0y} e^{i\phi_{y}} \\ 0 \end{pmatrix} e^{i(kz - \omega t)}.$$

Linear polarized in the x direction (typically called "horizontal"):
$$|H\rangle = \begin{pmatrix}1\\0\end{pmatrix}$$

Linear polarized in the y direction (typically called "vertical"):
$$|V\rangle = \begin{pmatrix}0\\1\end{pmatrix}$$

Linear polarized at 45° from the x axis (typically called "diagonal" L+45):
$$|D\rangle = \frac{1}{\sqrt{2}} \left(|H\rangle +|V\rangle \right) = \frac{1}{\sqrt{2}} \begin{pmatrix}1\\1\end{pmatrix}$$

Linear polarized at −45° from the x axis (typically called "anti-diagonal" L−45):
$$|A\rangle =\frac{1}{\sqrt{2}} \left(|H\rangle - |V\rangle\right) = \frac{1}{\sqrt{2}} \begin{pmatrix}1\\-1\end{pmatrix}$$

Right-hand circular polarized (typically called "RCP" or "RHCP"):
$$|R\rangle =\frac{1}{\sqrt{2}} \left(|H\rangle - i|V\rangle\right) = \frac{1}{\sqrt{2}} \begin{pmatrix}1\\-i\end{pmatrix}$$

Left-hand circular polarized (typically called "LCP" or "LHCP"):
$$|L\rangle =\frac{1}{\sqrt{2}} \left(|H\rangle + i|V\rangle\right) = \frac{1}{\sqrt{2}} \begin{pmatrix}1\\+i\end{pmatrix}$$

Arbitrary direction:
$$\mathbf{E}(\mathbf{r}, t) =
\mathbf{R} \begin{pmatrix} E_{0x} e^{i\phi_{x}} \\ E_{0y} e^{i\phi_{y}} \\ 0 \end{pmatrix} e^{i (\mathbf{R} \mathbf{k}) \cdot \mathbf{r} }e^{-i \omega t}.$$

"""
struct PlaneWave
    k::Float64
    e::SVector{2,ComplexF64}
    θ::Float64
    ϕ::Float64
end

function field(E₀, k , r)
    return E₀ * exp(-im * dot(k, r))
end

function field(E::PlaneWave, r)
    R = RotZY(E.θ, E.ϕ)
    E₀ = R[:,1:2] * E.e
    kvec = R[:,3] * E.k
    return field(E₀, kvec, r)
end

function field(E::PlaneWave, r::Vector{SVector{3, Float64}})
    out = Array{ComplexF64}(undef, (3, length(r)))
    for i in eachindex(r)
        out[:,i] = field(E, r[i])
    end
    return out
end





# """
# Incident wave:
# \$\$
# E_{inc,j} = E_0 \exp(-i k r_j - i \omega t)
# $$
# """
#
# https://discourse.julialang.org/t/struct-of-arrays-soa-vs-array-of-structs-aos/30015/16
#
# # polarization:
# # left and right hand circular polarized VS the orthonormal Cartesian vectors ex and ey
# struct PlaneWave <: Field
#     E₀::SVector{3,Float64}
#     k::SVector{3,Float64}
# end
#
# function field(r, k, E₀)
#     return E₀ * exp(-im * dot(k, r))
# end
#
# function field(r, E::PlaneWave)
#     return field(r, E.E₀, E.k)
# end
#
# function field(r::Vector{SVector{3, Float64}}, e::PlaneWave)
#     out = Array{ComplexF64}(undef, (3, length(r)))
#     for i in eachindex(r)
#         out[:,i] = field(r[i], e)
#     end
#     return out
# end
#
# # E_inc(r, E::PlaneWave) = E.E0 * exp(-im * E.k)
# # E_inc_td(t, r, E::PlaneWave, omega) = field(E) * exp(im * omega * t)
# #
# # k = [0 0 1]
# # # u = [0,0,1] polarisibility
# # E_inc(r, E_0, k) = E_0 * exp(-im * dot(k,r))
# #
#
#
# # function E_inc(E0, kvec, r)
# #     Ei = zeros(ComplexF64, 3, length(r))
# #
# #     for (i, ri) in enumerate(r)
# #         Ei[:,i] = E0 .* exp.(im * dot(kvec, ri))
# #     end
# #
# #     return reshape(Ei, :)
# # end




