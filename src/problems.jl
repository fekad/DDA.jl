# 1. load or create the coordinates of the dipoles,
# 2. load or assign the polarizability αj to each dipole,
# 3. calculated the incident field Einc,j at each dipole,
# 4. assemble the interaction matrix A and
# 5. solve for P in the system of linear equations
#

# doi:10.1016/j.jqsrt.2011.03.012
# Knowing Pj, other quantities such as
# - the scattered field,
# - dipole force,
# - Poynting vector,
# - extinction,
# - absorption and
# - scattering cross sections,
# - phase function,
# - Mueller matrix etc.
# can be calculated.



# TODO: allpha is a 3x3 (symmetric) tensor by defualt )
# calc_Ajj(alph) = Diagonal(1 / alph * I, 3)


# """
#     FrequencySimulation([particles::AbstractParticles=[],]
#                         source::AbstractSource)
# Build a FrequencySimulation. If particles are not provided, an empty array is used.
# After building, you can [`run`](@ref) the simulation to get a [`FrequencySimulationResult`](@ref).
# """
# mutable struct FrequencySimulation{Dim,P<:PhysicalMedium} <: Simulation{Dim}
#     "Vector of particles, can be of different types."
#     particles::AbstractParticles
#     "RegularSource wave, where source.medium is the background medium of the simulation."
#     source::AbstractSource{P}
# end

# High level
# mutable struct DipoleProblem <: AbstractProblem
#     "Vector of scatterers, can be of different types."
#     scatterers::Vector{S} where S<:AbstractScatterer
#     "single or as many as scatterers"
#     grid
#     "RegularSource wave, where source.medium is the background medium of the simulation."
#     Einc
# end

# low-level structure (with high level constructors)
struct DipoleProblem <: AbstractProblem
    k
    E0
    dipoles
    alphas
    Eincs
end


# struct DipoleProblem{Dim, T<:AbstractFloat} <: AbstractProblem
#     k::T
#     dipoles::Vector{V} where V<:StaticVector{Dim,T}
#     alphas::Vector{Complex{T}}
#     Eincs::Vector{C} where C<:StaticVector{Dim, <:Complex{T}}
# end


struct DipoleProblemHigh{Dim,T<:AbstractFloat} <: AbstractProblem
    scatterers::Vector{AbstractScatterer} # dipoles + alphas
    Eincs::AbstractIncidentField
end


struct GridProblem2 <: AbstractProblem
    grid
    inds
    alphas
    Eincs
    k
end

struct GridProblem <: AbstractProblem
    grid
    scatterers
    Einc
end


# TODO: move this to solvers + for s in scatterers (merging)
function discretize(p::GridProblem)
    # 1. create the coordinates of the dipoles,
    occ = discretize(p.grid, p.scatterers.target)

    # 2. assign the polarizability αj to each dipole
    alphas = polarisbility(p.scatterers.model, p)

    coords = p.grid[occ]
    return coords, occ, alphas
end



function discretize(p::GridProblem, s::Scatterer)
    # 1. create the coordinates of the dipoles,
    occ = discretize(p.grid, s.target)

    # 2. assign the polarizability αj to each dipole
    alphas = polarisbility(s.model, p)

    return occ, alphas
end



# TODO: fix spacing
function polarisbility(m::LDRModel, prob::GridProblem)
    ε = m.ε
    d = prob.grid.spacing[1]
    kvec, E₀  = prob.Einc.kvec, prob.Einc.E₀
    LDR(ε, d, kvec, E₀)
end



@doc raw"""
The electric field from a radiating electric dipole:
$$
E = \frac{1}{4 \pi \varepsilon_0} \left\{
    k^2(\hat{r} \times p) \times \hat{r} \frac{e^{ikr}}{r} +
    \left[3\hat{r}(\hat{r}\cdot p)-p\right] \left( \frac{1}{r^3}-\frac{ik}{r^2} \right) e^{ikr}
\right\}
$$

$A_{jk}$ is the tensor that represents the interaction between a receiving dipole at $r_j$ and the radiating dipole
at $r_k$.
$$
A_{jk} = \frac{\exp(i k r_{jk})}{r_{jk}} \left[
    k^2(\hat{r}_{jk} \hat{r}_{jk} - 11) + \frac{i k r_{jk} - 1}{r_{jk}^2} (3 \hat{r}_{jk} \hat{r}_{jk} - 11)
\right]
$$
"""
function calc_Ajk(k,  rj, rk)
    r_jk = rj - rk

    r = norm(r_jk)
    rn = r_jk / r

    A_jk = exp(im * k * r) / r * (k^2 * (rn * rn' - I) + (im * k * r - 1) / r^2 * (3 * rn * rn' - I))
    return A_jk
end

calc_Ajj(alph::T) where T<:Number = Diagonal(1 / alph * I, 3)


function interactions(k, r, alph)
    N = length(r)
    out = zeros(ComplexF64, 3, N, 3, N)

    for i in 1:N
        # out[:,i,:,i] = 1/alph[i] * I
        out[:,i,:,i] = calc_Ajj(alph[i])
    end

    for i in 2:N
        for j in 1:i - 1
            out[:,i,:,j] = calc_Ajk(k, r[i], r[j])
        end
    end

    # DDA.calc_Ajk(k, r[2], r[1])
    return Symmetric(reshape(out, 3N, 3N), :L)
end



# TODO: adding parametic types like here: https://github.com/JuliaMatrices/ToeplitzMatrices.jl/blob/master/src/ToeplitzMatrices.jl
# TODO: improve structure: Array of symmetric tensors
"""
    TensorConvolution(...)

A preplanned, circular convolution operator on an M x N matrix of data of type T

# Fields
- `Ĝ`: DFT coefficients of the convolution kernel
- `F`: plan for FFT operator
- `inds`: Cartesian indecies of each dipole
- `alphas`: Polirisability of each dipole
- `tmp`: scratch space to store the zero-pad the input matrix

# Constructors:
- `TensorConvolution(grid::CartesianGrid, occ, k_norm, alphas)
"""
struct TensorConvolution{P<:AbstractFFTs.Plan{ComplexF64}}
    Ĝ::Array{ComplexF64,4}
    F::P
    inds::Vector{CartesianIndex{3}}
    alphas::ComplexF64
    tmp::Array{ComplexF64,4}
end

Base.eltype(::TensorConvolution) = ComplexF64
Base.size(A::TensorConvolution) = (3 * length(A.inds), 3 * length(A.inds))
Base.size(A::TensorConvolution, i::Int) = size(A)[i]


function TensorConvolution(grid::CartesianGrid, occ, k_norm, alphas)

    Nx, Ny, Nz = size(grid)
    Ĝ = Array{ComplexF64}(undef, 6, 2Nx, 2Ny, 2Nz)

    for i = 1:Nx, j = 1:Ny, k = 1:Nz
        Rij = grid[i, j, k]
        A = DDA.calc_Ajk(k_norm, Rij)
        Ĝ[1, i, j, k] = A[1, 1] # xx
        Ĝ[2, i, j, k] = A[1, 2] # xy
        Ĝ[3, i, j, k] = A[1, 3] # xz
        Ĝ[4, i, j, k] = A[2, 2] # yy
        Ĝ[5, i, j, k] = A[2, 3] # yz
        Ĝ[6, i, j, k] = A[3, 3] # zz
    end
    Ĝ[:, 1, 1, 1] .= zero(ComplexF64)

    extend!(Ĝ, Nx, Ny, Nz)
    fft!(Ĝ, 2:4)

    tmp = zeros(ComplexF64, 3, 2Nx, 2Ny, 2Nz)
    F = FFTW.plan_fft!(tmp, 2:4, flags = FFTW.ESTIMATE)

    inds = findall(occ)
    TensorConvolution{typeof(F)}(Ĝ, F, inds, alphas, tmp)
end


A::TensorConvolution * x = mul!(similar(x), A, x)

# TODO: full grid input and output
function LinearAlgebra.mul!(y::StridedVector, A::TensorConvolution, x::StridedVector)
    x = reinterpret(SVector{3,ComplexF64}, x)
    y = reinterpret(SVector{3,ComplexF64}, y)

    fill!(A.tmp, 0)

    for i in eachindex(x)
        A.tmp[:, A.inds[i]] = x[i]
    end

    A.F * A.tmp # FFT

    # !! Matrix vector multiplication in the Fourier space
    _, Nx, Ny, Nz = size(A.Ĝ)
    for k = 1:Nz, j = 1:Ny, i = 1:Nx
        E = A.tmp[:, i, j, k]
        A.tmp[1, i, j, k] = A.Ĝ[1, i, j, k] * E[1] + A.Ĝ[2, i, j, k] * E[2] + A.Ĝ[3, i, j, k] * E[3] # x component
        A.tmp[2, i, j, k] = A.Ĝ[2, i, j, k] * E[1] + A.Ĝ[4, i, j, k] * E[2] + A.Ĝ[5, i, j, k] * E[3] # y component
        A.tmp[3, i, j, k] = A.Ĝ[3, i, j, k] * E[1] + A.Ĝ[5, i, j, k] * E[2] + A.Ĝ[6, i, j, k] * E[3] # z component
    end

    A.F \ A.tmp # inverse FFT

    for i in eachindex(y)
        y[i] = A.tmp[:, A.inds[i]] .+ 1 / A.alphas * x[i]
    end

    return y

end



function extend!(Ĝ::Array{ComplexF64,4}, Nx, Ny, Nz)

    # Periodic extension
    g = [+1 -1 -1 +1 +1 +1  # x
        +1 -1 +1 +1 -1 +1  # y
        +1 +1 -1 +1 -1 +1  # z
        +1 +1 -1 +1 -1 +1  # xy
        +1 -1 +1 +1 -1 +1  # xz
        +1 -1 -1 +1 +1 +1]  # yz

    Ĝ[:, Nx+1, :, :] .= zero(ComplexF64)
    Ĝ[:, :, Ny+1, :] .= zero(ComplexF64)
    Ĝ[:, :, :, Nz+1] .= zero(ComplexF64)

    for i = 1:6
        Ĝ[i, Nx+2:2Nx, 1:Ny, 1:Nz] = Ĝ[i, Nx:-1:2, 1:Ny, 1:Nz] * g[1, i]       # x
        Ĝ[i, 1:Nx, Ny+2:2Ny, 1:Nz] = Ĝ[i, 1:Nx, Ny:-1:2, 1:Nz] * g[2, i]       # y
        Ĝ[i, 1:Nx, 1:Ny, Nz+2:2Nz] = Ĝ[i, 1:Nx, 1:Ny, Nz:-1:2] * g[3, i]       # z
        Ĝ[i, Nx+2:2Nx, Ny+2:2Ny, 1:Nz] = Ĝ[i, Nx:-1:2, Ny:-1:2, 1:Nz] * g[4, i] # xy
        Ĝ[i, Nx+2:2Nx, 1:Ny, Nz+2:2Nz] = Ĝ[i, Nx:-1:2, 1:Ny, Nz:-1:2] * g[5, i] # xz
        Ĝ[i, 1:Nx, Ny+2:2Ny, Nz+2:2Nz] = Ĝ[i, 1:Nx, Ny:-1:2, Nz:-1:2] * g[6, i] # yz
    end
    Ĝ[:, Nx+2:2Nx, Ny+2:2Ny, Nz+2:2Nz] = Ĝ[:, Nx:-1:2, Ny:-1:2, Nz:-1:2]   # xyz

end



