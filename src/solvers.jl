abstract type AbstractMethod end

struct Direct <: AbstractMethod end
struct BiCGStabl <: AbstractMethod end
# CG
# MINRES
# GMRES
# IDRs


abstract type AbstractSolution end

C_sca(sol::AbstractSolution) = C_ext(sol) - C_abs(sol)

# TODO: storing the alg and its parameters
# TODO: storing the prepared (factorised) interaction matrix (like in EELS calculations)
# TODO: do not trim A*E because it is useful for the nearfield calculations

struct GridSolution{T} <: AbstractSolution
    P::Vector{SVector{3, Complex{T}}}
    prob::GridProblem{T}
    alg::AbstractMethod
end

C_abs(sol::GridSolution) = C_abs(sol.prob.k, sol.prob.E0, sol.P, sol.prob.alphas)
C_ext(sol::GridSolution) = C_ext(sol.prob.k, sol.prob.E0, sol.prob.Eincs, sol.P)
# C_abs(sol::GridSolution) = C_abs(norm(sol.prob.k), norm(sol.prob.E0), sol.P, sol.prob.alphas)
# C_ext(sol::GridSolution)= C_ext(norm(sol.prob.k), norm(sol.prob.E0), sol.prob.Eincs, sol.P)


# # TODO: move this to solvers + for s in scatterers (merging)
# function discretize(p::GridProblem)
#     # 1. create the coordinates of the dipoles,
#     occ = discretize(p.grid, p.scatterers.target)
#
#     # 2. assign the polarizability αj to each dipole
#     alphas = polarisbility(p.scatterers.model, p)
#
#     coords = p.grid[occ]
#     return coords, occ, alphas
# end
#
#
#
# function discretize(p::GridProblem, s::Scatterer)
#     # 1. create the coordinates of the dipoles,
#     occ = discretize(p.grid, s.target)
#
#     # 2. assign the polarizability αj to each dipole
#     alphas = polarisbility(s.model, p)
#
#     return occ, alphas
# end
#

# TODO: Warnings
# k = 2π/λ
# |m|kd = 1
# x = ka =2πa/λ


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
    alphas::Vector{ComplexF64}
    tmp::Array{ComplexF64,4}
end

Base.eltype(::TensorConvolution) = ComplexF64
Base.size(A::TensorConvolution) = (3 * length(A.inds), 3 * length(A.inds))
Base.size(A::TensorConvolution, i::Int) = size(A)[i]


function TensorConvolution(grid::CartesianGrid, inds, k_norm, alphas)

    Nx, Ny, Nz = size(grid)
    Ĝ = Array{ComplexF64}(undef, 6, 2Nx, 2Ny, 2Nz)

    for i = 1:Nx, j = 1:Ny, k = 1:Nz
        Rij = [i-1, j-1, k-1] *  grid.spacing
        A = calc_Ajk(k_norm, Rij)
        Ĝ[1, i, j, k] = A[1, 1] # xx
        Ĝ[2, i, j, k] = A[1, 2] # xy
        Ĝ[3, i, j, k] = A[1, 3] # xz
        Ĝ[4, i, j, k] = A[2, 2] # yy
        Ĝ[5, i, j, k] = A[2, 3] # yz
        Ĝ[6, i, j, k] = A[3, 3] # zz
    end
    Ĝ[:, 1, 1, 1] .= zero(ComplexF64)

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

    fft!(Ĝ, 2:4)

    tmp = zeros(ComplexF64, 3, 2Nx, 2Ny, 2Nz)
    F = FFTW.plan_fft!(tmp, 2:4, flags = FFTW.ESTIMATE)

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
        y[i] = A.tmp[:, A.inds[i]] .+ inv(A.alphas[i]) * x[i]
    end

    return y

end


function solve(p::GridProblem, alg::BiCGStabl; kwargs...)

    # 5. solve for P in the system of linear equations
    T = Float64
    V = SVector{3,Complex{Float64}}

    P = similar(p.alphas, V)
    fill!(P, zero(V))

    return solve!(P, p, alg; kwargs...)

end


function solve!(P::Vector{SVector{3,Complex{Float64}}}, p::GridProblem, alg::BiCGStabl;
    reltol = 1e-3, verbose = true, kwargs...)

    # # 1. create the coordinates of the dipoles
    # # 2. assign the polarizability αj to each dipole
    # # occ = DDA.discretize(p.grid, p.scatterer.target)
    # # coords = p.grid[occ]
    # # alphas = polarisbility(p.scatterer.model, p)
    # # coords, occ, alphas = DDA.discretize(p.grid, p.scatterer, p.Einc)
    #
    # coords, occ, alphas = discretize(p)
    #
    #
    # # 3. calculated the incident field Einc, at each dipole,
    # # Einc = field(p.Einc, p.grid[occ])
    #
    # Einc = similar(coords, SVector{3,Complex{Float64}})
    # for i = eachindex(coords)
    #     Einc[i] = field(p.Einc, coords[i])
    # end

    T = Float64
    # V = SVector{3,Complex{T}}
    #
    # P = similar(p.alphas, V)
    # fill!(P, zero(V))

    # 4. assemble the interaction matrix A and
    A = TensorConvolution(p.grid, p.inds, p.k, p.alphas)

    # 5. solve for P in the system of linear equations
    bicgstabl!(reinterpret(Complex{T}, P), A, reinterpret(Complex{T}, p.Eincs); reltol = reltol, verbose = verbose, kwargs...)

    return GridSolution(P, p, alg)
end




struct DipoleSolution <: AbstractSolution
    P
    prob
    alg
end

C_abs(sol::DipoleSolution) = C_abs(sol.prob.k, sol.prob.E0, sol.P, sol.prob.alphas)
C_ext(sol::DipoleSolution) = C_ext(sol.prob.k, sol.prob.E0, sol.prob.Eincs, sol.P)


function interactions(k, r, alphas)
    N = length(r)
    out = zeros(ComplexF64, 3, N, 3, N)

    for i in 1:N
        out[:, i, :, i] = Diagonal(inv(alphas[i]) * I, 3)
    end

    for i in 2:N
        for j in 1:i-1
            out[:, i, :, j] = calc_Ajk(k, r[i], r[j])
        end
    end
    return Symmetric(reshape(out, 3N, 3N), :L)
end


function solve(p::DipoleProblem, alg::BiCGStabl; reltol = 1e-3, verbose = true, kwargs...)

    # # 1. create the coordinates of the dipoles
    # # 2. assign the polarizability αj to each dipole
    # # occ = DDA.discretize(p.grid, p.scatterer.target)
    # # coords = p.grid[occ]
    # # alphas = polarisbility(p.scatterer.model, p)
    # # coords, occ, alphas = DDA.discretize(p.grid, p.scatterer, p.Einc)
    #
    # coords, occ, alphas = discretize(p)
    #
    #
    # # 3. calculated the incident field Einc, at each dipole,
    # # Einc = field(p.Einc, p.grid[occ])
    #
    # Einc = similar(coords, SVector{3,Complex{Float64}})
    # for i = eachindex(coords)
    #     Einc[i] = field(p.Einc, coords[i])
    # end
    #
    #
    # # 4. assemble the interaction matrix A and
    #
    # k = wavenumber(p.Einc)
    # A_conv = TensorConvolution(p.grid, occ, k, alphas)

    A = interactions(p.k, p.dipoles, p.alphas)


    # 5. solve for P in the system of linear equations

    T = Float64
    V = SVector{3,Complex{T}}

    P = similar(p.dipoles, V)
    fill!(P, zero(V))

    bicgstabl!(reinterpret(Complex{T}, P), A, reinterpret(Complex{T}, p.Eincs); reltol = reltol, verbose = verbose, kwargs...)

    return DipoleSolution(P, p, alg)
end


