using LinearAlgebra
using FFTW
using AbstractFFTs
using IterativeSolvers
using StaticArrays
using DDA
using BenchmarkTools

# TODO: return value: P(3, Nx, Ny, Nz) or P(3, N_occupied)

k_norm = 1
dx = 1
alpha = Inf
alpha = 2.1

Nx = 16;
Ny = 16;
Nz = 16;

grid = CartesianGrid(Nx, Ny, Nz)


occ = rand(Bool, Nx, Ny, Nz)
X = rand(SVector{3,ComplexF64}, sum(occ))


# 1. calc A (Nx, Ny, Nz, 3, 3) or (Nx, Ny, Nz, 6) becuse of symmetry
# 2. add padding (2Nx, 2Ny, 2Nz, 6)
# 3. fft (2Nx, 2Ny, 2Nz, 6)
# 3. construct (2Nx, 2Ny, 2Nz, 3) vector from dipoles and indecies
# 4. plan for fft (ifft) of (2Nx, 2Ny, 2Nz, 3)
# 5. multiplication
# 6. ifft (2Nx, 2Ny, 2Nz, 3)
# 7. set zero for vacuum
# 8. polarisibality


import Base: *



"""
    TensorConvolution(...)

A preplanned, circular convolution operator on an M × N matrix of data of type T

# Fields
- `Ĝ`: DFT coefficients of the convolution kernel
- `F`: plan for FFT operator
- `paddedSpace`: scratch space to zero-pad the input matrix
- `Â`: scratch space to store the DFT coefficients of the zero-padded input matrix

# Constructors:
- `TensorConvolution(G::Matrix{T})
"""
struct TensorConvolution{P<:AbstractFFTs.Plan{ComplexF64}}
    Ĝ::Array{ComplexF64, 4}
    F::P
    # inds::Array{Bool, 3}
    inds::Vector{CartesianIndex{3}}
    alphas::ComplexF64
    tmp::Array{ComplexF64, 4}
end


function TensorConvolution(grid::CartesianGrid, occ, k_norm, alphas)

    G = Array{ComplexF64}(undef, 6, Nx, Ny, Nz)

    for i=1:Nx, j=1:Ny, k=1:Nz
        Rij = grid[i,j,k]
        A = DDA.calc_Ajk(k_norm, Rij)
        G[1,i,j,k] = A[1, 1] # xx
        G[2,i,j,k] = A[1, 2] # xy
        G[3,i,j,k] = A[1, 3] # xz
        G[4,i,j,k] = A[2, 2] # yy
        G[5,i,j,k] = A[2, 3] # yz
        G[6,i,j,k] = A[3, 3] # zz
        # G[i] = SymmetricTensor{2,3,ComplexF64}((A[1, 1], A[1, 2], A[1, 3], A[2, 2], A[2, 3], A[3, 3]))
        # G[:,i] = (A[1, 1], A[1, 2], A[1, 3], A[2, 2], A[2, 3], A[3, 3])
    end

    G[:, 1, 1, 1] .= zero(ComplexF64)
    # G[1, 1, 1] = zero(SymmetricTensor{2,3,ComplexF64})

    TensorConvolution(G::Array{ComplexF64, 4}, occ::Array{Bool, 3}, alphas)

end

function TensorConvolution(G::Array{ComplexF64, 4}, occ::Array{Bool, 3}, alphas)

    _, Nx, Ny, Nz = size(G)

    # Ĝ = Array{SymmetricTensor{2,3,ComplexF64,6}}(undef,  2Nx, 2Ny, 2Nz)
    Ĝ = Array{ComplexF64}(undef, 6, 2Nx, 2Ny, 2Nz)

    # Ĝ[1:Nx, 1:Ny, 1:Nz, :] = G
    copyto!(Ĝ, CartesianIndices(G), G, CartesianIndices(G))

    extend!(Ĝ)

    fft!(Ĝ, 2:4)

    # Ê = Array{SVector{3,ComplexF64}}(undef,  2Nx, 2Ny, 2Nz)
    # Ê = Array{ComplexF64}(undef, 3, 2Nx, 2Ny, 2Nz)
    tmp = zeros(ComplexF64, 3, 2Nx, 2Ny, 2Nz)

    F = FFTW.plan_fft!(tmp, 2:4, flags = FFTW.ESTIMATE)

    inds = findall(occ)

    TensorConvolution{typeof(F)}(Ĝ, F, inds, alphas, tmp)
end

Base.eltype(::TensorConvolution) = ComplexF64
Base.size(A::TensorConvolution) = (length(A.inds),length(A.inds))
Base.size(A::TensorConvolution,i::Int) = size(A)[i]

# TODO: copyto!
# TODO: function mul!(y::StridedVector, A::TensorConvolution, x::StridedVector, α::Number, β::Number)
using TimerOutputs
const to = TimerOutput()

function LinearAlgebra.mul!(y::StridedVector, A::TensorConvolution, x::StridedVector)

    @timeit to "copy1" begin
    fill!(A.tmp, 0)

    for i in eachindex(x)
        A.tmp[:, A.inds[i]] = x[i]
    end
    end  # timeit

    @timeit to "fft" begin
    A.F * A.tmp # FFT
    end  # timeit

    @timeit to "matvec" begin
    # !! Matrix vector multiplication in the Fourier space

    _, Nx, Ny, Nz = size(A.Ĝ)
    for k=1:Nz, j=1:Ny, i=1:Nx
        E = A.tmp[:,i,j,k]
        A.tmp[1,i,j,k] = A.Ĝ[1,i,j,k] * E[1] + A.Ĝ[2,i,j,k] * E[2] + A.Ĝ[3,i,j,k] * E[3] # x component
        A.tmp[2,i,j,k] = A.Ĝ[2,i,j,k] * E[1] + A.Ĝ[4,i,j,k] * E[2] + A.Ĝ[5,i,j,k] * E[3] # y component
        A.tmp[3,i,j,k] = A.Ĝ[3,i,j,k] * E[1] + A.Ĝ[5,i,j,k] * E[2] + A.Ĝ[6,i,j,k] * E[3] # z component
    end
    end  # timeit

    @timeit to "ifft" begin
    inv(A.F) * A.tmp # inverse FFT
    end  # timeit

    @timeit to "copy2" begin

    for i in eachindex(y)
        y[i] = A.tmp[:, A.inds[i]] .+ 1/A.alphas * x[i]
    end
    end  # timeit
    return y

end

A::TensorConvolution * x = mul!(similar(x), A, x)



function extend!(Ĝ::Array{ComplexF64, 4})

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

    # Ĝ[Nx+1, :, :] = zero(SymmetricTensor{2,3,ComplexF64})
    # Ĝ[:, Ny+1, :] = zero(SymmetricTensor{2,3,ComplexF64})
    # Ĝ[:, :, Nz+1] = zero(SymmetricTensor{2,3,ComplexF64})

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


function test_fft_multiply(grid, X, k_norm, alpha, occ)
    A = TensorConvolution(grid, occ, k_norm, alpha)
    return A * X
end

x_fft = test_fft_multiply(grid, X, k_norm, alpha, occ)

function test_direct_multiply(grid, X, k_norm, alpha, occ)
    r =  grid[occ]
    A = DDA.interactions(k_norm, r,  ones(length(r)) * alpha)
    # A[diagind(A)] .= 0

    m = A * reinterpret(ComplexF64, X)
    return reshape(m, 3,:)
end


x_direct = test_direct_multiply(grid, X, k_norm, alpha, occ)

reinterpret(reshape, ComplexF64, x_fft) ≈ x_direct





A = TensorConvolution(grid, occ, k_norm, alpha)
x_fft  = similar(X)
mul!(x_fft, A, X)
@btime  mul!(x_fft, A, X)



nothing








#
# # function test_fft_multiply(grid, X, k_norm, alpha, inds)
#     Nx, Ny, Nz = size(grid)
#
#     # Ĝ = Array{SymmetricTensor{2,3,ComplexF64,6}}(undef,  2Nx, 2Ny, 2Nz)
#     Ĝ = Array{ComplexF64}(undef, 6, 2Nx, 2Ny, 2Nz)
#     # Ê = Array{SVector{3,ComplexF64}}(undef,  2Nx, 2Ny, 2Nz)
#     # Ê = Array{ComplexF64}(undef, 3, 2Nx, 2Ny, 2Nz)
#     Ê = zeros(ComplexF64, 3, 2Nx, 2Ny, 2Nz)
#
#     for i=1:Nx, j=1:Ny, k=1:Nz
#         Rij = grid[i,j,k]
#         A = DDA.calc_Ajk(k_norm, Rij)
#         Ĝ[1,i,j,k] = A[1, 1] # xx
#         Ĝ[2,i,j,k] = A[1, 2] # xy
#         Ĝ[3,i,j,k] = A[1, 3] # xz
#         Ĝ[4,i,j,k] = A[2, 2] # yy
#         Ĝ[5,i,j,k] = A[2, 3] # yz
#         Ĝ[6,i,j,k] = A[3, 3] # zz
#         # Ĝ[i] = SymmetricTensor{2,3,ComplexF64}((A[1, 1], A[1, 2], A[1, 3], A[2, 2], A[2, 3], A[3, 3]))
#         # Ĝ[:,i] = (A[1, 1], A[1, 2], A[1, 3], A[2, 2], A[2, 3], A[3, 3])
#     end
#
#     Ĝ[:, 1, 1, 1] .= zero(ComplexF64)
#     # Ĝ[1, 1, 1] = zero(SymmetricTensor{2,3,ComplexF64})
#
#     # Periodic extension
#     g = [+1 -1 -1 +1 +1 +1  # x
#         +1 -1 +1 +1 -1 +1  # y
#         +1 +1 -1 +1 -1 +1  # z
#         +1 +1 -1 +1 -1 +1  # xy
#         +1 -1 +1 +1 -1 +1  # xz
#         +1 -1 -1 +1 +1 +1]  # yz
#
#     Ĝ[:, Nx+1, :, :] .= zero(ComplexF64)
#     Ĝ[:, :, Ny+1, :] .= zero(ComplexF64)
#     Ĝ[:, :, :, Nz+1] .= zero(ComplexF64)
#     # Ĝ[Nx+1, :, :] = zero(SymmetricTensor{2,3,ComplexF64})
#     # Ĝ[:, Ny+1, :] = zero(SymmetricTensor{2,3,ComplexF64})
#     # Ĝ[:, :, Nz+1] = zero(SymmetricTensor{2,3,ComplexF64})
#
#     for i = 1:6
#         Ĝ[i, Nx+2:2Nx, 1:Ny, 1:Nz] = Ĝ[i, Nx:-1:2, 1:Ny, 1:Nz] * g[1, i]       # x
#         Ĝ[i, 1:Nx, Ny+2:2Ny, 1:Nz] = Ĝ[i, 1:Nx, Ny:-1:2, 1:Nz] * g[2, i]       # y
#         Ĝ[i, 1:Nx, 1:Ny, Nz+2:2Nz] = Ĝ[i, 1:Nx, 1:Ny, Nz:-1:2] * g[3, i]       # z
#         Ĝ[i, Nx+2:2Nx, Ny+2:2Ny, 1:Nz] = Ĝ[i, Nx:-1:2, Ny:-1:2, 1:Nz] * g[4, i] # xy
#         Ĝ[i, Nx+2:2Nx, 1:Ny, Nz+2:2Nz] = Ĝ[i, Nx:-1:2, 1:Ny, Nz:-1:2] * g[5, i] # xz
#         Ĝ[i, 1:Nx, Ny+2:2Ny, Nz+2:2Nz] = Ĝ[i, 1:Nx, Ny:-1:2, Nz:-1:2] * g[6, i] # yz
#     end
#     Ĝ[:, Nx+2:2Nx, Ny+2:2Ny, Nz+2:2Nz] = Ĝ[:, Nx:-1:2, Ny:-1:2, Nz:-1:2]   # xyz
#
#     fft!(Ĝ, 2:4)
#
#     fftw_flags = FFTW.ESTIMATE
#     F = FFTW.plan_fft!(Ê, 2:4, flags = fftw_flags)
#
#     cinds = findall(inds)
#
#     return do_multiply(F, Ĝ, Ê, X, alpha, cinds, Nx, Ny, Nz)
#
# # end
# using TimerOutputs
# const to = TimerOutput()
#
# function do_multiply(F, Ĝ, Ê, X, alpha, cinds, Nx, Ny, Nz)
#
#     @timeit to "nest 1" begin
#         @timeit to "copy1" begin
#             fill(Ê,0)
#             for i = 1:length(cinds)
#                 Ê[:, cinds[i]] = X[i]
#             end
#         end
#         @timeit to "fft" begin
#
#             F * Ê
#         end
#         @timeit to "mul" begin
#
#             # !! Matrix vector multiplication
#             for k=1:2Nz, j=1:2Ny, i=1:2Nx
#                 tmp = Ê[:,i,j,k]
#
#                 # !! -> x component
#                 Ê[1,i,j,k] = Ĝ[1,i,j,k] * tmp[1] + Ĝ[2,i,j,k] * tmp[2] + Ĝ[3,i,j,k] * tmp[3]
#                 # !! -> y component
#                 Ê[2,i,j,k] = Ĝ[2,i,j,k] * tmp[1] + Ĝ[4,i,j,k] * tmp[2] + Ĝ[5,i,j,k] * tmp[3]
#                 # !! -> z component
#                 Ê[3,i,j,k] = Ĝ[3,i,j,k] * tmp[1] + Ĝ[5,i,j,k] * tmp[2] + Ĝ[6,i,j,k] * tmp[3]
#             end
#         end
#         @timeit to "ifft" begin
#
#             inv(F) * Ê
#         end
#         @timeit to "alpha" begin
#
#             for i = 1:length(cinds)
#                 Ê[:, cinds[i]] += 1/alpha * X[i]
#             end
#         end
#         @timeit to "copy2" begin
#
#            out =  Ê[:,cinds]
#         end
#     end
#     return out
# end