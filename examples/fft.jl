using Revise
using LinearAlgebra
using FFTW
# using AbstractFFTs
using IterativeSolvers

k = 1
dx = 1
alpha = 1

Nx = 4; Ny = 4; Nz = 4


# dipoles = zeros(Float64, Nx÷2, Ny÷2, Nz÷2, 3)
# for k in 1:Nz÷2, j in 1:Ny÷2,i in 1:Nx÷2
#     # @show [i,j,k].-1
#     dipoles[i,j,k,:] = ([i,j,k].-1)*dx
# end
#

dipoles = Vector{Vector{Float64}}()
for k in 1:Nz÷2, j in 1:Ny÷2,i in 1:Nx÷2
    # @show [i,j,k].-1
    position = ([i,j,k].-1)*dx
    push!(dipoles, position )
end

# A = DDA.interaction(k,dipoles, alpha*ones(length(dipoles)))

# Figuring out why it is Topleitz????
# AA = zeros(ComplexF64, Nx, Ny, Nz, Nx, Ny, Nz, 3, 3)
# for k in 1:Nz, j in 1:Ny,i in 1:Nx, kk in 1:Nz, kj in 1:Ny,ki in 1:Nx
#     Rij = ([i,j,k].-[ki,kj,kk])*dx
#     if i==ki && j==kj && k==kk
#         @show i,j,k
#         continue
#     end
#     # AA[i,j,k,ki,kj,kk,:,:] =calc_Ajk(k,Rij)
#     AA[i,ki,j,kj,k,kk,:,:] =calc_Ajk(k,Rij)
# end
#
# AA[:,2,3,:,1,1,3,3] # Block-Toeplitz ???

# Figuring out why it is Topleitz????
# dipoles = Vector{Vector{Float64}}()
# for k in 1:Nz, j in 1:Ny,i in 1:Nx
#     # @show [i,j,k].-1
#     position = ([i,j,k].-1)*dx
#     push!(dipoles, position )
# end
#
# A = zeros(ComplexF64, length(dipoles), length(dipoles), 3,3)
# for i in 1:length(dipoles)
#     for j in 1:length(dipoles)
#         # @show i, j
#         # @show k, r[i], r[j]
#         # @show DDA.calc_Ajk(k, r[i], r[j])
#         A[i,j,:,:] = calc_Ajk(k, dipoles[i], dipoles[j])
#     end
# end

# Figuring out why it is Topleitz????
# A = interaction(k,dipoles, alpha*ones(length(dipoles)))




import Base: *


# struct CircularConvolution{M, N, T, K, KI}
#     Ĝ::Matrix{ComplexF64}
#     F::K
#     F⁻¹::KI
#     nthreads::Int
#
#     paddedSpace::Matrix{T}
#     Â::Matrix{ComplexF64}
# end
#
# struct ToeplitzFactorization{T<:Number,A<:AbstractToeplitz{T},S<:Number,P<:Plan{S}} <: Factorization{T}
#     vcvr_dft::Vector{S}
#     tmp::Vector{S}
#     dft::P
# end





function TensorConvolution2(G::AbstractMatrix)

    M, N = size(G)

    #paddedSpace = Matrix{Float64}(undef, 2M-1, 2N-1)
    paddedSpace = Matrix{dtype}(undef, 2M, 2N)

    if dtype == ComplexF64
      F = FFTW.plan_fft(paddedSpace, flags = fftw_flags)
    else
      F = FFTW.plan_rfft(paddedSpace, flags = fftw_flags)
    end


    mirror!(paddedSpace, G)
    Ĝ = F * paddedSpace

    Â = similar(Ĝ)
    #F⁻¹ = FFTW.plan_irfft(Â, 2M - 1, flags = fftw_flags)
    if dtype == ComplexF64
      F⁻¹ = FFTW.plan_ifft(Â, flags = fftw_flags)
    else
      F⁻¹ = FFTW.plan_irfft(Â, 2M, flags = fftw_flags)
    end

    CircularConvolution{M, N, dtype, typeof(F), typeof(F⁻¹)}(Ĝ, F, F⁻¹, nthreads, paddedSpace, Â)
end


function Base.show(io::IO, c::CircularConvolution{M, N, T}) where {M, N, T}
    print(io, "Circular convolution on a $M × $N matrix of data type $T")
end

function CircularConvolution(G::AbstractMatrix{T},fftw_flags = FFTW.ESTIMATE; dtype = Float64, nthreads = length(Sys.cpu_info())) where {T}
    FFTW.set_num_threads(nthreads)

    M, N = size(G)
    #paddedSpace = Matrix{Float64}(undef, 2M-1, 2N-1)
    paddedSpace = Matrix{dtype}(undef, 2M, 2N)

    if dtype == ComplexF64
      F = FFTW.plan_fft(paddedSpace, flags = fftw_flags)
    else
      F = FFTW.plan_rfft(paddedSpace, flags = fftw_flags)
    end


    mirror!(paddedSpace, G)
    Ĝ = F * paddedSpace

    Â = similar(Ĝ)
    #F⁻¹ = FFTW.plan_irfft(Â, 2M - 1, flags = fftw_flags)
    if dtype == ComplexF64
      F⁻¹ = FFTW.plan_ifft(Â, flags = fftw_flags)
    else
      F⁻¹ = FFTW.plan_irfft(Â, 2M, flags = fftw_flags)
    end

    CircularConvolution{M, N, dtype, typeof(F), typeof(F⁻¹)}(Ĝ, F, F⁻¹, nthreads, paddedSpace, Â)
end

function mul!(out, C::CircularConvolution{M, N, T}, B) where {M, N, T}
    FFTW.set_num_threads(C.nthreads)

    MB, NB = size(B)
    #@assert size(out) == size(B) == (M, N)
    @assert size(out) == (MB, NB)

    inds = CartesianIndices((MB,NB))
    fill!(C.paddedSpace, 0)
    copyto!(C.paddedSpace, inds, B, inds)
    mul!(C.Â, C.F, C.paddedSpace)

    C.Â .*= C.Ĝ

    mul!(C.paddedSpace, C.F⁻¹, C.Â)

    #copyto!(out, inds, C.paddedSpace, CartesianIndices((M:2M-1,N:2N-1)))
    copyto!(out, inds, C.paddedSpace, CartesianIndices((M+1:M+MB,N+1:N+NB)))

end

C::CircularConvolution * B = mul!(similar(B), C, B)

function mirror!(A)
    Nr, Nc = size(a)

    A .= 0
    A[2:Nr, 2:Nc] .= a[Nr:-1:2, Nc:-1:2]
    A[2:Nr, Nc+1:end] .= a[Nr:-1:2, 1:Nc]
    A[Nr+1:end, 2:Nc] .= a[1:Nr, Nc:-1:2]
    A[Nr+1:end, Nc+1:end] .= a[1:Nr, 1:Nc]
    A
end

function mirror(a::AbstractArray{T,2}) where {T}
    Nr, Nc = size(a)
    mirror!(zeros(T, 2Nr-1, 2Nc-1), a)
end






# 1. calc A (Nx, Ny, Nz, 3, 3) or (Nx, Ny, Nz, 6) becuse of symmetry
# 2. add padding (2Nx, 2Ny, 2Nz, 6)
# 3. fft (2Nx, 2Ny, 2Nz, 6)
# 3. construct (2Nx, 2Ny, 2Nz, 3) vector from dipoles and indecies
# 4. plan for fft (ifft) of (2Nx, 2Ny, 2Nz, 3)


# //////////////////////////////////////////////////////////
#  building "toeplitz matix"
# Fill the interaction tensor A(i,j,k,1:6)
# A = zeros(ComplexF64, Nx, Ny, Nz, 6)
# A = Array{ComplexF64}(undef, Nx, Ny, Nz, 6)


"""
    CircularConvolution{M, N, T}

A preplanned, circular convolution operator on an M × N matrix of data of type T

# Fields
- `Ĝ`: DFT coefficients of the convolution kernel
- `F`: preplanned rFFT operator
- `F⁻¹`: preplanned irFFT operator
- `paddedSpace`: scratch space to zero-pad the input matrix
- `Â`: scratch space to store the DFT coefficients of the zero-padded input matrix

# Constructors:
- `CircularConvolution(G::Matrix{T})
"""
struct TensorConvolution
    Ĝ::Matrix{ComplexF64}
    F::AbstractFFTs.Plan{ComplexF64}
    inds_vacuum::Vector{CartesianIndex}
    inds::Vector{CartesianIndex}
    alphas
    Xtmp::Matrix{ComplexF64}
end







function TensorConvolution(g, inds, k, alphas)

    Nx, Ny, Nz = size(g)

    # Ĝ = zeros(ComplexF64, 2Nx, 2Ny, 2Nz, 6)
    Ĝ = Array{ComplexF64}(undef, 2Nx, 2Ny, 2Nz, 6)
    Ê = zeros(ComplexF64, 2Nx, 2Ny, 2Nz, 3)

    for i in Base.
        Rij = g[i]
        A = DDA.calc_Ajk(k, Rij)
        Ĝ[i,j,k,1] = A[1,1] # xx
        Ĝ[i,j,k,2] = A[1,2] # xy
        Ĝ[i,j,k,3] = A[1,3] # xz
        Ĝ[i,j,k,4] = A[2,2] # yy
        Ĝ[i,j,k,5] = A[2,3] # yz
        Ĝ[i,j,k,6] = A[3,3] # zz

    end

    periodic_extension!(Ĝ)

    F = FFTW.plan_fft(view(Ĝ, :,:,:,1), flags = fftw_flags)

    for i in 1:6
        fft!(view(Ĝ, :,:,:,i))
    end


    CircularConvolution{M, N, dtype, typeof(F), typeof(F⁻¹)}(Ĝ, F)
end

function periodic_extension!(G::AbstractMatrix)
    # Periodic extension
    g = [+1 -1 -1 +1 +1 +1;  # x
        +1 -1 +1 +1 -1 +1;  # y
        +1 +1 -1 +1 -1 +1;  # z
        +1 +1 -1 +1 -1 +1;  # xy
        +1 -1 +1 +1 -1 +1;  # xz
        +1 -1 -1 +1 +1 +1]  # yz

    # G[1:Nx, 1:Ny, 1:Nz, :] = A
    for i = 1:6
        G[Nx+2:2Nx, 1:Ny, 1:Nz, i] = G[Nx:-1:2, 1:Ny, 1:Nz, i] * g[1, i]       # x
        G[1:Nx, Ny+2:2Ny, 1:Nz, i] = G[1:Nx, Ny:-1:2, 1:Nz, i] * g[2, i]       # y
        G[1:Nx, 1:Ny, Nz+2:2Nz, i] = G[1:Nx, 1:Ny, Nz:-1:2, i] * g[3, i]       # z
        G[Nx+2:2Nx, Ny+2:2Ny, 1:Nz, i] = G[Nx:-1:2,Ny:-1:2, 1:Nz, i] * g[4, i] # xy
        G[Nx+2:2Nx, 1:Ny, Nz+2:2Nz, i] = G[Nx:-1:2, 1:Ny,Nz:-1:2, i] * g[5, i] # xz
        G[1:Nx, Ny+2:2Ny, Nz+2:2Nz, i] = G[1:Nx,Ny:-1:2, Nz:-1:2, i] * g[6, i] # yz
        G[Nx+2:2Nx, Ny+2:2Ny, Nz+2:2Nz, i] = G[Nx:-1:2, Ny:-1:2, Nz:-1:2, i]   # xyz
    end
end


































# Allocate vectors for the matrice vector operator
X3D = zeros(ComplexF64, 2Nx, 2Ny, 2Nz, 3)

Xcomp = zeros(ComplexF64, 2Nx, 2Ny, 2Nz)
Ycomp = zeros(ComplexF64, 2Nx, 2Ny, 2Nz)
Zcomp = zeros(ComplexF64, 2Nx, 2Ny, 2Nz)




# //////////////////////////////////////////////////////////
#  building "toeplitz matix"
# Fill the interaction tensor A(i,j,k,1:6)
A = zeros(ComplexF64, Nx, Ny, Nz, 6)
# A = Array{ComplexF64}(undef, Nx, Ny, Nz, 6)


for i = 1:Nx, j = 1:Ny, k=1:Nz
    Rij = ([i,j,k].-1)*dx
    if i==1 && j==1 && k==1
        @show i,j,k
        continue
    end
    g = DDA.calc_Ajk(k, Rij)
    A[i,j,k,1] = g[1,1] # xx
    A[i,j,k,2] = g[1,2] # xy
    A[i,j,k,3] = g[1,3] # xz
    A[i,j,k,4] = g[2,2] # yy
    A[i,j,k,5] = g[2,3] # yz
    A[i,j,k,6] = g[3,3] # zz

end


# Periodic extension
g = [+1 -1 -1 +1 +1 +1;  # x
     +1 -1 +1 +1 -1 +1;  # y
     +1 +1 -1 +1 -1 +1;  # z
     +1 +1 -1 +1 -1 +1;  # xy
     +1 -1 +1 +1 -1 +1;  # xz
     +1 -1 -1 +1 +1 +1]  # yz

G = zeros(ComplexF64, 2Nx, 2Ny, 2Nz, 6)
G[1:Nx, 1:Ny, 1:Nz, :] = A
for i = 1:6
    G[Nx+2:2Nx, 1:Ny, 1:Nz, i] = G[Nx:-1:2, 1:Ny, 1:Nz, i] * g[1, i]       # x
    G[1:Nx, Ny+2:2Ny, 1:Nz, i] = G[1:Nx, Ny:-1:2, 1:Nz, i] * g[2, i]       # y
    G[1:Nx, 1:Ny, Nz+2:2Nz, i] = G[1:Nx, 1:Ny, Nz:-1:2, i] * g[3, i]       # z
    G[Nx+2:2Nx, Ny+2:2Ny, 1:Nz, i] = G[Nx:-1:2,Ny:-1:2, 1:Nz, i] * g[4, i] # xy
    G[Nx+2:2Nx, 1:Ny, Nz+2:2Nz, i] = G[Nx:-1:2, 1:Ny,Nz:-1:2, i] * g[5, i] # xz
    G[1:Nx, Ny+2:2Ny, Nz+2:2Nz, i] = G[1:Nx,Ny:-1:2, Nz:-1:2, i] * g[6, i] # yz
    G[Nx+2:2Nx, Ny+2:2Ny, Nz+2:2Nz, i] = G[Nx:-1:2, Ny:-1:2, Nz:-1:2, i]   # xyz
end


Ĝ=similar(G)
for i in 1:6
    Ĝ[:,:,:,i] = fft(G[:,:,:,i])
end

for i in 1:6
    # doesn't work: fft!(G[:,:,:,i])
    fft!(view(G, :,:,:,i))
end






# import Base: convert, *, \, getindex, print_matrix, size, Matrix, +, -, copy, similar, sqrt, copyto!,
#     adjoint, transpose
# import LinearAlgebra: Cholesky, DimensionMismatch, cholesky, cholesky!, eigvals, inv, ldiv!,
#     mul!, pinv, rmul!, tril, triu
#
# using LinearAlgebra: LinearAlgebra, Adjoint, Factorization, factorize


import Base: convert, *, \, getindex, print_matrix, size, Matrix, +, -, copy, similar, sqrt, copyto!,
    adjoint, transpose
import LinearAlgebra: Cholesky, DimensionMismatch, cholesky, cholesky!, eigvals, inv, ldiv!,
    mul!, pinv, rmul!, tril, triu

using LinearAlgebra: LinearAlgebra, Adjoint, Factorization, factorize

using AbstractFFTs
using AbstractFFTs: Plan



# For matrix-free types of A the following interface is expected to be defined:
# - A*v computes the matrix-vector product on a v::AbstractVector;
# - mul!(y, A, v) computes the matrix-vector product on a v::AbstractVector in-place;
# - eltype(A) returns the element type implicit in the equivalent matrix representation of A;
# - size(A, d) returns the nominal dimensions along the dth axis in the equivalent matrix representation of A.







# gg = copy(G[:,:,:, 1])
# @show fft(gg) ≈ Ĝ[:,:,:,1]


# function test(k,Rij)
#     cik0 = 1im * k
#     idr = OffsetArray(zeros(6), 0:5)
#     constvar = zeros(ComplexF64, 4)
#
#     idr[2] = 1.0/dot(Rij,Rij)  # [idr[2]] = m**(-2)
#     idr[1] = sqrt(idr[2])
#     idr[3] = idr[2]*idr[1]
#     idr[4] = idr[3]*idr[1]
#     idr[5] = idr[4]*idr[1]
#     idr[0] = 1.0/idr[1]                #! [idr[0]] = m
#
#     constvar[1] = (k^2*idr[3] + 3.0*(cik0*idr[4] - idr[5]))
#     constvar[2] = idr[3] - cik0*idr[2] - k^2*idr[1]
#     constvar[3] = exp(cik0*idr[0])           #! without '-'
#     constvar[4] = constvar[1]*constvar[3]
#
#
#     Tij= zeros(ComplexF64,3,3)
#     Tij[1,1] = (Rij[1]*Rij[1]*constvar[1] + constvar[2])*constvar[3]
#     Tij[2,1] = Rij[2]*Rij[1]*constvar[4]
#     Tij[3,1] = Rij[3]*Rij[1]*constvar[4]
#     Tij[1,2] = Tij[2,1]
#     Tij[2,2] = (Rij[2]*Rij[2]*constvar[1] + constvar[2])*constvar[3]
#     Tij[3,2] = Rij[3]*Rij[2]*constvar[4]
#     Tij[1,3] = Tij[3,1]
#     Tij[2,3] = Tij[3,2]
#     Tij[3,3] = (Rij[3]*Rij[3]*constvar[1] + constvar[2])*constvar[3]
#
#     return Tij
# end
#
# k = rand()
# Rij=rand(3)
#
# Tij = test(k, Rij)
# @show Tij ≈ calc_Ajk(k,Rij)



#
# # A working demonstraion of toeplitz matrices
# Nx, Ny, Nz = 6, 6, 6
#
# A = zeros( Nx, Ny, Nz, Nx, Ny, Nz)
# for i = 1:Nx, ii = 1:Nx, j = 1:Ny, jj = 1:Ny, k = 1:Nz,  kk = 1:Nz
#
#     A[i, j, k, ii, jj, kk] = round(norm( ([i, j, k] .- 1) .- ([ii, jj, kk] .- 1)), sigdigits=2)
#
# end
#
# A[:,1,1,:,1,1]
# A[1,:,1,1,:,1]
# A[1,1,:,1,1,:]
# A[1,4,:,2,3,:]
