using Revise

using LinearAlgebra
using FFTW
# using AbstractFFTs
using IterativeSolvers
using DDA

# module ToeplitzMatrix
# using StatsBase



# Stored data:
# Nx, Ny, Nz; N = Nx*Ny*Nz
# - A[Nx, Ny, Nz, 6] (A[Nx, Ny, Nz, 3,3])
# - G[Nx+1, 2Ny+1, 2Nz+1, 6] (G[2Nx, 2Ny, 2Nz, 6] or G[2Nx, 2Ny, 2Nz, 3,3])
# - x[Nx, Ny, Nz, 3]
# - X[2Nx, 2Ny, 2Nz, 3] fourier space
# - Occ[Nx, Ny, Nz]


g = CubicGrid(10,10,10)

target = 1

# For matrix-free types of A the following interface is expected to be defined:
# - A*v computes the matrix-vector product on a v::AbstractVector;
# - mul!(y, A, v) computes the matrix-vector product on a v::AbstractVector in-place;
# - eltype(A) returns the element type implicit in the equivalent matrix representation of A;
# - size(A, d) returns the nominal dimensions along the dth axis in the equivalent matrix representation of A.




import Base: convert, *, \, getindex, print_matrix, size, Matrix, +, -, copy, similar, sqrt, copyto!,
    adjoint, transpose, eltype
import LinearAlgebra: Cholesky, DimensionMismatch, cholesky, cholesky!, eigvals, inv, ldiv!,
    mul!, pinv, rmul!, tril, triu

using LinearAlgebra: LinearAlgebra, Adjoint, Factorization, factorize

using AbstractFFTs
using AbstractFFTs: Plan








struct InteractionTensorFactorization
    grid
    kvec
end

Base


"""
    InteractionTensor

A Toeplitz matrix.
"""
struct InteractionTensor{T<:Number}
    vc::Vector{T}
    vr::Vector{T}

    # function InteractionTensor{T}(vc::Vector{T}, vr::Vector{T}) where {T<:Number}
    #     if first(vc) != first(vr)
    #         error("First element of the vectors must be the same")
    #     end
    #     return new{T}(vc, vr)
    # end
end


# Size of a general InteractionTensor matrix
function size(A::InteractionTensor, dim::Int)
    if dim == 1
        return length(A.vc)
    elseif dim == 2
        return length(A.vr)
    elseif dim > 2
        return 1
    else
        error("arraysize: dimension out of range")
    end
end

function getindex(A::InteractionTensor, i::Integer)
    return A[mod(i - 1, size(A, 1)) + 1, div(i - 1, size(A, 1)) + 1]
end
# Retrieve an entry
function getindex(A::InteractionTensor, i::Integer, j::Integer)
    m = size(A,1)
    n = size(A,2)
    if i > m || j > n
        error("BoundsError()")
    end

    if i >= j
        return A.vc[i - j + 1]
    else
        return A.vr[1 - i + j]
    end
end







# Application of a general InteractionTensor matrix to a general matrix
function mul!(
    C::StridedMatrix, A::InteractionTensor, B::StridedMatrix, α::Number, β::Number
)
    return mul!(C, factorize(A), B, α, β)
end


# Fast application of a general InteractionTensor matrix to a column vector via FFT
function mul!(
    y::StridedVector, A::InteractionTensor, x::StridedVector, α::Number, β::Number
)
    m, n = size(A)
    if length(y) != m
        throw(DimensionMismatch(
            "first dimension of A, $(m), does not match length of y, $(length(y))"
        ))
    end
    if length(x) != n
        throw(DimensionMismatch(
            "second dimension of A, $(n), does not match length of x, $(length(x))"
        ))
    end



    mul!(y, factorize(A), x, α, β)

    return y
end






function LinearAlgebra.factorize(A::InteractionTensor)
    T = eltype(A)
    m, n = size(A)
    S = promote_type(T, Complex{Float32})
    tmp = Vector{S}(undef, m + n - 1)
    copyto!(tmp, A.vc)
    copyto!(tmp, m + 1, Iterators.reverse(A.vr), 1, n - 1)
    dft = plan_fft!(tmp)
    return ToeplitzFactorization{T,typeof(A),S,typeof(dft)}(dft * tmp, similar(tmp), dft)
end













"""
    ToeplitzFactorization

Factorization of a InteractionTensor matrix using FFT.
"""
struct ToeplitzFactorization{T<:Number,A<:InteractionTensor{T},S<:Number,P<:Plan{S}} <: Factorization{T}
    vcvr_dft::Vector{S}
    tmp::Vector{S}
    dft::P
end


function mul!(
    y::StridedVector, A::ToeplitzFactorization, x::StridedVector, α::Number, β::Number
)
    n = length(x)
    m = length(y)
    vcvr_dft = A.vcvr_dft
    N = length(vcvr_dft)
    if m > N || n > N
        throw(DimensionMismatch(
            "InteractionTensor factorization does not match size of input and output vector"
        ))
    end

    T = Base.promote_eltype(y, A, x, α, β)
    tmp = A.tmp
    dft = A.dft
    @inbounds begin
        for i in 1:n
            tmp[i] = x[i]
        end
        for i in (n+1):N
            tmp[i] = 0
        end
        mul!(tmp, dft, tmp)
        for i in 1:N
            tmp[i] *= vcvr_dft[i]
        end
        dft \ tmp
        if iszero(β)
            for i in 1:m
                y[i] = α * maybereal(T, tmp[i])
            end
        else
            for i in 1:m
                y[i] = muladd(α, maybereal(T, tmp[i]), β * y[i])
            end
        end
    end

    return y
end

function mul!(
    C::StridedMatrix, A::ToeplitzFactorization, B::StridedMatrix, α::Number, β::Number
)
    l = size(B, 2)
    if size(C, 2) != l
        throw(DimensionMismatch("input and output matrices must have same number of columns"))
    end
    for j = 1:l
        mul!(view(C, :, j), A, view(B, :, j), α, β)
    end
    return C
end












# Circulant
"""
    Circulant

A circulant matrix.
"""
struct Circulant{T<:Number} <: InteractionTensor{T}
    vc::Vector{T}
end

"""
    Circulant(vc::AbstractVector{<:Number})

Create a circulant matrix from its first column `vc`.
"""
Circulant(vc::AbstractVector) = Circulant{eltype(vc)}(vc)

"""
    Circulant(A::AbstractMatrix)

Create a circulant matrix from the first column of matrix `A`.
"""
Circulant(A::AbstractMatrix) = Circulant{eltype(A)}(A)
Circulant{T}(A::AbstractMatrix) where {T<:Number} = Circulant{T}(A[:,1])

const CirculantFactorization{T<:Number} = ToeplitzFactorization{T,Circulant{T}}
function LinearAlgebra.factorize(C::Circulant)
    T = eltype(C)
    vc = C.vc
    S = promote_type(T, Complex{Float32})
    tmp = Vector{S}(undef, length(vc))
    copyto!(tmp, vc)
    dft = plan_fft!(tmp)
    return ToeplitzFactorization{T,typeof(C),S,typeof(dft)}(dft * tmp, similar(tmp), dft)
end

convert(::Type{InteractionTensor{T}}, A::Circulant) where {T} = convert(Circulant{T}, A)
convert(::Type{Circulant{T}}, A::Circulant) where {T} = Circulant(convert(Vector{T}, A.vc))


function size(C::Circulant, dim::Integer)
    if 1 <= dim <= 2
        return length(C.vc)
    else
        error("arraysize: dimension out of range")
    end
end

function getindex(C::Circulant, i::Integer, j::Integer)
    n = size(C, 1)
    if i > n || j > n
        error("BoundsError()")
    end
    return C.vc[mod(i - j, length(C.vc)) + 1]
end






