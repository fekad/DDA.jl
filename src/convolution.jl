# TODO: adding parametic types like here: https://github.com/JuliaMatrices/ToeplitzMatrices.jl/blob/master/src/ToeplitzMatrices.jl
# TODO: improve structure: Array of symmetric tensors
"""
    TensorConvolution(...)

A preplanned, circular convolution operator on an M × N matrix of data of type T

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





struct TensorConvolutionFull{P<:AbstractFFTs.Plan{ComplexF64}}
    Ĝ::Array{ComplexF64,4}
    F::P
    inds::Vector{CartesianIndex{3}}
    alphas::ComplexF64
    tmp::Array{ComplexF64,4}
end

Base.eltype(::TensorConvolutionFull) = ComplexF64
Base.size(A::TensorConvolutionFull) = (3 * length(A.inds), 3 * length(A.inds))
Base.size(A::TensorConvolutionFull, i::Int) = size(A)[i]


function TensorConvolutionFull(grid::CartesianGrid, occ, k_norm, alphas)

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
    G[:, 1, 1, 1] .= zero(ComplexF64)

    extend!(Ĝ, Nx, Ny, Nz)
    fft!(Ĝ, 2:4)

    tmp = zeros(ComplexF64, 3, 2Nx, 2Ny, 2Nz)
    F = FFTW.plan_fft!(tmp, 2:4, flags = FFTW.ESTIMATE)

    inds = findall(occ)
    TensorConvolutionFull{typeof(F)}(Ĝ, F, inds, alphas, tmp)
end


A::TensorConvolutionFull * x = mul!(similar(x), A, x)

function LinearAlgebra.mul!(y::StridedVector, A::TensorConvolutionFull, x::StridedVector)
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


#
# function extend!(Ĝ::Array{ComplexF64,4}, Nx, Ny, Nz)
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
# end
#

# Future improvements:

# using Tensors
# using Tensorial

# function TensorConvolution(grid::CartesianGrid, occ, k_norm, alphas)
#     T = SymmetricTensor{2,3,ComplexF64} # T((1,2,3,4,5,6)) # using Tenors
#     T = SymmetricSecondOrderTensor{3,ComplexF64} # T(1,2,3,4,5,6) # using Tensorial
#
#     Nx, Ny, Nz = size(grid)
#     Ĝ = Array{T}(undef, 2Nx, 2Ny, 2Nz)  # BlockArrays ArraysofStruct # StaticArrays
#     for i in eachindex(Ĝ)
#         r_ij = grid[i]
#         Ĝ[i] = DDA.calc_Ajk(k_norm, r_ij)
#     end
#
#     # TODO: merging this two together and doing fft by slices (1-6) and trimming the data afterwards
#     extend!(Ĝ, Nx, Ny, Nz)
#     fft!(Ĝ, 2:4)  # This doesn't work on arrays of tensors
#
#     tmp = zeros(ComplexF64, 2Nx, 2Ny, 2Nz)
#     F = FFTW.plan_fft!(tmp, 2:4) # reinterpret??
#
#     inds = findall(occ)
#     TensorConvolution2{typeof(F)}(Ĝ, F, inds, alphas, tmp)
# end

# function LinearAlgebra.mul!(y::StridedVector, A::TensorConvolution2, x::StridedVector)
#     x = reinterpret(SVector{3,ComplexF64}, x)
#     y = reinterpret(SVector{3,ComplexF64}, y)
#
#     tmp = padded(A.inds, x)
#
#     A.F * tmp # FFT
#
#     # !! Matrix vector multiplication in the Fourier space
#     for i in eachindex(A.Ĝ)
#         tmp[i] = A.Ĝ[i] * tmp[i]
#     end
#
#     A.F \ tmp # inverse FFT
#
#     # TODO: alpha must be a tensor
#     for i in eachindex(y)
#         y[i] = tmp[A.inds[i]] .+ inv(A.alpha[i]) * x[i]
#     end
#
#     return y
#
# end

# function extend!(Ĝ::Array{...}, Nx, Ny, Nz)
#
#     T = eltype(Ĝ)
#
#
#     for k = 1:Nz, j = 1:Ny
#         Ĝ[Nx+1, j, k] = zero(T)
#     end
#
#     for k = 1:Nz, i = 1:Nx-1
#         Ĝ[i, Ny+1, k] = zero(T)
#     end
#
#     for j = 1:Ny, i = 1:Nx-1
#         Ĝ[i, j, Nz+1] = zero(T)
#     end
#
#
#     sym = [-1 1 1]   # x
#     g = T(sym'*sym)
#     for k = 1:Nz, j = 1:Ny, i = 1:Nx-1
#         Ĝ[Nx+1+i, j, k] = Ĝ[Nx-i, j, k] * g
#     end
#
#     sym = [1 -1 1]   # y
#     g = T(sym'*sym)
#     for k = 1:Nz, j = 1:Ny-1, i = 1:Nx
#         Ĝ[i, Ny+1+j, k] = Ĝ[i, Ny-j, k] * g
#     end
#
#     sym = [1 1 -1]   # z
#     g = T(sym'*sym)
#     for k = 1:Nz-1, j = 1:Ny, i = 1:Nx
#         Ĝ[i, j, Nz+1+k] = Ĝ[i, 1:Ny, Nz-k] * g
#     end
#
#     sym = [-1 -1 1]   # xy
#     g = T(sym'*sym)
#     for k = 1:Nz, j = 1:Ny-1, i = 1:Nx-1
#         Ĝ[Nx+1+i, Ny+1+j, k] = Ĝ[Nx-i, Ny-j, k] * g
#     end
#
#     sym = [-1 1 -1]   # xz
#     g = T(sym'*sym)
#     for k = 1:Nz-1, j = 1:Ny, i = 1:Nx-1
#         Ĝ[Nx+1+i, j, Nz+1+k] = Ĝ[Nx-i, j, Nz-k] * g
#     end
#
#     sym = [1 -1 1]   # yz
#     g = T(sym'*sym)
#     for k = 1:Nz-1, j = 1:Ny-1, i = 1:Nx
#         Ĝ[Nx+i, Ny+1+j, Nz+1+k] = Ĝ[i, Ny-j, Nz-k] * g
#     end
#
#     sym = [1 1 1]   # xyz
#     g = T(sym'*sym)
#     for k = 1:Nz-1, j = 1:Ny, i = 1:Nx
#         Ĝ[Nx+1+i, Ny+1+j, Nz+1+k] = Ĝ[Nx-i, Ny-j, Nz-k]
#     end
#
# end
