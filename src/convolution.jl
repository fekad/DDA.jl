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
    G[:, 1, 1, 1] .= zero(ComplexF64)

    extend!(Ĝ, Nx, Ny, Nz)
    fft!(Ĝ, 2:4)

    tmp = zeros(ComplexF64, 3, 2Nx, 2Ny, 2Nz)
    F = FFTW.plan_fft!(tmp, 2:4, flags = FFTW.ESTIMATE)

    inds = findall(occ)
    TensorConvolution{typeof(F)}(Ĝ, F, inds, alphas, tmp)
end


A::TensorConvolution * x = mul!(similar(x), A, x)

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

    inv(A.F) * A.tmp # inverse FFT

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
