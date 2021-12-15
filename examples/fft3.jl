using LinearAlgebra
using FFTW
# using AbstractFFTs
using IterativeSolvers

using DDA

# TODO: return value: P(3, Nx, Ny, Nz) or P(3, N_occupied)

k_norm = 1
dx = 1
alpha = Inf
alpha = 2.1

Nx = 32;
Ny = 32;
Nz = 32;

grid = CartesianGrid(Nx, Ny, Nz)


inds = rand(Bool, Nx, Ny, Nz)
X = rand(SVector{3,ComplexF64}, sum(inds))


# 1. calc A (Nx, Ny, Nz, 3, 3) or (Nx, Ny, Nz, 6) becuse of symmetry
# 2. add padding (2Nx, 2Ny, 2Nz, 6)
# 3. fft (2Nx, 2Ny, 2Nz, 6)
# 3. construct (2Nx, 2Ny, 2Nz, 3) vector from dipoles and indecies
# 4. plan for fft (ifft) of (2Nx, 2Ny, 2Nz, 3)
# 5. multiplication
# 6. ifft (2Nx, 2Ny, 2Nz, 3)
# 7. set zero for vacuum
# 8. polarisibality

# function test_fft_multiply(grid, X, k_norm, alpha, inds)
    Nx, Ny, Nz = size(grid)

    # Ĝ = Array{SymmetricTensor{2,3,ComplexF64,6}}(undef,  2Nx, 2Ny, 2Nz)
    Ĝ = Array{ComplexF64}(undef, 6, 2Nx, 2Ny, 2Nz)
    # Ê = Array{SVector{3,ComplexF64}}(undef,  2Nx, 2Ny, 2Nz)
    # Ê = Array{ComplexF64}(undef, 3, 2Nx, 2Ny, 2Nz)
    Ê = zeros(ComplexF64, 3, 2Nx, 2Ny, 2Nz)

    for i=1:Nx, j=1:Ny, k=1:Nz
        Rij = grid[i,j,k]
        A = DDA.calc_Ajk(k_norm, Rij)
        Ĝ[1,i,j,k] = A[1, 1] # xx
        Ĝ[2,i,j,k] = A[1, 2] # xy
        Ĝ[3,i,j,k] = A[1, 3] # xz
        Ĝ[4,i,j,k] = A[2, 2] # yy
        Ĝ[5,i,j,k] = A[2, 3] # yz
        Ĝ[6,i,j,k] = A[3, 3] # zz
        # Ĝ[i] = SymmetricTensor{2,3,ComplexF64}((A[1, 1], A[1, 2], A[1, 3], A[2, 2], A[2, 3], A[3, 3]))
        # Ĝ[:,i] = (A[1, 1], A[1, 2], A[1, 3], A[2, 2], A[2, 3], A[3, 3])
    end

    Ĝ[:, 1, 1, 1] .= zero(ComplexF64)
    # Ĝ[1, 1, 1] = zero(SymmetricTensor{2,3,ComplexF64})

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

    fft!(Ĝ, 2:4)

    fftw_flags = FFTW.ESTIMATE
    F = FFTW.plan_fft!(Ê, 2:4, flags = fftw_flags)

    cinds = findall(inds)

    return do_multiply(F, Ĝ, Ê, X, alpha, cinds, Nx, Ny, Nz)

# end
# using TimerOutputs
# const to = TimerOutput()

# function do_multiply(F, Ĝ, Ê, X, alpha, cinds, Nx, Ny, Nz)

    # @timeit to "nest 1" begin
        # @timeit to "copy1" begin
            fill(Ê,0)
            for i = 1:length(cinds)
                Ê[:, cinds[i]] = X[i]
            end
        # end
        # @timeit to "fft" begin

            F * Ê
        # end
        # @timeit to "mul" begin

            # !! Matrix vector multiplication
            for k=1:2Nz, j=1:2Ny, i=1:2Nx
                tmp = Ê[:,i,j,k]

                # !! -> x component
                Ê[1,i,j,k] = Ĝ[1,i,j,k] * tmp[1] + Ĝ[2,i,j,k] * tmp[2] + Ĝ[3,i,j,k] * tmp[3]
                # !! -> y component
                Ê[2,i,j,k] = Ĝ[2,i,j,k] * tmp[1] + Ĝ[4,i,j,k] * tmp[2] + Ĝ[5,i,j,k] * tmp[3]
                # !! -> z component
                Ê[3,i,j,k] = Ĝ[3,i,j,k] * tmp[1] + Ĝ[5,i,j,k] * tmp[2] + Ĝ[6,i,j,k] * tmp[3]
            end
        # end
        # @timeit to "ifft" begin

            inv(F) * Ê
        # end
        # @timeit to "alpha" begin

            for i = 1:length(cinds)
                Ê[:, cinds[i]] += 1/alpha * X[i]
            end
        # end
        # @timeit to "copy2" begin

           out =  Ê[:,cinds]
        # end
    # end
    # return out
# end


function test_direct_multiply(grid, X, k_norm, alpha,inds)
    r =  grid[inds]
    A = DDA.interactions(k_norm, r,  ones(length(r)) * alpha)
    # A[diagind(A)] .= 0

    m = A * reinterpret(ComplexF64, X)
    return reshape(m, 3,:)
end

test_direct_multiply(grid, X, k_norm, alpha,inds) ≈ test_fft_multiply(grid, X, k_norm, alpha,inds)

