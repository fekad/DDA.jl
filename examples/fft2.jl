using LinearAlgebra
using FFTW
# using AbstractFFTs
using IterativeSolvers

using DDA


k = 1
dx = 1
alpha = Inf
alpha = 2.1

Nx = 8;
Ny = 8;
Nz = 8;

grid = CartesianGrid(Nx, Ny, Nz)



inds = rand(Bool, Nx, Ny, Nz)
r =  grid[inds]
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

function test_fft_multiply(grid, X, k_norm, alpha, inds)


    Ĝ = Array{ComplexF64}(undef, 2Nx, 2Ny, 2Nz, 6)
    Ê = zeros(ComplexF64, 2Nx, 2Ny, 2Nz, 3)

    for i in eachindex(grid)
        Rij = grid[i]
        A = DDA.calc_Ajk(k_norm, Rij)
        Ĝ[i, 1] = A[1, 1] # xx
        Ĝ[i, 2] = A[1, 2] # xy
        Ĝ[i, 3] = A[1, 3] # xz
        Ĝ[i, 4] = A[2, 2] # yy
        Ĝ[i, 5] = A[2, 3] # yz
        Ĝ[i, 6] = A[3, 3] # zz

    end

    Ĝ[1, 1, 1, :] .= zero(ComplexF64)

    # Periodic extension
    g = [+1 -1 -1 +1 +1 +1  # x
        +1 -1 +1 +1 -1 +1  # y
        +1 +1 -1 +1 -1 +1  # z
        +1 +1 -1 +1 -1 +1  # xy
        +1 -1 +1 +1 -1 +1  # xz
        +1 -1 -1 +1 +1 +1]  # yz

    Ĝ[Nx+1, :, :, :] .= zero(ComplexF64)
    Ĝ[:, Ny+1, :, :] .= zero(ComplexF64)
    Ĝ[:, :, Nz+1, :] .= zero(ComplexF64)

    for i = 1:6
        Ĝ[Nx+2:2Nx, 1:Ny, 1:Nz, i] = Ĝ[Nx:-1:2, 1:Ny, 1:Nz, i] * g[1, i]       # x
        Ĝ[1:Nx, Ny+2:2Ny, 1:Nz, i] = Ĝ[1:Nx, Ny:-1:2, 1:Nz, i] * g[2, i]       # y
        Ĝ[1:Nx, 1:Ny, Nz+2:2Nz, i] = Ĝ[1:Nx, 1:Ny, Nz:-1:2, i] * g[3, i]       # z
        Ĝ[Nx+2:2Nx, Ny+2:2Ny, 1:Nz, i] = Ĝ[Nx:-1:2, Ny:-1:2, 1:Nz, i] * g[4, i] # xy
        Ĝ[Nx+2:2Nx, 1:Ny, Nz+2:2Nz, i] = Ĝ[Nx:-1:2, 1:Ny, Nz:-1:2, i] * g[5, i] # xz
        Ĝ[1:Nx, Ny+2:2Ny, Nz+2:2Nz, i] = Ĝ[1:Nx, Ny:-1:2, Nz:-1:2, i] * g[6, i] # yz
        Ĝ[Nx+2:2Nx, Ny+2:2Ny, Nz+2:2Nz, i] = Ĝ[Nx:-1:2, Ny:-1:2, Nz:-1:2, i]   # xyz
    end


    # g1 = copy(Ĝ)
    # g2 = copy(Ĝ)
    #
    #
    #
    # fftw_flags = FFTW.ESTIMATE
    # FG = FFTW.plan_fft!(view(Ĝ, :, :, :, 1), flags = fftw_flags)
    # #
    # # for i = 1:6
    # #     F * view(Ĝ, :, :, :, i)
    # # end
    fft!(Ĝ, 1:3)


    Ê = zeros(ComplexF64, 2Nx, 2Ny, 2Nz, 3)

    fftw_flags = FFTW.ESTIMATE
    F = FFTW.plan_fft!(Ê, 1:3, flags = fftw_flags)

    cinds = findall(inds)

    return do_multiply(F, Ĝ, Ê, X, alpha, cinds, Nx, Ny, Nz)

end

function do_multiply(F, Ĝ, Ê, X, alpha, cinds, Nx, Ny, Nz)

    for i = 1:length(cinds)
        Ê[cinds[i], :] = X[i]
    end


    # for i = 1:3
    #     F * view(Ê, :, :, :, i)
    # end
    F * Ê

    Etmp = similar(Ê)

#
#     # !! Matrix vector multiplication
#     # !! -> x component
#     # !! -> y component
#     # !! -> z component
#
#     Etmp[:, :, :, 1] = @. Ĝ[:, :, :, 1] * Ê[:, :, :, 1] + Ĝ[:, :, :, 2] * Ê[:, :, :, 2] + Ĝ[:, :, :, 3] * Ê[:, :, :, 3]
#     Etmp[:, :, :, 2] = @. Ĝ[:, :, :, 2] * Ê[:, :, :, 1] + Ĝ[:, :, :, 4] * Ê[:, :, :, 2] + Ĝ[:, :, :, 5] * Ê[:, :, :, 3]
#     Etmp[:, :, :, 3] = @. Ĝ[:, :, :, 3] * Ê[:, :, :, 1] + Ĝ[:, :, :, 5] * Ê[:, :, :, 2] + Ĝ[:, :, :, 6] * Ê[:, :, :, 3]
#
#     # !! X, Y and Z components stored in X_3D
#     Ê[:] = Etmp[:]


    # !! Matrix vector multiplication
    for k=1:2Nz, j=1:2Ny, i=1:2Nx
        tmp = Ê[i,j,k,:]

        # !! -> x component
        Ê[i,j,k,1] = Ĝ[i,j,k,1] * tmp[1] + Ĝ[i,j,k,2] * tmp[2] + Ĝ[i,j,k,3] * tmp[3]
        # !! -> y component
        Ê[i,j,k,2] = Ĝ[i,j,k,2] * tmp[1] + Ĝ[i,j,k,4] * tmp[2] + Ĝ[i,j,k,5] * tmp[3]
        # !! -> z component
        Ê[i,j,k,3] = Ĝ[i,j,k,3] * tmp[1] + Ĝ[i,j,k,5] * tmp[2] + Ĝ[i,j,k,6] * tmp[3]
    end


    # inds = CartesianIndices(size(Ĝ)[1:3])
    # for i in inds
    #
    # end
    #

    # for i = 1:3
    #     inv(FG) * view(Ê, :, :, :, i)
    # end

    inv(F) * Ê


    for i = 1:length(cinds)
        Ê[cinds[i], :] += 1/alpha * X[i]
    end



    return Ê[cinds, :]
end

function test_direct_multiply(grid, X, k_norm, alpha,inds)
    r =  grid[inds]


    A = DDA.interactions(k_norm, r,  ones(length(r)) * alpha)
    # A[diagind(A)] .= 0

    m = A * reinterpret(ComplexF64, X)
    return transpose(reshape(m, 3,:))
end

test_direct_multiply(grid, X, k_norm, alpha,inds) ≈ test_fft_multiply(grid, X, k_norm, alpha,inds)
