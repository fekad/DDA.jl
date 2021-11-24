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

A = interaction(k,dipoles, alpha*ones(length(dipoles)))

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



# Allocate vectors for the matrice vector operator
X3D = zeros(ComplexF64, 2Nx, 2Ny, 2Nz, 3)

Xcomp = zeros(ComplexF64, 2Nx, 2Ny, 2Nz)
Ycomp = zeros(ComplexF64, 2Nx, 2Ny, 2Nz)
Zcomp = zeros(ComplexF64, 2Nx, 2Ny, 2Nz)




# //////////////////////////////////////////////////////////
#  building "toeplitz matix"
# Fill the interaction tensor A(i,j,k,1:6)
A = zeros(ComplexF64, Nx, Ny, Nz, 6)

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
