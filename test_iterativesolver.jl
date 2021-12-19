using IterativeSolvers
using Quaternions
import Random: AbstractRNG, SamplerType

# rand not implemented in Quaternions, needed to run, but ignore this
Base.rand(r::AbstractRNG, ::SamplerType{Quaternion{T}}) where {T<:Real} =
    Quaternion(rand(r,T), rand(r,T), rand(r,T), rand(r,T))

q = quat(1.0, 1.0, 1.0, 1.0)
A = [q q'; q' q]
b = [quat(1.0), quat(1.0)]
bicgstabl(A, b)



using IterativeSolvers
using StaticArrays

A = rand(SMatrix{3, 3, Float64, 9}, 2, 2)
b = rand(SVector{3, Float64}, 2)

cg(A, b)
bicgstabl(A, b)

Ar = reshape(permutedims(reshape(reinterpret(reshape, Float64, A), 3,3,2,2),(1,3,2,4)),6,6)
br = reinterpret(Float64, b)

yr = cg(Ar,br, verbose=true)
yr = bicgstabl(Ar, br; verbose=true, abstol=1e-8, reltol=1e-8)
Ar*yr

using Tensors

b = [rand(Vec{3, Float64}), rand(Vec{3, Float64})]
A = reshape(Array([rand(Tensor{2,3,Float64}), rand(Tensor{2,3,Float64}), rand(Tensor{2,3,Float64}), rand(Tensor{2,3,Float64})]), 2,2)

cg(A, b)
bicgstabl(A, b)









using BenchmarkTools
# using AbstractFFTs
using FFTW

N = 250
G1 = rand(N, N, N, 3)
F1 = FFTW.plan_fft(G1, 1:3, flags = FFTW.ESTIMATE)

@btime $F1*$G1

G11 = rand(N, N, N)
F11 = FFTW.plan_fft(G11, flags = FFTW.ESTIMATE)

@btime $F11*$G11

G2 = rand(3, N, N, N)
F2 = FFTW.plan_fft(G2, 2:4, flags = FFTW.ESTIMATE)

@btime $F2*$G2

# it is about 10-20% percent difference






using Einsum

## If the destination array is preallocated, then use =:

A = ones(5, 6, 7) # preallocated space
X = randn(5, 2)
Y = randn(6, 2)
Z = randn(7, 2)

# Store the result in A, overwriting as necessary:
@einsum A[i, j, k] = X[i, r] * Y[j, r] * Z[k, r]

## If destination is not preallocated, then use := to automatically create a new array for the result:

X = randn(5, 2)
Y = randn(6, 2)
Z = randn(7, 2)

# Create new array B with appropriate dimensions:
@einsum B[i, j, k] := X[i, r] * Y[j, r] * Z[k, r]




N = 100
A = copy(reshape(1.:(3^2 * N^2), 3, 3, N, N)) # preallocated space
B = copy(reshape(1.:3N, 3,  N))

X = zeros(size(B))

fun1(A,B) = @einsum X[i, k] := A[i, j, k, l] * B[j, l]
fun2(A,B,N) = reshape(reshape(permutedims(A,(1,3,2,4)), 3N, 3N) * reshape(B, 3N), 3, N)

Ar = copy(reshape(permutedims(A,(1,3,2,4)), 3N, 3N))
@btime Br = copy(reshape(B, 3N))
fun3(A,B,N) = reshape(A*B, 3, N)
fun4(A) = copy(A)
@assert fun1(A,B) == fun2(A,B,N)
@assert fun1(A,B) == fun3(Ar,Br,N)

@btime fun1($A, $B)
@btime fun2($A, $B, $N)
@btime fun3($Ar, $Br, $N)
@btime fun4(A);

# 112.026 μs (1 allocation: 2.50 KiB)
# 130.619 μs (10 allocations: 706.02 KiB)
# 27.120 μs (3 allocations: 2.59 KiB)
# 40.250 μs (2 allocations: 703.20 KiB)


using HybridArrays, StaticArrays
using DifferentialEquations
using BenchmarkTools, Random

N = 50
Dim = 2

# just some random dynamics to show speed benefit of static arrays
function f_ode(du, u, p, t)
    N = size(u,2)
    du[:] .= 0.
    for i = 1:N
        for j = 1:i-1
            v = u[:,i] - u[:,j]
            dist_squared = sum( x -> x^2, v)
            du[:,i] -= v * dist_squared
            du[:,j] += v * dist_squared
        end
    end
end

# same as f_ode but for vector of svectors
function f_ode_vec_svec(du, u, p, t)
    N = size(u,1)
    du[:] .-= du[:]  #(sorry, I don't know how to zero a vector of svectors fast.)
    for i = 1:N
        for j = 1:i-1
            v = u[i] - u[j]
            dist_squared = sum( x -> x^2, v)
            du[i] -= v * dist_squared
            du[j] += v * dist_squared
        end
    end
end


Random.seed!(1)
z0 = rand(Dim,N)
z0_hybrid = HybridArray{Tuple{Dim,StaticArrays.Dynamic()}}(z0)

Random.seed!(1)
z0_vec_svec = rand(SVector{Dim,Float64},N)


prob          = ODEProblem(f_ode,          z0,          (0., 10.))
prob_hybrid   = ODEProblem(f_ode,          z0_hybrid,   (0., 10.))
prob_vec_svec = ODEProblem(f_ode_vec_svec, z0_vec_svec, (0., 10.))

sol          = @btime solve(prob,          Euler(), dt=0.01) seconds=2
# 682.401 ms (33153208 allocations: 1.65 GiB)
sol_hybrid   = @btime solve(prob_hybrid,   Euler(), dt=0.01) seconds=2
# 7.337 ms (14070 allocations: 3.04 MiB)
sol_vec_svec = @btime solve(prob_vec_svec, Euler(), dt=0.01) seconds=2
# 5.019 ms (20074 allocations: 4.78 MiB)


sol2          = solve(prob,          Rosenbrock23())
# works
sol2_hybrid   = solve(prob_hybrid,   Rosenbrock23())
# MethodError: copyto!(::DiffEqBase.DiffEqBC{Array{Float64,2}}, ::Base.Broadcast.Broadcasted{HybridArrays.HybridArrayStyle{2},Tuple{Base.OneTo{Int64},Base.OneTo{Int64}},typeof(muladd),Tuple{Base.Broadcast.Broadcasted{HybridArrays.HybridArrayStyle{2},Nothing,typeof(DiffEqBase.ODE_DEFAULT_NORM),Tuple{HybridArray{Tuple{2,StaticArrays.Dynamic()},Float64,2,2,Array{Float64,2}},Float64}},Float64,Float64}})
# is ambiguous.
sol2_vec_svec = solve(prob_vec_svec, Rosenbrock23())
# ArgumentError: Cannot create a dual over scalar type SArray{Tuple{2},Float64,1,2}.
#  If the type behaves as a scalar, define FowardDiff.can_dual.

# using Plots
# plot( sol[1,1,:], sol[2,1,:])
# plot( sol_hybrid[1,1,:], sol_hybrid[2,1,:])
# plot( getindex.(sol_vec_svec[1,:],1) , getindex.(sol_vec_svec[1,:],2))
# # plots look all the same
#
# plot( sol )
# # works
# plot( sol_hybrid )
# # MethodError: no method matching u_n(...)
# # Note: Here a ::CartesianIndex{2}  does not match the expected ::Int64
# plot( sol_vec_svec )
# # MethodError: no method matching one(::Type{SArray{Tuple{2},Float64,1,2}})


#
#
sol   = @btime solve(prob,   ImplicitEuler(), dt=0.01); seconds=2
#   760.786 ms (5048795 allocations: 534.85 MiB)
# 2
#
sol_hybrid   = @btime solve(prob_hybrid,   ImplicitEuler(), dt=0.01); seconds=2
#   27.109 ms (10635 allocations: 976.16 KiB)
# 2
#
sol_hybrid   = @btime solve(prob_hybrid,   Rosenbrock23(), dt=0.01); seconds=2
#   38.912 ms (5940 allocations: 622.34 KiB)
# 2
#
sol   = @btime solve(prob,   Rosenbrock23(), dt=0.01); seconds=2
#   884.613 ms (3313352 allocations: 723.02 MiB)
# 2
#
#
