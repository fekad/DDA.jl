using StructArrays
using StaticArrays
using LinearAlgebra
using BenchmarkTools


struct Point3DAoS
  x::Float64
  y::Float64
  z::Float64
end

struct Point3Dv1 <: FieldVector{3, Float64}
    x::Float64
    y::Float64
    z::Float64
end

const Point3Dv2 = SVector{3,Float64}

struct Point3DSoA
    x::Vector{Float64}
    y::Vector{Float64}
    z::Vector{Float64}
end


N = 1_000_000

arrAoSv1 = [Point3DAoS(rand(),rand(),rand()) for i in 1:N]
arrAoSv2 = [Point3Dv1(rand(),rand(),rand()) for i in 1:N]
arrAoSv3 = [Point3Dv2(rand(),rand(),rand()) for i in 1:N]

arrSoAv1 = Point3DSoA(rand(N),rand(N),rand(N))
arrSoAv2 = StructArray(Point3DAoS(rand(),rand(),rand()) for i in 1:N)
arrSoAv3 = rand(3, N)


function normAoSv1(x)
    out = 0.0
    # @simd for i in 1:length(x)
    for i in 1:length(x)
        @inbounds out += sqrt(x[i].x^2 + x[i].y^2+ x[i].z^2)
    end
    out
end

function normAoSv2(x)
    out = 0.0
    # @simd for i in 1:length(x)
    for i in 1:length(x)
        # @inbounds out += sqrt(x[i][1]^2 + x[i][2]^2+ x[i][3]^2)
        @inbounds out += norm(x[i])
    end
    out
end

# normAoSv2(x) = sum(norm(x[i]) for i in 1:length(x))

function normSoAv1(x)
    out = 0.0
    # @simd for i in 1:length(x.x)
    for i in 1:length(x.x)
        @inbounds out += sqrt(x.x[i]^2 + x.y[i]^2+ x.z[i]^2)
    end
    out
end

normSoAv2(x) = sum(norm(x[i]) for i in 1:length(x))
normSoAv3(x) = sum(sqrt.(sum(abs2, x, dims=1)))

function normSoAv4(x)
    out = 0.0
    for i in eachcol(x)
        # @inbounds out += sqrt(i[1]^2 + i[2]^2+ i[3]^2)
        out += norm(i)
    end
    out
end

function normSoAv5(x)
    out = 0.0
    for i in eachcol(x)
        out += sqrt(sum(abs2, i))
    end
    out
end

# 1.536 ms (0 allocations: 0 bytes)
@btime normAoSv1($arrAoSv1);
# 1.569 ms (0 allocations: 0 bytes)
@btime normAoSv1($arrAoSv2);
# 1.515 ms (0 allocations: 0 bytes)
@btime normAoSv2($arrAoSv2);
# 1.515 ms (0 allocations: 0 bytes)
@btime normAoSv2($arrAoSv3);

# 970.068 μs (0 allocations: 0 bytes)
@btime normSoAv1($arrSoAv1);
# 968.515 μs (0 allocations: 0 bytes)
@btime normSoAv1($arrSoAv2);

#  3.641 ms (0 allocations: 0 bytes)
@btime normSoAv2($arrSoAv3);
# 5.850 ms (4 allocations: 15.26 MiB)
@btime normSoAv3($arrSoAv3);
# 1.977 ms (0 allocations: 0 bytes)
@btime normSoAv4($arrSoAv3);
# 4.383 ms (0 allocations: 0 bytes)
@btime normSoAv5($arrSoAv3);

# without @simd











using GeometryBasics
using LinearAlgebra
using BenchmarkTools


N = 1_000_000
arr = [Point3(rand(),rand(),rand()) for i in 1:N]

function mynorm(x::AbstractVector{P}) where P<:Point3{T<:Real} where T
    out = zero(T)
    @simd for i in 1:length(x)
        @inbounds out += norm(x[i])
    end
    return out
end

function mynorm3(x::Vector{<:Point3})
    out = 0.0
    for i in 1:length(x)
        out += norm(x[i])
    end
    return out
end
mynorm3(arr)
@code_warntype mynorm(arr)
@code_warntype mynorm3(arr)

@btime mynorm($arr);
@btime mynorm3($arr);
@btime norm($arr);




# syntax sugar

a = [x * 2 for x=1:100]

a = map(x -> x * 2, 1:100)

a = map(1:100) do x
    x * 2
end
