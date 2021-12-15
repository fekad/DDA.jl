# see https://discourse.julialang.org/t/vector-of-upper-triangle/7764

function vec_triu(M::AbstractMatrix{T}) where T
    m, n = size(M)
    m == n || throw(error("not square"))
    l = n*(n+1) ÷ 2
    v = Vector{T}(undef, l)
    k = 0
    for i in 1:n
        v[k .+ (1:i)] .= M[1:i, i]
        k += i
    end
    v
end

function vec_triu_collect(M::AbstractMatrix{T}) where T
    m, n = size(M)
    m == n || throw(error("not square"))
    [M[i, j] for j in axes(M, 1) for i in 1:j]
end

function vec_triu_inbounds(M::AbstractMatrix{T}) where T
    m, n = size(M)
    m == n || throw(error("not square"))
    l = n*(n+1) ÷ 2
    v = Vector{T}(undef, l)
    k = 0
    @inbounds for i in 1:n
        v[k .+ (1:i)] .= M[1:i, i]
        k += i
    end
    v
end

function vec_triu_loop(M::AbstractMatrix{T}) where T
    m, n = size(M)
    m == n || throw(error("not square"))
    l = n*(n+1) ÷ 2
    v = Vector{T}(undef, l)
    k = 0
    @inbounds for i in 1:n
        for j in 1:i
            v[k + j] = M[j, i]
        end
        k += i
    end
    v
end


x = Float64.(reshape(1:9, 3, :))

using BenchmarkTools: @benchmark

@benchmark vec_triu($x)
@benchmark vec_triu_collect($x)
@benchmark vec_triu_inbounds($x)
@benchmark vec_triu_loop($x)