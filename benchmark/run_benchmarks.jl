# https://github.com/Wikunia/ConstraintSolver.jl/blob/master/benchmark/run_benchmarks.jl

using BenchmarkTools


function f(x)
    while true
        x+=1
    end
end


using BenchmarkPlots, StatsPlots

b = @benchmarkable sin(10) samples=1000 seconds=2

t=run(b)
plot(t, yscale=:log10)


# Execution time
# number of dipoles VS sec
#
r = Int[]
for i = 2:100
    @show i^3
    append!(r, i^3)
end




using BenchmarkTools
using Plots


function get_primes(to)
    possible_primes = collect(2:to)
    # iterate over the primes and remove multiples
    i = 1
    while i <= length(possible_primes)
        prime = possible_primes[i]
        possible_primes = filter(n->n <= prime || n % prime != 0, possible_primes)
        i += 1
    end
    return possible_primes
end

# create a list of numbers from 2:100
possible_primes = get_primes(100)

xs = Float64[]
ys = Float64[]

for n in range(100, stop=5000, length=20)
    println("n: $n")
    t = @benchmark get_primes($n) seconds=1
    push!(xs, n)
    # convert from nano seconds to seconds
    push!(ys, minimum(t).time / 10^9)
end

plot(xs, ys,
    legend=:topleft,
    label="Calculating prime numbers",
    xaxis="# of dipoles in",
    yaxis="Time [seconds]"
)

