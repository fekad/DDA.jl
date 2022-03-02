
# CG
# MINRES
# GMRES
# IDRs

struct Direct <: AbstractMethod end
struct BiCGStabl <: AbstractMethod end
struct BiCGStablFFT <: AbstractMethod end

C_sca(sol::AbstractSolution) = C_ext(sol) - C_abs(sol)

# TODO: storing the alg and its parameters
# TODO: storing the prepared (factorised) interaction matrix
# TODO: do not trim A*E because it is useful for the nearfield calculations

struct GridSolution <: AbstractSolution
    P
    alphas
    Einc
    prob
    alg
end
C_abs(sol::GridSolution) = C_abs(norm(sol.prob.Einc.kvec), norm(sol.prob.Einc.E₀), sol.P, sol.alphas)
C_ext(sol::GridSolution)= C_ext(norm(sol.prob.Einc.kvec), norm(sol.prob.Einc.E₀), sol.Eincs, sol.P)



function solve(p::GridProblem, alg::BiCGStablFFT;
    reltol=1e-3, verbose=true, kwargs...)

    # 1. create the coordinates of the dipoles
    # 2. assign the polarizability αj to each dipole
    # occ = DDA.discretize(p.grid, p.scatterer.target)
    # coords = p.grid[occ]
    # alphas = polarisbility(p.scatterer.model, p)
    # coords, occ, alphas = DDA.discretize(p.grid, p.scatterer, p.Einc)

    coords, occ, alphas = discretize(p)


    # 3. calculated the incident field Einc, at each dipole,
    # Einc = field(p.Einc, p.grid[occ])

    Einc = similar(coords, SVector{3,Complex{Float64}})
    for i = eachindex(coords)
        Einc[i] = field(p.Einc, coords[i])
    end


    # 4. assemble the interaction matrix A and

    k = wavenumber(p.Einc)
    A_conv = TensorConvolution(p.grid, occ, k, alphas)


    # 5. solve for P in the system of linear equations

    P = similar(coords, SVector{3,Complex{Float64}})
    fill!(P, zero(SVector{3,Complex{Float64}}))

    bicgstabl!(reinterpret(ComplexF64, P), A_conv, reinterpret(ComplexF64, Einc); reltol=reltol, verbose=verbose)

    return GridSolution(P, alphas, Einc, p, alg)
end



struct DipoleSolution <: AbstractSolution
    P
    prob
    alg
end

C_abs(sol::DipoleSolution) = C_abs(sol.prob.k, sol.prob.E0, sol.P, sol.prob.alphas)
C_ext(sol::DipoleSolution) = C_ext(sol.prob.k, sol.prob.E0, sol.prob.Eincs, sol.P)

function solve(p::DipoleProblem, alg::BiCGStabl; reltol=1e-3, verbose=true, kwargs...)

    # # 1. create the coordinates of the dipoles
    # # 2. assign the polarizability αj to each dipole
    # # occ = DDA.discretize(p.grid, p.scatterer.target)
    # # coords = p.grid[occ]
    # # alphas = polarisbility(p.scatterer.model, p)
    # # coords, occ, alphas = DDA.discretize(p.grid, p.scatterer, p.Einc)
    #
    # coords, occ, alphas = discretize(p)
    #
    #
    # # 3. calculated the incident field Einc, at each dipole,
    # # Einc = field(p.Einc, p.grid[occ])
    #
    # Einc = similar(coords, SVector{3,Complex{Float64}})
    # for i = eachindex(coords)
    #     Einc[i] = field(p.Einc, coords[i])
    # end
    #
    #
    # # 4. assemble the interaction matrix A and
    #
    # k = wavenumber(p.Einc)
    # A_conv = TensorConvolution(p.grid, occ, k, alphas)

    A = DDA.interactions(p.k, p.dipoles, p.alphas)


    # 5. solve for P in the system of linear equations

    T = Float64
    V = SVector{3,Complex{T}}

    P = similar(p.dipoles, V)
    fill!(P, zero(V))

    bicgstabl!(reinterpret(Complex{T}, P), A, reinterpret(Complex{T}, p.Eincs); reltol=reltol, verbose=verbose, kwargs...)

    return DipoleSolution(P, p, alg)
end


