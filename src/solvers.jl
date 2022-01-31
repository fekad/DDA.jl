
# CG
# MINRES
# GMRES
# IDRs

struct Direct <: AbstractMethod end
struct BiCGStabl <: AbstractMethod end
struct BiCGStablFFT <: AbstractMethod end


# TODO: storing the alg and its parameters
# TODO: storing the prepared (factorised) interaction matrix
# TODO: do not trim A*E because it is useful for the nearfield calculations

struct GridSolution <: AbstractSolution
    P
    alphas
    prob
    alg
end


# TODO: make it unitless (x=k*a) (Works only with spheres)
function solve(p::GridProblem, alg::BiCGStablFFT;
    reltol=1e-3, verbose=true, kwargs...)
    # @show p kwargs

#     # 1. create the coordinates of the dipoles,
#     # coords, inds = DDA.dipoles(p.grid, p.scatterer.target)
#     occ = DDA.discretize(p.grid, p.scatterer.target)
#
#     # 2. assign the polarizability αj to each dipole
#     alphas = polarisbility(p.scatterer.model, p)


    occ, alphas = DDA.discretize(p, p.scatterer)
    # coords, alphas, occ = DDA.discretize(p.scatterer, p.grid, p.Einc)


    # 3. calculated the incident field Einc, at each dipole,
    Einc = field(p.Einc, p.grid[occ])


    # 4. assemble the interaction matrix A and
    k = wavenumber(p.Einc)
    A_conv = TensorConvolution(p.grid, occ, k, alphas)

    # 5. solve for P in the system of linear equations

    Ei = reinterpret(ComplexF64, Einc)
    P = bicgstabl(A_conv, Ei; reltol=reltol, verbose=verbose)

    P = reinterpret(SVector{3,ComplexF64}, P)

    return GridSolution(P, alphas, p, alg)
end


