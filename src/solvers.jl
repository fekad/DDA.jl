
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
    Einc
    prob
    alg
end


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


