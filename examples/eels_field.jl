using SpecialFunctions: besselk
import PhysicalConstants.CODATA2018: ħ, ε_0, e as q

K₀(z) = besselk(0, z)
K₁(z) = besselk(1, z)

struct EELSField
    k::Float64 # wavenumber
    v::Float64 # charge velocity
end


function field(E::EELSField, r)
    r_q = [0, 0]
    d = r[1:2] - r_q
    dnorm = norm(d)

    A = q.val / (2π * ε_0.val) * E.k / E.v * exp(1im * E.k * r[3])

    return [
        A * d[1] / dnorm * K₁(E.k * dnorm),
        A * d[2] / dnorm * K₁(E.k * dnorm),
        - A * 1im *  K₀(E.k * dnorm)
    ]
end

# struct Planewave
#     k::Float64 # wavenumber
#     e::SVector{2, Float64} # polaristaion
# end

# # TODO Translation and rotation of the input positions
# struct Rotated
#     E::Field
#     θ::Float64
#     ϕ::Float64
# end
#
# struct TranslatedRotated
#     E::Field
#     θ::Float64
#     ϕ::Float64
#     r::SVector{3, Float64}
# end
using LinearAlgebra
using PlotlyJS

E = EELSField(1,.5)

N = 300
r = Diagonal([2π,2π,4π]) * rand(3, N) .- π

E = similar(r)
for (r, E) in zip(eachcol(r), eachcol(E))
    # @show field(Einc, r)
    E[:] = real.(field(Einc, r))
end

# amplitude.(Einc, collect(eachcol(r)))

PlotlyJS.plot(
    cone(
        z=r[3,:],
        x=r[1,:],
        y=r[2,:],
        u=E[1,:],
        v=E[2,:],
        w=E[3,:],
        sizemode="relative",
        sizeref=1,
        anchor="tail",
        colorscale=colors.Blues_8
    ),
    Layout(
        scene=attr(layoutdomain_x=[0,1])
    )
)