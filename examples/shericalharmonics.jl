using SpecialFunctions
using Plots
using SphericalHarmonics

plotlyjs()


# Grids of polar and azimuthal angles
n = 100
theta_range = LinRange(0, π, n)
phi_range = LinRange(0, 2π, 2n)

# x_grid = [x for x = theta_range for y = theta_range]
# a_grid = [y for x = theta_range for y = as]

# xyz = [[sin(t) * sin(p), sin(t) * cos(p), cos(t)] for t =theta_range, p=phi_range]
# Y = [SphericalHarmonics.sphericalharmonic(t, p, l = 3, m = 0, SHType=SphericalHarmonics.RealHarmonics()) for t =theta_range, p=phi_range]

Yx = zeros(length(theta_range), length(phi_range))
Yy = zeros(length(theta_range), length(phi_range))
Yz = zeros(length(theta_range), length(phi_range))
for i = eachindex(theta_range), j = eachindex(phi_range)
    # @show i, j
    l, m = 3,1

    y = abs(SphericalHarmonics.sphericalharmonic(theta_range[i], phi_range[j], l=l, m=m, SHType=SphericalHarmonics.RealHarmonics()))
    # xyz = [sin(theta_range[i]) * sin(phi_range[j]), sin(theta_range[i]) * cos(phi_range[j]), cos(theta_range[i])]
    # Y[i,j, :] = y * xyz
    Yx[i,j] = y * sin(theta_range[i]) * sin(phi_range[j])
    Yy[i,j] = y * sin(theta_range[i]) * cos(phi_range[j])
    Yz[i,j] = y * cos(theta_range[i])
end



surface(Yx, Yy, Yz)


SphericalHarmonics.sphericalharmonic(0, 0, l = 3, m = 0, SHType=SphericalHarmonics.RealHarmonics())
SphericalHarmonics.sphericalharmonic(π/3, π/3, l = 3, m = 0, SHType=SphericalHarmonics.RealHarmonics())

Yxyz = xyz .* Y
