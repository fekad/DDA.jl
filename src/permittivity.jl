# using Plots
# plotlyjs()

#=[md]
Material models
# Dielectrics and Conductors: This includes perfect (lossless) dielectrics, conductive (lossy) materials, and perfect electric conductors.
# $$ \varepsilon (\omega) = \varepsilon $$

# Permittivity
# Electric susceptibility
=#




abstract type AbstractPermittivity end
abstract type ConstantPermittivity <: AbstractPermittivity end
abstract type TablePermittivity <: AbstractPermittivity end


# # only real values
# struct DielectricConstant <: ConstantPermittivity
#     eps::Float64
# end
# permittivity(e::DielectricConstant) = e.eps
# permittivity(e::DielectricConstant, lambda) = permittivity(e)
#
# # complex (dielectric + absorption)
# struct PermittivityConstant <: ConstantPermittivity
#     eps::ComplexF64
# end
# permittivity(e::PermittivityConstant) = e.eps
# permittivity(e::PermittivityConstant, lambda) = permittivity(e)


# complex (real or complex)
struct Permittivity{T<:Number} <: ConstantPermittivity
    eps::T
end

permittivity(e::Permittivity) = e.eps
permittivity(e::Permittivity, lambda) = permittivity(e)

struct DrudeModel

end

# Interpolated frequency dependent permitivirty
struct PermittivityTable <: TablePermittivity
    lambda::Vector{Float64}
    eps::Vector{ComplexF64}
    itp_real::Spline1D
    itp_imag::Spline1D
    function PermittivityTable(lambda::Vector{Float64}, eps::Vector{ComplexF64})

        inds = sortperm(lambda)
        # itp_real = interpolate((lambda[inds],), real(eps[inds]), Gridded(Linear()))
        # itp_imag = interpolate((lambda[inds],), imag(eps[inds]), Gridded(Linear()))

        itp_real = Spline1D(lambda[inds], real(eps[inds]))
        itp_imag = Spline1D(lambda[inds], imag(eps[inds]))

        new(lambda, eps, itp_real, itp_imag)
    end
end

permittivity(e::PermittivityTable, lambda) = e.itp_real(lambda) + e.itp_imag(lambda) * im

# TODO: proper broadcasting
reflective_index(e::AbstractPermittivity, lambda) = sqrt(permittivity(e, lambda))
