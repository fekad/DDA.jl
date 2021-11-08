# https://qutip.org/docs/latest/guide/guide-bloch.html
"""
	(type) Spherical
Representation of a complex value on the Riemann sphere.
"""
struct Spherical{T<:AbstractFloat} <: Number
	lat::T
	lon::T
end
#
# constructors
Spherical{T}(z::Spherical{T}) where {T<:AbstractFloat} = z

latitude(z::Number) = π/2 - 2*acot(abs(z))
function Spherical{T}(z::Number) where T<:AbstractFloat
	θ,ϕ = latitude(z),angle(z)
	Spherical{T}(convert(T,θ),convert(T,ϕ))
end

# Constructors without subtype
"""
	Spherical(latitude,azimuth)
Construct a spherical representation with given `latitude` in [-π/2,π/2] and `azimuth`.
"""
function Spherical(θ::Real,ϕ::Real)
	θ,ϕ = promote(float(θ),float(ϕ))
	Spherical{typeof(θ)}(θ,ϕ)
end
"""
Spherical(z)
Construct a spherical representation of the value `z`.
"""
Spherical(z::Number) = Spherical(latitude(z),angle(z))

# # one and zero
# one(::Type{Spherical{T}}) where T<:Real = Spherical{T}(zero(T),zero(T))
# one(::Type{Spherical}) = one(Spherical{Float64})
# zero(::Type{Spherical{T}}) where T<:Real = Spherical{T}(T(-π/2),zero(T))
# zero(::Type{Spherical}) = zero(Spherical{Float64})
#
# # conversion into standard complex
# function Complex(z::Spherical{S}) where S<:Real
# 	if iszero(z)
# 		zero(Complex{S})
# 	else
# 		cot(π/4-z.lat/2) * exp(complex(zero(z.lon),z.lon))
# 	end
# end
#
# """
# 	S2coord(u::Spherical)
# Convert the spherical value to a 3-vector of coordinates on the unit sphere.
# """
# S2coord(u::Spherical) = [cos(u.lat)*[cos(u.lon),sin(u.lon)];sin(u.lat)]
#
# # basic arithmetic
# function +(u::Spherical,v::Spherical)
# 	if isinf(u)
# 		isinf(v) ? NaN : u
# 	elseif isinf(v)
# 		v
# 	else
# 		Spherical(Complex(u)+Complex(v))  # faster way?
# 	end
# end
# -(u::Spherical) = Spherical(u.lat,cleanangle(u.lon+π))
# -(u::Spherical,v::Spherical) = u + (-v)
# *(u::Spherical,v::Spherical) = Spherical(Polar(u)*Polar(v))   # faster way?
# inv(u::Spherical) = Spherical(-u.lat,cleanangle(-u.lon))
# /(u::Spherical,v::Spherical) = u*inv(v)
#
# # common complex overloads
# angle(u::Spherical) = cleanangle(u.lon)
# function abs(z::Spherical{T}) where T
# 	if iszero(z)
# 		zero(T)
# 	elseif isinf(z)
# 		T(Inf)
# 	else
# 		cot(π/4-z.lat/2)
# 	end
# end
# abs2(u::Spherical) = abs(u)^2
# real(u::Spherical) = abs(u)*cos(u.lon)
# imag(u::Spherical) = abs(u)*sin(u.lon)
# conj(u::Spherical) = Spherical(u.lat,-u.lon)
# sign(u::Spherical) = Spherical(zero(u.lat),u.lon)
#
# # numerical comparisons
# iszero(u::Spherical) = u.lat == convert(typeof(u.lat),-π/2)
# isinf(u::Spherical) = u.lat == convert(typeof(u.lat),π/2)
# isfinite(u::Spherical) = ~isinf(u)
# isapprox(u::Spherical,v::Spherical;args...) = S2coord(u) ≈ S2coord(v)
#
# # pretty output
# show(io::IO,z::Spherical) = print(io,"(latitude = $(z.lat/π)⋅π, angle = $(z.lon/π)⋅π)")
# show(io::IO,::MIME"text/plain",z::Spherical) = print(io,"Complex Spherical: ",z)
#

using RecipesBase

# @recipe function f(z::Array{Polar{T}}) where T
#     projection --> :polar
#     angle.(z),abs.(z)
# end

@userplot BlochSphere

@recipe function f(s::BlochSphere; sphere=true)
    # @show s.args
    # @show length(s.args)
    # @show sphere
    z, = s.args

    delete!(plotattributes,:sphere)
    # markersize --> 1
    aspect_ratio --> 1
    xlims --> (-1, 1)
    ylims --> (-1, 1)
    zlims --> (-1, 1)
    framestyle := :none
    grid := false
    # axis := nothing
    ticks := false
    hover --> :none


    # θ ∈ [0, π],  ϕ ∈ [0, 2π]
    inner(θ, ϕ) = cos(ϕ) * sin(θ), sin(ϕ) * sin(θ), cos(θ)

    function latcurves(θ)
        nvertices = 100
        ϕ = LinRange(0, 2π, nvertices)
        return [inner(θ, ϕ) for θ in θ for θ in θ, ϕ in ϕ ]
    end

    function loncurves(ϕ)
        nvertices = 100
        θ = LinRange(0, π, nvertices)
        return [inner(θ, ϕ) for θ in θ for θ in θ, ϕ in ϕ ]
    end

    n = 4

    for θ = LinRange(0, π, 2n+1)
        @series begin
            group := 2
            seriestype := :path3d
            color := :gray
            # linealpha := 0.5
            linestyle := :solid
            markershape := :none
            linewidth := 1.
            label := ""
            hover := :none
            latcurves(θ)
        end
    end


    for ϕ = LinRange(0, 2π, 2n+1)
        @series begin
            group := 2
            seriestype := :path3d
            color := :gray
            # linealpha := 0.5
            linestyle := :solid
            markershape := :none
            linewidth := 1.
            label := ""
            hover := :none
            loncurves(ϕ)
        end
    end


    # add axes
    axis = [-1, 1]
    other = [0, 0]
    @series begin
        # series_annotations:=["|V>","|H>"] # TODO: no 3d support
        hover := ["|V>", "|H>"]
        label := "|V>, |H>"
        # label := ""
        color := :red
        linewidth := 4.
        group := 2
        hover := :none
        axis, other, other
    end

    @series begin
        hover := ["|A>", "|D>"]
        label := "|A>, |D>"
        # label := ""
        color := :green
        linewidth := 4.
        group := 2
        other, axis, other
    end

    @series begin
        hover := ["|L>", "|R>"]
        label := "|L>, |R>"
        # label := ""
        color := :blue
        linewidth := 4.
        group := 2
        seriestype:=:path3d
        arrow:=true
        # seriestype:=:scatter3d
        # mode:="markers+text"
        # text := ["|L>", "|R>"]
        # textposition:="top center"
        # textfont_family="Raleway, sans-serif"
        other, other, axis
    end


    # # Labels
    # text!(scene, "|R>", textsize=30, position = Point3f0(0, 0,  1.4), align = (:center, :center))
    # text!(scene, "|L>", textsize=30, position = Point3f0(0, 0, -1.4), align = (:center, :center))
    #
    # text!(scene, "|H>", textsize=30, position = Point3f0( 1.4, 0, 0), align = (:center, :center))
    # text!(scene, "|V>", textsize=30, position = Point3f0(-1.4, 0, 0), align = (:center, :center))
    #
    # text!(scene, "|D>", textsize=30, position = Point3f0(0, 1.4, 0), align = (:center, :center))
    # text!(scene, "|A>", textsize=30, position = Point3f0(0,-1.4, 0), align = (:center, :center))

    @series begin
        # arrow :=true
        # hover := :none
        x = [ cos(z.lat)*cos(z.lon) for z in z ]
        y = [ cos(z.lat)*sin(z.lon) for z in z ]
        z = [ sin(z.lat) for z in z ]
        x, y, z
    end

end

using Plots
plotlyjs()
# gr()
# z = Spherical(1.,1.)
# # plot(z,sphere=true)
# plot([Spherical(0.,0.),z],sphere=(8,8))


zl = collect(LinRange(50-50im,-50+50im,601));
# plot(Spherical.(zc/2),l=3,leg=false)  # plot on the Riemann sphere
# plot(Spherical.(-1 .+ zl), l=3, sphere=true, label="jhge")
f = blochsphere(Spherical.(-1 .+ zl); l=3, sphere=false, label="jhge")
