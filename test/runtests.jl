using SafeTestsets, Test

@safetestset "discretisaition" begin
    include("test_grid.jl")
end

@safetestset "Incindent field" begin
    include("test_incident_field.jl")
end

@safetestset "Polarizability model" begin
    include("test_polarizability.jl")
end

@safetestset "dipole" begin
    include("test_dipole.jl")
end

