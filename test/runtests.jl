using SafeTestsets, Test

@safetestset "Polarizability model" begin
    include("test_polarizability.jl")
end

@safetestset "Incindent field" begin
    include("test_incident_field.jl")
end

@safetestset "dipole" begin
    include("test_dipole.jl")
end

@safetestset "permitivity" begin
    include("test_permitivity.jl")
end