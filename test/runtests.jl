using Test
using SafeTestsets

@safetestset "Legendre tests" begin include("legendre.jl") end
@safetestset "Bessel tests" begin include("bessel.jl") end
@safetestset "Mie series tests" begin include("scattering.jl") end