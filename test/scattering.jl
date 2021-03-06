using Test
using StaticArrays
using MieSeries

# cross product
mycross(a,b) = SVector(a[2]*b[3]-a[3]*b[2],
                       a[3]*b[1]-a[1]*b[3],
                       a[1]*b[2]-a[2]*b[1])
# norm
mynorm(a) = (b=abs.(a); sqrt(b[1]^2+b[2]^2+b[3]^2))

@testset "Test boundary condition" begin
    TOL = 1e-7     # tolerance
    npoints = 10   # points per dimension in mesh
    n_terms = 80   # number of terms in Mie series
    for a in [1, 5]
        # spherical mesh
        θrange = range(0, 2π, length=npoints)[1:end-1]
        ϕrange = range(0, π, length=npoints)[2:end-1]
        normals = [SVector(sin(ϕ)*cos(θ), 
                        sin(ϕ)*sin(θ), 
                        cos(ϕ)) for ϕ in ϕrange for θ in θrange]
        mesh = a.*normals
        for k in [1,5,10]
            Es = mieseries(mesh; k, a, n_terms)  # scattered field
            Ei = incident_planewave.(mesh; k)     # incident field
            Et = Es + Ei     # total field
            # test boundary condition 
            # n × Eᵗ = 0  on Γ
            for i in eachindex(Es)
                n = normals[i]
                @test mynorm(mycross(n, Et[i])) < TOL
            end
            # test conjugation
            @test conj(Es) == mieseries(mesh; k, a, n_terms, expjω=true)
            @test conj(Ei) == incident_planewave.(mesh; k, expjω=true)
        end
    end
end

@testset "Far field: monostatic RCS" begin
    TOL = 1e-3     # tolerance
    n_terms = 80   # number of terms in Mie series
    a = 0.1        # sphere radius

    obs_point = SVector(1e-7, 1e-7, -1.)          # observation point
    radius_in_λ = range(0.01, 1.6, length=50)
    λ_list = a ./ radius_in_λ                     # wavelengths
    for λ in λ_list
        k = 2π/λ       # wavenumber
        E_far = MieSeries.mieseries(obs_point; k, a, n_terms, farfield=true)
        E_far_norm2 = sum(abs2.(E_far))
        σ_mono = 4*π*E_far_norm2
        @test abs(σ_mono-MieSeries.sphere_monostatic_rcs(;k, a, n_terms)) < TOL
    end
end

@testset "Mie series debug mode" begin
    a = 1          # sphere radius
    n_terms = 80   # number of terms in Mie series
    k = 5.5        # wavenumber
    point = SVector(3.0, -5.0, 6.0)

    Es = mieseries(point; k, a, n_terms)
    Ex_terms, Ey_terms, Ez_terms = MieSeries.mieseries_debug(point; k, a, n_terms)
    @test length(Ex_terms) == n_terms
    @test length(Ey_terms) == n_terms
    @test length(Ez_terms) == n_terms
    @test sum(Ex_terms) ≈ Es[1]
    @test sum(Ey_terms) ≈ Es[2]
    @test sum(Ez_terms) ≈ Es[3]
end
