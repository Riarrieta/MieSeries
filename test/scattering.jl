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
