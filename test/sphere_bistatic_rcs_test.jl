using StaticArrays
using MieSeries
using Plots

n_terms = 80   # number of terms in Mie series
f = 300e6      # frequency
λ = 3e8 / f    # wavelength
k = 2π / λ     # wavenumber
a = λ          # sphere radius

## Vertical x-polarized
TOL = 1e-5   # tolerance
θlist = range(0+TOL, π-TOL, length=100)
mesh = [SVector(sin(θ), 0., cos(θ)) for θ in θlist]

E_far = MieSeries.mieseries(mesh; k, a, n_terms, farfield=true)
E_far_norm2 = [sum(abs2.(E)) for E in E_far]
σ = 4*π*E_far_norm2
σ_dBsm = 10*log10.(σ/(π*a^2))

# plot
plot(θlist, σ_dBsm, legend=false)
title!("Bistatic RCS vertical plane x-polarized")
xlabel!("θ (rad)")
ylabel!("bistatic RCS (dBsm)")

## Vertical y-polarized
θlist = range(0+TOL, π-TOL, length=100)
mesh = [SVector(0., -sin(θ), cos(θ)) for θ in θlist]

E_far = MieSeries.mieseries(mesh; k, a, n_terms, farfield=true)
E_far_norm2 = [sum(abs2.(E)) for E in E_far]
σ = 4*π*E_far_norm2
σ_dBsm = 10*log10.(σ/(π*a^2))

# plot
plot(θlist, σ_dBsm, legend=false)
title!("Bistatic RCS vertical plane y-polarized")
xlabel!("θ (rad)")
ylabel!("bistatic RCS (dBsm)")

## Horizontal
ϕlist = range(0, π, length=100)
mesh = [SVector(cos(ϕ), sin(ϕ), 0.) for ϕ in ϕlist]

E_far = MieSeries.mieseries(mesh; k, a, n_terms, farfield=true)
E_far_norm2 = [sum(abs2.(E)) for E in E_far]
σ = 4*π*E_far_norm2
σ_dBsm = 10*log10.(σ/(π*a^2))

# plot
plot(ϕlist, σ_dBsm, legend=false)
title!("Bistatic RCS horizontal plane")
xlabel!("ϕ (rad)")
ylabel!("bistatic RCS (dBsm)")

