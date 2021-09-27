using StaticArrays
using MieSeries
using Plots

a = 5          # sphere radius
n_terms = 80   # number of terms in Mie series
k = 5          # wavenumber
point = SVector(3.0, -5.0, 6.0)

Ex_terms, Ey_terms, Ez_terms = MieSeries.mieseries_debug(point; k, a, n_terms)
Ex_abs, Ey_abs, Ez_abs = abs.(Ex_terms), abs.(Ey_terms), abs.(Ez_terms)

plot(Ex_abs, yscale=:log10, label="Ex")
plot!(Ey_abs, label="Ey")
plot!(Ez_abs, label="Ez")
