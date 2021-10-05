using StaticArrays
using MieSeries
using Plots

n_terms = 80   # number of terms in Mie series
a = 2          # sphere radius

# normalized monostatic RCS (σ/(πa^2))
obs_point = SVector(1e-7, 1e-7, -1.)          # observation point
radius_in_λ = range(0.01, 1.6, length=100)
λ_list = a ./ radius_in_λ                     # wavelengths
σ_mono_list = Float64[]
for λ in λ_list
    k = 2π/λ       # wavenumber
    E_far = MieSeries.mieseries(obs_point; k, a, n_terms, farfield=true)
    E_far_norm2 = sum(abs2.(E_far))
    σ_mono = 4*π*E_far_norm2
    push!(σ_mono_list, σ_mono)
end
σ_mono_list /= π*a^2;  # normalize 

## plot
plot(radius_in_λ, σ_mono_list, yscale=:log10, legend=false)
hline!([1], linestyle=:dash)
xlabel!("Sphere radius in wavelengths")
ylabel!("Normalized monostatic RCS")
