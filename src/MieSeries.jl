module MieSeries
using StaticArrays
import SpecialFunctions
import Legendre

export mieseries, incident_planewave

include("legendre.jl")
include("bessel.jl")
include("scattering.jl")

end # module

