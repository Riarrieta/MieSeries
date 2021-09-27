"""
    incident_planewave(mesh::Vector{SVector{3,Float64}}; k, expjω=false)

Evaluates the electric field of a x-polarized plane wave traveling along 
the z-direction. This field impinging on a PEC sphere give rise to the scattered field
computed by [`mieseries`](@ref).
# Arguments
- `point::SVector{3,Float64}`: point in space where the electric field is evaluated.
- `k`: wavenumber.
- `expjω`: set to `true` (`false`) if `exp(jω)` (`exp(-jω)`) time-harmonic dependence is used.
"""
function incident_planewave(point::SVector{3,Float64}; k, expjω=false)
    _, _, z = point
    if expjω
        return SVector{3,ComplexF64}(exp(-im*k*z), 0, 0)
    else
        return SVector{3,ComplexF64}(exp(im*k*z), 0, 0)
    end
end


"""
    mieseries(mesh::Vector{SVector{3,Float64}}; k, a=1, n_terms=20, expjω=false)
    mieseries(point::SVector{3,Float64}; k, a=1, n_terms=20, expjω=false)

Evaluates the scattered field produced by an incident plane wave ([`incident_planewave`](@ref)) 
on a PEC sphere. This scattered field is known as Mie series.
Refer to: C. A. Balanis. Advanced Engineering Electromatnetics (2013), Chapter 11.8. 
# Arguments
- `mesh::Vector{SVector{3,Float64}}`: list of points in space where the electric field is evaluated.
- `k`: wavenumber.
- `a`: radius of the PEC sphere.
- `n_terms`: the number of terms to compute in the Mie series.
- `expjω`: set to `true` (`false`) if `exp(jω)` (`exp(-jω)`) time-harmonic dependence is used.
"""
function mieseries(mesh::Vector{SVector{3,Float64}}; k, a=1, n_terms=20, expjω=false)
    isreal(k) || @error "Mie series is not working for complex wavenumbers"
    # compute Ricatti-Hankel functions, second kind 
    rh, rh_d, _ = riccatihankel2_and_derivatives(1:n_terms, k*a)
    # compute coefficients aₙ, bₙ, cₙ
    an = [(1.0im)^(-n)*(2n + 1)/n/(n + 1) for n in 1:n_terms]
    bn = -an.*real(rh_d)./rh_d
    cn = -an.*real(rh)./rh
    # compute scattered field
    n_points = length(mesh)
    Es = Vector{SVector{3,ComplexF64}}(undef,n_points)
    for i in 1:n_points
        Es[i] = _mieseries(mesh[i], bn, cn; k, n_terms)
    end
    # conjugate result if exp(-jω) time-harmonic dependence is used
    if !expjω
        Es = conj(Es)
    end
    return Es
end
function mieseries(point::SVector{3,Float64}; k, a=1, n_terms=20, expjω=false)
    Es_list = mieseries([point]; k, a, n_terms, expjω)
    return first(Es_list)
end

function _mieseries(point::SVector{3,Float64}, bn, cn; k, n_terms, _debug::Union{Val{false}, Val{true}}=Val(false))
    x, y, z = point
    # convert to spherical coordinates
    r = sqrt(x^2 + y^2 + z^2)
    kr = k*r
    cos_theta = z / r
    sin_theta = sqrt(1 - cos_theta^2)
    cos_phi = x / (r*sin_theta)
    sin_phi = y / (r*sin_theta)
    # compute Ricatti-Hankel functions, second kind 
    rh, rh_d, rh_dd = riccatihankel2_and_derivatives(1:n_terms, kr)
    # compute associated Legendre functions, m=1, for l=1,..,n_term
    lf, lf_d = MieSeries.associated_legendre(cos_theta; l_max=n_terms, m=1)
    # compute scattered field in spherical coordinates
    # return all the terms of the series if !(_debug === Val(false))
    func = (_debug===Val(false)) ? sum : identity
    Es_r = -im*cos_phi*func(@. bn*(rh_dd + rh)*lf)
    Es_theta = cos_phi/kr*func(@. im*bn*rh_d*sin_theta*lf_d - cn*rh*lf/sin_theta)
    Es_phi = sin_phi/kr*func(@. im*bn*rh_d*lf/sin_theta - cn*rh*sin_theta*lf_d)
    # convert to cartesian coordinates
    Es_x = sin_theta*cos_phi*Es_r + cos_theta*cos_phi*Es_theta - sin_phi*Es_phi
    Es_y = sin_theta*sin_phi*Es_r + cos_theta*sin_phi*Es_theta + cos_phi*Es_phi
    Es_z = cos_theta*Es_r - sin_theta*Es_theta
    return Es_x, Es_y, Es_z
end

"""
    mieseries_debug(point::SVector{3,Float64}; k, a=1, n_terms=20, expjω=false) 

Same as ([`mieseries`](@ref)), but returns all the terms of the series, for each cartesian component.
"""
function mieseries_debug(point::SVector{3,Float64}; k, a=1, n_terms=20, expjω=false)
    isreal(k) || @error "Mie series is not working for complex wavenumbers"
    # compute Ricatti-Hankel functions, second kind 
    rh, rh_d, _ = riccatihankel2_and_derivatives(1:n_terms, k*a)
    # compute coefficients aₙ, bₙ, cₙ
    an = [(1.0im)^(-n)*(2n + 1)/n/(n + 1) for n in 1:n_terms]
    bn = -an.*real(rh_d)./rh_d
    cn = -an.*real(rh)./rh
    # compute scattered field
    Es_x, Es_y, Es_z = _mieseries(point, bn, cn; k, n_terms, _debug=Val(true))
    # conjugate result if exp(-jω) time-harmonic dependence is used
    if !expjω
        Es_x, Es_y, Es_z = conj.((Es_x,Es_y,Es_z))
    end
    return Es_x, Es_y, Es_z
end
