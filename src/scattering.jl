
function mieseries(mesh::Vector{SVector{3,Float64}}; k, a=1, n_terms=20, expjω=false)
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

function _mieseries(point::SVector{3,Float64}, bn, cn; k, n_terms)
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
    Es_r = -im*cos_phi*sum(@. bn*(rh_dd + rh)*lf)
    Es_theta = cos_phi/kr*sum(@. im*bn*rh_d*sin_theta*lf_d - cn*rh*lf/sin_theta)
    Es_phi = sin_phi/kr*sum(@. im*bn*rh_d*lf/sin_theta - cn*rh*sin_theta*lf_d)
    # convert to cartesian coordinates
    Es_x = sin_theta*cos_phi*Es_r + cos_theta*cos_phi*Es_theta - sin_phi*Es_phi
    Es_y = sin_theta*sin_phi*Es_r + cos_theta*sin_phi*Es_theta + cos_phi*Es_phi
    Es_z = cos_theta*Es_r - sin_theta*Es_theta
    return SVector{3,ComplexF64}(Es_x, Es_y, Es_z)
end

function incident_planewave(mesh::Vector{SVector{3,Float64}}; k, expjω=false)
    n_points = length(mesh)
    planewave = Vector{SVector{3,ComplexF64}}(undef, n_points)
    for i in 1:n_points
        _, _, z = mesh[i]
        planewave[i] = SVector{3,ComplexF64}(exp(-im*k*z), 0, 0)
    end
    # conjugate result if exp(-jω) time-harmonic dependence is used
    if !expjω
        planewave = conj(planewave)
    end
    return planewave
end