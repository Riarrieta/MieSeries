
"""
    sphericalhankel2(nu, x)

Spherical Hankel function of the second kind, order `nu`
and argument `x`.
"""
function sphericalhankel2(nu, x)
    jn = SpecialFunctions.sphericalbesselj(nu, x)
    yn = SpecialFunctions.sphericalbessely(nu, x)
    return jn - im*yn
end

"""
    sphericalhankel2_and_derivatives(nu::UnitRange, x)

Evaluates the spherical Hankel function of the second kind, its first and 
second derivative, for orders `nu::UnitRange` and argument `x`.
"""
function sphericalhankel2_and_derivatives(nu::UnitRange, x)
    start = nu.start
    stop = nu.stop
    len = stop-start+1   # number of terms to compute
    sph = sphericalhankel2.(start:stop+2, x)  # compute two more terms
    # compute first derivative using recurrence
    # h_n'(x) = -h_{n+1}(x) + n/x*h_n(z)
    sph_d = similar(sph, len+1)   # compute one more term
    n = stop+1   # order
    for index in len+1:-1:1
        sph_d[index] = -sph[index+1] + n/x*sph[index]
        n -= 1
    end
    # compute second derivative using recurrence
    # h_n''(x) = n / x^2 * (x * h_n'(x) - h_n(x)) - h_{n+1}'(x)
    sph_dd = similar(sph, len)
    n = stop   # order
    for index in len:-1:1
        sph_dd[index] = n / x^2 * (x*sph_d[index] - sph[index]) - sph_d[index+1]
        n -= 1
    end
    return sph[1:len], sph_d[1:len], sph_dd
end

"""
    riccatihankel2_and_derivatives(nu::UnitRange, x)

Evaluates the spherical Riccati-Bessel Hankel function of the second kind, its first and 
second derivative, for orders `nu::UnitRange` and argument `x`.
"""
function riccatihankel2_and_derivatives(nu::UnitRange, x)
    sph, sph_d, sph_dd = sphericalhankel2_and_derivatives(nu, x)
    rh = x.*sph
    rh_d = sph + x.*sph_d
    rh_dd = 2*sph_d + x.*sph_dd
    return rh, rh_d, rh_dd
end