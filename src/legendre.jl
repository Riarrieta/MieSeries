
"""
    associated_legendre(x; l_max, m)

Returns the associated Legendre polynomials `Pₗᵐ(x)` and their
first derivative `(d/dx)Pₗᵐ(x)`, for `l=1,...,l_max`. The polynomials
are unnormalized and include the Condon–Shortley phase.
"""
function associated_legendre(x; l_max, m)
    abs(x) < 1 || @error "Argument 'x' must be between -1 and 1."
    plm = Legendre.Plm(0:l_max, m, x)
    # compute first derivate using recurrence
    # (x^2 - 1)(d/dx)p_nm(x) = n * x * p_nm(x) - (n + m) * p_{n-1, m}(x)
    plm_d = similar(plm)
    for n in 1:l_max
        n_index = n + 1
        plm_d[n_index] = (n*x*plm[n_index] - (n+m)*plm[n_index-1]) / (x^2 - 1)
    end
    # return values that correspond to l=1,...,l_max
    result_plm = plm[2:end]
    result_plm_d = plm_d[2:end]
    @assert length(result_plm) == length(result_plm_d) == l_max
    return result_plm, result_plm_d
end