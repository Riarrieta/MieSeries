using Test
using MieSeries

@testset "Test Hankel" begin
    x = 4
    nu = 0:4
    sh, sh_d, sh_dd = MieSeries.sphericalhankel2_and_derivatives(nu, x)
    @test sh[1] ≈ -0.18920062382698205-0.16341090521590299im  # order 0
    @test sh[2] ≈ 0.11611074925915747-0.2300533501309578im
    @test sh[3] ≈ 0.27628368577135015-0.009129107382315343im
    @test sh[4] ≈ 0.22924385795503022+0.21864196590306362im
    @test sh[5] ≈ 0.12489306564995283+0.39175254771267665im
    @test sh_d[1] ≈ -0.11611074925915747+0.2300533501309578im  # order 0
    @test sh_d[2] ≈ -0.2472559984565608-0.04838423015042409im
    @test sh_d[3] ≈ -0.09110201506935513-0.22320651959422128im
    @test sh_d[4] ≈ 0.047039827816319935-0.22777107328537896im
    @test sh_d[5] ≈ 0.0731275258925892-0.2710487187377822im
    @test sh_dd[1] ≈ 0.247255998456560+0.0483842301504240im  # order 0
    @test sh_dd[2] ≈ 0.022031093626517+0.2254887964398001im
    @test sh_dd[3] ≈ -0.127126296072416+0.1173089519110577im
    @test sh_dd[4] ≈ -0.080830878396917+0.0592250451669235im
    @test sh_dd[5] ≈ -0.005340496533806+0.2334624962970602im

    for (index,i) in enumerate(nu)
        @test MieSeries.sphericalhankel2_and_derivatives(i, x) == (sh[index], sh_d[index], sh_dd[index])
    end
    @test MieSeries.sphericalhankel1_and_derivatives(nu, x) == conj.((sh,sh_d,sh_dd))
    
    sh2, sh_d2, sh_dd2 = MieSeries.sphericalhankel2_and_derivatives(1:4, x)
    @test sh2 == sh[2:end]
    @test sh_d2 == sh_d[2:end]
    @test sh_dd2 == sh_dd[2:end]

    # Spherical Bessel differential equation
    for n in nu
        index = n+1
        @test sh_dd[index] ≈ -(2*x*sh_d[index] + (x^2 - n*(n+1))*sh[index])/x^2
    end
end

@testset "Test Riccati-Bessel Hankel" begin
    x = 3.3
    nu = 0:2
    rh, rh_d, rh_dd = MieSeries.riccatihankel2_and_derivatives(nu, x)
    @test rh[1] ≈ -0.15774569-0.98747977im  # order 0
    @test rh[2] ≈ 0.93967804-0.45698199im
    @test rh[3] ≈ 1.01199846+0.5720416im
    @test rh_d[1] ≈ -0.98747977+0.15774569im  # order 0
    @test rh_d[2] ≈ -0.44249662-0.84900038im
    @test rh_d[3] ≈ 0.32634564-0.80367387im
    # Riccati-Bessel differential equation
    for n in nu
        index = n+1
        @test rh_dd[index] ≈ -(1-n*(n+1)/x^2)*rh[index]
    end
end