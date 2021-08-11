using Test
using MieSeries

@testset "Test 1" begin
    x = 0.3
    l_max = 3
    m = 0
    p, dp = MieSeries.associated_legendre(x, l_max, m)
    @test p[1] ≈ 0.3
    @test p[2] ≈ -0.365
    @test p[3] ≈ -0.3825
    @test dp[1] ≈ 1
    @test dp[2] ≈ 0.9
    @test dp[3] ≈ -0.825
end

@testset "Test 2" begin
    x = -0.3
    l_max = 3
    m = 0
    p, dp = MieSeries.associated_legendre(x, l_max, m)
    @test p[1] ≈ -0.3
    @test p[2] ≈ -0.365
    @test p[3] ≈ 0.3825
    @test dp[1] ≈ 1
    @test dp[2] ≈ -0.9
    @test dp[3] ≈ -0.825
end

@testset "Test 3" begin
    x = 0.3
    l_max = 3
    m = 1
    p, dp = MieSeries.associated_legendre(x, l_max, m)
    @test p[1] ≈ -0.953939201416945649152621586023226540254623425250545753908
    @test p[2] ≈ -0.858545281275251084237359427420903886229161082725491178517
    @test p[3] ≈ 0.7869998411689801605509128084691618957100643258317002469742
    @test dp[1] ≈ 0.3144854510165755
    @test dp[2] ≈ -2.578780698335919
    @test dp[3] ≈ -4.55217690346493
end

@testset "Test 4" begin
    x = -0.3
    l_max = 3
    m = 1
    p, dp = MieSeries.associated_legendre(x, l_max, m)
    @test p[1] ≈ -0.953939201416945649152621586023226540254623425250545753908
    @test p[2] ≈ 0.858545281275251084237359427420903886229161082725491178517
    @test p[3] ≈ 0.7869998411689801605509128084691618957100643258317002469742
    @test dp[1] ≈ -0.3144854510165755
    @test dp[2] ≈ -2.578780698335919
    @test dp[3] ≈ 4.55217690346493
end

@testset "Test 5" begin
    x = 0.8
    l_max = 50
    m = 1
    p, dp = MieSeries.associated_legendre(x, l_max, m)
    @test p[48] ≈ 5.96348895739654
    @test p[49] ≈ 2.40302362973485973633385371736914932219523828838
    @test p[50] ≈ -2.201121967241330605081500626464831230994512944496
    @test dp[48] ≈ 328.9545358407495 
    @test dp[49] ≈ 566.5997821783876 
    @test dp[50] ≈ 584.9974550170323 
end