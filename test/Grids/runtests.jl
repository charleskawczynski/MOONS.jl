using Test
using MOONS.Grids
const G = Grids

@testset "Coordinates init" begin
  for FT in (Float32, Float64)
    a = FT(0.0)
    b = FT(1.0)
    n = 10
    c = Coordinates(a,b,n)
    @test c.c.Δh[1] ≈ FT(1/10)
    @test c.n.Δh[1] ≈ FT(1/10)
    @test h_extreme(c.n, Min()) ≈ FT(0)
    @test h_extreme(c.n, Max()) ≈ FT(1)
    sprint(show, c)
    @test c.c == coords(c, Center1D())
    @test c.n == coords(c, Vertex1D())
    # FT==Float64 && @show c
  end
end

@testset "Helper funcs" begin
    @test n_hat(Min()) == -1
    @test binary(Min()) == 0
    @test ghost_vec(Min()) == (1, 0, 0)
    @test ghost_dual(Min()) == (1, 0)

    @test n_hat(Max()) == 1
    @test binary(Max()) == 1
    @test ghost_vec(Max()) == (0, 0, 1)
    @test ghost_dual(Max()) == (0, 1)
end

@testset "Stretch funcs" begin
    FT = Float64
    Re = FT(100)
    Ha = FT(200)
    hmin = FT(0)
    hmax = FT(1)
    δ = hmax*99/100
    @test all(G.ReynoldsBL(Re, FT[0, 0, 0], FT[1, 1, 1]) .≈ [1.0540925533894598, 1.0540925533894598, 1.0540925533894598])
    @test all(G.Re_Ha_BL(Re, Ha, FT[0, 0, 0], FT[1, 1, 1]) .≈ [1.002509414234171, 1.002509414234171, 1.002509414234171])
    @test G.Re_Ha_BL_1D(Re, Ha, hmin, hmax) ≈ 1.002509414234171
    @test G.Reynolds_BL_1D(Re, hmin, hmax) ≈ 1.0540925533894598
    @test G.Hartmann_BL_1D(Ha, hmin, hmax) ≈ 1.002509414234171
    @test G.RobertsBL1D(δ, hmax) ≈ FT(1000)
    @test G.β_Δh_big(hmin, hmax, 100, FT(0.02)) ≈ 1.0443565558255983
    @test G.β_Δh_both(hmin, hmax, 100, FT(0.001)) ≈ 1.0211975072167816
    @test G.β_Δh_small(hmin, hmax, 100, FT(0.001)) ≈ 1.0218301816802464

    Δh = FT(0.1)
    @test all(uniform_left(hmin, Δh, 3) .≈ [-0.3, -0.2, -0.1, 0.0])
    @test all(uniform_right(hmin, Δh, 3) .≈ [0.0, 0.1, 0.2, 0.3])
end

@testset "Grid init" begin
  for FT in (Float32, Float64)
    a = FT(0.0)
    b = FT(1.0)
    n = 10
    β = FT(1.1)
    c1 = Coordinates(a,b,n; warpfun=Roberts_both, args=(β,))
    c2 = Coordinates(a,b,n; warpfun=Roberts_both, args=(β,))
    c3 = Coordinates(a,b,n; warpfun=Roberts_both, args=(β,))
    grid = Grid((c1,c2,c3))

    yc = FT(0.5)
    τ = FT(10)
    c1 = Coordinates(a,b,n; warpfun=cluster, args=(yc,τ))
    c2 = Coordinates(a,b,n; warpfun=cluster, args=(yc,τ))
    c3 = Coordinates(a,b,n; warpfun=cluster, args=(yc,τ))
    grid = Grid((c1,c2,c3))

    c1 = Coordinates(a,b,n; warpfun=Roberts_right, args=(β,))
    c2 = Coordinates(a,b,n; warpfun=Roberts_left, args=(β,))
    c3 = Coordinates(a,b,n; warpfun=Roberts_both, args=(β,))
    grid = Grid((c1,c2,c3))
    sprint(show, grid)
    # FT==Float64 && @show grid

    @test over_points(c1, All(), Vertex1D()) == 1:13
    @test over_points(c1, Interior(), Vertex1D()) == 3:11
    @test over_points(c1, Boundary(), Vertex1D()) == 2:12
    @test over_points(c1, Ghost(), Vertex1D()) == (1:1, 13:13)

    @test over_points(c1, All(), Center1D()) == 1:12
    @test over_points(c1, Interior(), Center1D()) == 2:11
    @test over_points(c1, Ghost(), Center1D()) == (1:1, 12:12)
  end
end

@testset "Stretched parameter matching" begin
  for FT in (Float32, Float64)
    a = FT(0.0)
    b = FT(1.0)
    n = 10
    β = FT(1.1)
    c = Coordinates(a,b,n; warpfun=Roberts_both, args=(β,))
    β = G.β_Δh_small(a, b, n, c.n.Δh[1])
    @test β ≈ 1.1036862256343674
  end
end
