using Test
using UnPack
using MOONS.Grids
using MOONS.Fields
using MOONS.FFTs

@testset "FFTs" begin

    # Solve ∇² u = f
    # where f = -2*(3*π)^2*cos(3*π*x)*cos(3*π*y)
    #       u = cos(3*π*x)*cos(3*π*y))

    FT = Float64
    n = 128
    cx = Coordinates(FT(0),FT(1),n)
    cy = Coordinates(FT(0),FT(1),n)
    cz = Coordinates(FT(0),FT(1),2)
    g = Grid((cx,cy,cz))

    f_exact = CellCenter(g)
    u_exact = CellCenter(g)
    f = CellCenter(g)
    u = CellCenter(g)
    u_err = CellCenter(g)

    @inbounds for local_grid in GridField(g,f_exact)
        @unpack ijk,x,y = local_grid
        u_exact[ijk] = cos(3*π*x)*cos(3*π*y)
        f_exact[ijk] = -2*(3*π)^2*cos(3*π*x)*cos(3*π*y)
    end

    ∇²!(f, u_exact, g)
    extrap!(f, g)
    solve!(u, f_exact, g, 3)

    @test max(average(u, 3)...) ≈ 1.0003012584260522
    @test max(abs.(diff(u; dims=3))...) ≈ 0

    u_err .= abs.(u_exact - u);
    assign!(u_err, 0, Ghost())
    @test max(u_err...) ≈ 0.0004518490779501505
end
