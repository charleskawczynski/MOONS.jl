using Test
using UnPack
using MOONS.Grids
using MOONS.Fields

@testset "Fields" begin
    FT = Float64
    a, b = FT(0.0), FT(1.0)
    n = 10
    β = FT(1.1)
    c1 = Coordinates(a,b,n; warpfun=Roberts_both, args=(β,))
    c2 = Coordinates(a,b,n; warpfun=Roberts_both, args=(β,))
    c3 = Coordinates(a,b,n; warpfun=Roberts_both, args=(β,))
    grid = Grid((c1,c2,c3))

    cent = CellCenter(grid)
    @test size(cent.data) == (n+2,n+2,n+2)

    corn = CellCorner(grid)
    @test size(corn.data) == (n+3,n+3,n+3)

    face = CellFace(grid, 1)
    @test size(face.data) == (n+3,n+2,n+2)
    face = CellFace(grid, 2)
    @test size(face.data) == (n+2,n+3,n+2)
    face = CellFace(grid, 3)
    @test size(face.data) == (n+2,n+2,n+3)

    edge = CellEdge(grid, 1)
    @test size(edge.data) == (n+2,n+3,n+3)
    edge = CellEdge(grid, 2)
    @test size(edge.data) == (n+3,n+2,n+3)
    edge = CellEdge(grid, 3)
    @test size(edge.data) == (n+3,n+3,n+2)

end

@testset "Fields iterator" begin

    FT = Float64
    a, b = FT(0.0), FT(1.0)
    n = 3
    β = FT(1.01)
    grid = Cube(a, b, n; warpfun=Roberts_both, args=(β,))

    # CellCenter field
    f = CellCenter(grid)
    fi = CellCenter(grid)
    @inbounds for local_grid in GridField(grid,fi)
        @unpack ijk,x,y,z = local_grid
        fi[ijk] = sin(3*π*x)* sin(3*π*y)* sin(3*π*z)
    end
    @inbounds begin
    for i in 1:grid.c[1].c.s, j in 1:grid.c[2].c.s, k in 1:grid.c[3].c.s
        x = grid.c[1].c.h[i]
        y = grid.c[2].c.h[j]
        z = grid.c[3].c.h[k]
        f[i,j,k] = sin(3*π*x)* sin(3*π*y)* sin(3*π*z)
    end; end
    @test all(fi .≈ f)

    # CellCorner field
    f = CellCorner(grid)
    fi = CellCorner(grid)
    @inbounds for local_grid in GridField(grid,fi)
        @unpack ijk,x,y,z = local_grid
        fi[ijk] = sin(3*π*x)* sin(3*π*y)* sin(3*π*z)
    end
    @inbounds begin
    for i in 1:grid.c[1].n.s, j in 1:grid.c[2].n.s, k in 1:grid.c[3].n.s
        x = grid.c[1].n.h[i]
        y = grid.c[2].n.h[j]
        z = grid.c[3].n.h[k]
        f[i,j,k] = sin(3*π*x)* sin(3*π*y)* sin(3*π*z)
    end; end
    @test all(fi .≈ f)

    # CellFace field
    for dir in 1:3
        f = CellFace(grid, dir)
        fi = CellFace(grid, dir)
        @inbounds for local_grid in GridField(grid,fi)
            @unpack ijk,x,y,z = local_grid
            fi[ijk] = sin(3*π*x)* sin(3*π*y)* sin(3*π*z)
        end
        nx = dir == 1 ? grid.c[1].n.s : grid.c[1].c.s
        ny = dir == 2 ? grid.c[2].n.s : grid.c[2].c.s
        nz = dir == 3 ? grid.c[3].n.s : grid.c[3].c.s
        @inbounds begin
        for i in 1:nx, j in 1:ny, k in 1:nz
            x = dir == 1 ? grid.c[1].n.h[i] : grid.c[1].c.h[i]
            y = dir == 2 ? grid.c[2].n.h[j] : grid.c[2].c.h[j]
            z = dir == 3 ? grid.c[3].n.h[k] : grid.c[3].c.h[k]
            f[i,j,k] = sin(3*π*x)* sin(3*π*y)* sin(3*π*z)
        end; end
        @test all(fi .≈ f)
    end

    # CellEdge field
    for dir in 1:3
        f = CellEdge(grid, dir)
        fi = CellEdge(grid, dir)
        @inbounds for local_grid in GridField(grid,fi)
            @unpack ijk,x,y,z = local_grid
            fi[ijk] = sin(3*π*x)* sin(3*π*y)* sin(3*π*z)
        end
        nx = dir == 1 ? grid.c[1].c.s : grid.c[1].n.s
        ny = dir == 2 ? grid.c[2].c.s : grid.c[2].n.s
        nz = dir == 3 ? grid.c[3].c.s : grid.c[3].n.s
        @inbounds begin
        for i in 1:nx, j in 1:ny, k in 1:nz
            x = dir == 1 ? grid.c[1].c.h[i] : grid.c[1].n.h[i]
            y = dir == 2 ? grid.c[2].c.h[j] : grid.c[2].n.h[j]
            z = dir == 3 ? grid.c[3].c.h[k] : grid.c[3].n.h[k]
            f[i,j,k] = sin(3*π*x)* sin(3*π*y)* sin(3*π*z)
        end; end
        @test all(fi .≈ f)
    end

end
