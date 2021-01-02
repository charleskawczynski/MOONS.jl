using Test
using LinearAlgebra
using UnPack
using MOONS.Grids
using MOONS.Fields

#####
##### CellCenter <-> CellFace
#####

function interp_simple!(cent::CellCenter, face::CellFace{1}, grid::Grid)
    for k in over_points(grid.c[3], All(), Center1D())
        for j in over_points(grid.c[2], All(), Center1D())
            for i in over_points(grid.c[1], All(), Center1D())
                D = 1//2
                U = 1//2
                cent[i, j, k] = D*face[i, j, k] + U*face[i+1, j, k]
            end
        end
    end
end

function interp_simple!(face::CellFace{1}, cent::CellCenter, grid::Grid)
    face .= 0 # boundaries are not populated.
    for k in over_points(grid.c[3], All(), Center1D())
        for j in over_points(grid.c[2], All(), Center1D())
            for i in over_points(grid.c[1], Interior(), Vertex1D())
                D = grid.c[1].interp.diag[i]
                U = grid.c[1].interp.upper[i]
                face[i, j, k] = D*cent[i, j, k] + U*cent[i+1, j, k]
            end
        end
    end
end

@testset "Interpolations - Simple implementation Center <-> Face" begin
    FT = Float64
    n, a, b, β = 10, FT(0.0), FT(1.0), FT(1.01)
    grid = Cube(a, b, n; warpfun=Roberts_both, args=(β,))

    # Face -> Center
    face = CellFace(grid, 1)
    cent_correct = CellCenter(grid)
    cent = CellCenter(grid)
    @inbounds for local_grid in GridField(grid,face)
        @unpack ijk,x,y,z = local_grid
        face[ijk] = sin(3*π*x)* sin(3*π*y)* sin(3*π*z)
    end
    interp_simple!(cent_correct, face, grid)
    interp!(cent, face, grid)
    @test all(cent .≈ cent_correct)

    # Center -> Face
    cent = CellCenter(grid)
    face_correct = CellFace(grid, 1)
    face = CellFace(grid, 1)
    @inbounds for local_grid in GridField(grid,cent)
        @unpack ijk,x,y,z = local_grid
        cent[ijk] = sin(3*π*x)* sin(3*π*y)* sin(3*π*z)
    end
    interp_simple!(face_correct, cent, grid)
    interp!(face, cent, grid)
    @test all(face .≈ face_correct)
end

#####
##### CellCorner <-> CellEdge
#####

function interp_simple!(corn::CellCorner, edge::CellEdge{1}, grid::Grid)
    corn .= 0 # boundaries are not populated.
    for k in over_points(grid.c[3], All(), Vertex1D())
        for j in over_points(grid.c[2], All(), Vertex1D())
            for i in over_points(grid.c[1], Interior(), Vertex1D())
                D = grid.c[1].interp.diag[i]
                U = grid.c[1].interp.upper[i]
                corn[i, j, k] = D*edge[i, j, k] + U*edge[i+1, j, k]
            end
        end
    end
end

function interp_simple!(edge::CellEdge{1}, corn::CellCorner, grid::Grid)
    for k in over_points(grid.c[3], All(), Vertex1D())
        for j in over_points(grid.c[2], All(), Vertex1D())
            for i in over_points(grid.c[1], All(), Center1D())
                D = 1//2
                U = 1//2
                edge[i, j, k] = D*corn[i, j, k] + U*corn[i+1, j, k]
            end
        end
    end
end

@testset "Interpolations - Simple implementation Corner <-> Edge" begin
    FT = Float64
    n, a, b, β = 10, FT(0.0), FT(1.0), FT(1.01)
    grid = Cube(a, b, n; warpfun=Roberts_both, args=(β,))

    # Edge -> Corner
    edge = CellEdge(grid, 1)
    corn_correct = CellCorner(grid)
    corn = CellCorner(grid)
    @inbounds for local_grid in GridField(grid,edge)
        @unpack ijk,x,y,z = local_grid
        edge[ijk] = sin(3*π*x)* sin(3*π*y)* sin(3*π*z)
    end
    interp_simple!(corn_correct, edge, grid)
    interp!(corn, edge, grid)
    @test all(corn .≈ corn_correct)

    # corner -> edge
    corn = CellCorner(grid)
    edge_correct = CellEdge(grid, 1)
    edge = CellEdge(grid, 1)
    @inbounds for local_grid in GridField(grid,corn)
        @unpack ijk,x,y,z = local_grid
        corn[ijk] = sin(3*π*x)* sin(3*π*y)* sin(3*π*z)
    end
    interp_simple!(edge_correct, corn, grid)
    interp!(edge, corn, grid)
    @test all(edge .≈ edge_correct)
end

#####
##### CellFace <-> CellEdge
#####

function interp_simple!(face::CellFace{1}, edge::CellEdge{2}, grid::Grid)
    for k in over_points(grid.c[3], All(), Center1D())
        for j in over_points(grid.c[2], All(), Center1D())
            for i in over_points(grid.c[1], All(), Vertex1D())
                D = 1//2
                U = 1//2
                face[i, j, k] = D*edge[i, j, k] + U*edge[i, j, k+1]
            end
        end
    end
end

function interp_simple!(edge::CellEdge{2}, face::CellFace{1}, grid::Grid)
    edge .= 0 # boundaries are not populated.
    for k in over_points(grid.c[3], Interior(), Vertex1D())
        for j in over_points(grid.c[2], All(), Center1D())
            for i in over_points(grid.c[1], All(), Vertex1D())
                D = grid.c[3].interp.diag[k]
                U = grid.c[3].interp.upper[k]
                edge[i, j, k] = D*face[i, j, k] + U*face[i, j, k+1]
            end
        end
    end
end

@testset "Interpolations - Simple implementation Face <-> Edge" begin
    FT = Float64
    n, a, b, β = 10, FT(0.0), FT(1.0), FT(1.01)
    grid = Cube(a, b, n; warpfun=Roberts_both, args=(β,))

    # Edge -> Face
    edge = CellEdge(grid, 2)
    face_correct = CellFace(grid, 1)
    face = CellFace(grid, 1)
    @inbounds for local_grid in GridField(grid,edge)
        @unpack ijk,x,y,z = local_grid
        edge[ijk] = sin(3*π*x)* sin(3*π*y)* sin(3*π*z)
    end
    interp_simple!(face_correct, edge, grid)
    interp!(face, edge, grid)
    @test all(face .≈ face_correct)

    # Face -> edge
    face = CellFace(grid, 1)
    edge_correct = CellEdge(grid, 2)
    edge = CellEdge(grid, 2)
    @inbounds for local_grid in GridField(grid,face)
        @unpack ijk,x,y,z = local_grid
        face[ijk] = sin(3*π*x)* sin(3*π*y)* sin(3*π*z)
    end
    interp_simple!(edge_correct, face, grid)
    interp!(edge, face, grid)
    @test all(edge .≈ edge_correct)
end
