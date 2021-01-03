using Test
using UnPack
using MOONS.Grids
using MOONS.Fields

#####
##### Self-extrapolation
#####

extrap_simple!(f::AbstractField, grid::Grid, dim::Int) =
    extrap_simple!(f, grid, Val(dim))

function extrap_simple!(f::AbstractField, grid::Grid, ::Val{1})
    f[ 1 ,:,:] .= 2*f[  2  ,:,:] .- f[  3  ,:,:]
    f[end,:,:] .= 2*f[end-1,:,:] .- f[end-2,:,:]
end
function extrap_simple!(f::AbstractField, grid::Grid, ::Val{2})
    f[:, 1 ,:] .= 2*f[:,  2  ,:] .- f[:,  3  ,:]
    f[:,end,:] .= 2*f[:,end-1,:] .- f[:,end-2,:]
end
function extrap_simple!(f::AbstractField, grid::Grid, ::Val{3})
    f[:,:, 1 ] .= 2*f[:,:,  2  ] .- f[:,:,  3  ]
    f[:,:,end] .= 2*f[:,:,end-1] .- f[:,:,end-2]
end

@testset "Extrapolations - self extrap" begin
    FT = Float64
    n, a, b, β = 10, FT(0.0), FT(1.0), FT(1.01)
    grid = Cube(a, b, n; warpfun=Roberts_both, args=(β,))

    cent = CellCenter(grid)
    corn = CellCorner(grid)
    face = CellFaces(grid)
    edge = CellEdges(grid)
    cent_simple = CellCenter(grid)
    corn_simple = CellCorner(grid)
    face_simple = CellFaces(grid)
    edge_simple = CellEdges(grid)
    @inbounds for local_grid in GridField(grid, cent)
        @unpack ijk,x,y,z = local_grid
        cent[ijk] = cos(3*π*x)* cos(3*π*y)* cos(3*π*z)
        cent_simple[ijk] = cos(3*π*x)* cos(3*π*y)* cos(3*π*z)
    end
    @inbounds for local_grid in GridField(grid, corn)
        @unpack ijk,x,y,z = local_grid
        corn[ijk] = cos(3*π*x)* cos(3*π*y)* cos(3*π*z)
        corn_simple[ijk] = cos(3*π*x)* cos(3*π*y)* cos(3*π*z)
    end
    for dir in 1:3
        @inbounds for local_grid in GridField(grid, face[dir])
            @unpack ijk,x,y,z = local_grid
            face[dir][ijk] = cos(3*π*x)* cos(3*π*y)* cos(3*π*z)
            face_simple[dir][ijk] = cos(3*π*x)* cos(3*π*y)* cos(3*π*z)
        end
        @inbounds for local_grid in GridField(grid, edge[dir])
            @unpack ijk,x,y,z = local_grid
            edge[dir][ijk] = cos(3*π*x)* cos(3*π*y)* cos(3*π*z)
            edge_simple[dir][ijk] = cos(3*π*x)* cos(3*π*y)* cos(3*π*z)
        end
    end
    for dir in 1:3
        extrap!(face[1], grid, dir)
        extrap!(face[2], grid, dir)
        extrap!(face[3], grid, dir)
        extrap!(edge[1], grid, dir)
        extrap!(edge[2], grid, dir)
        extrap!(edge[3], grid, dir)
        extrap!(cent, grid, dir)
        extrap!(corn, grid, dir)
        extrap_simple!(face_simple[1], grid, dir)
        extrap_simple!(face_simple[2], grid, dir)
        extrap_simple!(face_simple[3], grid, dir)
        extrap_simple!(edge_simple[1], grid, dir)
        extrap_simple!(edge_simple[2], grid, dir)
        extrap_simple!(edge_simple[3], grid, dir)
        extrap_simple!(cent_simple, grid, dir)
        extrap_simple!(corn_simple, grid, dir)
        @test all(face_simple[1] .≈ face[1])
        @test all(face_simple[2] .≈ face[2])
        @test all(face_simple[3] .≈ face[3])
        @test all(edge_simple[1] .≈ edge[1])
        @test all(edge_simple[2] .≈ edge[2])
        @test all(edge_simple[3] .≈ edge[3])
        @test all(cent_simple .≈ cent)
        @test all(corn_simple .≈ corn)
    end
end

#####
##### Staggered 2-field extrapolation
#####

function extrap_x_simple!(f::AbstractField, g::AbstractField,
    grid::Grid, f_dl::Center1D, g_dl::Vertex1D)
    f[ 1 ,:,:] .= 2*g[  2  ,:,:] .- f[  2  ,:,:]
    f[end,:,:] .= 2*g[end-1,:,:] .- f[end-1,:,:]
end
function extrap_x_simple!(f::AbstractField, g::AbstractField,
    grid::Grid, f_dl::Vertex1D, g_dl::Center1D)
    f[ 1 ,:,:] .= 2*g[ 1 ,:,:] .- f[  2  ,:,:]
    f[end,:,:] .= 2*g[end,:,:] .- f[end-1,:,:]
end

@testset "Extrapolations - Staggered 2-field" begin
    FT = Float64
    n, a, b, β = 10, FT(0.0), FT(1.0), FT(1.01)
    grid = Cube(a, b, n; warpfun=Roberts_both, args=(β,))

    cent = CellCenter(grid)
    corn = CellCorner(grid)
    face = CellFace(grid, 1)
    edge = CellEdge(grid, 1)
    cent_simple = CellCenter(grid)
    corn_simple = CellCorner(grid)
    face_simple = CellFace(grid, 1)
    edge_simple = CellEdge(grid, 1)
    cent_out = CellCenter(grid)
    corn_out = CellCorner(grid)
    face_out = CellFace(grid, 1)
    edge_out = CellEdge(grid, 1)
    @inbounds for local_grid in GridField(grid, cent)
        @unpack ijk,x,y,z = local_grid
        cent[ijk] = cos(3*π*x)* cos(3*π*y)* cos(3*π*z)
        cent_simple[ijk] = cos(3*π*x)* cos(3*π*y)* cos(3*π*z)
        cent_out[ijk] = cos(3*π*x)* cos(3*π*y)* cos(3*π*z)
    end
    @inbounds for local_grid in GridField(grid, corn)
        @unpack ijk,x,y,z = local_grid
        corn[ijk] = cos(3*π*x)* cos(3*π*y)* cos(3*π*z)
        corn_simple[ijk] = cos(3*π*x)* cos(3*π*y)* cos(3*π*z)
        corn_out[ijk] = cos(3*π*x)* cos(3*π*y)* cos(3*π*z)
    end
    @inbounds for local_grid in GridField(grid, face)
        @unpack ijk,x,y,z = local_grid
        face[ijk] = cos(3*π*x)* cos(3*π*y)* cos(3*π*z)
        face_simple[ijk] = cos(3*π*x)* cos(3*π*y)* cos(3*π*z)
        face_out[ijk] = cos(3*π*x)* cos(3*π*y)* cos(3*π*z)
    end
    @inbounds for local_grid in GridField(grid, edge)
        @unpack ijk,x,y,z = local_grid
        edge[ijk] = cos(3*π*x)* cos(3*π*y)* cos(3*π*z)
        edge_simple[ijk] = cos(3*π*x)* cos(3*π*y)* cos(3*π*z)
        edge_out[ijk] = cos(3*π*x)* cos(3*π*y)* cos(3*π*z)
    end

    extrap_x_simple!(cent_simple, face, grid, Center1D(), Vertex1D())
    extrap_x_simple!(face_simple, cent, grid, Vertex1D(), Center1D())
    extrap_x_simple!(corn_simple, edge, grid, Vertex1D(), Center1D())
    extrap_x_simple!(edge_simple, corn, grid, Center1D(), Vertex1D())

    extrap!(cent_out, face, grid, 1)
    extrap!(face_out, cent, grid, 1)
    extrap!(corn_out, edge, grid, 1)
    extrap!(edge_out, corn, grid, 1)

    @test all(cent_out .≈ cent_simple)
    @test all(face_out .≈ face_simple)
    @test all(corn_out .≈ corn_simple)
    @test all(edge_out .≈ edge_simple)
end


@testset "Extrapolations - test throws" begin
    FT = Float64
    n, a, b = 10, FT(0.0), FT(1.0)
    grid = Cube(a, b, n)
    face = CellFaces(grid)
    edge = CellEdges(grid)
    for dir in 1:3
        @test_throws MethodError extrap!(face[dir], edge[dir], grid, dir)
        @test_throws MethodError extrap!(edge[dir], face[dir], grid, dir)
    end
end

@testset "Extrapolations - Staggered 2-field convergence" begin
    FT = Float64
    a, b = FT(0.0), FT(1.0)
    nr = 2 .^ (3,4,5,6)
    Δh = zeros(length(nr))
    err_c2v = zeros(length(nr))
    err_v2c = zeros(length(nr))
    for (k,n) in enumerate(2 .^ (3,4,5,6))
        grid = Cube(a, b, n)
        Δh[k] = grid.c[1].n.Δh[1]
        cent = CellCenter(grid)
        corn = CellCorner(grid)
        face = CellFace(grid, 1)
        edge = CellEdge(grid, 1)
        cent_exact = CellCenter(grid)
        corn_exact = CellCorner(grid)
        face_exact = CellFace(grid, 1)
        edge_exact = CellEdge(grid, 1)
        @inbounds for local_grid in GridField(grid,face)
            @unpack ijk,x,y,z = local_grid
            face[ijk] = sin(3*π*x)* sin(3*π*y)* sin(3*π*z)
            face_exact[ijk] = sin(3*π*x)* sin(3*π*y)* sin(3*π*z)
        end
        @inbounds for local_grid in GridField(grid,cent)
            @unpack ijk,x,y,z = local_grid
            cent[ijk] = sin(3*π*x)* sin(3*π*y)* sin(3*π*z)
            cent_exact[ijk] = sin(3*π*x)* sin(3*π*y)* sin(3*π*z)
        end
        extrap!(cent, face, grid, 1)
        extrap!(face, cent, grid, 1)
        err_c2v[k] = norm(cent .- cent_exact) / length(cent)
        err_v2c[k] = norm(face .- face_exact) / length(face)
    end
    conv_c2v = convergence_rate(err_c2v, Δh)
    conv_v2c = convergence_rate(err_v2c, Δh)
    @test conv_c2v[1] ≈ 4.093446045776282
    @test conv_c2v[2] ≈ 2.6665448437905694
    @test conv_c2v[3] ≈ 3.1919109381432618

    @test conv_v2c[1] ≈ 4.5665684649141545
    @test conv_v2c[2] ≈ 4.719101505009674
    @test conv_v2c[3] ≈ 4.848188053541373
end
