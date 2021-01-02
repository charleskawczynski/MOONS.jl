using Test
using LinearAlgebra
using UnPack
using MOONS.Grids
using MOONS.Fields

"""
    convergence_rate(err, Δh)

Estimate convergence rate given
vectors `err` and `Δh`
```
err = C Δh^p+ H.O.T
err_k ≈ C Δh_k^p
err_k/err_m ≈ Δh_k^p/Δh_m^p
log(err_k/err_m) ≈ log((Δh_k/Δh_m)^p)
log(err_k/err_m) ≈ p*log(Δh_k/Δh_m)
log(err_k/err_m)/log(Δh_k/Δh_m) ≈ p
```
"""
convergence_rate(err, Δh) =
    [log(err[i]/err[i-1])/log(Δh[i]/Δh[i-1]) for i in 2:length(Δh)]

@testset "Face -> Center" begin
    FT = Float64
    a, b = FT(0.0), FT(1.0)
    nr = 2 .^ (3,4,5,6)
    err, Δh = zeros(length(nr)), zeros(length(nr))
    for (k,n) in enumerate(2 .^ (3,4,5,6))
        grid = Cube(a, b, n; warpfun=Roberts_both, args=(FT(1.01),))
        Δh[k] = grid.c[1].n.Δh[1]
        face = CellFace(grid, 1)
        cent_exact = CellCenter(grid)
        cent = CellCenter(grid)
        @inbounds for local_grid in GridField(grid,face)
            @unpack ijk,x,y,z = local_grid
            face[ijk] = sin(3*π*x)* sin(3*π*y)* sin(3*π*z)
        end
        @inbounds for local_grid in GridField(grid,cent_exact)
            @unpack ijk,x,y,z = local_grid
            cent_exact[ijk] = sin(3*π*x)* sin(3*π*y)* sin(3*π*z)
        end
        interp!(cent, face, grid)
        err[k] = norm(cent .- cent_exact) / length(cent)
    end
    conv = convergence_rate(err, Δh)
    @test conv[1] ≈ 2.4741587422769133
    @test conv[2] ≈ 2.691392600417213
    @test conv[3] ≈ 3.0586931843562843
end

@testset "Center -> Face" begin
    FT = Float64
    a, b = FT(0.0), FT(1.0)
    nr = 2 .^ (3,4,5,6)
    err, Δh = zeros(length(nr)), zeros(length(nr))
    for (k,n) in enumerate(2 .^ (3,4,5,6))
        grid = Cube(a, b, n; warpfun=Roberts_both, args=(FT(1.01),))
        Δh[k] = grid.c[1].n.Δh[1]
        cent = CellCenter(grid)
        face_exact = CellFace(grid, 1)
        face = CellFace(grid, 1)
        @inbounds for local_grid in GridField(grid,cent)
            @unpack ijk,x,y,z = local_grid
            cent[ijk] = sin(3*π*x)* sin(3*π*y)* sin(3*π*z)
        end
        @inbounds for local_grid in GridField(grid,face_exact)
            @unpack ijk,x,y,z = local_grid
            face_exact[ijk] = sin(3*π*x)* sin(3*π*y)* sin(3*π*z)
        end
        face_exact[1,:,:] .= 0 # interpolation does not populate this
        face_exact[end,:,:] .= 0 # interpolation does not populate this
        interp!(face, cent, grid)
        err[k] = norm(face .- face_exact) / length(face)
    end
    conv = convergence_rate(err, Δh)
    @test conv[1] ≈ 1.6467834910768682
    @test conv[2] ≈ 2.63750140470393
    @test conv[3] ≈ 3.031320955254459
end

@testset "Corner -> Edge" begin
    FT = Float64
    a, b = FT(0.0), FT(1.0)
    nr = 2 .^ (3,4,5,6)
    err, Δh = zeros(length(nr)), zeros(length(nr))
    for (k,n) in enumerate(2 .^ (3,4,5,6))
        grid = Cube(a, b, n; warpfun=Roberts_both, args=(FT(1.01),))
        Δh[k] = grid.c[1].n.Δh[1]
        corn = CellCorner(grid)
        edge_exact = CellEdge(grid, 1)
        edge = CellEdge(grid, 1)
        @inbounds for local_grid in GridField(grid,corn)
            @unpack ijk,x,y,z = local_grid
            corn[ijk] = sin(3*π*x)* sin(3*π*y)* sin(3*π*z)
        end
        @inbounds for local_grid in GridField(grid,edge_exact)
            @unpack ijk,x,y,z = local_grid
            edge_exact[ijk] = sin(3*π*x)* sin(3*π*y)* sin(3*π*z)
        end
        interp!(edge, corn, grid)
        err[k] = norm(edge .- edge_exact) / length(edge)
    end
    conv = convergence_rate(err, Δh)
    @test conv[1] ≈ 2.026575204117072
    @test conv[2] ≈ 2.6180379005577397
    @test conv[3] ≈ 3.0183074193213058
end

@testset "Edge -> Corner" begin
    FT = Float64
    a, b = FT(0.0), FT(1.0)
    nr = 2 .^ (3,4,5,6)
    err, Δh = zeros(length(nr)), zeros(length(nr))
    for (k,n) in enumerate(2 .^ (3,4,5,6))
        grid = Cube(a, b, n; warpfun=Roberts_both, args=(FT(1.01),))
        Δh[k] = grid.c[1].n.Δh[1]
        edge = CellEdge(grid, 1)
        corn_exact = CellCorner(grid)
        corn = CellCorner(grid)
        @inbounds for local_grid in GridField(grid,edge)
            @unpack ijk,x,y,z = local_grid
            edge[ijk] = sin(3*π*x)* sin(3*π*y)* sin(3*π*z)
        end
        @inbounds for local_grid in GridField(grid,corn_exact)
            @unpack ijk,x,y,z = local_grid
            corn_exact[ijk] = sin(3*π*x)* sin(3*π*y)* sin(3*π*z)
        end
        corn_exact[1,:,:] .= 0 # interpolation does not populate this
        corn_exact[end,:,:] .= 0 # interpolation does not populate this
        interp!(corn, edge, grid)
        err[k] = norm(corn .- corn_exact) / length(corn)
    end
    conv = convergence_rate(err, Δh)
    @test conv[1] ≈ 1.1991999529170272
    @test conv[2] ≈ 2.5641467048444557
    @test conv[3] ≈ 2.9909351902194867
end

@testset "Edge -> Face" begin
    FT = Float64
    a, b = FT(0.0), FT(1.0)
    nr = 2 .^ (3,4,5,6)
    err, Δh = zeros(length(nr)), zeros(length(nr))
    for (k,n) in enumerate(2 .^ (3,4,5,6))
        grid = Cube(a, b, n; warpfun=Roberts_both, args=(FT(1.01),))
        Δh[k] = grid.c[1].n.Δh[1]
        edge = CellEdge(grid, 2)
        face_exact = CellFace(grid, 1)
        face = CellFace(grid, 1)
        @inbounds for local_grid in GridField(grid,edge)
            @unpack ijk,x,y,z = local_grid
            edge[ijk] = sin(3*π*x)* sin(3*π*y)* sin(3*π*z)
        end
        @inbounds for local_grid in GridField(grid,face_exact)
            @unpack ijk,x,y,z = local_grid
            face_exact[ijk] = sin(3*π*x)* sin(3*π*y)* sin(3*π*z)
        end
        interp!(face, edge, grid)
        err[k] = norm(face .- face_exact) / length(face)
    end
    conv = convergence_rate(err, Δh)
    @test conv[1] ≈ 2.2503669731969924
    @test conv[2] ≈ 2.6547152504874756
    @test conv[3] ≈ 3.038500301838792
end

@testset "Face -> Edge" begin
    FT = Float64
    a, b = FT(0.0), FT(1.0)
    nr = 2 .^ (3,4,5,6)
    err, Δh = zeros(length(nr)), zeros(length(nr))
    for (k,n) in enumerate(2 .^ (3,4,5,6))
        grid = Cube(a, b, n; warpfun=Roberts_both, args=(FT(1.01),))
        Δh[k] = grid.c[1].n.Δh[1]
        face = CellFace(grid, 1)
        edge_exact = CellEdge(grid, 2)
        edge = CellEdge(grid, 2)
        @inbounds for local_grid in GridField(grid,face)
            @unpack ijk,x,y,z = local_grid
            face[ijk] = sin(3*π*x)* sin(3*π*y)* sin(3*π*z)
        end
        @inbounds for local_grid in GridField(grid,edge_exact)
            @unpack ijk,x,y,z = local_grid
            edge_exact[ijk] = sin(3*π*x)* sin(3*π*y)* sin(3*π*z)
        end
        interp!(edge, face, grid)
        edge_exact[:,:,1] .= 0 # interpolation does not populate this
        edge_exact[:,:,end] .= 0 # interpolation does not populate this
        err[k] = norm(edge .- edge_exact) / length(edge)
    end
    conv = convergence_rate(err, Δh)
    @test conv[1] ≈ 1.4229917219969475
    @test conv[2] ≈ 2.600824054774193
    @test conv[3] ≈ 3.011128072736973
end

