using Test
using LinearAlgebra
using UnPack
using MOONS.Grids
using MOONS.Fields

convergence_rate(err, Δh) =
    [log(err[i]/err[i-1])/log(Δh[i]/Δh[i-1]) for i in 2:length(Δh)]

@testset "Discrete operators - ∇(f_corn)=(∇f)_edge" begin
    FT = Float64
    a, b = FT(0.0), FT(1.0)
    nr = 2 .^ (3,4,5,6)
    Δh = zeros(length(nr))
    errx = zeros(length(nr))
    erry = zeros(length(nr))
    errz = zeros(length(nr))
    for (k,n) in enumerate(2 .^ (3,4,5,6))
        grid = Cube(a, b, n; warpfun=Roberts_both, args=(FT(1.01),))
        Δh[k] = grid.c[1].n.Δh[1]

        f = CellCorner(grid)
        ∇f = CellEdges(grid)
        ∇f_exact = CellEdges(grid)

        @inbounds for local_grid in GridField(grid,f)
            @unpack ijk,x,y,z = local_grid
            f[ijk] = sin(3*π*x)*sin(3*π*y)*sin(3*π*z)
        end
        @inbounds for local_grid in GridField(grid,∇f_exact[1])
            @unpack ijk,x,y,z = local_grid
            ∇f_exact[1][ijk] = 3*π*cos(3*π*x)*sin(3*π*y)*sin(3*π*z)
        end
        @inbounds for local_grid in GridField(grid,∇f_exact[2])
            @unpack ijk,x,y,z = local_grid
            ∇f_exact[2][ijk] = 3*π*sin(3*π*x)*cos(3*π*y)*sin(3*π*z)
        end
        @inbounds for local_grid in GridField(grid,∇f_exact[3])
            @unpack ijk,x,y,z = local_grid
            ∇f_exact[3][ijk] = 3*π*sin(3*π*x)*sin(3*π*y)*cos(3*π*z)
        end
        ∇!(∇f, f, grid)
        errx[k] = norm(∇f_exact[1] .- ∇f[1]) / length(∇f[1])
        erry[k] = norm(∇f_exact[2] .- ∇f[2]) / length(∇f[2])
        errz[k] = norm(∇f_exact[3] .- ∇f[3]) / length(∇f[3])
    end

    @test errx[1] ≈ 0.010636252411867275
    @test errx[2] ≈ 0.0008811926803278781
    @test errx[3] ≈ 9.917217369054536e-5
    @test errx[4] ≈ 9.925963845483006e-6
    @test erry[1] ≈ 0.010636252411867273
    @test erry[2] ≈ 0.0008811926803278777
    @test erry[3] ≈ 9.917217369054523e-5
    @test erry[4] ≈ 9.925963845482962e-6
    @test errz[1] ≈ 0.010636252411867275
    @test errz[2] ≈ 0.0008811926803278781
    @test errz[3] ≈ 9.917217369054519e-5
    @test errz[4] ≈ 9.925963845482955e-6

    convx = convergence_rate(errx, Δh)
    convy = convergence_rate(erry, Δh)
    convz = convergence_rate(errz, Δh)

    @test convx[1] ≈ 2.328658024470567
    @test convx[2] ≈ 2.5112742673767268
    @test convx[3] ≈ 2.9570095611906884
    @test convy[1] ≈ 2.328658024470567
    @test convy[2] ≈ 2.511274267376728
    @test convy[3] ≈ 2.9570095611906924
    @test convz[1] ≈ 2.328658024470567
    @test convz[2] ≈ 2.511274267376729
    @test convz[3] ≈ 2.9570095611906932
end

@testset "Discrete operators - ∇(f_cent)=(∇f)_face" begin
    FT = Float64
    a, b = FT(0.0), FT(1.0)
    nr = 2 .^ (3,4,5,6)
    Δh = zeros(length(nr))
    errx = zeros(length(nr))
    erry = zeros(length(nr))
    errz = zeros(length(nr))
    for (k,n) in enumerate(2 .^ (3,4,5,6))
        grid = Cube(a, b, n; warpfun=Roberts_both, args=(FT(1.01),))
        Δh[k] = grid.c[1].n.Δh[1]

        f = CellCenter(grid)
        ∇f = CellFaces(grid)
        ∇f_exact = CellFaces(grid)

        @inbounds for local_grid in GridField(grid,f)
            @unpack ijk,x,y,z = local_grid
            f[ijk] = sin(3*π*x)*sin(3*π*y)*sin(3*π*z)
        end
        @inbounds for local_grid in GridField(grid,∇f_exact[1])
            @unpack ijk,x,y,z = local_grid
            ∇f_exact[1][ijk] = 3*π*cos(3*π*x)*sin(3*π*y)*sin(3*π*z)
        end
        @inbounds for local_grid in GridField(grid,∇f_exact[2])
            @unpack ijk,x,y,z = local_grid
            ∇f_exact[2][ijk] = 3*π*sin(3*π*x)*cos(3*π*y)*sin(3*π*z)
        end
        @inbounds for local_grid in GridField(grid,∇f_exact[3])
            @unpack ijk,x,y,z = local_grid
            ∇f_exact[3][ijk] = 3*π*sin(3*π*x)*sin(3*π*y)*cos(3*π*z)
        end
        ∇!(∇f, f, grid)

        # Ghost points are not populated:
        ∇f[1][1,:,:] .= 0;   ∇f_exact[1][1,:,:] .= 0
        ∇f[1][end,:,:] .= 0; ∇f_exact[1][end,:,:] .= 0
        ∇f[2][:,1,:] .= 0;   ∇f_exact[2][:,1,:] .= 0
        ∇f[2][:,end,:] .= 0; ∇f_exact[2][:,end,:] .= 0
        ∇f[3][:,:,1] .= 0;   ∇f_exact[3][:,:,1] .= 0
        ∇f[3][:,:,end] .= 0; ∇f_exact[3][:,:,end] .= 0

        errx[k] = norm(∇f_exact[1] .- ∇f[1]) / length(∇f[1])
        erry[k] = norm(∇f_exact[2] .- ∇f[2]) / length(∇f[2])
        errz[k] = norm(∇f_exact[3] .- ∇f[3]) / length(∇f[3])
    end
    @test errx[1] ≈ 0.006504291238410591
    @test errx[2] ≈ 0.0015861146672118875
    @test errx[3] ≈ 0.0001769587725572569
    @test errx[4] ≈ 1.7467942413498958e-5
    @test erry[1] ≈ 0.006504291238410591
    @test erry[2] ≈ 0.0015861146672118875
    @test erry[3] ≈ 0.00017695877255725676
    @test erry[4] ≈ 1.7467942413498904e-5
    @test errz[1] ≈ 0.006504291238410591
    @test errz[2] ≈ 0.0015861146672118875
    @test errz[3] ≈ 0.0001769587725572567
    @test errz[4] ≈ 1.7467942413498917e-5

    convx = convergence_rate(errx, Δh)
    convy = convergence_rate(erry, Δh)
    convz = convergence_rate(errz, Δh)

    @test convx[1] ≈ 1.3193403930281424
    @test convx[2] ≈ 2.5212841871228213
    @test convx[3] ≈ 2.974797468220822
    @test convy[1] ≈ 1.3193403930281424
    @test convy[2] ≈ 2.521284187122822
    @test convy[3] ≈ 2.9747974682208245
    @test convz[1] ≈ 1.3193403930281424
    @test convz[2] ≈ 2.521284187122822
    @test convz[3] ≈ 2.9747974682208236
end

@testset "Discrete operators - ∇(f_face)=(∇f)_edge" begin
    FT = Float64
    a, b = FT(0.0), FT(1.0)
    nr = 2 .^ (3,4,5,6)
    Δh = zeros(length(nr))
    err = zeros(length(nr))
    for (k,n) in enumerate(2 .^ (3,4,5,6))
        grid = Cube(a, b, n; warpfun=Roberts_both, args=(FT(1.01),))
        Δh[k] = grid.c[1].n.Δh[1]

        f = CellFace(grid, 1)
        ∇f = CellEdge(grid, 2)
        ∇f_exact = CellEdge(grid, 2)

        @inbounds for local_grid in GridField(grid,f)
            @unpack ijk,x,y,z = local_grid
            f[ijk] = sin(3*π*x)*sin(3*π*y)*sin(3*π*z)
        end
        @inbounds for local_grid in GridField(grid,∇f_exact)
            @unpack ijk,x,y,z = local_grid
            ∇f_exact[ijk] = 3*π*sin(3*π*x)*sin(3*π*y)*cos(3*π*z)
        end
        ∇!(∇f, f, grid, 3)

        # Ghost points are not populated:
        ∇f[:,:,1] .= 0;   ∇f_exact[:,:,1] .= 0
        ∇f[:,:,end] .= 0; ∇f_exact[:,:,end] .= 0

        err[k] = norm(∇f_exact .- ∇f) / length(∇f)
    end
    @test err[1] ≈ 0.007431514201977017
    @test err[2] ≈ 0.0014919756305751361
    @test err[3] ≈ 0.00017156454627127846
    @test err[4] ≈ 1.7198655291695e-5

    conv = convergence_rate(err, Δh)

    @test conv[1] ≈ 1.5011397571937848
    @test conv[2] ≈ 2.486531977140251
    @test conv[3] ≈ 2.9549859273624386
end

function ∇²_simple!(∇²f::T, f::T, grid::Grid) where {T<:AbstractField}
    ∇f = dual_fields(grid, f)
    ∇²f_x = deepcopy(f)
    ∇²f_y = deepcopy(f)
    ∇²f_z = deepcopy(f)
    ∇!(∇f[1], f, grid)
    ∇!(∇f[2], f, grid)
    ∇!(∇f[3], f, grid)
    ∇!(∇²f_x, ∇f[1], grid)
    ∇!(∇²f_y, ∇f[2], grid)
    ∇!(∇²f_z, ∇f[3], grid)
    ∇²f .= ∇²f_x .+ ∇²f_y .+ ∇²f_z
    assign!(∇²f, 0, Ghost())
end

@testset "Discrete operators - ∇²f test against simple" begin
    FT = Float64
    a, b = FT(0.0), FT(1.0)
    n = 10
    grid = Cube(a, b, n; warpfun=Roberts_both, args=(FT(1.01),))
    f = CellCenter(grid)
    ∇²f = CellCenter(grid)
    ∇²f_simple = CellCenter(grid)
    @inbounds for local_grid in GridField(grid,f)
        @unpack ijk,x,y,z = local_grid
        f[ijk] = sin(3*π*x)*sin(3*π*y)*sin(3*π*z)
    end
    ∇²_simple!(∇²f_simple, f, grid)
    ∇²!(∇²f, f, grid)
    @test all(∇²f .≈ ∇²f_simple)

    f = CellCorner(grid)
    ∇²f = CellCorner(grid)
    ∇²f_simple = CellCorner(grid)
    @inbounds for local_grid in GridField(grid,f)
        @unpack ijk,x,y,z = local_grid
        f[ijk] = sin(3*π*x)*sin(3*π*y)*sin(3*π*z)
    end
    ∇²_simple!(∇²f_simple, f, grid)
    ∇²!(∇²f, f, grid)

    # We get hit by precision error a bit with ∇², so need to relax approx:
    @test all(isapprox.(∇²f, ∇²f_simple; atol=10^4*eps()))
end
