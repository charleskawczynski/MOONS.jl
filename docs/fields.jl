using Plots, UnPack
using MOONS.Grids
using MOONS.Fields
FT = Float64
a, b, n = FT(0.0), FT(1.0), 50
β = FT(1.01)
g = Cube(a, b, n; warpfun=Roberts_both, args=(β,))

cent = CellCenter(g)

@inbounds for local_grid in GridField(g,cent)
    @unpack ijk,x,y,z = local_grid
    cent[ijk] = sin(3*π*x)*sin(3*π*y)*sin(3*π*z)
end

@unpack xn, yn, zn, xc, yc, zc = all_coords(g)

plot(yc,zc, midslice(cent, 1), seriestype=:wireframe)
savefig("wireframe.svg")
contourf(yc,zc, midslice(cent, 1), c = :viridis)
savefig("contourf.svg")

# Extrap

FT = Float64
a, b, n = FT(0.0), FT(1.0), 20
g = Cube(a, b, n)
@unpack xn, yn, zn, xc, yc, zc = all_coords(g)

cent = CellCenter(g)
face = CellFace(g, 1)
corn = CellCorner(g)
edge = CellEdge(g, 1)
@inbounds for local_grid in GridField(g,cent)
    cent[local_grid.ijk] = cos(3*π*local_grid.x)
end
@inbounds for local_grid in GridField(g,face)
    face[local_grid.ijk] = cos(3*π*local_grid.x)
end
@inbounds for local_grid in GridField(g,corn)
    corn[local_grid.ijk] = cos(3*π*local_grid.x)
end
@inbounds for local_grid in GridField(g,edge)
    edge[local_grid.ijk] = cos(3*π*local_grid.x)
end
extrap!(face, cent, g, 1)
extrap!(edge, corn, g, 1)
plot(xn[1:4], midline(face, 1)[1:4], label="face (extrapolated)")
plot!(xc[1:4], midline(cent, 1)[1:4], label="cent")
savefig("extrap_c2v.svg")
plot(xc[1:4], midline(edge, 1)[1:4], label="edge (extrapolated)")
plot!(xn[1:4], midline(corn, 1)[1:4], label="corn")
savefig("extrap_v2c.svg")

FT = Float64
ax, bx = FT(0.0), FT(1.0)
ay, by = FT(0.0), FT(2.0)
n = 2 ^ 5
g = Grid(
    (Coordinates(ax,bx,n),
     Coordinates(ay,by,n),
     Coordinates(-FT(0.5),FT(0.5), 1))
    )

f = CellCorner(g)
∇f = CellEdges(g)
∇f_exact = CellEdges(g)

@inbounds for local_grid in GridField(g,f)
    @unpack ijk,x,y = local_grid
    f[ijk] = sin(3*π*x)*sin(3*π*y)
end
@inbounds for local_grid in GridField(g,∇f_exact[1])
    @unpack ijk,x,y = local_grid
    ∇f_exact[1][ijk] = 3*π*cos(3*π*x)*sin(3*π*y)
end
@inbounds for local_grid in GridField(g,∇f_exact[2])
    @unpack ijk,x,y = local_grid
    ∇f_exact[2][ijk] = 3*π*sin(3*π*x)*cos(3*π*y)
end
∇!(∇f, f, g)

@unpack xc, yc, zc, xn, yn, zn = all_coords(g)

contourf(xn,yn, midslice(f, 3)', c = :viridis, xlabel="x", ylabel="y")
savefig("f.svg")

contourf(xc,yn, midslice(∇f[1], 3)', c = :viridis, xlabel="x", ylabel="y")
savefig("∇f_1.svg")
contourf(xn,yc, midslice(∇f[2], 3)', c = :viridis, xlabel="x", ylabel="y")
savefig("∇f_2.svg")

contourf(xc,yn, midslice(∇f_exact[1], 3)', c = :viridis, xlabel="x", ylabel="y")
savefig("∇f_1_exact.svg")
contourf(xn,yc, midslice(∇f_exact[2], 3)', c = :viridis, xlabel="x", ylabel="y")
savefig("∇f_2_exact.svg")
