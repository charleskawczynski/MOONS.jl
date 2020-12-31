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
    cent[ijk] = sin(3*π*x)* sin(3*π*y)* sin(3*π*z)
end

# Avoid `Plots.coords` scope collision
x = Grids.coords(g.c[1], Center1D()).h
y = Grids.coords(g.c[2], Center1D()).h
z = Grids.coords(g.c[3], Center1D()).h

plot(y,z, midslice(cent, 1), seriestype=:wireframe) # nice option
savefig("wireframe.svg")
contourf(y,z, midslice(cent, 1), c = :viridis)
savefig("contourf.svg")
