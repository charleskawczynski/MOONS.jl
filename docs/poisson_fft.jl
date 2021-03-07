using Plots, UnPack, FFTW
using MOONS.Grids
using MOONS.Fields
using MOONS.FFTs

FT = Float64
n = 128
cx = Coordinates(FT(0),FT(1),n)
cy = Coordinates(FT(0),FT(1),n)
cz = Coordinates(FT(0),FT(1),2)
g = Grid((cx,cy,cz))

f_exact = CellCenter(g)
f_recreate = CellCenter(g)
u_exact = CellCenter(g)
f = CellCenter(g)
u = CellCenter(g)

@inbounds for local_grid in GridField(g,f_exact)
    @unpack ijk,x,y = local_grid
    u_exact[ijk] = cos(3*π*x)*cos(3*π*y)
    f_exact[ijk] = -2*(3*π)^2*cos(3*π*x)*cos(3*π*y)
end
∇²!(f, u_exact, g)
extrap!(f, g)
solve!(u, f_exact, g, 3)

∇²!(f_recreate, u, g)

@unpack xc, yc = all_coords(g)

contourf(xc,yc, midslice(f, 3); c = :viridis)
savefig("f_discrete.svg")

contourf(xc,yc, midslice(f_exact, 3); c = :viridis)
savefig("f_exact.svg")

u_ave = average(u, 3);
@show max(u_ave...)
max(abs.(diff(u; dims=3))...) ≈ 0 || @warn "solution is not 2D"

contourf(xc,yc, u_ave; c = :viridis)
savefig("u_discrete.svg")

contourf(xc,yc, midslice(u_exact, 3); c = :viridis, clims=(-1,1))
savefig("u_exact.svg")

u_err = CellCenter(g)
u_err .= abs.(u_exact - u);
assign!(u_err, 0, Ghost())
@show max(u_err...)

contourf(xc,yc, midslice(u_err, 3); c = :viridis)
savefig("u_err.svg")
