using Plots

abstract type AbstractGridPlot end
struct Isosurfaces3D <: AbstractGridPlot end
struct Streamlines3D <: AbstractGridPlot end
struct GridVis{slice} <: AbstractGridPlot end
struct Contour2D{slice} <: AbstractGridPlot end
struct Profile1D{axis} <: AbstractGridPlot end

plot_grid(g::Grid, filename::AbstractString, grid_viz::GridVis{slice}) where {F, slice} = nothing
plot_func(f::F, g::Grid, filename::AbstractString, ::GridVis{slice}) where {F, slice} = nothing

function plot_grid(g::Grid, filename::AbstractString, grid_viz::GridVis{slice}) where {F, slice}
  return plot_func((x,y) -> 1, g, filename, grid_viz, Plots.wireframe)
end

function plot_func(f::F,
                   g::Grid,
                   filename::AbstractString,
                   ::GridVis{slice},
                   viz_type::G=contourf) where {G, F, slice}
    i_orth = orthogs(Val(slice))
    _x = g.c[i_orth[1]]
    _y = g.c[i_orth[2]]
    x = _x.n
    y = _y.n
    data = [f(x.h[i], y.h[i]) for i in over_points(_x, Interior(), Vertex1D()),
                                  j in over_points(_y, Interior(), Vertex1D())]
    X = [x.h[i] for i in over_points(_x, Interior(), Vertex1D())]
    Y = [y.h[j] for j in over_points(_y, Interior(), Vertex1D())]
    viz_type(X,Y,data, camera=(90, 90), size=(2*320,2*300), aspect_ratio=:equal)
    viz_type(X,Y,data)
    savefig(filename)
end
