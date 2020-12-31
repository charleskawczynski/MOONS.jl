module Grids

using DocStringExtensions
export Grid, Cube, SquarePlane

include("directional_funcs.jl")
include("stretch_funcs.jl")
include("stretch_params.jl")
include("stretch_param_match.jl")
include("coordinates.jl")

"""
    Grid

An N-dimensional Cartesian, structured, computational grid.
"""
struct Grid{FT, T, N, CE, CO}
  c::NTuple{N,T}
  cent::CE
  corn::CO
end
function Grid(c::NTuple{N,T}) where {N,T}
    cent = CenterSpace(c)
    corn = VertexSpace(c)
    Grid{eltype(c[1].c.h),T,N,typeof(cent), typeof(corn)}(c, cent, corn)
end

abstract type Dimension end
struct XDim <: Dimension end
struct YDim <: Dimension end
struct ZDim <: Dimension end

struct Space{X,Y,Z}
    x::X
    y::Y
    z::Z
end
CenterSpace(c::Tuple) = Space(
        Space1D{XDim}(c[1].c.h),
        Space1D{YDim}(c[2].c.h),
        Space1D{ZDim}(c[3].c.h)
    )
VertexSpace(c::Tuple) = Space(
        Space1D{XDim}(c[1].n.h),
        Space1D{YDim}(c[2].n.h),
        Space1D{ZDim}(c[3].n.h)
    )

struct Space1D{dim,H}
    h::H
end
Space1D{dim}(h::H) where {H,dim} = Space1D{dim,H}(h)

# Allow for field-like indexing:
Base.@propagate_inbounds Base.getindex(c::Space1D{XDim}, i,j,k) = Base.getindex(c.h, i)
Base.@propagate_inbounds Base.getindex(c::Space1D{YDim}, i,j,k) = Base.getindex(c.h, j)
Base.@propagate_inbounds Base.getindex(c::Space1D{ZDim}, i,j,k) = Base.getindex(c.h, k)

function Cube(a::FT, b::FT, n::Int; kwargs...) where {FT}
    c = Coordinates(a,b,n; kwargs...)
    return Grid((c,c,c))
end

function SquarePlane(a::FT, b::FT, n::Int; kwargs...) where {FT}
    c_plane = Coordinates(FT(-0.5),FT(0.5),1)
    c = Coordinates(a,b,n; kwargs...)
    return Grid((c,c,c_plane))
end

function Base.size(grid::Grid, t::Tuple)
    @assert length(t) == 3
    return (length(grid.c[k], dl) for (k,dl) in enumerate(t))
end

function Base.show(io::IO, g::Grid{FT,T,N}) where {FT,T,N}
  println(io, "---------------------- Grid")
  for i in 1:N
    println(io, g.c[i])
  end
  println(io, "---------------------- end grid")
end

end