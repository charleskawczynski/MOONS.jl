module Grids

using DocStringExtensions
export Grid

include("directional_funcs.jl")
include("stretch_funcs.jl")
include("stretch_params.jl")
include("stretch_param_match.jl")
include("coordinates.jl")

"""
    Grid

An N-dimensional Cartesian, structured, computational grid.
"""
struct Grid{T, N}
  c::NTuple{N,T}
end

function Base.show(io::IO, g::Grid{T,N}) where {T,N}
  println(io, "---------------------- Grid")
  for i in 1:N
    println(io, g.c[i])
  end
  println(io, "---------------------- end grid")
end

end