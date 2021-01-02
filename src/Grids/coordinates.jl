#### Coordinates

export Coordinates, coords
export h_extreme, boundary_pts, over_points
export DataLocation, Center1D, Vertex1D
export CellBoundary, NodeBoundary
export All, Ghost, Boundary, Interior

"""
    DataLocation

Subtypes are used for dispatching
on the data location
"""
abstract type DataLocation end

""" 1D cell center """
struct Center1D <: DataLocation end

""" 1D cell vertex """
struct Vertex1D <: DataLocation end

abstract type CollocatedCoordinates{DL,FT} end

struct Coordinates1D{DL,FT,VFT,IT} <: CollocatedCoordinates{DL,FT}
  s::IT
  h::VFT
  Δh::VFT
end

function Base.show(io::IO, c::CollocatedCoordinates)
  print(io, "\n")
  println(io, "  s     = ",c.s)
  println(io, "  h     = ",c.h)
  println(io, "  Δh    = ",c.Δh)
end

struct SparseMat{D,U}
  diag::D
  upper::U
end

"""
    Coordinates

Cell node (`n`) and center (`c`)
1D coordinates.
"""
struct Coordinates{N, C, I, D}
  n::N
  c::C
  interp::I
  deriv::D
end

coords(c::Coordinates, ::Center1D) = c.c
coords(c::Coordinates, ::Vertex1D) = c.n

Base.length(c::Coordinates, ::Center1D) = c.c.s
Base.length(c::Coordinates, ::Vertex1D) = c.n.s

include("extremum_funcs.jl")

abstract type BoundaryPoints end

""" Cell center boundary type, containing ghost and first interior indices """
struct CellBoundary{IT<:Integer} <: BoundaryPoints
  "ghost index"
  g::IT
  "interior index"
  i::IT
end

""" Node boundary type, containing ghost, boundary, and first interior indices """
struct NodeBoundary{IT<:Integer} <: BoundaryPoints
  "ghost index"
  g::IT
  "boundary index"
  b::IT
  "interior index"
  i::IT
end

boundary_pts(::Coordinates1D{Center1D}, ::Min) = CellBoundary(1, 2)
boundary_pts(::Coordinates1D{Vertex1D}, ::Min) = NodeBoundary(1, 2, 3)

boundary_pts(c::Coordinates1D{Center1D}, ::Max) = CellBoundary(c.s, c.s-1)
boundary_pts(n::Coordinates1D{Vertex1D}, ::Max) = NodeBoundary(n.s, n.s-1, n.s-2)

h_extreme(n::Coordinates1D{Vertex1D}, e::Extremum) = n.h[boundary_pts(n,e).b]

abstract type DomainDecomp end
""" Type for dispatching on all points """
struct All <: DomainDecomp end
""" Type for dispatching on boundary points """
struct Boundary <: DomainDecomp end
""" Type for dispatching on interior points """
struct Interior <: DomainDecomp end
""" Type for dispatching on ghost points """
struct Ghost <: DomainDecomp end

function over_points(c::Coordinates, ::All     , dl::DataLocation)
  c = coords(c, dl)
  boundary_pts(c,Min()).g:boundary_pts(c,Max()).g
end
function over_points(c::Coordinates, ::Boundary, dl::DataLocation)
  c = coords(c, dl)
  boundary_pts(c,Min()).b:boundary_pts(c,Max()).b
end
function over_points(c::Coordinates, ::Interior, dl::DataLocation)
  c = coords(c, dl)
  boundary_pts(c,Min()).i:boundary_pts(c,Max()).i
end
function over_points(c::Coordinates, ::Ghost   , dl::DataLocation)
  c = coords(c, dl)
  b_pts_zmin = boundary_pts(c,Min())
  b_pts_zmax = boundary_pts(c,Max())
  return (b_pts_zmin.g:b_pts_zmin.g,
          b_pts_zmax.g:b_pts_zmax.g)
end

function pad!(h::Vector{FT}) where FT
  pushfirst!(h, h[1] - (h[2]-h[1]))
  push!(h, h[end] + (h[end]-h[end-1]))
end

function Coordinates(a::FT, b::FT, n::IT; warpfun::F=nothing, args=nothing) where {FT, IT<:Integer, F}
  sn = n+3

  if !(warpfun==nothing)
    hn = warpfun(a, b, n, args...)
  else
    hn = uniform(a,b,n)
  end
  pad!(hn)
  @assert length(hn) == n+3
  Δhn = FT[hn[i+1]-hn[i] for i in 1:sn-1]
  cn = Coordinates1D{Vertex1D,FT,typeof(hn),typeof(sn)}(sn, hn, Δhn)

  sc = n+2
  hc = FT[hn[i]+Δhn[i]/2 for i in 1:sc]
  @assert length(hc) == n+2
  Δhc = FT[hc[i+1]-hc[i] for i in 1:sc-1]
  cc = Coordinates1D{Center1D,FT,typeof(hc),typeof(sc)}(sc, hc, Δhc)

  if !(all(isfinite.(hn)) && all(isfinite.(hc)))
    throw(ArgumentError("Coordinates are not finite."))
  end

  interp = SparseMat(interp_diag(hn, hc, sc), interp_upper_diag(hn, hc, sc))
  deriv = (
    c2v=SparseMat(deriv_diag_C2V(Δhc, sc), deriv_upper_diag_C2V(Δhc, sc)),
    v2c=SparseMat(deriv_diag_V2C(Δhn, sn), deriv_upper_diag_V2C(Δhn, sn)),
    )
  return Coordinates(cn, cc, interp, deriv)
end

interp_diag(hn, hc, sc) =
    [(hn[i+1] - hc[i])/(hc[i+1] - hc[i]) for i in 1:sc-1]

interp_upper_diag(hn, hc, sc) = 1 .- interp_diag(hn, hc, sc)

deriv_diag_C2V(Δh, sc) = [-1/Δh[i] for i in 1:sc-1]
deriv_upper_diag_C2V(Δh, sc) = [1/Δh[i] for i in 1:sc-1]

deriv_diag_V2C(Δh, sc) = [-1/Δh[i] for i in 1:sc-1]
deriv_upper_diag_V2C(Δh, sc) = [1/Δh[i] for i in 1:sc-1]

