module Fields

import LinearAlgebra
using UnPack
using ..Grids

using Base.Broadcast: Broadcasted, BroadcastStyle, ArrayStyle

export AbstractField, CellCenter, CellCorner, CellFace, CellEdge
export CellFaces, CellEdges, CellCenters, CellCorners
export dual_fields
export average, midslice, midline
export primary_along, dual_along
export GridField

abstract type FieldLocations end
struct Center <: FieldLocations end
struct Face <: FieldLocations end
struct Edge <: FieldLocations end
struct Corner <: FieldLocations end

abstract type GridType end
struct Primary <: GridType end
struct Dual <: GridType end

"""
    get_dl

A 3-element tuple of `DataLocation`'s,
based on the grid type and direction.
"""
function get_dl end

get_dl(::Primary) = (Center1D(),Center1D(),Center1D())
get_dl(::Dual) = (Vertex1D(),Vertex1D(),Vertex1D())
get_dl(::Primary, dir::Val{1}) = (Vertex1D(),Center1D(),Center1D())
get_dl(::Primary, dir::Val{2}) = (Center1D(),Vertex1D(),Center1D())
get_dl(::Primary, dir::Val{3}) = (Center1D(),Center1D(),Vertex1D())
get_dl(::Dual, dir::Val{1}) = (Center1D(),Vertex1D(),Vertex1D())
get_dl(::Dual, dir::Val{2}) = (Vertex1D(),Center1D(),Vertex1D())
get_dl(::Dual, dir::Val{3}) = (Vertex1D(),Vertex1D(),Center1D())
get_dl(gt::GridType, dir::Int) = get_dl(gt, Val(dir))

abstract type AbstractField{FT,NDIMS} <: AbstractArray{FT, NDIMS} end

Base.size(f::AbstractField) = Base.size(f.data)
Base.size(f::AbstractField, dim) = Base.size(f.data, dim)
Base.length(f::AbstractField) = Base.length(f.data)
Base.eltype(f::AbstractField) = Base.eltype(f.data)
Base.ndims(::Type{T}) where {NDIMS,FT,T<:AbstractField{FT,NDIMS}} = NDIMS
Base.iterate(f::AbstractField, state=1) =
    state > length(f) ? nothing : (f[state], state+1)
Base.lastindex(f::AbstractField, dim) = size(f, dim)

Base.copyto!(dest::Array, f::AbstractField) = copyto!(dest, f.data)
Base.BroadcastStyle(::Type{<:AbstractField}) = ArrayStyle{AbstractField}()

find_field(bc::Broadcasted) = find_field(bc.args)
find_field(args::Tuple) = find_field(find_field(args[1]), Base.tail(args))
find_field(x) = x
find_field(a::AbstractField, rest) = a
find_field(::Any, rest) = find_field(rest)

Base.similar(bc::Broadcasted{ArrayStyle{AbstractField}}, ::Type{FT}) where {FT} =
    similar(find_field(bc), FT)

function Base.copyto!(dest::AbstractField, src::AbstractField)
    copyto!(dest.data, src.data)
    dest
end

transform_array(f::AbstractField) = f.data
transform_array(x) = x

@inline function Base.copyto!(dest::AbstractField, bc::Broadcasted{Nothing})
    # check for the case a .= b, where b is an array
    if bc.f === identity && bc.args isa Tuple{AbstractArray}
        copyto!(dest.data, bc.args[1])
    else
        copyto!(dest.data, Broadcasted(bc.f, transform_array.(bc.args), bc.axes))
    end
    dest
end


LinearAlgebra.norm(f::AbstractField) = LinearAlgebra.norm(f.data)
Base.@propagate_inbounds Base.getindex(f::AbstractField, key...) =
    Base.getindex(f.data, key...)
Base.@propagate_inbounds Base.setindex!(f::AbstractField, v, key...) =
    Base.setindex!(f.data, v, key...)

struct CellCenter{A,FT,NDIMS} <: AbstractField{FT,NDIMS}
    data::A
    function CellCenter(grid::Grid{FT}) where {FT}
        s = size(grid, get_dl(Primary()))
        data = zeros(FT, s...)
        new{typeof(data),FT,length(s)}(data)
    end
end
CellCenters(grid::Grid) = (CellCenter(grid),CellCenter(grid),CellCenter(grid))

struct CellCorner{A,FT,NDIMS} <: AbstractField{FT,NDIMS}
    data::A
    function CellCorner(grid::Grid{FT}) where {FT}
        s = size(grid, get_dl(Dual()))
        data = zeros(FT, s...)
        new{typeof(data),FT,length(s)}(data)
    end
end
CellCorners(grid::Grid) = (CellCorner(grid),CellCorner(grid),CellCorner(grid))

struct CellFace{dir,A,FT,NDIMS} <: AbstractField{FT,NDIMS}
    data::A
    function CellFace(grid::Grid{FT}, dir::Int) where {FT}
        @assert 1 <= dir <= 3
        s = size(grid, get_dl(Primary(), dir))
        data = zeros(FT, s...)
        new{dir,typeof(data),FT,length(s)}(data)
    end
end
CellFaces(grid::Grid) = (CellFace(grid, 1),CellFace(grid, 2),CellFace(grid, 3))

struct CellEdge{dir,A,FT,NDIMS} <: AbstractField{FT,NDIMS}
    data::A
    function CellEdge(grid::Grid{FT}, dir::Int) where {FT}
        @assert 1 <= dir <= 3
        s = size(grid, get_dl(Dual(), dir))
        data = zeros(FT, s...)
        new{dir,typeof(data),FT,length(s)}(data)
    end
end
CellEdges(grid::Grid) = (CellEdge(grid, 1),CellEdge(grid, 2),CellEdge(grid, 3))

direction(f::CellCenter) = nothing
direction(f::CellCorner) = nothing
direction(f::CellFace{dir}) where {dir} = dir
direction(f::CellEdge{dir}) where {dir} = dir

midslice(f, dir::Int) = midslice(f, Val(dir))
midslice(f, ::Val{1}) = f[round(Int,size(f, 1)/2),:,:]
midslice(f, ::Val{2}) = f[:,round(Int,size(f, 2)/2),:]
midslice(f, ::Val{3}) = f[:,:,round(Int,size(f, 3)/2)]

average(f, dir::Int) = average(f, Val(dir))
average(f, ::Val{1}) = dropdims(sum(f; dims=1), dims=1) ./ size(f, 1)
average(f, ::Val{2}) = dropdims(sum(f; dims=2), dims=2) ./ size(f, 2)
average(f, ::Val{3}) = dropdims(sum(f; dims=3), dims=3) ./ size(f, 3)

midline(f, dir::Int) = midline(f, Val(dir))
midline(f, ::Val{1}) = f[:,round(Int,size(f, 2)/2),round(Int,size(f, 3)/2)]
midline(f, ::Val{2}) = f[round(Int,size(f, 1)/2),:,round(Int,size(f, 3)/2)]
midline(f, ::Val{3}) = f[round(Int,size(f, 1)/2),round(Int,size(f, 2)/2),:]

struct GridField{G,F,CI}
    grid::G
    f::F
    ci::CI
    function GridField(grid::Grid, f::AbstractField)
        ci = CartesianIndices(size(f))
        new{typeof(grid), typeof(f),typeof(ci)}(grid, f, ci)
    end
end

get_dl(::CellCenter, _) = get_dl(Primary())
get_dl(::CellCorner, _) = get_dl(Dual())
get_dl(::CellFace, ::Val{N}) where {N} = get_dl(Primary(), Val(N))
get_dl(::CellEdge, ::Val{N}) where {N} = get_dl(Dual(), Val(N))
get_dl(f::AbstractField, dir::Int) = get_dl(f, Val(dir))

primary_along(f::AbstractField) = (primary_along(f, 1),primary_along(f, 2),primary_along(f, 3))
primary_along(f::AbstractField, dir::Int) = primary_along(f, Val(dir))
primary_along(f::CellCenter, ::Val{dim}) where {dim} = Center1D()
primary_along(f::CellCorner, ::Val{dim}) where {dim} = Vertex1D()
primary_along(f::CellFace{dim}, ::Val{dim}) where {dim} = Vertex1D()
primary_along(f::CellFace{fdim}, ::Val{dim}) where {dim,fdim} = Center1D()
primary_along(f::CellEdge{dim}, ::Val{dim}) where {dim} = Center1D()
primary_along(f::CellEdge{fdim}, ::Val{dim}) where {dim,fdim} = Vertex1D()

dual_along(f::AbstractField) = (dual_along(f, 1),dual_along(f, 2),dual_along(f, 3))
dual_along(f::AbstractField, dir::Int) = dual_along(f, Val(dir))
dual_along(f::CellCenter, ::Val{dim}) where {dim} = Vertex1D()
dual_along(f::CellCorner, ::Val{dim}) where {dim} = Center1D()
dual_along(f::CellFace{dim}, ::Val{dim}) where {dim} = Center1D()
dual_along(f::CellFace{fdim}, ::Val{dim}) where {dim,fdim} = Vertex1D()
dual_along(f::CellEdge{dim}, ::Val{dim}) where {dim} = Vertex1D()
dual_along(f::CellEdge{fdim}, ::Val{dim}) where {dim,fdim} = Center1D()

space(grid::Grid, ::Center1D) = grid.cent
space(grid::Grid, ::Vertex1D) = grid.corn

dual_fields(grid::Grid, f::CellCenter) = CellFaces(grid)
dual_fields(grid::Grid, f::CellCorner) = CellEdges(grid)

function Base.iterate(gf::GridField{G,F}, state=1) where {G,F<:AbstractField}
    if state > length(gf)
        return nothing
    else
        ijk = gf.ci[state]
        dl = get_dl(gf.f, direction(gf.f))
        x = space(gf.grid, dl[1]).x[Tuple(ijk)...]
        y = space(gf.grid, dl[2]).y[Tuple(ijk)...]
        z = space(gf.grid, dl[3]).z[Tuple(ijk)...]
        return ((;x,y,z,ijk), state+1)
    end
end
Base.length(gf::GridField) = length(gf.f)

function sweep(A, dim, iR = 1:size(A, dim))
    Rpre = CartesianIndices(size(A)[1:dim-1])
    Rpost = CartesianIndices(size(A)[dim+1:end])
    return Iterators.product(Rpre, iR, Rpost)
end

# TODO: improve implementation
export assign!
function assign!(f::AbstractField, val, ::Ghost)
    f[1,:,:] .= val; f[end,:,:] .= val
    f[:,1,:] .= val; f[:,end,:] .= val
    f[:,:,1] .= val; f[:,:,end] .= val
end

include("interpolations_base.jl")
include("interpolations_md.jl")
include("extrapolations.jl")
include("discrete_ops.jl")

end