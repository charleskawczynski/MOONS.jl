module Fields

using ..Grids

export AbstractField, CellCenter, CellCorner, CellFace, CellEdge
export midslice
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

abstract type AbstractField end

Base.size(f::AbstractField) = Base.size(f.data)
Base.size(f::AbstractField, dim) = Base.size(f.data, dim)
Base.length(f::AbstractField) = Base.length(f.data)
Base.iterate(f::AbstractField, state=1) =
    state > length(f) ? nothing : (f[state], state+1)

Base.@propagate_inbounds Base.getindex(f::AbstractField, key...) =
    Base.getindex(f.data, key...)
Base.@propagate_inbounds Base.setindex!(f::AbstractField, v, key...) =
    Base.setindex!(f.data, v, key...)

struct CellCenter{A} <: AbstractField
    data::A
    function CellCenter(grid::Grid{FT}) where {FT}
        data = zeros(FT, size(grid, get_dl(Primary()))...)
        new{typeof(data)}(data)
    end
end

struct CellCorner{A} <: AbstractField
    data::A
    function CellCorner(grid::Grid{FT}) where {FT}
        data = zeros(FT, size(grid, get_dl(Dual()))...)
        new{typeof(data)}(data)
    end
end

struct CellFace{A,dir} <: AbstractField
    data::A
    function CellFace(grid::Grid{FT}, dir::Int) where {FT}
        @assert 1 <= dir <= 3
        data = zeros(FT, size(grid, get_dl(Primary(), dir))...)
        new{typeof(data),dir}(data)
    end
end

struct CellEdge{A,dir} <: AbstractField
    data::A
    function CellEdge(grid::Grid{FT}, dir::Int) where {FT}
        @assert 1 <= dir <= 3
        data = zeros(FT, size(grid, get_dl(Dual(), dir))...)
        new{typeof(data),dir}(data)
    end
end

direction(f::CellCenter) = nothing
direction(f::CellCorner) = nothing
direction(f::CellFace{A,dir}) where {A,dir} = dir
direction(f::CellEdge{A,dir}) where {A,dir} = dir

midslice(f::AbstractField, dir::Int) = midslice(f, Val(dir))
midslice(f::AbstractField, ::Val{1}) = f[round(Int,size(f, 1)/2),:,:]
midslice(f::AbstractField, ::Val{2}) = f[:,round(Int,size(f, 2)/2),:]
midslice(f::AbstractField, ::Val{3}) = f[:,:,round(Int,size(f, 2)/2)]

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

space(grid::Grid, ::Center1D) = grid.cent
space(grid::Grid, ::Vertex1D) = grid.corn

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

end