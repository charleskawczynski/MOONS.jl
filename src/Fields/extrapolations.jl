##### Extrapolations

export extrap!
export pre_post_colons

function extrap! end

#####
##### Self-extrapolation
#####

pre_post_colons(f, dim) = (;
    Ipre = ntuple(i->Colon(), Val(length(size(f)[1:dim-1]))),
    Ipost = ntuple(i->Colon(), Val(length(size(f)[dim+1:end])))
)

""" Simple linear extrapolation `f_g = 2f_b - f_i` """
extrap!(f::AbstractField, grid::Grid, dir::Int) =
    extrap!(f, grid, Val(dir))

function extrap!(f::AbstractField, grid::Grid)
    extrap!(f, grid, 1)
    extrap!(f, grid, 2)
    extrap!(f, grid, 3)
end

function extrap!(f::AbstractField, grid::Grid, ::Val{dim}) where {dim}
    @unpack Ipre, Ipost = pre_post_colons(f, dim)
    f[Ipre..., 1 , Ipost...] .= 2*f[Ipre...,  2  , Ipost...] .- f[Ipre...,  3  , Ipost...]
    f[Ipre...,end, Ipost...] .= 2*f[Ipre...,end-1, Ipost...] .- f[Ipre...,end-2, Ipost...]
    nothing
end

""" Simple linear extrapolation with two staggered fields: `f_g = 2g_b - f_i` """
extrap!(f::AbstractField, g::AbstractField, grid::Grid, dir::Int) =
    extrap!(f, g, grid, Val(dir))

extrap!(face::CellFace{dim}, cent::CellCenter, grid::Grid, ::Val{dim}) where {dim} =
    extrap!(face, cent, grid, Val(dim), Vertex1D(), Center1D())
extrap!(cent::CellCenter, face::CellFace{dim}, grid::Grid, ::Val{dim}) where {dim} =
    extrap!(cent, face, grid, Val(dim), Center1D(), Vertex1D())
extrap!(corn::CellCorner, edge::CellEdge{dim}, grid::Grid, ::Val{dim}) where {dim} =
    extrap!(corn, edge, grid, Val(dim), Vertex1D(), Center1D())
extrap!(edge::CellEdge{dim}, corn::CellCorner, grid::Grid, ::Val{dim}) where {dim} =
    extrap!(edge, corn, grid, Val(dim), Center1D(), Vertex1D())

extrap!(edge::CellEdge{dir}, face::CellFace{dir}, grid::Grid, ::Val{dir}) where {dir} =
    throw(MethodError("Cannot extrap CellEdge{dir} with CellFace{dir}"))
extrap!(face::CellFace{dir}, edge::CellEdge{dir}, grid::Grid, ::Val{dir}) where {dir} =
    throw(MethodError("Cannot extrap CellFace{dir} with CellEdge{dir}"))

extrap!(edge::CellEdge{edgedir}, face::CellFace{facedir}, grid::Grid, ::Val{dim}) where {dim,edgedir,facedir} =
    extrap!(edge, face, grid, Val(dim), Vertex1D(), Center1D())

extrap!(face::CellFace{facedir}, edge::CellEdge{edgedir}, grid::Grid, ::Val{dim}) where {dim,edgedir,facedir} =
    extrap!(face, edge, grid, Val(dim), Center1D(), Vertex1D())

function extrap!(f, g, grid::Grid, ::Val{dim}, ::Center1D, ::Vertex1D) where {dim}
    @unpack Ipre, Ipost = pre_post_colons(f, dim)
    f[Ipre..., 1 , Ipost...] .= 2*g[Ipre...,  2  , Ipost...] .- f[Ipre...,  2  , Ipost...]
    f[Ipre...,end, Ipost...] .= 2*g[Ipre...,end-1, Ipost...] .- f[Ipre...,end-1, Ipost...]
    nothing
end

function extrap!(f, g, grid::Grid, ::Val{dim}, ::Vertex1D, ::Center1D) where {dim}
    @unpack Ipre, Ipost = pre_post_colons(f, dim)
    f[Ipre..., 1 , Ipost...] .= 2*g[Ipre..., 1 , Ipost...] .- f[Ipre...,  2  , Ipost...]
    f[Ipre...,end, Ipost...] .= 2*g[Ipre...,end, Ipost...] .- f[Ipre...,end-1, Ipost...]
    nothing
end
