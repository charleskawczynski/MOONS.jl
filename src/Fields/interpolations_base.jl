##### Base interpolations (do not require intermediate / nested computations)

export interp!

function interp! end

interp!(A::T, B::T, grid::Grid) where {T <: AbstractField} = (A .= B)

#####
##### CellCenter <-> CellFace
#####

function interp!(cent::CellCenter, face::CellFace{dim}, grid::Grid) where {dim}
    FT = eltype(cent)
    @inbounds for (Ipre, i, Ipost) in sweep(cent, dim)
        cent[Ipre, i, Ipost] = FT(0.5)*face[Ipre, i, Ipost] +
                               FT(0.5)*face[Ipre, i+1, Ipost]
    end
    nothing
end

function interp!(face::CellFace{dim}, cent::CellCenter, grid::Grid) where {dim}
    iR = over_points(grid.c[dim], Interior(), Vertex1D())
    @inbounds for (Ipre, i, Ipost) in sweep(cent, dim, iR)
        D = grid.c[dim].interp.diag[i]
        U = grid.c[dim].interp.upper[i]
        face[Ipre, i, Ipost] = D*cent[Ipre, i-1, Ipost] + U*cent[Ipre, i, Ipost]
    end
    nothing
end

#####
##### CellCorner <-> CellEdge
#####

function interp!(corn::CellCorner, edge::CellEdge{dim}, grid::Grid) where {dim}
    iR = over_points(grid.c[dim], Interior(), Vertex1D())
    @inbounds for (Ipre, i, Ipost) in sweep(corn, dim, iR)
        D = grid.c[dim].interp.diag[i]
        U = grid.c[dim].interp.upper[i]
        corn[Ipre, i, Ipost] = D*edge[Ipre, i-1, Ipost] + U*edge[Ipre, i, Ipost]
    end
    nothing
end

function interp!(edge::CellEdge{dim}, corn::CellCorner, grid::Grid) where {dim}
    FT = eltype(corn)
    @inbounds for (Ipre, i, Ipost) in sweep(edge, dim)
        edge[Ipre, i, Ipost] = FT(0.5)*corn[Ipre, i, Ipost] +
                               FT(0.5)*corn[Ipre, i+1, Ipost]
    end
    nothing
end

#####
##### CellFace{facedim} <-> CellEdge{edgedim}, facedim ≠ edgedim
#####

function interp!(face::CellFace{facedim}, edge::CellEdge{edgedim}, grid::Grid) where {facedim,edgedim}
    FT = eltype(edge)
    interp_dim = orthog(facedim, edgedim)
    @inbounds for (Ipre, i, Ipost) in sweep(face, interp_dim)
        face[Ipre, i, Ipost] = FT(0.5)*edge[Ipre, i, Ipost] +
                               FT(0.5)*edge[Ipre, i+1, Ipost]
    end
    nothing
end

function interp!(edge::CellEdge{edgedim}, face::CellFace{facedim}, grid::Grid) where {facedim,edgedim}
    interp_dim = orthog(facedim, edgedim)
    iR = over_points(grid.c[interp_dim], Interior(), Vertex1D())
    @inbounds for (Ipre, i, Ipost) in sweep(edge, interp_dim, iR)
        D = grid.c[interp_dim].interp.diag[i]
        U = grid.c[interp_dim].interp.upper[i]
        edge[Ipre, i, Ipost] = D*face[Ipre, i-1, Ipost] + U*face[Ipre, i, Ipost]
    end
    nothing
end
