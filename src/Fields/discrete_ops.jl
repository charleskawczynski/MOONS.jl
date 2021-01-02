##### Discrete operators

export ∇!
export ∇²!

#####
##### cent <-> face
#####

function ∇!(face::CellFace{∇_dir}, cent::CellCenter, grid::Grid) where {∇_dir}
    iR = over_points(grid.c[∇_dir], Boundary(), Vertex1D())
    for (Ipre, i, Ipost) in sweep(cent, ∇_dir, iR)
        D = grid.c[∇_dir].deriv.c2v.diag[i-1]
        U = grid.c[∇_dir].deriv.c2v.upper[i-1]
        face[Ipre, i, Ipost] = D*cent[Ipre, i-1, Ipost] + U*cent[Ipre, i, Ipost]
    end
    nothing
end

function ∇!(cent::CellCenter, face::CellFace{∇_dir}, grid::Grid) where {∇_dir}
    for (Ipre, i, Ipost) in sweep(cent, ∇_dir)
        D = grid.c[∇_dir].deriv.v2c.diag[i]
        U = grid.c[∇_dir].deriv.v2c.upper[i]
        cent[Ipre, i, Ipost] = D*face[Ipre, i, Ipost] + U*face[Ipre, i+1, Ipost]
    end
    nothing
end

#####
##### corn <-> edge
#####

function ∇!(edge::CellEdge{∇_dir}, corn::CellCorner, grid::Grid) where {∇_dir}
    for (Ipre, i, Ipost) in sweep(edge, ∇_dir)
        D = grid.c[∇_dir].deriv.v2c.diag[i]
        U = grid.c[∇_dir].deriv.v2c.upper[i]
        edge[Ipre, i, Ipost] = D*corn[Ipre, i, Ipost] + U*corn[Ipre, i+1, Ipost]
    end
    nothing
end

function ∇!(corn::CellCorner, edge::CellEdge{∇_dir}, grid::Grid) where {∇_dir}
    iR = over_points(grid.c[∇_dir], Boundary(), Vertex1D())
    for (Ipre, i, Ipost) in sweep(edge, ∇_dir, iR)
        D = grid.c[∇_dir].deriv.c2v.diag[i-1]
        U = grid.c[∇_dir].deriv.c2v.upper[i-1]
        corn[Ipre, i, Ipost] = D*edge[Ipre, i-1, Ipost] + U*edge[Ipre, i, Ipost]
    end
    nothing
end

#####
##### face <-> edge
#####

∇!(edge::CellEdge, face::CellFace, grid::Grid, ∇_dir::Int) =
    ∇!(edge, face, grid, Val(∇_dir))

∇!(edge::CellEdge{dir}, face::CellFace{dir}, grid::Grid, ::Val{dir}) where {dir} =
    MethodError("∇(CellFace{dir}) does not live on CellEdge{dir}!")

function ∇!(edge::CellEdge{edgedir}, face::CellFace{facedir}, grid::Grid, ::Val{∇_dir}) where {edgedir,facedir,∇_dir}
    iR = over_points(grid.c[∇_dir], Boundary(), Vertex1D())
    for (Ipre, i, Ipost) in sweep(edge, ∇_dir, iR)
        D = grid.c[∇_dir].deriv.c2v.diag[i-1]
        U = grid.c[∇_dir].deriv.c2v.upper[i-1]
        edge[Ipre, i, Ipost] = D*face[Ipre, i-1, Ipost] + U*face[Ipre, i, Ipost]
    end
    nothing
end

#####
##### Generic forwarding
#####

VecTup(T) = Tuple{<:T,<:T,<:T}

function ∇!(∇f::VecTup(T∇F), f::TF, grid::Grid) where {T∇F,TF}
    ∇!(∇f[1], f, grid)
    ∇!(∇f[2], f, grid)
    ∇!(∇f[3], f, grid)
end

#####
##### ∇²!
#####

real_domain(::Center1D) = Interior()
real_domain(::Vertex1D) = Boundary()
is_cent(::Vertex1D) = 0
is_cent(::Center1D) = 1

# TODO: this can be much better optimized.
function ∇²!(∇²f::T, f::T, grid::Grid) where {T<:AbstractField}
    p_along = primary_along(f)
    d_along = dual_along(f)
    i_cent = is_cent.(d_along)

    @inbounds for k in over_points(grid.c[3], real_domain(p_along[3]), p_along[3])
    @inbounds for j in over_points(grid.c[2], real_domain(p_along[2]), p_along[2])
    @inbounds for i in over_points(grid.c[1], real_domain(p_along[1]), p_along[1])
        Δx⁺ = coords(grid.c[1], p_along[1]).Δh[i]
        Δy⁺ = coords(grid.c[2], p_along[2]).Δh[j]
        Δz⁺ = coords(grid.c[3], p_along[3]).Δh[k]
        Δx⁻ = coords(grid.c[1], p_along[1]).Δh[i-1]
        Δy⁻ = coords(grid.c[2], p_along[2]).Δh[j-1]
        Δz⁻ = coords(grid.c[3], p_along[3]).Δh[k-1]
        Δx = coords(grid.c[1], d_along[1]).Δh[i-i_cent[1]]
        Δy = coords(grid.c[2], d_along[2]).Δh[j-i_cent[2]]
        Δz = coords(grid.c[3], d_along[3]).Δh[k-i_cent[3]]
        f_ijk = f[i,j,k]
        ∇fx⁺ = (f[i+1,j,k]-f_ijk)/Δx⁺
        ∇fy⁺ = (f[i,j+1,k]-f_ijk)/Δy⁺
        ∇fz⁺ = (f[i,j,k+1]-f_ijk)/Δz⁺
        ∇fx⁻ = (f_ijk-f[i-1,j,k])/Δx⁻
        ∇fy⁻ = (f_ijk-f[i,j-1,k])/Δy⁻
        ∇fz⁻ = (f_ijk-f[i,j,k-1])/Δz⁻
        ∇²f_x = (∇fx⁺-∇fx⁻)/Δx
        ∇²f_y = (∇fy⁺-∇fy⁻)/Δy
        ∇²f_z = (∇fz⁺-∇fz⁻)/Δz
        ∇²f[i,j,k] = ∇²f_x+∇²f_y+∇²f_z
    end
    end
    end
    nothing
end
