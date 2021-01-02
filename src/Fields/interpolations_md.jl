##### Multi-dimensional interpolations

# TODO: Add implementations that perform interps with a single pass and no alloc

#####
##### CellCenter <-> CellEdge
#####

function interp!(cent::CellCenter, edge::CellEdge{dim}, grid::Grid) where {dim}
    face_orth = CellFace(grid, first(orthogs(dim)))
    interp!(face_orth, edge, grid)
    interp!(cent, face_orth, grid)
end

function interp!(edge::CellEdge{dim}, cent::CellCenter, grid::Grid) where {dim}
    face_orth = CellFace(grid, first(orthogs(dim)))
    interp!(face_orth, cent, grid)
    interp!(edge, face_orth, grid)
end

#####
##### CellCorner <-> CellFace
#####

function interp!(corn::CellCorner, face::CellFace{dim}, grid::Grid) where {dim}
    edge_orth = CellEdge(grid, first(orthogs(dim)))
    interp!(edge_orth, face, grid)
    interp!(corn, edge_orth, grid)
end

function interp!(face::CellFace{dim}, corn::CellCorner, grid::Grid) where {dim}
    edge_orth = CellEdge(grid, first(orthogs(dim)))
    interp!(edge_orth, corn, grid)
    interp!(face, edge_orth, grid)
end

#####
##### CellFace <-> CellEdge
#####

function interp!(face::CellFace{dim}, edge::CellEdge{dim}, grid::Grid) where {dim}
    face_orth = CellFace(grid, first(orthogs(dim)))
    cent = CellCenter(grid)
    interp!(face_orth, edge, grid)
    interp!(cent, face_orth, grid)
    interp!(face, cent, grid)
end

function interp!(edge::CellEdge{dim}, face::CellFace{dim}, grid::Grid) where {dim}
    cent = CellCenter(grid)
    face_orth = CellFace(grid, first(orthogs(dim)))
    interp!(cent, face, grid)
    interp!(face_orth, cent, grid)
    interp!(edge, face_orth, grid)
end

#####
##### CellCenter <-> CellCorner
#####

function interp!(cent::CellCenter, corn::CellCorner, grid::Grid)
    edge = CellEdge(grid, 1)
    interp!(edge, corn, grid)
    interp!(cent, edge, grid)
end

function interp!(corn::CellCorner, cent::CellCenter, grid::Grid)
    edge = CellEdge(grid, 1)
    interp!(edge, cent, grid)
    interp!(corn, edge, grid)
end
