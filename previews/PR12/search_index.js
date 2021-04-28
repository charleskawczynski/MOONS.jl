var documenterSearchIndex = {"docs":
[{"location":"api/#API","page":"API","title":"API","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"CurrentModule = MOONS","category":"page"},{"location":"api/","page":"API","title":"API","text":"Grids.DataLocation\nGrids.Center1D\nGrids.Vertex1D\nGrids.Coordinates\nGrids.ghost_dual\nGrids.binary\nGrids.n_hat\nGrids.ghost_vec\nGrids.Boundary\nGrids.Interior\nGrids.All\nGrids.NodeBoundary\nGrids.CellBoundary\nGrids.Ghost\nGrids.Grid\nGrids.Max\nGrids.Min","category":"page"},{"location":"api/#MOONS.Grids.DataLocation","page":"API","title":"MOONS.Grids.DataLocation","text":"DataLocation\n\nSubtypes are used for dispatching on the data location\n\n\n\n\n\n","category":"type"},{"location":"api/#MOONS.Grids.Center1D","page":"API","title":"MOONS.Grids.Center1D","text":"1D cell center \n\n\n\n\n\n","category":"type"},{"location":"api/#MOONS.Grids.Vertex1D","page":"API","title":"MOONS.Grids.Vertex1D","text":"1D cell vertex \n\n\n\n\n\n","category":"type"},{"location":"api/#MOONS.Grids.Coordinates","page":"API","title":"MOONS.Grids.Coordinates","text":"Coordinates\n\nCell node (n) and center (c) 1D coordinates.\n\n\n\n\n\n","category":"type"},{"location":"api/#MOONS.Grids.ghost_dual","page":"API","title":"MOONS.Grids.ghost_dual","text":"ghost_dual(::Extremum)\n\nA 2-element vector near the boundary with 1's on ghost cells and 0's on interior cells\n\n\n\n\n\n","category":"function"},{"location":"api/#MOONS.Grids.binary","page":"API","title":"MOONS.Grids.binary","text":"binary(::Extremum)\n\nReturns 0 for Min and 1 for Max\n\n\n\n\n\n","category":"function"},{"location":"api/#MOONS.Grids.n_hat","page":"API","title":"MOONS.Grids.n_hat","text":"n_hat(::Extremum)\n\nThe outward normal vector to the boundary\n\n\n\n\n\n","category":"function"},{"location":"api/#MOONS.Grids.ghost_vec","page":"API","title":"MOONS.Grids.ghost_vec","text":"ghost_vec(::Extremum)\n\nA 3-element vector near the boundary with 1's on ghost cells and 0's on interior cells\n\n\n\n\n\n","category":"function"},{"location":"api/#MOONS.Grids.Boundary","page":"API","title":"MOONS.Grids.Boundary","text":"Type for dispatching on boundary points \n\n\n\n\n\n","category":"type"},{"location":"api/#MOONS.Grids.Interior","page":"API","title":"MOONS.Grids.Interior","text":"Type for dispatching on interior points \n\n\n\n\n\n","category":"type"},{"location":"api/#MOONS.Grids.All","page":"API","title":"MOONS.Grids.All","text":"Type for dispatching on all points \n\n\n\n\n\n","category":"type"},{"location":"api/#MOONS.Grids.NodeBoundary","page":"API","title":"MOONS.Grids.NodeBoundary","text":"Node boundary type, containing ghost, boundary, and first interior indicies \n\n\n\n\n\n","category":"type"},{"location":"api/#MOONS.Grids.CellBoundary","page":"API","title":"MOONS.Grids.CellBoundary","text":"Cell center boundary type, containing ghost and first interior indicies \n\n\n\n\n\n","category":"type"},{"location":"api/#MOONS.Grids.Ghost","page":"API","title":"MOONS.Grids.Ghost","text":"Type for dispatching on ghost points \n\n\n\n\n\n","category":"type"},{"location":"api/#MOONS.Grids.Grid","page":"API","title":"MOONS.Grids.Grid","text":"Grid\n\nAn N-dimensional Cartesian, structured, computational grid.\n\n\n\n\n\n","category":"type"},{"location":"api/#MOONS.Grids.Max","page":"API","title":"MOONS.Grids.Max","text":"Max\n\nA type for dispatching on the maximum of the domain along a particular direction\n\n\n\n\n\n","category":"type"},{"location":"api/#MOONS.Grids.Min","page":"API","title":"MOONS.Grids.Min","text":"Min\n\nA type for dispatching on the minimum of the domain along a particular direction\n\n\n\n\n\n","category":"type"},{"location":"grids/#Grids","page":"Grids","title":"Grids","text":"","category":"section"},{"location":"grids/","page":"Grids","title":"Grids","text":"CurrentModule = MOONS","category":"page"},{"location":"grids/","page":"Grids","title":"Grids","text":"using MOONS.Grids\ninclude(joinpath(\"..\", \"plothelpers.jl\"))\ninclude(joinpath(\"..\", \"grids.jl\"))\nnothing","category":"page"},{"location":"grids/","page":"Grids","title":"Grids","text":"(Image: ) (Image: ) (Image: )","category":"page"},{"location":"#MOONS.jl","page":"Home","title":"MOONS.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"MOONS.jl is a pure julia implementation of MOONS.","category":"page"}]
}
