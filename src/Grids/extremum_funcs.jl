#### ExtremumFuncs

export Extremum, Min, Max
export n_hat, binary, ghost_vec, ghost_dual

abstract type Extremum end

"""
    Min

A type for dispatching on the minmimum of
the domain along a particular direction
"""
struct Min <: Extremum end

"""
    Max

A type for dispatching on the maximum of
the domain along a particular direction
"""
struct Max <: Extremum end

"""
    n_hat(::Extremum)

The outward normal vector to the boundary
"""
n_hat(::Min) = -1
n_hat(::Max) = 1

"""
    binary(::Extremum)

Returns 0 for `Min` and 1 for `Max`
"""
binary(::Min) = 0
binary(::Max) = 1

"""
    ghost_vec(::Extremum)

A 3-element vector near the boundary with 1's
on ghost cells and 0's on interior cells
"""
ghost_vec(::Min) = (1, 0, 0)
ghost_vec(::Max) = (0, 0, 1)

"""
    ghost_dual(::Extremum)

A 2-element vector near the boundary with 1's
on ghost cells and 0's on interior cells
"""
ghost_dual(::Min) = (1, 0)
ghost_dual(::Max) = (0, 1)
