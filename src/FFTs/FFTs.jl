"""
    FFTs (Fast Fourier Transforms)

The Fast Fourier Transform of the scalar field, `f`,
w.r.t direction `dir` (1,2,3) which corresponds to (x,y,z).

# Implementation

```julia
fft!(ω,f,dir,pad)
```

# Input
```
f    = f(x,y,z)
dir  = direction along which to take the FT (1,2,3)
pad  = (1,0) = (exclude,include) boundary calc along FT direction
       |0000000|     |-------|
       |-------|  ,  |-------| Look at fft for implementation details
       |0000000|     |-------|
```

!!! note
    The index range of the incoming scalar
    field is assumed to begin at one.

# References
 - https://jakevdp.github.io/blog/2013/08/28/understanding-the-fft/
"""
module FFTs

using FFTW: fft
using ..Grids: orthogs, Vertex1D, Center1D, Coordinates
export solve!

abstract type FFTType end; export FFTType
struct Radix2Type <: FFTType end; export Radix2Type
struct FullType <: FFTType end; export FullType

include("dct.jl")
include("idct.jl")
include("wrappers3D.jl")

n_points_real(c::Coordinates, dl::Center1D) = c.c.s - 2
n_points_real(c::Coordinates, dl::Vertex1D) = c.n.s - 2

function solve!(u, f_in, grid, dir)
    FT = eltype(u)
    f = similar(f_in)
    f .= f_in
    i_orth = orthogs(dir)

    fv = @view f[2:end-1,2:end-1,:]
    dct3D!(fv, 3)
    fv = permutedims(fv, (2,1,3))
    dct3D!(fv, 3)
    fv = permutedims(fv, (2,1,3))
    f[2:end-1,2:end-1,:] .= fv

    Nx, Δx = n_points_real(grid.c[1], Center1D()), first(grid.c[1].c.Δh)
    Ny, Δy = n_points_real(grid.c[2], Center1D()), first(grid.c[2].c.Δh)

    @assert Δx ≈ Δy # Assumes (see poisson2 docs in solution) Δx == Δy

    ΔxΔy = Δx*Δy
    Δx⁻² = 1/Δx^2
    Δy⁻² = 1/Δy^2

    coeff = similar(u)
    for k in 1:size(u, 3)
        for j in 2:size(u, 2)-1
            for i in 2:size(u, 1)-1
                cosi = cos(π*FT(i-2)/FT(Nx))
                cosj = cos(π*FT(j-2)/FT(Ny))
                if !((i==2) && (j==2))
                    λ = cosi + cosj - 2
                    coeff[i,j,k] = 0.5*ΔxΔy/λ
                else
                    coeff[i,j,k] = 1
                end
            end
        end
    end

    u[2:end-1,2:end-1,:] .= f[2:end-1,2:end-1,:] .* coeff[2:end-1,2:end-1,:]
    u[2,2,:] .= 0 # "Pin down pressure" in Fourier space
    uv = @view u[2:end-1, 2:end-1,:]

    idct3D!(uv, 3)
    uv = permutedims(uv, (2,1,3))
    idct3D!(uv, 3)
    uv = permutedims(uv, (2,1,3))
    u[2:end-1, 2:end-1,:] .= uv
end

end # module
