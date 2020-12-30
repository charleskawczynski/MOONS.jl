
export uniform,
       linspace,
       uniform_left,
       uniform_right,
       Roberts_left,
       Roberts_right,
       Roberts_both,
       cluster

function uniform_direction(hstart::FT, Δh::FT, N::IT, dir::IT) where {FT, IT}
  @assert (N>0)
  hn = [hstart+FT(dir)*FT(i-1)*Δh for i in 1:N+1]
  hn = dir>0 ? hn : reverse(hn)
  return hn
end

function assert_monotonically_increasing!(h::VFT, hmin::FT, hmax::FT, β::FT, caller::AbstractString) where {FT,VFT<:AbstractVector{FT}}
  s = length(h)
  for i in 1:s-2
    R = (h[i+2]-h[i+1])/(h[i+1]-h[i])
    if R < FT(1)
      println("Error: coordinates not monotonically increasing in ", caller)
      println("hmin, hmax, β = ", hmin, hmax, β)
      println("R = ", R)
      println("h = ", h)
      throw(MethodError("coordinates not monotonically increasing"))
    end
  end
end

function insist_monotonically_decreasing!(h::VFT, hmin::FT, hmax::FT, β::FT, caller::AbstractString) where {FT,VFT<:AbstractVector{FT}}
  s = length(h)
  for i in 1:s-2
    R = (h[i+2]-h[i+1])/(h[i+1]-h[i])
    if R > FT(1)
     println("Error: coordinates not monotonically decreasing in ", caller)
     println("in coordinate_distribution_funcs.jl")
     println("hmin, hmax, β = ", hmin, hmax, β)
     println("R = ", R)
     println("h = ", h)
     throw(MethodError("coordinates not monotonically increasing"))
    end
  end
end

"""
    transformation1(hmin::FT, hmax::FT, N::FT, β::FT) where FT

The coordinates and differences of a Robert's stretching
function as described in section 5.6 (page 333) of
Computational Fluid Mechanics and Heat Transfer,
2nd edition, J. Tannehill et al. (Transformation 1)

# Input
```
    hmin     = minimum value
    hmax     = maximum value
    N        = N segments of Δh
    β     = 1.0 <= β <= infinity = stretching factor
             = larger -> less stretching
  y=0                         y=h
                               |
   |-|--|---|-------|----------|
   |------> y
```

!!! note
    I have abused notation a bit to provide consistent notation
    with the reference as well as generalize the returned grid
    so that it need not start at y=0.
"""
function transformation1(hmin::FT, hmax::FT, N::IT, β::FT) where {FT, IT}
  a = FT(1); b = β
  g = (b+FT(1))/(b-FT(1))
  Δh = (hmax - hmin)/FT(N)
  hnbar = FT[hmin+FT(i-1)*Δh for i in 1:N+1]
  hnbar = (hnbar .- hmin)./(hmax-hmin)
  hn = [((b+a)-(b-a)*g^(a-hnbar[i]))/(g^(a-hnbar[i])+a) for i in 1:N+1]
  assert_monotonically_increasing!(hn, hmin, hmax, β, "transformation1")
  hn = hmin .+ (hmax - hmin).*hn
  return hn
end

"""
    transformation2

Coordinates and differences of a Robert's
stretching function as described in section 5.6 (page 333) of
Computational Fluid Mechanics and Heat Transfer,
2nd edition, J. Tannehill et al.

# Input
```
hmin     = minimum value
hmax     = maximum value
N        = N segments of Δh
α        = 0      stretching at y=h only
α        = 0.5    stretching at y=0 and y=hmax
β        = 1.0 <= β <= infinity = stretching factor
         = larger -> less stretching
```

Here is a picture illustration for α = 0:
```
                             y=h
                              |
  |----------|-------|---|--|-|
  |------> y
```

Note that this must be used in reverse for the lid driven
cavity geometry for the 'front' and 'back' walls.

!!! note
    I have abused notation a bit to provide consistent notation
    with the reference as well as generalize the returned grid
    so that it need not start at y=0.
"""
function transformation2(hmin::FT, hmax::FT, N::IT, α::FT, β::FT) where {FT, IT}
  a = α; b = β
  g = (b+FT(1))/(b-FT(1))
  Δh = (hmax - hmin)/FT(N)
  hnbar = FT[hmin+FT(i-1)*Δh for i in 1:N+1]
  hnbar = (hnbar .- hmin)./(hmax-hmin)
  hn = [ ((b+FT(2)*a)*g^((hnbar[i]-a)/(FT(1)-a))-
          b+FT(2)*a)/((FT(2)*a+FT(1))*(FT(1)+
          g^((hnbar[i]-a)/(FT(1)-a)))) for i in 1:N+1]
  if abs(α)<eps(FT)
   insist_monotonically_decreasing!(hn, hmin, hmax, β, "transformation2")
  end
  hn .= hmin .+ (hmax - hmin).*hn
  return hn
end

"""
    transformation3

This function returns the coordinates and differences of a Robert's
stretching function as described in section 5.6 (page 333) of
Computational Fluid Mechanics and Heat Transfer,
2nd edition, J. Tannehill et al. (Transformation 3)

# Input
```
hmin     = minimum value
hmax     = maximum value
N        = N segments of Δh
τ        = 0 <= τ <= infinity = stretching factor

τ        = 0            --> no stretching
τ        = large values --> strong stretching
```

Here is a picture illustration for α = 0:

```
                                 y=yc                        y=h
                                  |
      |----------|-------|---|--|-|-|--|---|-------|----------|
      |------> y
```

Note that this must be used in reverse for the lid driven
cavity geometry for the 'front' and 'back' walls.

!!! note
    I have abused notation a bit to provide consistent notation
    with the reference as well as generalize the returned grid
    so that it need not start at y=0.
"""
function transformation3(hmin::FT, hmax::FT, N::IT, yc::FT, τ::FT) where {FT, IT}
  a = FT(1); c = FT(2)
  Δh = (hmax - hmin)/FT(N)
  hnbar = [hmin+FT(i-1)*Δh for i in 1:N+1]
  hnbar = (hnbar - hmin)/(hmax-hmin)
  d = yc/(hmax-hmin); e = FT(exp(τ))
  B = a/(c*τ)*log((a+(e-a)*d)/(a+(a/e-a)*d))
  hn = [d*(a+FT(sinh(τ*(hnbar(i)-B)))/FT(sinh(τ*B))) for i in 1:N+1]
  hn .= hmin + (hmax - hmin).*hn
  return hn
end

uniform_left(hstart::FT, Δh::FT, N::IT) where {FT, IT} = uniform_direction(hstart, Δh, N,  IT(-1))
uniform_right(hstart::FT, Δh::FT, N::IT) where {FT, IT} = uniform_direction(hstart, Δh, N, IT( 1))

uniform(hmin::FT, hmax::FT, N::IT) where {FT, IT} = [hmin+FT(i-1)*(hmax - hmin)/FT(N) for i in 1:N+1]
linspace(hmin::FT, hmax::FT, N::IT) where {FT, IT} = [hmin+FT(i-1)*(hmax - hmin)/FT(N) for i in 1:N+1]
Roberts_left(hmin::FT, hmax::FT, N::IT, β::FT) where {FT, IT} = transformation1(hmin, hmax, N, β)
Roberts_right(hmin::FT, hmax::FT, N::IT, β::FT) where {FT, IT} = transformation2(hmin, hmax, N, FT(0), β)
Roberts_both(hmin::FT, hmax::FT, N::IT, β::FT) where {FT, IT} = transformation2(hmin, hmax, N, FT(1/2), β)
cluster(hmin::FT, hmax::FT, N::IT, yc::FT, tau::FT) where {FT, IT} = transformation3(hmin, hmax, N, yc, tau)

