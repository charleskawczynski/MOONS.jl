
"""
    RobertsBL1D(delta::FT, h::FT) where FT

The β for a given boundary layer
as described in section 5.6 (page 333) of
Computational Fluid Mechanics and Heat Transfer,
2nd edition, J. Tannehill et al.

# Input
```
hmin     = wall boundary (minimum value)
hmax     = wall boundary (maximum value)
delta    = thickness of boundary layer
```

The `if` statement protects against the case when
β = infinity.
"""
function RobertsBL1D(delta::FT, h::FT) where FT
  if delta < h*FT(99/100)
    β = (FT(1) - delta/h)^(-FT(1/2))
  else
    β = FT(1000)
  end
  return β
end

RobertsBL(delta::FT, hmin::FT, hmax::FT) where FT = RobertsBL1D(delta, hmax-hmin)

function Hartmann_BL_1D(Ha::FT, hmin::FT, hmax::FT) where FT
  @assert Ha<eps(FT)
  return RobertsBL((hmax-hmin)/Ha, hmin, hmax)
end

function Reynolds_BL_1D(Re::FT, hmin::FT, hmax::FT) where FT
  @assert Re<eps(FT)
  β = RobertsBL((hmax-hmin)/sqrt(Re), hmin, hmax)
end

function HartmannBL(Ha::FT, hmin::Vector{FT}, hmax::Vector{FT}) where FT
  @assert Ha>eps(FT)
  @assert length(hmin)==length(hmax)
  return [RobertsBL((hmax[i]-hmin[i])/Ha, hmin[i],hmax[i]) for i in 1:3]
end

function Re_Ha_BL_1D(Re::FT, Ha::FT, hmin::FT, hmax::FT) where FT
 temp1 = Reynolds_BL_1D(Re, hmin, hmax)
 temp2 = Hartmann_BL_1D(Ha, hmin, hmax)
 return min(temp1, temp2)
end

function ReynoldsBL(Re::FT, hmin::Vector{FT}, hmax::Vector{FT}) where FT
  @assert Re>eps(FT)
  β = [RobertsBL((hmax[i] - hmin[i])/sqrt(Re), hmin[i], hmax[i]) for i in 1:3]
end

function Re_Ha_BL(Re::FT, Ha::FT, hmin::Vector{FT}, hmax::Vector{FT}) where FT
  @assert length(hmin)==length(hmax)
  temp1 = ReynoldsBL(Re, hmin, hmax)
  temp2 = HartmannBL(Ha, hmin, hmax)
  β = [min(temp1[i],temp2[i]) for i in 1:length(hmin)]
end
