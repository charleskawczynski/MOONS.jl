using ForwardDiff
using RootSolvers

function T2_h_i(β, α::FT, i::IT, N::IT) where {FT, IT}
  a = α; b = β
  g = (b+FT(1))/(b-FT(1))
  t = (FT(i)/FT(N) - a)/(FT(1)-a)
  f = ((b+FT(2)*a)*g^t - b + FT(2)*a)/((FT(2)*a+FT(1))*(FT(1)+g^t))
  return f
end

T2_root_Δh_near_hmax(β, hmin::FT,hmax::FT,α::FT,N::IT,Δh::FT) where {FT,IT} =
  T2_h_i(β,α,N,N) - T2_h_i(β,α,N-1,N) - Δh/(hmax-hmin)

T2_root_Δh_near_hmin(β, hmin::FT,hmax::FT,α::FT,N::IT,Δh::FT) where {FT,IT} =
  T2_h_i(β,α,1,N) - T2_h_i(β,α,0,N) - Δh/(hmax-hmin)

"""
    β_Δh_small

Compute `β` such that `N` grid points
between `hmin` and `hmax` stretch so
that the smallest `Δh` matches the one
given.

!!! warn
    If `N` is small, `β->∞`, and the
    given `Δh` may not be satisfied.
    `N` can be increased to fix this.
"""
function β_Δh_small(hmin::FT, hmax::FT, N::IT, Δh::FT) where {FT, IT}
  tol = eps(FT)
  β_initial = FT(1) + tol # Initial guess
  root_closure(β) = T2_root_Δh_near_hmax(β, hmin, hmax, FT(0), N, Δh)
  args = (NewtonsMethodAD(β_initial), CompactSolution(), SolutionTolerance(eps(FT)))
  sol = find_zero(β -> root_closure(β), args...)
  return sol.root
end

"""
    β_Δh_big

Compute `β` such that `N` grid points
between `hmin` and `hmax` stretch so
that the largest `Δh` matches the one
given.

!!! warn
    If `N` is small, `β->∞`, and the
    given `Δh` may not be satisfied.
    `N` can be increased to fix this.
"""
function β_Δh_big(hmin::FT, hmax::FT, N::IT, Δh::FT) where {FT, IT}
  tol = eps(FT)
  β_initial = FT(1) + tol # Initial guess
  root_closure(β) = T2_root_Δh_near_hmin(β, hmin, hmax, FT(0), N, Δh)
  args = (NewtonsMethodAD(β_initial), CompactSolution(), SolutionTolerance(eps(FT)))
  sol = find_zero(β -> root_closure(β), args...)
  return sol.root
end

"""
    β_Δh_both

Compute `β` such that `N` grid points
between `hmin` and `hmax` stretch so
that the smallest `Δh` matches the one
given.

!!! warn
    If `N` is small, `β->∞`, and the
    given `Δh` may not be satisfied.
    `N` can be increased to fix this.
"""
function β_Δh_both(hmin::FT, hmax::FT, N::IT, Δh::FT) where {FT, IT}
  tol = eps(FT)
  β_initial = FT(1) + tol # Initial guess
  root_closure(β) = T2_root_Δh_near_hmax(β, hmin, hmax, FT(1/2), N, Δh)
  args = (NewtonsMethodAD(β_initial), CompactSolution(), SolutionTolerance(eps(FT)))
  sol = find_zero(β -> root_closure(β), args...)
  return sol.root
end
