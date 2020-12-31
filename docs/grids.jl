using Plots
using MOONS.Grids
FT = Float64
a, b, n = FT(0.0), FT(1.0), 10
n_βs = 4
R = 1 .+ (10 .^ range(FT(-3), stop=FT(1), length=n_βs))
y_buff = 2

p = plot()
for (y,β) in enumerate(R)
  c = Coordinates(a,b,n; warpfun=Roberts_both, args=(β,))
  plot!(c.n.h, zeros(length(c.n.h)) .+ y,
    markershape=:circle, label="β=$β",
    ylims=(0,n_βs+y_buff), yaxis=nothing)
  plot!(layout=(n_βs,1), title="Roberts stretching (both)")
end
savefig("Roberts_both.svg")

p = plot()
for (y,β) in enumerate(R)
  c = Coordinates(a,b,n; warpfun=Roberts_left, args=(β,))
  plot!(c.n.h, zeros(length(c.n.h)) .+ y,
    markershape=:circle, label="β=$β",
    ylims=(0,n_βs+y_buff), yaxis=nothing)
  plot!(layout=(n_βs,1), title="Roberts stretching (left)")
end
savefig("Roberts_left.svg")

p = plot()
for (y,β) in enumerate(R)
  c = Coordinates(a,b,n; warpfun=Roberts_right, args=(β,))
  plot!(c.n.h, zeros(length(c.n.h)) .+ y,
    markershape=:circle, label="β=$β",
    ylims=(0,n_βs+y_buff), yaxis=nothing)
  plot!(layout=(n_βs,1), title="Roberts stretching (right)")
end
savefig("Roberts_right.svg")
