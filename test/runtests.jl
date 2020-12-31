using Test

for submodule in [
                  "Grids",
                  "Fields",
                 ]

  println("Testing $submodule")
  t = @elapsed include(joinpath(submodule,"runtests.jl"))
  println("Completed tests for $submodule, $(round(Int, t)) seconds elapsed")
end
