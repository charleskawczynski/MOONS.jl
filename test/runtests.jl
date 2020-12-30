using Test

for submodule in [
                  "Grids",
                 ]

  println("Testing $submodule")
  t = @elapsed include(joinpath(submodule,"runtests.jl"))
  println("Completed tests for $submodule, $(round(Int, t)) seconds elapsed")
end
