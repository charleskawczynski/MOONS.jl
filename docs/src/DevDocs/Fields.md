# Fields

```@meta
CurrentModule = MOONS
```

```@example
include(joinpath("..", "..", "fields.jl"))
nothing
```
![](wireframe.svg)
![](contourf.svg)

## Extrapolations

![](extrap_c2v.svg)
![](extrap_v2c.svg)

## Gradients

### f(x,y,z)

![](f.svg)

### ∇f(x,y,z)

| Gradient            | Exact                |  Approx          |
:--------------------:|:--------------------:|:-----------------:
∇f[1]                | ![](∇f_1_exact.svg)  |   ![](∇f_1.svg)
∇f[2]                | ![](∇f_2_exact.svg)  |   ![](∇f_2.svg)

