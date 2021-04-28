# Poisson using FFT

```@meta
CurrentModule = MOONS
```

```@example
include(joinpath("..", "..", "poisson_fft.jl"))
nothing
```

| Function            | Exact                |  Approx               |
|--------------------:|:--------------------:|:----------------------|
`f`                   | ![](f_exact.svg)     |   ![](f_discrete.svg)
`u`                   | ![](u_exact.svg)     |   ![](u_discrete.svg)

![](u_err.svg)

