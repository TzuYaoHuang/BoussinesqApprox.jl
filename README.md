# BoussinesqApprox

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://TzuYaoHuang.github.io/BoussinesqApprox.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://TzuYaoHuang.github.io/BoussinesqApprox.jl/dev/)
[![Build Status](https://github.com/TzuYaoHuang/BoussinesqApprox.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/TzuYaoHuang/BoussinesqApprox.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/TzuYaoHuang/BoussinesqApprox.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/TzuYaoHuang/BoussinesqApprox.jl)

![Julia flow](assets/HotColdCylinder_temp.gif)

This repository extends the [WaterLily.jl](https://github.com/WaterLily-jl/WaterLily.jl) flow solver to simulate thermal flows under the [Boussinesq approximation](https://en.wikipedia.org/wiki/Boussinesq_approximation_(buoyancy)). Immersed solid boundaries are implemented using a method similar to the first-order [Boundary Data Immersion Method](https://doi.org/10.1016/j.jcp.2011.04.022).

### WaterLily.jl Simulations with BoussinesqApprox.jl

Enabling the thermal simulation requires only loading an extra package and few addition to the `Simulation` constructor!
```julia
using WaterLily, BoussinesqApprox
```
The simulation struct is identical to WaterLily's, with a few extra parameters: `flow_ctor=ThermalFlow` activates the thermal simulation, `θ0` sets the initial temperature field, `θb` sets the body temperature field, `κ` is the thermal diffusivity, and `α` is the thermal expansion coefficient.
```julia
ThermalSimulation(NN, (0,0), radius; U, Δt=0.05, ν, α, κ, g, perdir, θ0, θb, body)
```
You can step or visualize the simulation exactly as you would a standard WaterLily `Simulation`:
```julia
sim_step!(sim, t_end; remeasure::Bool)
```
See the `examples` folder for more examples of how to run simulations. The hot/cold cylinder simulation shown in the GIF above can be found in [examples/HotColdCylinder.jl](examples/HotColdCylinder.jl).