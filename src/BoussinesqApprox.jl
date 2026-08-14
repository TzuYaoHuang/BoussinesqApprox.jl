module BoussinesqApprox

using WaterLily
import WaterLily: δ, @loop, CI, CIj, slice, inside_u,BC!

include("ThermalFlow.jl")
export ThermalFlow

include("util.jl")

function ThermalSimulation(args...; θ0=nothing, θb=nothing, α=0.001, κ=0.1, kwargs...)
    Simulation(args...; 
        flow_ctor=(flargs...;flkwargs...) -> ThermalFlow(
            flargs...; θ0, θb, α, κ, flkwargs...
        ), kwargs...
    )
end
export ThermalSimulation

end
