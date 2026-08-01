module BoussinesqApprox

using WaterLily
import WaterLily: δ, @loop, CIj, slice, inside_u,BC!

include("Thermal.jl")
export ThermalFlow

include("util.jl")

end
