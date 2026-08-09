module BoussinesqApprox

using WaterLily
import WaterLily: δ, @loop, CI, CIj, slice, inside_u,BC!

include("ThermalFlow.jl")
export ThermalFlow

include("util.jl")

"""
    check_sfn(f, N, T, nargs)

Check a scalar field function `f` has `nargs` arguments (`f(x)` or `f(x,t)`) and returns a `T`.
Scalar counterpart to [`check_fn`](@ref); `apply!` won't catch a wrong arity itself.
"""
check_sfn(f,N,T,nargs) = nothing # fallback: `nothing`, arrays, numbers are all fine
function check_sfn(f::Function,N,T,nargs)
    @assert first(methods(f)).nargs==nargs+1 "$f signature needs $nargs argument(s)" # +1 for the function object itself
    @assert typeof(f(sargs(Val{nargs}(),N,T)...))==T "$f is not type stable" # catches a Float64 literal in a Float32 field
end
sargs(::Val{1},N,T) = (zeros(SVector{N,T}),) # f(x)
sargs(::Val{2},N,T) = (zeros(SVector{N,T}),zero(T)) # f(x,t)
function ThermalSimulation(args...; θ0=nothing, θb=nothing, α=0.001, κ=0.1, kwargs...)
    # check function type stability
    check_sfn(θ0,D,T,1)          # θ0(x)
    check_sfn(θb,D,T,1)          # θb(x)
    Simulation(args...; 
        flow_ctor=(flargs...;flkwargs...) -> ThermalFlow(
            flargs...; θ0, θb, α, κ, flkwargs...
        ), kwargs...
    )
end
export ThermalSimulation

end
