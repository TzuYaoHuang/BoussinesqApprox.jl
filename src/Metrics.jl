using StaticArrays
import WaterLily: nds
import LinearAlgebra:⋅

heat_flux(sim) = heat_flux(sim.flow,sim.body)
heat_flux(flow,body) = heat_flux(flow.θ,flow.κ,flow.Ψ,body,time(flow))
function heat_flux(θ,κ,dq,body,t=0)
    Tθ = eltype(θ); To = promote_type(Float64,Tθ)
    dq .= zero(Tu)
    @loop dq[I] .= -𝛁θ(I,θ) ⋅ nds(body,loc(0,I,Tu),t) over I ∈ inside(θ)
    sum(To,dq)*κ
end
@inline 𝛁θ(I::CartesianIndex{D},θ) where D = SVector{D}((θ[I+δ(i,I)]-θ[I-δ(i,I)])/2 for i∈1:D)