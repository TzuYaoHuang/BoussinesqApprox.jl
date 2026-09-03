using StaticArrays
import WaterLily: nds
import LinearAlgebra:⋅

"""
    heat_flux(sim)

Total diffusive heat flux leaving the immersed body, `κ*∮-∇θ·n̂ dS`, computed
via the BDIM-masked surface normal `nds` (see `WaterLily.pressure_force`).
"""
heat_flux(sim) = heat_flux(sim.flow,sim.body)
heat_flux(flow,body) = heat_flux(flow.θ,flow.κ,flow.Ψ,body,WaterLily.time(flow))
function heat_flux(θ,κ,dq,body,t=0)
    Tθ = eltype(θ); To = promote_type(Float64,Tθ)
    dq .= zero(Tθ)
    @loop dq[I] = -𝛁θ(I,θ) ⋅ nds(body,loc(0,I,Tθ),t) over I ∈ inside(θ)
    sum(To,dq)*κ
end
@inline 𝛁θ(I::CartesianIndex{D},θ) where D = SVector{D}((θ[I+δ(i,I)]-θ[I-δ(i,I)])/2 for i∈1:D)