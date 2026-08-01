
import WaterLily: ∂,ϕ,ϕu,ϕuL,ϕuR,ϕuP,mom_predict!,mom_correct!,conv_diff!,accelerate!, measure!

struct ThermalFlow{D, T, Sf<:AbstractArray{T}, Vf<:AbstractArray{T}, Tf<:AbstractArray{T}, Lf} <: AbstractFlow{D,T}
    flow :: Flow{D,T}
    # Thermal fields
    θ :: Sf   # Temprature field  (θ = T-T₀)
    θ⁰:: Sf   # previous temprature
    Ψ :: Sf   # force vector

    # Body fields
    Λ :: Sf   # Body temperature field
    ξ₀:: Sf   # zeroth moment on the collocated grid

    # Properties
    α :: T    # coefficient of thermal expansion
    κ :: T    # temperature diffusivity
    function ThermalFlow(N::NTuple{D}, uBC; θ0=nothing, θb=nothing, α=0.001, κ=0.1, kwargs...) where D
        flow = Flow(N,uBC; kwargs...)
        θ = zero(flow.σ); θ⁰= zero(flow.σ); Λ = zero(flow.σ)
        ξ₀ = zero(flow.σ); fill!(ξ₀,1)
        if isa(θ0, Function)
            apply!(θ0, θ)
            BC!(θ,flow.perdir)
        end
        if isa(θb, Function)
            apply!(θb, Λ)
            BC!(Λ,flow.perdir)
        end
        Ψ = zero(flow.σ)

        new{D,eltype(flow.p),typeof(flow.p),typeof(flow.u),typeof(flow.μ₁),typeof(flow.λ)}(flow,θ,θ⁰,Λ,ξ₀,Ψ,α,κ)
    end
end
Base.getproperty(f::ThermalFlow, s::Symbol) = s in propertynames(f) ? getfield(f, s) : getfield(f.flow, s)
Base.setproperty!(f::ThermalFlow, s::Symbol, x) = s in propertynames(f) ? setproperty!(f,s,x) : setproperty!(f.flow,s,x)

# ∂ₜθ = -∂ⱼuⱼθ + κ∂ⱼⱼθ
function conv_diff!(r::AbstractArray{T,n},θ,u,Φ,λ::F; κ=0.1,perdir=()) where {T,n,F}
    fill!(r,0)
    N = size(r)
    for j∈1:n
        # if it is periodic direction
        tagper = (j in perdir)
        # treatment for bottom boundary with BCs
        lowerBoundary!(r,θ,u,Φ,κ,j,N,λ,Val{tagper}())
        # inner cells
        @loop (Φ[I] = ϕu(j,I,θ,u[I,j],λ) - κ*∂(j,I,θ);
               r[I] += Φ[I]) over I ∈ inside_u(N,j)
        @loop r[I-δ(j,I)] -= Φ[I] over I ∈ inside_u(N,j)
        # treatment for upper boundary with BCs
        upperBoundary!(r,θ,u,Φ,κ,j,N,λ,Val{tagper}())
    end
end

# Neumann BC Building block
lowerBoundary!(r,θ,u,Φ,κ,j,N,λ,::Val{false}) = @loop r[I] += ϕuL(j,I,θ,u[I,j],λ) - κ*∂(j,I,θ) over I ∈ slice(N,2,j,2)
upperBoundary!(r,θ,u,Φ,κ,j,N,λ,::Val{false}) = @loop r[I-δ(j,I)] += -ϕuR(j,I,θ,u[I,j],λ) + κ*∂(j,I,θ) over I ∈ slice(N,N[j],j,2)

# Periodic BC Building block
lowerBoundary!(r,θ,u,Φ,κ,j,N,λ,::Val{true}) = @loop (
    Φ[I] = ϕuP(j,CIj(j,I,N[j]-2),I,θ,u[I,j],λ) -κ*∂(j,I,θ); r[I] += Φ[I]) over I ∈ slice(N,2,j,2)
upperBoundary!(r,θ,u,Φ,κ,j,N,λ,::Val{true}) = @loop r[I-δ(j,I)] -= Φ[CIj(j,I,2)] over I ∈ slice(N,N[j],j,2)


accelerate!(r,t,::Nothing,::Union{Nothing,Tuple},θ,α) = nothing
accelerate!(r,t,f::Function,θ,α) = @loop r[Ii] += f(last(Ii),loc(Ii,eltype(r)),t)*(1-α*ϕ(last(Ii),CI(Base.front(Ii)),θ)) over Ii ∈ CartesianIndices(r)
accelerate!(r,t,g::Function,::Union{Nothing,Tuple},θ,α) = accelerate!(r,t,g,θ,α)
accelerate!(r,t,::Nothing,U::Function,θ,α) = accelerate!(r,t,(i,x,t)->derivative(τ->U(i,x,τ),t),θ,α)
accelerate!(r,t,g::Function,U::Function,θ,α) = accelerate!(r,t,(i,x,t)->g(i,x,t)+derivative(τ->U(i,x,τ),t),θ,α)

function BDIMΘ!(a::ThermalFlow{T},w=1) where T # include 0.5
    wT = T(w)
    dt = a.Δt[end]
    @loop a.Ψ[I] = a.θ⁰[I] + dt*a.Ψ[I] - a.Λ[I] over I in CartesianIndices(a.Ψ)
    @loop a.θ[I] =(a.θ[I]  + a.Λ[I]  + a.ξ₀[I]*a.Ψ[I])*wT over I in CartesianIndices(a.Ψ)
end

function mom_predict!(a::ThermalFlow, t₀, t₁; udf=nothing, kwargs...)
    a.θ⁰ .= a.θ; fill!(a.θ, 0)
    conv_diff!(a.f,a.u⁰,a.σ,a.λ;ν=a.ν,perdir=a.perdir)
    udf!(a,udf,a.u⁰,t₀; kwargs...) # advect with u⁰ (a.u is zeroed by scale_u!)
    # Thermal flow
    accelerate!(a.f,t₀,a.g,a.uBC,a.θ⁰,a.α)
    conv_diff!(a.Ψ,a.θ⁰,a.u⁰,a.σ,a.λ;κ=a.κ,perdir=a.perdir); BDIMΘ!(a); BC!(a.θ,a.perdir)

    BDIM!(a); BC!(a.u,a.uBC,a.exitBC,a.perdir,t₁) # BC MUST be at t₁
    a.exitBC && exitBC!(a.u,a.u⁰,a.Δt[end]) # convective exit
end

function mom_correct!(a::ThermalFlow, t; udf=nothing, kwargs...)
    conv_diff!(a.f,a.u,a.σ,a.λ;ν=a.ν,perdir=a.perdir)
    udf!(a,udf,a.u,t; kwargs...) # advect with projected a.u
    # Thermal flow
    accelerate!(a.f,t,a.g,a.uBC,a.θ,a.α)
    conv_diff!(a.Ψ,a.θ,a.u,a.σ,a.λ;κ=a.κ,perdir=a.perdir); BDIMΘ!(a); BC!(a.θ,a.perdir)

    BDIM!(a); scale_u!(a,0.5); BC!(a.u,a.uBC,a.exitBC,a.perdir,t)
end

function CFL(a::ThermalFlow;Δt_max=10)
    @inside a.σ[I] = flux_out(I,a.u)
    min(Δt_max,inv(maximum(a.σ)+5*(a.ν+a.κ)))
end

using StaticArrays
function measure!(a::ThermalFlow{N,T},body::AbstractBody;t=zero(T),ϵ=1) where {N,T}
    a.V .= zero(T); a.μ₀ .= one(T); a.μ₁ .= zero(T); d²=T(2+ϵ)^2
    measure_sdf!(a.σ, body, t; fastd²=d²) # measure separately to allow specialization
    @fastmath @inline function fill!(μ₀,μ₁,ξ₀,V,d,I)
        if d[I]^2<d²
            for i ∈ 1:N
                dᵢ,nᵢ,Vᵢ = measure(body,loc(i,I,T),t,fastd²=d²)
                dᵢ = abs(dᵢ) ≤ 0.5 ? dᵢ : copysign(dᵢ,d[I]) # enforce sign consistency
                V[I,i] = Vᵢ[i]
                μ₀[I,i] = WaterLily.μ₀(dᵢ,ϵ)
                for j ∈ 1:N
                    μ₁[I,i,j] = WaterLily.μ₁(dᵢ,ϵ)*nᵢ[j]
                end
            end
            ξ₀[I] = WaterLily.μ₀(d[I],ϵ)
        elseif d[I]<zero(T)
            for i ∈ 1:N
                μ₀[I,i] = zero(T)
            end
            ξ₀[I] = zero(T)
        end
    end
    @loop fill!(a.μ₀,a.μ₁,a.ξ₀,a.V,a.σ,I) over I ∈ inside(a.p)
    BC!(a.μ₀,zeros(SVector{N,T}),false,a.perdir) # BC on μ₀, don't fill normal component yet
    BC!(a.ξ₀,a.perdir)
    BC!(a.V ,zeros(SVector{N,T}),a.exitBC,a.perdir)
end