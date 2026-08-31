
import WaterLily: ∂,ϕ,ϕu,ϕuL,ϕuR,ϕuP,mom_predict!,mom_correct!,conv_diff!,accelerate!, measure!,udf!,
                  BDIM!,scale_u!,exitBC!,flux_out

using StaticArrays

"""
    check_sfn(f, N, T, nargs)

Check a scalar field function `f` has `nargs` arguments (`f(x)` or `f(x,t)`) and returns a `T`.
Scalar counterpart to `WaterLily.check_fn`; `apply!` won't catch a wrong arity itself.
"""
check_sfn(f,N,T,nargs) = nothing # fallback: `nothing`, arrays, numbers are all fine
function check_sfn(f::Function,N,T,nargs)
    @assert first(methods(f)).nargs==nargs+1 "$f signature needs $nargs argument(s)" # +1 for the function object itself
    @assert typeof(f(sargs(Val{nargs}(),N,T)...))==T "$f is not type stable" # catches a Float64 literal in a Float32 field
end
sargs(::Val{1},N,T) = (zeros(SVector{N,T}),) # f(x)
sargs(::Val{2},N,T) = (zeros(SVector{N,T}),zero(T)) # f(x,t)

struct ThermalFlow{D,T,Sf<:AbstractArray{T},Vf<:AbstractArray{T},F<:Flow{D,T}} <: AbstractFlow{D,T}
    flow :: F
    # Thermal fields
    θ :: Sf   # Temprature field  (θ = T-T₀)
    θ⁰:: Sf   # previous temprature
    Ψ :: Sf   # force vector

    # Body fields
    Λ :: Sf   # Body temperature field
    ξ₀:: Sf   # zeroth moment on the collocated grid
    ξ₁:: Vf   # first  moment on the staggered  grid

    # Properties
    α :: T    # coefficient of thermal expansion
    κ :: T    # temperature diffusivity
    function ThermalFlow(N::NTuple{D}, uBC; θ0=nothing, θb=nothing, α=0.001, κ=0.1, kwargs...) where D
        flow = Flow(N,uBC; kwargs...)
        T = eltype(flow.p)
        # check function type stability
        check_sfn(θ0,D,T,1)          # θ0(x)
        check_sfn(θb,D,T,2)          # θb(x,t)
        θ = zero(flow.σ); θ⁰= zero(flow.σ); Λ = zero(flow.σ)
        ξ₀ = zero(flow.σ); fill!(ξ₀,1)
        ξ₁ = zero(flow.u); fill!(ξ₁,0)
        zeroT=zero(T)
        isa(θ0,Function) && (apply!(θ0,θ); BC!(θ,flow.perdir))
        isa(θb,Function) && (apply!((x)->θb(x,zeroT),Λ); BC!(Λ,flow.perdir))
        Ψ = zero(flow.σ)

        new{D,T,typeof(flow.p),typeof(flow.u),typeof(flow)}(flow,θ,θ⁰,Ψ,Λ,ξ₀,ξ₁,α,κ)
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

function BDIMΘ!(a::ThermalFlow{n,T},w=1) where {n,T} # include 0.5
    wT = T(w)
    dt = a.Δt[end]
    @loop a.Ψ[I] = a.θ⁰[I] + dt*a.Ψ[I] - a.Λ[I] over I in CartesianIndices(a.Ψ)
    @loop a.θ[I] += ξddn(I,a.ξ₁,a.Ψ) + a.Λ[I] + a.ξ₀[I]*a.Ψ[I] over I in inside(a.Ψ)
    a.θ .*= wT
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
    conv_diff!(a.Ψ,a.θ,a.u,a.σ,a.λ;κ=a.κ,perdir=a.perdir); BDIMΘ!(a,0.5); BC!(a.θ,a.perdir)

    BDIM!(a); scale_u!(a,0.5); BC!(a.u,a.uBC,a.exitBC,a.perdir,t)
end

function CFL(a::ThermalFlow;Δt_max=10)
    @inside a.σ[I] = flux_out(I,a.u)
    min(Δt_max,inv(maximum(a.σ)+5*(a.ν+a.κ)))
end

function measure!(a::ThermalFlow{N,T},body::AbstractBody;t=zero(T),ϵ=1) where {N,T}
    @invoke measure!(a::AbstractFlow,body;t,ϵ)
    a.ξ₀ .= one(T); a.ξ₁ .= zero(T); d²=T(2+ϵ)^2
    @fastmath @inline function fill!(ξ₀,ξ₁,dd,I)
        d = dd[I]
        if d^2<d²
            ξ₀[I] = WaterLily.μ₀(d,ϵ)
            _,nᵢ,_ = measure(body,loc(0,I,T),t,fastd²=d²)
            for i ∈ 1:N
                ξ₁[I,i] = WaterLily.μ₁(d,ϵ)*nᵢ[i]
            end
        elseif d<zero(T)
            ξ₀[I] = zero(T)
        end
    end
    @loop fill!(a.ξ₀,a.ξ₁,a.σ,I) over I ∈ inside(a.p)
    BC!(a.ξ₀,a.perdir)
end

@fastmath @inline function ξddn(I::CartesianIndex{n},μ,f) where n
    s = zero(eltype(f))
    for j ∈ 1:n
        s+= @inbounds μ[I,j]*(f[I+δ(j,I)]-f[I-δ(j,I)])
    end
    return s/2
end