
import WaterLily: ∂,ϕu,ϕuL,ϕuR,ϕuP

struct Thermal{D, T, Sf<:AbstractArray{T}, Lf}
    # Thermal fields
    Θ :: Sf   # Temprature field  (θ = T-T₀)
    θ⁰:: Sf   # previous temprature
    f :: Sf   # force vector    

    # Properties
    α :: T    # coefficient of thermal expansion
    κ :: T    # temperature diffusivity

    perdir :: NTuple # tuple of periodic direction
    λ :: Lf # convective scheme λ(u,c,d); a type parameter so conv_diff!'s kernel specializes on it
    function Thermal()
    end
end

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