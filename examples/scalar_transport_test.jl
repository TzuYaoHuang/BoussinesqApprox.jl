#=
Validate the scalar transport/diffusion kernel `conv_diff!` (added to WaterLily's
generic `conv_diff!` in src/Thermal.jl) against the analytical solution of the
1D advection-diffusion equation

    ∂ₜθ + U∂ₓθ = κ∂ₓₓθ,   θ(x,0) = sin(2πk x/Lx)

    ⟹  θ(x,t) = sin(2πk(x-Ut)/Lx) · exp(-κ(2πk/Lx)² t)

`conv_diff!` itself only knows about a D-dimensional grid, so the 1D problem is
embedded in a thin, doubly-periodic 2D domain: the initial condition is uniform
along y, so no y-gradients ever develop and the test isolates the x-direction
discretization.

Run with:
    julia --project=examples examples/scalar_transport_test.jl
=#
using WaterLily, BoussinesqApprox
using Printf, Test, Plots

# The scalar-transport method lives on the same generic function as WaterLily's
# momentum conv_diff! (Thermal.jl extends it via `import`), and isn't exported.
const conv_diff! = WaterLily.conv_diff!
# BoussinesqApprox's scalar BC! (util.jl) is distinct from -- and not exported
# alongside -- WaterLily's vector-field BC!, so it must be called qualified.
const scalarBC! = BoussinesqApprox.BC!

"""
    setup(Nx, Ny, U, k)

Build ghost-padded scalar (`θ`,`r`,`Φ`) and velocity (`u`) arrays on an
`Nx×Ny` periodic grid with a `sin(2πk x/Nx)` initial condition and uniform
velocity `(U,0)`.
"""
function setup(Nx, Ny, U, k)
    N = (Nx+2, Ny+2)
    θ = zeros(Float64, N); r = zero(θ); Φ = zero(θ)
    u = zeros(Float64, N..., 2); u[:,:,1] .= U
    K = 2π*k/Nx
    apply!(x -> sin(K*x[1]), θ)
    scalarBC!(θ, (1,2))
    θ, r, Φ, u, K
end

"""
    advect_diffuse!(θ, r, u, Φ, λ, tend; κ, cfl=0.2)

Explicit-Euler time integration of `∂ₜθ = -∂ⱼ(uⱼθ) + κ∂ⱼⱼθ` up to time `tend`,
using `conv_diff!` for the spatial operator. Step size is set from the
diffusive and advective stability limits (grid spacing is 1 in WaterLily units).
"""
function advect_diffuse!(θ, r, u, Φ, λ, tend; κ, cfl=0.2)
    U = maximum(abs, @view u[:,:,1])
    dt = min(κ > 0 ? cfl/κ : Inf, U > 0 ? cfl/U : Inf)
    nsteps = max(ceil(Int, tend/dt), 1)
    dt = tend/nsteps
    for _ ∈ 1:nsteps
        conv_diff!(r, θ, u, Φ, λ; κ, perdir=(1,2))
        θ .+= dt .* r
        scalarBC!(θ, (1,2))
    end
    nsteps
end

analytical(x, t, K, U, κ) = sin(K*(x-U*t))*exp(-κ*K^2*t)

"""
    xslice(θ, Nx)

Return the interior-cell x-locations and the corresponding θ profile along an
interior y-row (any row works since the field is uniform in y).
"""
function xslice(θ, Nx)
    xs = [loc(0,CartesianIndex(i,3))[1] for i ∈ 2:Nx+1]
    xs, θ[2:Nx+1, 3]
end

l2(e)  = sqrt(sum(abs2, e)/length(e))
linf(e) = maximum(abs, e)

@testset "conv_diff! vs analytical advection-diffusion" begin

    # -----------------------------------------------------------------------
    # Test 1: pure diffusion (U=0) -- grid convergence of the diffusive stencil.
    # With no advection, ϕu(...)=u*λ(...)=0 regardless of λ, so this isolates
    # the κ∂ⱼⱼθ discretization.
    # -----------------------------------------------------------------------
    println("Test 1: pure diffusion, grid convergence")
    @printf("%6s  %12s  %12s\n", "Nx", "L2 err", "L∞ err")
    κ, k, tend = 0.02, 1, 3.0
    errs = Float64[]
    for Nx ∈ (32,64,128,256)
        θ,r,Φ,u,K = setup(Nx, 4, 0.0, k)
        advect_diffuse!(θ,r,u,Φ,WaterLily.cds,tend; κ)
        xs, θnum = xslice(θ, Nx)
        e = θnum .- analytical.(xs, tend, K, 0.0, κ)
        push!(errs, l2(e))
        @printf("%6d  %12.3e  %12.3e\n", Nx, l2(e), linf(e))
    end
    rates = log2.(errs[1:end-1]./errs[2:end])
    println("Observed L2 convergence rates: ", round.(rates,digits=2))
    @test errs[end] < 1e-3
    @test all(rates .> 1.8) # finite-volume Laplacian stencil is at least 2nd-order accurate

    # -----------------------------------------------------------------------
    # Test 2: advection-diffusion of a translating, decaying sine wave (U≠0),
    # using the central scheme `cds` so the linear PDE is matched by a linear
    # discretization (no TVD-limiter clipping near extrema).
    # -----------------------------------------------------------------------
    println("\nTest 2: advection-diffusion (U=1, cds scheme)")
    Nx, Ny, U, κ2, k2, tend2 = 128, 4, 1.0, 0.02, 1, 2.0
    θ,r,Φ,u,K = setup(Nx, Ny, U, k2)
    nsteps = advect_diffuse!(θ,r,u,Φ,WaterLily.cds,tend2; κ=κ2)
    xs, θnum = xslice(θ, Nx)
    θex = analytical.(xs, tend2, K, U, κ2)
    e = θnum .- θex
    @printf("nsteps=%d  L2 err=%.3e  L∞ err=%.3e\n", nsteps, l2(e), linf(e))
    @test l2(e) < 5e-3

    # -----------------------------------------------------------------------
    # Test 3: same case with the default flux-limited scheme (`quick`), as used
    # by ThermalFlow in practice. The nonlinear limiter clips to first-order
    # near the sine's extrema, so accuracy is lower than Test 2 -- checked with
    # a looser tolerance rather than asserted to machine-level accuracy.
    # -----------------------------------------------------------------------
    println("\nTest 3: advection-diffusion (U=1, default quick scheme)")
    θq,rq,Φq,uq,Kq = setup(Nx, Ny, U, k2)
    advect_diffuse!(θq,rq,uq,Φq,WaterLily.quick,tend2; κ=κ2)
    _, θnum_q = xslice(θq, Nx)
    eq = θnum_q .- θex
    @printf("L2 err=%.3e  L∞ err=%.3e\n", l2(eq), linf(eq))
    @test l2(eq) < 5e-2

    global xs, θex, θnum, θnum_q, tend2 # used below for the plot
end

# ---------------------------------------------------------------------------
# Visual comparison for the advection-diffusion case (cds vs quick vs exact)
# ---------------------------------------------------------------------------
plt = plot(xs, θex; label="analytical", lw=2, color=:black,
           xlabel="x", ylabel="θ", title="Advection-diffusion of θ(x,0)=sin(2πx/Lx), t=$tend2")
plot!(plt, xs, θnum; label="conv_diff! (cds)", ls=:dash, lw=2)
plot!(plt, xs, θnum_q; label="conv_diff! (quick)", ls=:dot, lw=2)
savefig(plt, joinpath(@__DIR__, "scalar_transport_test.png"))
println("\nSaved comparison plot to examples/scalar_transport_test.png")
