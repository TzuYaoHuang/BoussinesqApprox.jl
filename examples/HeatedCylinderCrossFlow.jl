using WaterLily, BoussinesqApprox
using Statistics: mean
using Plots

# Isothermal cylinder held at ΔT above the ambient, in a uniform cross flow.
# No buoyancy (g=nothing) -- this is a forced-convection validation case.
function heated_circle(N=64;Re=100,Pr=0.7,U=1,mem=Array,T=Float32)
    NN = (2N,N)
    radius = N/20 |> T
    center = N÷2
    body = AutoBody((x,t)->√sum(abs2, x .- center) - radius)
    oneT = one(T)

    ν = U*radius/Re |> T
    κ = ν/Pr |> T

    ThermalSimulation(NN, (U,0), radius; U, ν, κ, θb=(x,t)->oneT, body, mem, T)
end

# Churchill-Bernstein correlation for forced convection across a cylinder (valid for Re*Pr>0.2)
Nu_churchill_bernstein(Re,Pr) = 0.3 + 0.62*√Re*Pr^(1/3)/(1+(0.4/Pr)^(2/3))^(1/4)*(1+(Re/282000)^(5/8))^(4/5)

# Advance the sim in `step`-sized chunks up to `dur`, recording the total heat flux and
# animating the temperature field; discard the initial transient before averaging
function mean_Nu(sim;dur=100,step=0.2,transient=50,video=joinpath(@__DIR__,"HeatedCylinderCrossFlow_temp.gif"))
    t₀ = sim_time(sim)
    ts = Float64[]; Nu = Float64[]
    anim = @animate for tᵢ in range(t₀,t₀+dur;step)
        sim_step!(sim,tᵢ;remeasure=false)
        push!(ts, tᵢ)
        push!(Nu, heat_flux(sim)/(π*sim.flow.κ)) # ΔT=1, D=2radius, Nu=Q/(π*κ*ΔT)
        flood(sim.flow.θ[inside(sim.flow.p)]; clims=(-0.01,1), cfill=:thermal, hidedecorations=true)
        body_plot!(sim)
        println("tU/L=",round(tᵢ,digits=3)," Nu=",round(Nu[end],digits=3))
    end
    gif(anim,video)
    keep = findall(>=(transient), ts)
    ts, Nu, mean(Nu[keep])
end

using CUDA
Re,Pr = 500,0.71
sim = heated_circle(512;Re,Pr,mem=CuArray)

ts, Nu, Nu_sim = mean_Nu(sim)
Nu_ref = Nu_churchill_bernstein(Re,Pr)
@info "Nu (simulation, time-averaged) = $(round(Nu_sim,digits=2))"
@info "Nu (Churchill-Bernstein)       = $(round(Nu_ref,digits=2))"
@info "relative error                 = $(round(100*(Nu_sim-Nu_ref)/Nu_ref,digits=1))%"

plot(ts,Nu,label="simulation",xlabel="tU/L",ylabel="Nu")
hline!([Nu_ref],label="Churchill-Bernstein",ls=:dash)
savefig(joinpath(@__DIR__,"HeatedCylinderCrossFlow_Nu.png"))
