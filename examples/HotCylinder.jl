using WaterLily, BoussinesqApprox
using Plots

function hot_circle(N;Ra=10^6,ε=0.001,Pr=10,mem=Array,T=Float32)
    NN = (N,N)
    radius = N÷8
    center = NN.÷2
    body = AutoBody((x,t)->√sum(abs2, x .- center) - radius)

    ΔT = T(1)
    α = ε/ΔT |> T
    g = T(1)

    νκ = g*α*ΔT*radius^3/Ra
    ν = sqrt(νκ*Pr)
    κ = sqrt(νκ/Pr)

    Uf = sqrt(g*α*ΔT*radius)
    U = Uf/√Pr |> T

    Simulation(NN, (0,0), radius; U, Δt=0.05, ν, κ, flow_ctor=ThermalFlow, g=(i,x,t)-> i==2 ? -g : zero(g), perdir=(1,),θ0=nothing,θb=(x)->ΔT, body)
end

sim = hot_circle(64)

# duration=8,step=0.08 covers plume formation through the steady rising plume reached around tU/L≈5
sim_gif!(sim; duration=16, step=0.08, field=sim->sim.flow.θ, clims=(0,1), cfill=:thermal,
              plotbody=true, fname=joinpath(@__DIR__,"HotCylinder.gif"))