using WaterLily, BoussinesqApprox
using Plots

function hot_circle(N;Ra=10^6,ε=0.001,Pr=10,mem=Array,T=Float32)
    NN = (2N,N)
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

    Simulation(NN, (0,0), radius; U, Δt=0.05, ν, κ, flow_ctor=ThermalFlow, g=(i,x,t)->g, perdir=(1,),θ0=nothing,θb=(x)->ΔT, body)
end

sim = hot_circle(64)

sim_step!(sim)

heatmap(sim.flow.θ',aspect_ratio=:equal)