using WaterLily, BoussinesqApprox
using Plots

function hot_cold_circle(N;Ra=10^6,ε=0.001,Pr=100,mem=Array,T=Float32)
    NN = (N,N)
    radius = N÷8
    body = AutoBody((x,t)->√sum(abs2, x .- (N÷3,N÷4)) - radius) + AutoBody((x,t)->√sum(abs2, x .- (2N÷3,3N÷4)) - radius)

    ΔT = T(1)
    α = ε/ΔT |> T
    g = T(1)

    νκ = g*α*ΔT*radius^3/Ra
    ν = sqrt(νκ*Pr)
    κ = sqrt(νκ/Pr)

    Uf = sqrt(g*α*ΔT*radius)
    U = Uf/√Pr |> T

    ThermalSimulation(NN, (0,0), radius; U, Δt=0.05, ν, α, κ, g=(i,x,t)-> i==2 ? -g : zero(g), perdir=(1,),θ0=nothing,θb=(x,t)-> -sign(x[2]-N÷2)*ΔT, body)
end

# using CUDA
sim = hot_cold_circle(128)#; mem=CuArray)

# duration=8,step=0.08 covers plume formation through the steady rising plume reached around tU/L≈5

# Visualize the temperature
sim_gif!(sim; duration=20, step=0.08, f=(dat,sim)->copyto!(dat,sim.flow.θ), clims=(-1,1), cfill=:seismic,
              plotbody=true, video=joinpath(@__DIR__,"HotColdCylinder_temp.gif"))

# Visualize the vorticity
# sim_gif!(sim; duration=40, step=0.08, clims=(-30,30), cfill=:seismic,
#               plotbody=true, video=joinpath(@__DIR__,"HotColdCylinder_vort.gif"))