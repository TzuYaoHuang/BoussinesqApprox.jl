using WaterLily, BoussinesqApprox
using LinearAlgebra: dot, norm
using Plots

function hot_cold_circle(N;Ra=10^6,ε=0.001,Pr=10,mem=Array,T=Float32)
    NN = (N,2N)
    radius = N÷8
    body = AutoBody((x,t)->√sum(abs2, x .- (N÷2,7N÷5)) - radius) + AutoBody((x,t)->√sum(abs2, x .- (N÷2,3N÷5)) - radius)

    ΔT = T(1)
    α = ε/ΔT |> T
    g = T(1)

    νκ = g*α*ΔT*radius^3/Ra
    ν = sqrt(νκ*Pr)
    κ = sqrt(νκ/Pr)

    Uf = sqrt(g*α*ΔT*radius)
    U = Uf/√Pr |> T

    ThermalSimulation(NN, (0,0), radius; U, Δt=0.05, ν, α, κ, g=(i,x,t)-> i==2 ? -g : zero(g), perdir=(1,),θ0=nothing,θb=(x,t)-> -sign(x[2]-N)*ΔT, body, mem)
end

# `body` is `AutoBody(cold) + AutoBody(hot)` = `SetBody(min,cold,hot)`. Recurse through the set-op
# tree and apply the force kernel to each leaf's own SDF, so each cylinder's contribution is
# isolated even though `sim.flow` was advanced with the two bodies combined.
each_force(f,flow,body) = f(flow,body)
each_force(f,flow,body::WaterLily.SetBody{typeof(min)}) = mapreduce(bod->each_force(f,flow,bod),hcat,(body.a,body.b))

uncentered_cor(x,y) = dot(x,y)/(norm(x)*norm(y))

# top cylinder is cold (body.a), bottom cylinder is hot (body.b) -- see θb's sign in hot_cold_circle
function fx_correlation(Ra,Pr;N=64,duration=20,n0=4,mem=Array,T=Float32)
    sim = hot_cold_circle(N;Ra,Pr,mem,T)
    t₀ = sim_time(sim)

    cold_x = Float64[]; hot_x = Float64[]
    t = sum(sim.flow.Δt[1:end-1])
    while t < (t₀+duration)*sim.L/sim.U
        mom_step!(sim.flow,sim.pois)

        # force *on* each cylinder = -(force the kernel measures on the fluid at its surface)
        Fp = -each_force(WaterLily.pressure_force,sim.flow,sim.body)
        Fν = -each_force(WaterLily.viscous_force,sim.flow,sim.body)
        push!(cold_x, Fp[1,1]+Fν[1,1])
        push!(hot_x,  Fp[1,2]+Fν[1,2])

        t += sim.flow.Δt[end]
    end

    keep = n0:length(hot_x) # discard the initial transient before the pressure solver has settled
    uncentered_cor(hot_x[keep], cold_x[keep])
end

# --- regime sweep ---
Ras = 10.0 .^ range(4,7,length=4)
Prs = 10.0 .^ range(-1,2,length=4)

C = Array{Float64}(undef,length(Ras),length(Prs))
for (j,Pr) in enumerate(Prs), (i,Ra) in enumerate(Ras)
    et = @elapsed C[i,j] = fx_correlation(Ra,Pr)
    @info "Ra=$Ra, Pr=$Pr -> corr=$(round(C[i,j],digits=3))  ($(round(et,digits=1))s)"
end

heatmap(log10.(Ras),log10.(Prs),C';
    xlabel="log10(Ra)", ylabel="log10(Pr)",
    colorbar_title="corr(hot Fx, cold Fx)", clims=(-1,1), c=:balance,
    title="Hot/cold cylinder Fx-correlation regime map")
savefig(joinpath(@__DIR__,"HotColdCylinder_regime_map.png"))
