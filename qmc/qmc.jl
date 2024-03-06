
julia --project=qmc qmc/qmc.jlusing BloqadeQMC
using Random

using Bloqade
using Yao: mat, ArrayReg
using Bloqade.BloqadeExpr.LinearAlgebra
using BloqadeQMC.Measurements
using BloqadeQMC.Measurements: value, uncertainty
using BloqadeQMC.BinningAnalysis
using Statistics

using Makie, CairoMakie

nsites = 9
atoms = generate_sites(ChainLattice(), nsites, scale = 5.72)
@info "generated atoms on a 1D chain: $atoms"

Ω = 2π * 4
Δ = 2π * 10
h = rydberg_h(atoms; Δ, Ω)
@info "generated Rydberg Hamiltonian: $h"

h_qmc = rydberg_qmc(h);
@info "generated QMC Hamiltonian: $h_qmc"

EQ_MCS = 100;
MCS = 100_000;

M = 50;
ts = BinaryThermalState(h_qmc, M);
@info "generated thermal state: $ts"

d = Diagnostics();

β = 0.5;

rng = MersenneTwister(3214);

[mc_step_beta!(rng, ts, h_qmc,β, d, eq=true) for i in 1:EQ_MCS] # equilibration phase

densities_QMC = zeros(nsites)
occs = zeros(MCS, nsites)

for i in 1:MCS # Monte Carlo Steps
    mc_step_beta!(rng, ts, h_qmc,β, d, eq=false) do lsize, ts, h_qmc
        SSE_slice = sample(h_qmc,ts, 1)
        occs[i, :] = ifelse.(SSE_slice .== true, 1.0, 0.0)
    end
end

for jj in 1:nsites
    densities_QMC[jj] = mean(occs[:,jj])
end

fig, = barplot(1:nsites, densities_QMC; xlabel = "Site number", ylabel = "Occupation density")
Makie.save("occupation_density.png", fig)
@info "saved occupation density plot to `occupation_density.png`"

nsites = 9;
atoms = generate_sites(ChainLattice(), nsites, scale = 5.72);

Δ_step = 15;
Δ = LinRange(-2π * 9, 2π * 9, Δ_step);

energy_QMC = []

for ii in 1:Δ_step
    h_ii = rydberg_h(atoms; Δ = Δ[ii], Ω)
    h_ii_qmc = rydberg_qmc(h_ii)
    ts_ii = BinaryThermalState(h_ii_qmc,M)
    d_ii = Diagnostics()

    [mc_step_beta!(rng, ts_ii, h_ii_qmc, β, d_ii, eq=true) for i in 1:EQ_MCS] #equilibration phase

    ns = zeros(MCS)

    for i in 1:MCS # Monte Carlo Steps
        ns[i] = mc_step_beta!(rng, ts_ii, h_ii_qmc, β, d_ii, eq=false)
    end

    energy(x) = -x / β + h_ii_qmc.energy_shift  # The energy shift here ensures that all matrix elements are non-negative. See Merali et al for details.
    BE = LogBinner(energy.(ns)) # Binning analysis
    τ_energy = tau(BE)
    ratio = 2 * τ_energy + 1
    energy_binned = measurement(mean(BE), std_error(BE)*sqrt(ratio))
    append!(energy_QMC, energy_binned)
end

energy_ED = zeros(Δ_step)

for ii in 1:Δ_step
    h_ii = rydberg_h(atoms; Δ = Δ[ii], Ω)
    h_m = Matrix(mat(h_ii))
    energies, vecs = LinearAlgebra.eigen(h_m)

    w = exp.(-β .* (energies .- energies[1]))
    energy_ED[ii] = sum(w .* energies) / sum(w)
end

fig, = Makie.scatter(Δ/2π, value.(energy_QMC); yerror=uncertainty.(energy_QMC), label="QMC", marker=:x)
Makie.save("order_param_QMC.png", fig)
@info "saved order parameter plot to `order_param_QMC.png`"

densities_QMC = []
order_param_QMC = []

for ii in 1:Δ_step
    h_ii = rydberg_h(atoms; Δ = Δ[ii], Ω)
    h_ii_qmc = rydberg_qmc(h_ii)
    ts_ii = BinaryThermalState(h_ii_qmc,M)
    d_ii = Diagnostics()

    [mc_step_beta!(rng, ts_ii, h_ii_qmc, β, d_ii, eq=true) for i in 1:EQ_MCS] #equilibration phase

    order_param = zeros(MCS)

    for i in 1:MCS # Monte Carlo Steps
        mc_step_beta!(rng, ts_ii, h_ii_qmc, β, d_ii, eq=false) do lsize, ts_ii, h_ii_qmc
            SSE_slice = sample(h_ii_qmc, ts_ii, 1) # occ = 0,1
            spin = 2 .* SSE_slice .- 1 # spin = -1, 1
            order_param[i] = abs(sum(spin[1:2:end]) - sum(spin[2:2:end]))/length(spin)
        end
    end

    BD = LogBinner(order_param)
    τ_energy = tau(BD)
    ratio = 2 * τ_energy + 1
    energy_binned = measurement(mean(BD), std_error(BD)*sqrt(ratio))
    append!(order_param_QMC, measurement(mean(BD), std_error(BD)*sqrt(ratio)) )
end

fig, = Makie.scatter(Δ/2π, value.(order_param_QMC); yerror=uncertainty.(order_param_QMC), label="", marker=:x)
Makie.save("order_param_QMC.png", fig)
@info "saved order parameter plot to `order_param_QMC.png`"
