import Pkg; 

Pkg.activate("../MorphoMolMonteCarlo/Project.toml")
Pkg.instantiate()

using MorphoMol

Pkg.activate("Project.toml")
ENV["PYCALL_JL_RUNTIME_PYTHON"] = Sys.which("python3")
using PyCall
using JLD2
using LinearAlgebra
using Rotations
println("Python exec path: ", pyimport("sys").executable)

py"""
import oineus as oin
import numpy as np
import torch
import diode

def calculate_total_persistence(points):
    points = np.asarray(points)
    simplices = diode.fill_alpha_shapes(points)
    fil = oin.Filtration_double([oin.Simplex_double(s[0], s[1]) for s in simplices])
    # no cohomology
    dualize = False
    # create VRU decomposition object, does not perform reduction yet
    dcmp = oin.Decomposition(fil, dualize)
    rp = oin.ReductionParams()
    rp.compute_u = rp.compute_v = True
    rp.n_threads = 1
    # perform reduction
    dcmp.reduce(rp)
    # now we can acess V, R and U
    # indices are sorted_ids of simplices == indices in fil.cells()
    V = dcmp.v_data
    simplices = fil.simplices()
    dgms = dcmp.diagram(fil, include_inf_points=False)
    dim=1
    dgm_dim = dcmp.diagram(fil, include_inf_points=False)[dim]
    return sum([e[1] - e[0] for e in dgm_dim])
"""

function calculate_total_persistence(points::Vector{Vector{Float64}})
    py"calculate_total_persistence"(points)
end


perturb_all(x, Σ) = x .+ (randn(length(x)) .* Σ)

function perturb_single_randomly_chosen(x, σ_r, σ_t)
    x_cand = deepcopy(x)
    i  = rand(0:(length(x)÷6)-1)
    x_cand[(i*6)+1:(i*6)+6] = x_cand[(i*6)+1:(i*6)+6] .+ randn(6) .* [σ_r, σ_r, σ_r, σ_t, σ_t, σ_t]
    x_cand
end

function solvation_free_energy_and_measures_in_bounds(x::Vector{Float64}, template_mol::Matrix{Float64}, radii::Vector{Float64}, rs::Float64, prefactors::AbstractVector, overlap_jump::Float64, overlap_slope::Float64, persistence_weight::Float64, bounds::Float64, delaunay_eps::Float64)
    if any(0.0 >= e || e >= bounds for e in x[4:6:end]) || any(0.0 >= e || e >= bounds for e in x[5:6:end]) || any(0.0 >= e || e >= bounds for e in x[6:6:end])
        return Inf, [Inf, Inf, Inf, Inf, Inf]
    end
    n_mol = length(x) ÷ 6
    n_atoms_per_mol = size(template_mol)[2]
    flat_realization = MorphoMol.Utilities.get_flat_realization(x, template_mol)
    tp = calculate_total_persistence(Vector{Vector{Float64}}(eachcol(reshape(flat_realization, (3, Int(length(flat_realization) / 3))))))
    measures = MorphoMol.Energies.get_geometric_measures_and_overlap_value(flat_realization, n_atoms_per_mol, radii, rs, overlap_jump, overlap_slope, delaunay_eps)
    measures = [measures; tp]
    # TODO: Slightly irritating that we have to set the prefactor for persistence here, while the prefactor for overlap penalty is overlap_slope and passed to measure calc
    sum(measures .* [prefactors; [1.0, persistence_weight]]), measures
end

function simulate(simulation_time_minutes = 1.0)
    template_file = "assets/input/jld2/2_6r7m_init/6r7m_protor_1.jld2"
    T = 2.5
    β = 1.0 / T

    rs = 1.4
    η = 0.3665
    pf = MorphoMol.Energies.get_prefactors(rs, η)
    overlap_slope = 0.85
    persistence_weight = -0.1
    bnds = 150.0
    delaunay_eps = 100.0

    σ_r = 0.15
    σ_t = 1.25

    @load template_file template_mol template_radii x_init
    n_atoms_per_mol = length(template_mol) ÷ 3
    n_mol = length(x_init) ÷ 6
    template_mol = reshape(template_mol,(3,n_atoms_per_mol))
    radii = vcat([template_radii for i in 1:n_mol]...);

    β = 1.0 / T
    pf = MorphoMol.Energies.get_prefactors(rs, η)
    Σ = vcat([[σ_r, σ_r, σ_r, σ_t, σ_t, σ_t] for _ in 1:n_mol]...)

    energy(x) = solvation_free_energy_and_measures_in_bounds(x, template_mol, radii, rs, pf, 0.0, overlap_slope, persistence_weight, bnds, delaunay_eps)
    perturbation(x) = perturb_single_randomly_chosen(x, σ_r, σ_t)
    #perturbation(x) = perturb_all(x, Σ)

    rwm = MorphoMol.Algorithms.RandomWalkMetropolis(energy, perturbation, β)

    input = MorphoMol.Algorithms.MorphometricSimulationInput(
        template_mol,
        template_radii,
        n_mol,
        σ_r,
        σ_t,
        rs,
        η,
        pf,
        0.0,
        overlap_slope,
        T,
        0.0,
        0
    )

    output = MorphoMol.Algorithms.SimulationOutput(
        Vector{Vector{Float64}}([]),
        Vector{Float64}([]),
        Vector{Vector{Float64}}([
            Vector{Float64}([]),
            Vector{Float64}([]),
            Vector{Float64}([]),
            Vector{Float64}([]),
            Vector{Float64}([]),
            Vector{Float64}([])
            ]),
        Vector{Float32}([])
    )

    MorphoMol.Algorithms.simulate!(rwm, output, deepcopy(x_init), simulation_time_minutes);

    output_directory = "assets/output/persistence_testerino/"
    name = "a"
    mkpath(output_directory)
    @save "$(output_directory)/$(name).jld2" input output
end