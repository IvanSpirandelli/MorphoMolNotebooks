module MorphoMolNotebooks

using Distances
using MorphoMol
using GeometryBasics
using LinearAlgebra
using JLD2
using Graphs
using Graphs.Experimental
using GLMakie
using Random
using Rotations
using StaticArrays

export visualize_energy_and_theta!, visualize_interface_persistence_measures!, visualize_alpha_shape_persistence_measures!, visualize_ma_measures!, visualize_interface_sequence!, visualize_configuration_sequence!
export configuration_to_poly, poly_to_configuration, STANDARD_COLORS
export get_mccs, minimal_mcc, get_energies_of_mccs
export is_configuration_equivalent_to_mcc, is_equivalent_configuration, is_dispersed
export calculate_initial_temperature, add_to_simulation_minima!, has_hard_sphere_overlap, get_random_configuration_without_overlap_in_bounds
export simulated_annealing_scan

include("visualization.jl")
include("configurations_as_polys.jl")
include("hs_cluster_utilities.jl")
include("sa_scan_setup.jl")
include("evaluation_utilities.jl")

STANDARD_COLORS = [(0.32, 0.06, 0.25), (0.98, 0.54, 0.14), (0.06, 0.3, 0.36), (0.34, 0.54, 0.98)]

end #module MorphoMolNotebooks