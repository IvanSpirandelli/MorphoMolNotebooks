function get_mccs(n::Int)
    @load "src/assets/mcc_configurations.jld2" mcc_configurations
    filter(x -> length(x.radii) == n, mcc_configurations)
end

function minimal_mcc(n::Int, rs::Float64, pf::Vector{Float64})
    all = get_energies_of_mccs(n, rs, pf)
    best = argmin(x -> x.E, all)
end

function get_energies_of_mccs(n::Int, rs::Float64, pf::Vector{Float64})
    mcc_configurations = get_mccs(n)
    [(name = cluster.name, centers = cluster.centers, radii = cluster.radii, E = sum(MorphoMol.Energies.get_geometric_measures(cluster.centers, cluster.radii, rs, 100.0) .* pf)) for cluster in filter(x -> length(x.radii) == n, mcc_configurations)]
end

are_balls_intersecting(c1, c2, r1, r2) = euclidean(c1, c2) <= r1 + r2

function is_dispersed(flat_centers, rs)
    centers = [e for e in eachcol(reshape(flat_centers, (3, length(flat_centers)÷3)))]
    n = length(centers)
    for i in 1:n
        for j in i+1:n
            if euclidean(centers[i], centers[j]) < 2.0*(1.0 + rs)
                return false
            end
        end
    end
    return true
end

function get_dispersed_energy(n, rs, pf)
    n * MorphoMol.Energies.solvation_free_energy([0.0, 0.0, 0.0], [1.0], rs, pf)
end

function construct_graph_from_configuration(flat_centers; threshold = 1.0)
    centers = [e for e in eachcol(reshape(flat_centers, (3, length(flat_centers)÷3)))]
    n = length(centers)
    edges = [(i,j) for i in 1:n for j in i:n if are_balls_intersecting(centers[i], centers[j], threshold, threshold) && i != j]
    g = Graph(n)
    for edge in edges
        add_edge!(g,edge[1], edge[2])
    end
    g
end

function is_equivalent_configuration(fc1, fc2; threshold = 1.0)
    g1 = construct_graph_from_configuration(fc1; threshold = threshold)
    g2 = construct_graph_from_configuration(fc2; threshold = threshold)
    has_isomorph(g1, g2)
end

function is_configuration_equivalent_to_mcc(configuration; threshold = 1.0)
    n = div(length(configuration), 3)
    any([is_equivalent_configuration(configuration, mcc.centers; threshold = threshold) for mcc in MorphoMolNotebooks.get_mccs(n)])
end

function calculate_initial_temperature(n_mol, rs, pf)
    E_disp = n_mol * sum(MorphoMol.Energies.get_geometric_measures([0.0, 0.0, 0.0], [1.0], rs, 100.0) .* pf)
    E_min = MorphoMolNotebooks.minimal_mcc(n_mol, rs, pf).E
    return 0.02 * (E_disp - E_min)
end

function add_to_simulation_minima!(
    c_sim::Vector{Float64}, 
    E_sim::Float64, 
    n::Int, 
    rs::Float64, 
    eta::Float64, 
    lowest_energy_simulated_states::Dict{@NamedTuple{n::Int64, rs::Float64, eta::Float64}, @NamedTuple{centers::Vector{Float64}, E::Float64}}, 
    lowest_energy_no_mcc_simulated_states::Dict{@NamedTuple{n::Int64, rs::Float64, eta::Float64}, @NamedTuple{centers::Vector{Float64}, E::Float64}}
    )
    if !is_configuration_equivalent_to_mcc(c_sim, threshold = 1.0 + rs / 2.0)
        if !haskey(lowest_energy_no_mcc_simulated_states, (n = n, rs = rs, eta = eta))
            lowest_energy_no_mcc_simulated_states[(n = n, rs = rs, eta = eta)] = (centers = c_sim, E = E_sim)
        elseif no_mcc_minima[(n = n, rs = rs, eta = eta)].E > E_sim
            lowest_energy_no_mcc_simulated_states[(n = n, rs = rs, eta = eta)] = (centers = c_sim, E = E_sim)
        end
        if !haskey(lowest_energy_simulated_states, (n = n, rs = rs, eta = eta))
            lowest_energy_simulated_states[(n = n, rs = rs, eta = eta)] = (centers = c_sim, E = E_sim)
        elseif lowest_energy_simulated_states[(n = n, rs = rs, eta = eta)].E > E_sim
            lowest_energy_simulated_states[(n = n, rs = rs, eta = eta)] = (centers = c_sim, E = E_sim)
        end
    else #Equivalent to mcc
        if !haskey(lowest_energy_simulated_states, (n = n, rs = rs, eta = eta))
            lowest_energy_simulated_states[(n = n, rs = rs, eta = eta)] = (centers = c_sim, E = E_sim)
        elseif lowest_energy_simulated_states[(n = n, rs = rs, eta = eta)].E > E_sim
            lowest_energy_simulated_states[(n = n, rs = rs, eta = eta)] = (centers = c_sim, E = E_sim)
        end
    end
end

function has_hard_sphere_overlap(x::Vector{Float64})
    centers = [e for e in eachcol(reshape(x, (3, length(x)÷3)))]
    n = length(centers)
    for i in 1:n
        for j in i+1:n
            if euclidean(centers[i], centers[j]) < 2.0
                return true
            end
        end
    end
    return false
end

function get_random_configuration_without_overlap_in_bounds(n, bounds; max_tries = 100)
    x = rand(n*3) .* bounds
    for i in 1:max_tries 
        if has_hard_sphere_overlap(x)
            x = rand(n*3) .* bounds
        else
            return x
        end
    end
    @assert false
end