function quadratic_additive(step::Int, iterations::Int, T_init::Float64, T_min::Float64)
    T_min + (T_init - T_min)*((iterations - step)/iterations)^2
end

function zig_zag(step::Int, iterations::Int, T_init::Float64, T_min::Float64, level::Vector{Float64})
    iteration_level = [Int(round(i * iterations / length(level))) for i in 0:length(level)]
    current = findfirst(x -> x >= step, iteration_level)
    if current == length(level)+1
        return T_min
    end
    quadratic_additive(step-iteration_level[current], Int(round(iterations / length(level))), T_init * level[current], T_min)
end

function solvation_free_energy_in_bounds(x::Vector{Float64}, radii::Vector{Float64}, rs::Float64, prefactors::AbstractVector, overlap_jump::Float64, overlap_slope::Float64, bounds::Float64, delaunay_eps::Float64)
    if any(0.0 >= e || e >= bounds for e in x[1:3:end]) || any(0.0 >= e || e >= bounds for e in x[2:3:end]) || any(0.0 >= e || e >= bounds for e in x[3:3:end])
        return Inf
    end
    measures = MorphoMol.Energies.get_geometric_measures_and_overlap_value(x, 1, radii, rs, overlap_jump, overlap_slope, delaunay_eps)
    sum(measures .* [prefactors; 1.0])
end

function perturb_single_randomly_chosen(x, σ_t)
    x_cand = deepcopy(x)
    i  = rand(0:(length(x)÷3)-1)
    x_cand[(i*3)+1:(i*3)+3] = x_cand[(i*3)+1:(i*3)+3] .+ (randn(3) .* σ_t)
    x_cand
end

function simulated_annealing_scan(;n = 8, rss = 0.025:0.025:1.5, etas = 0.025:0.025:0.4925, pf_id = "wb", temp_finder_target_rate = 0.5)
    lowest_energy_simulated_states = Dict{@NamedTuple{n::Int64, rs::Float64, eta::Float64}, @NamedTuple{centers::Vector{Float64}, E::Float64}}()
    lowest_energy_no_mcc_simulated_states = Dict{@NamedTuple{n::Int64, rs::Float64, eta::Float64}, @NamedTuple{centers::Vector{Float64}, E::Float64}}()
    for rs in rss, eta in etas
        println("rs: $rs, eta: $eta")
        template_centers = [0.0, 0.0, 0.0]
        bounds = n / 2.0 + 3.0
        x_init = get_random_configuration_without_overlap_in_bounds(n, bounds)
        radii = fill(1.0, n)

        pf = MorphoMol.Energies.get_prefactors(rs, eta, pf_id)
        σ_t = rs
        T_search = MorphoMolNotebooks.get_dispersed_energy(n, rs, pf) / 15.0
        T_init = temp_finder(
            n,
            rs,
            pf,
            σ_t,
            T_search, 
            50000, # Temp finder iterations
            temp_finder_target_rate # Target acceptance rate
        )

        T_end = 0.0
        iterations = 250000 * n
        
        energy(x) = solvation_free_energy_in_bounds(x, radii, rs, pf, 10.0^6, 10.0^6, bounds, 100.0)
        perturbation(x) = perturb_single_randomly_chosen(x, σ_t)
        temp_r(i) = zig_zag(i, iterations, T_init, T_end, [1.0, 0.5, 0.3, 0.3, 0.1, 0.5, 0.3, 0.1])

        output = Dict{String,Vector}(
            "states" => [],
            "Es" => [],
            "Ts" => []
        )

        sa = MorphoMol.Algorithms.SimulatedAnnealing(energy, perturbation, temp_r)
        x_end, E_end, alpha = MorphoMol.Algorithms.simulate!(sa, x_init, iterations, output);
        c_sim = output["states"][argmin(output["Es"])]
        E_sim = minimum(output["Es"])

        if n <= 8
            E_mcc = minimal_mcc(n, rs, pf).E
            mcc_name = minimal_mcc(n, rs, pf).name

            if E_mcc > E_sim && !occursin("dispersed", mcc_name)
                println("MCC: $(E_mcc), Sim: $(E_sim)")
                println("SUCCESS: Simulation minimum lower than lowest MCC")
            end

            add_to_simulation_minima!(c_sim, E_sim, n, rs, eta, lowest_energy_simulated_states, lowest_energy_no_mcc_simulated_states)
            println("T_init: $(T_init) | Alpha: $(alpha)")
            println("+-------------------+-------------------+")
        end
    end
    @save "assets/parameter_space_scans/$(pf_id)/$(pf_id)_scan_results_$(randstring(6)).jld2" lowest_energy_simulated_states lowest_energy_no_mcc_simulated_states
end

function temp_finder(
    n::Int,
    rs::Float64,
    pf::Vector{Float64},
    σ_t::Float64,
    T_search::Float64, 
    iterations::Int,
    target_acceptance_rate::Float64)

    @assert T_search > 0.0
    Es = generate_transitions(n, rs, pf, σ_t, T_search, iterations)
    calculate_T0(Es, target_acceptance_rate)
end

function generate_transitions(
    n::Int,
    rs::Float64,
    pf::Vector{Float64},
    σ_t::Float64, #Normal Deviation for perturbation moves
    T_search::Float64,
    iterations::Int,
    )
    
    template_centers = [0.0, 0.0, 0.0]
    bounds = n / 2.0 + 3.0
    x_init = get_random_configuration_without_overlap_in_bounds(n, bounds)
    radii = fill(1.0, n)

    energy(x) = solvation_free_energy_in_bounds(x, radii, rs, pf, 10.0^6, 10.0^6, bounds, 100.0)
    perturbation(x) = MorphoMolNotebooks.perturb_single_randomly_chosen(x, σ_t)
    temp_r(i) = T_search

    output = Dict{String,Vector}(
        "states" => [],
        "Es" => []
    )

    sa = MorphoMol.Algorithms.SimulatedAnnealing(energy, perturbation, temp_r)
    MorphoMol.Algorithms.simulate!(sa, x_init, iterations, output);
    output["Es"]
end