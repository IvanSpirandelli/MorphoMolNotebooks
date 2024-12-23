function get_theta(input, output)
    mol_type = input["mol_type"]
    if !(mol_type in keys(MorphoMol.TWOTMVSU_EXPERIMENTAL_ASSEMBLY))
        return Inf
    end
    if input["n_mol"] == 2
        exp_template_centers = MorphoMol.TWOTMVSU_EXPERIMENTAL_ASSEMBLY[mol_type]["template_centers"]
        exp_state = MorphoMol.TWOTMVSU_EXPERIMENTAL_ASSEMBLY[mol_type]["state"]
        return MorphoMol.average_offset_distance(exp_template_centers, input["template_centers"], exp_state, output["states"][argmin(output["Es"])])
    elseif input["n_mol"] == 3
        # Consecutive assembly
        R0 = RotMatrix(@SMatrix[1.000000  0.000000  0.000000; 0.000000  1.000000  0.000000; 0.000000  0.000000  1.000000])
        T0 = @SVector[0.00000, 0.00000, 0.00000]

        R1 = RotMatrix(@SMatrix[0.628642  0.777695  0.000000; -0.777695  0.628642  0.000000; 0.000000  0.000000  1.000000])
        T1 = @SVector[-69.28043, 195.91352, -49.35000]

        R2 = RotMatrix(@SMatrix[0.874450  0.485115  0.000000; -0.485115  0.874450  0.000000; 0.000000  0.000000  1.000000])
        T2 = @SVector[-61.30589, 104.11829, -47.94000]

        R3 = RotMatrix(@SMatrix[0.992567  0.121696  0.000000; -0.121696  0.992567  0.000000; 0.000000  0.000000  1.000000])
        T3 = @SVector[-19.48193, 22.01644, -46.53000]

        consecutive = Vector{Float64}([])
        for (R,T) in zip([log(R1), log(R2), log(R3)], [T1, T2, T3])
            consecutive = [consecutive; [R[3,2], R[1,3], R[2,1], T[1], T[2], T[3]]]
        end
        
        # Two top one bottom
        R1 = RotMatrix(@SMatrix[0.628642  0.777695  0.000000; -0.777695  0.628642  0.000000; 0.000000  0.000000 1.000000])
        T1 = @SVector[-69.28043, 195.91352, -49.35000]

        R2 = RotMatrix(@SMatrix[0.874450  0.485115  0.000000; -0.485115  0.874450  0.000000; 0.000000  0.000000 1.000000])
        T2 = @SVector[-61.30589, 104.11829, -47.94000]

        R3 = RotMatrix(@SMatrix[0.52145615 0.85327808 0.000000; -0.85327808 0.52145615 -0.000000; 0.000000 0.000000 1.000000])
        T3 = @SVector[-63.89221, 227.07555, -26.79]

        R4 = RotMatrix(@SMatrix[0.803441  0.595384  0.000000; -0.595384 0.803441 0.000000; 0.000000 0.000000 1.000000])
        T4 = @SVector[-67.99970, 135.02619, -25.38000]

        two_top_one_bottom = Vector{Float64}([])
        for (R,T) in zip([log(R1), log(R2), log(R4)], [T1, T2, T4])
            two_top_one_bottom = [two_top_one_bottom; [R[3,2], R[1,3], R[2,1], T[1], T[2], T[3]]]
        end

        # Two bottom one top
        R1 = RotMatrix(@SMatrix[0.628642  0.777695  0.000000; -0.777695  0.628642  0.000000; 0.000000  0.000000 1.000000])
        T1 = @SVector[-69.28043, 195.91352, -49.35000]

        R2 = RotMatrix(@SMatrix[0.874450  0.485115  0.000000; -0.485115  0.874450  0.000000; 0.000000  0.000000 1.000000])
        T2 = @SVector[-61.30589, 104.11829, -47.94000]

        R3 = RotMatrix(@SMatrix[0.52145615 0.85327808 0.000000; -0.85327808 0.52145615 -0.000000; 0.000000 0.000000 1.000000])
        T3 = @SVector[-63.89221, 227.07555, -26.79]

        R4 = RotMatrix(@SMatrix[0.803441  0.595384  0.000000; -0.595384 0.803441 0.000000; 0.000000 0.000000 1.000000])
        T4 = @SVector[-67.99970, 135.02619, -25.38000]

        two_bottom_one_top = Vector{Float64}([])
        for (R,T) in zip([log(R1), log(R3), log(R4)], [T1, T3, T4])
            two_bottom_one_top = [two_bottom_one_top; [R[3,2], R[1,3], R[2,1], T[1], T[2], T[3]]]
        end
        mindex = argmin(output["Es"])
        simulated_assembly_state = output["states"][mindex];

        exp_template_centers = MorphoMol.TWOTMVSU_EXPERIMENTAL_ASSEMBLY[mol_type]["template_centers"]
        sim_template_centers = input["template_centers"]
        permutations = [[1, 2, 3], [1, 3, 2], [2, 1, 3], [2, 3, 1], [3,1,2], [3,2,1]]
        min_theta = Inf
        for experimental_state in [consecutive, two_top_one_bottom, two_bottom_one_top]
            cand = minimum([MorphoMol.sum_of_permutation(sim_template_centers, exp_template_centers, simulated_assembly_state, experimental_state, [1, 2, 3], perm) for perm in permutations])
            if cand < min_theta
                min_theta = cand
            end
        end
        return min_theta
    end
    Inf
end

function get_low_energy_and_theta_states(folder::String, mol_type::String; E_cutoff_ratio = 0.95, theta_cutoff = 3.0)
    E_mins = Vector{Float64}([])
    thetas = []
    files = Vector{String}([])
    for file in readdir(folder)
        if split(file, ".")[end] == "jld2"
            try
                @load "$folder$file" input output
                theta = get_theta(input, output)
                println(theta)
                push!(files, file)
                push!(E_mins, minimum(output["Es"]))
                push!(thetas, theta)
            catch e
                println(e)
            end
        end
    end

    E_min = minimum(E_mins)
    E_max = maximum(E_mins)
    E_cutoff = E_max + E_cutoff_ratio * (E_min - E_max)
    low_energy_selection = [E < E_cutoff for E in E_mins]
    low_theta_selection =  [theta < theta_cutoff for theta in thetas]
    
    for file in files[low_energy_selection .| low_theta_selection]
        println(file)
    end
end