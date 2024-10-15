function interactive_simulation_analysis(input, output)
    realizations = [[e for e in eachcol(hcat(MorphoMol.Utilities.get_matrix_realization(state, input["template_centers"])...))] for state in output["states"]]
    realizations = [hcat(realizations[i]...) for i in 1:length(realizations)]
    dgms = output["PDGMs"]
    Es = output["Es"]
    mol_type = input["mol_type"]    
    exp_template_centers = MorphoMol.Utilities.TWOTMVSU_EXPERIMENTAL_ASSEMBLY[mol_type]["template_centers"]
    exp_state = MorphoMol.Utilities.TWOTMVSU_EXPERIMENTAL_ASSEMBLY[mol_type]["state"]
    thetas = [MorphoMol.Utilities.average_offset_distance(exp_template_centers, input["template_centers"], exp_state, state) for state in output["states"]]
    xs = [i for i in 1:length(dgms)];
    zero_ps = [MorphoMol.Energies.get_persistence(dgm[1], 0.1) for dgm in dgms]
    one_ps = [MorphoMol.Energies.get_persistence(dgm[2], -0.1) for dgm in dgms]
    two_ps = [MorphoMol.Energies.get_persistence(dgm[3], -0.1) for dgm in dgms]
    tps = [MorphoMol.Energies.get_total_persistence(dgm, [0.1, -0.1, -0.1, 0.0]) for dgm in dgms];

    f = Figure(fontsize = 12)
    #cm = :Accent_8
    cm = :Dark2_8
    cr = (1,8)

    #Slider
    sl_i = Slider(f[4, 1:4], range = 1:length(dgms), startvalue = 1)
    x = sl_i.value

    # Persistence Diagrams
    pd_ax_a = Axis(f[1, 1], title = "Persistence codim = 2", xticks = 0:25:100, yticks = 0:25:100)
    pd_ax_b = Axis(f[1, 2], title = "Persistence codim = 3", xticks = 0:25:100, yticks = 0:25:100)
    dgm_points_2 = @lift([Point2f(dgms[$x][2][i,1], dgms[$x][2][i,2]) for i in 1:size(dgms[$x][2])[1]])
    dgm_points_3 = @lift([Point2f(dgms[$x][3][i,1], dgms[$x][3][i,2]) for i in 1:size(dgms[$x][3])[1]])

    ms = 5
    scatter!(pd_ax_a, dgm_points_2, color = 2, colormap = cm, colorrange = cr, markersize = ms)
    lines!(pd_ax_a, [0, 100], [0, 100], color=:black)
    scatter!(pd_ax_b, dgm_points_3, color = 3, colormap = cm, colorrange = cr, markersize = ms)
    lines!(pd_ax_b, [0, 100], [0, 100], color=:black)

    #Energy and Theta
    theta_ax = Axis(f[2, 1], title = L"\Theta")
    E_ax = Axis(f[2, 2], title = L"F_{sol}")

    theta_mark = @lift(Point2f($x, $@lift(thetas[$x])))
    E_mark = @lift(Point2f($x, $@lift(Es[$x])))

    scatter!(theta_ax, xs, thetas, color = 5, colormap = cm, colorrange = cr, markersize = ms)
    scatter!(theta_ax, theta_mark, color=:black, markersize = 10)

    scatter!(E_ax, xs, Es, color = 6, colormap = cm, colorrange = cr, markersize = ms)
    scatter!(E_ax, E_mark, color=:black, markersize = 10)

    # Configuration
    conf_ax = LScene(f[1:2, 3:4])
    n_atoms = size(realizations[1])[2]
    n_mol = input["n_mol"]
    points = @lift([Point3f(realizations[$x][1,i], realizations[$x][2,i], realizations[$x][3,i]) for i in 1:n_atoms])
    colors = vcat([[j for _ in 1:1206] for j in 1:n_mol]...)
    conf_ms = 10

    scatter!(conf_ax, points, markersize = conf_ms, color = colors, colormap = :rainbow)

    #Persistence Measures
    zero_p_ax = Axis(f[3, 1], title = "0.1 * total 0 persistence")
    one_p_ax = Axis(f[3, 2], title = "-0.1 * total 1 persistence")
    two_p_ax = Axis(f[3, 3], title = "-0.1 * total 2 persistence")
    tp_ax = Axis(f[3, 4], title = "Total Persistences Summed")

    tp_mark = @lift(Point2f($x, $@lift(tps[$x])))

    zero_p_mark = @lift(Point2f($x, $@lift(zero_ps[$x])))
    one_p_mark = @lift(Point2f($x, $@lift(one_ps[$x])))
    two_p_mark = @lift(Point2f($x, $@lift(two_ps[$x])))

    scatter!(zero_p_ax, xs, zero_ps, color = 1, colormap = cm, colorrange = cr, markersize = ms)
    scatter!(zero_p_ax, zero_p_mark, color=:black, markersize = 10)

    scatter!(one_p_ax, xs, one_ps, color = 2, colormap = cm, colorrange = cr, markersize = ms)
    scatter!(one_p_ax, one_p_mark, color=:black, markersize = 10)

    scatter!(two_p_ax, xs, two_ps, color = 3, colormap = cm, colorrange = cr, markersize = ms)
    scatter!(two_p_ax, two_p_mark, color=:black, markersize = 10)

    scatter!(tp_ax, xs, tps, color = 4, colormap = cm, colorrange = cr, markersize = ms)
    scatter!(tp_ax, tp_mark, color = :black, markersize = 10)
    f
end