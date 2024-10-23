function visualize_persistence(input, output)
    realizations = [[e for e in eachcol(hcat(MorphoMol.Utilities.get_matrix_realization(state, input["template_centers"])...))] for state in output["states"]]
    realizations = [hcat(realizations[i]...) for i in 1:length(realizations)]
    dgms = output["PDGMs"]
    Es = output["Es"]
    mol_type = input["mol_type"]    
    exp_template_centers = MorphoMol.Utilities.TWOTMVSU_EXPERIMENTAL_ASSEMBLY[mol_type]["template_centers"]
    exp_state = MorphoMol.Utilities.TWOTMVSU_EXPERIMENTAL_ASSEMBLY[mol_type]["state"]
    thetas = [MorphoMol.Utilities.average_offset_distance(exp_template_centers, input["template_centers"], exp_state, state) for state in output["states"]]
    xs = [i for i in 1:length(dgms)];
    zero_ps = [MorphoMol.Energies.get_total_persistence(dgm[1], 0.1) for dgm in dgms]
    one_ps = [MorphoMol.Energies.get_total_persistence(dgm[2], -0.1) for dgm in dgms]
    two_ps = [MorphoMol.Energies.get_total_persistence(dgm[3], -0.1) for dgm in dgms]
    tps = [MorphoMol.Energies.get_total_persistence_summed(dgm, [0.1, -0.1, -0.1, 0.0]) for dgm in dgms];

    f = Figure(fontsize = 12)
    cm = :seaborn_muted
    cr = (1, 11)
    ms = 5
    title_fs = 15
    #Slider
    sl_i = Slider(f[4, 1:5], range = 1:length(dgms), startvalue = 1)
    x = sl_i.value

    # Persistence
    cm = :Spectral_4
    cr = (1, 4)

    pg = GridLayout(f[1:2, 3:5])
    Label(pg[0, 1:3], text = "Persistence", fontsize = title_fs)
    zero_p_ax = Axis(pg[2, 1], title = L"0.1 \text{total}(b_0)")
    one_p_ax = Axis(pg[1, 2], title = L"-0.1 \text{total}(b_1)")
    two_p_ax = Axis(pg[2, 2], title = L"-0.1 \text{total}(b_2)")
    tp_ax = Axis(pg[1, 1], title = L"Σ")

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

    pd_ax_a = Axis(pg[1, 3], title = L"b_1", xticks = 0:25:100, yticks = 0:25:100)
    pd_ax_b = Axis(pg[2, 3], title = L"b_2", xticks = 0:25:100, yticks = 0:25:100)
    dgm_points_1 = @lift([Point2f(dgms[$x][2][i,1], dgms[$x][2][i,2]) for i in 1:size(dgms[$x][2])[1]])
    dgm_points_2 = @lift([Point2f(dgms[$x][3][i,1], dgms[$x][3][i,2]) for i in 1:size(dgms[$x][3])[1]])


    scatter!(pd_ax_a, dgm_points_1, color = 2, colormap = cm, colorrange = cr, markersize = ms)
    lines!(pd_ax_a, [0, 100], [0, 100], color=:black)
    scatter!(pd_ax_b, dgm_points_2, color = 3, colormap = cm, colorrange = cr, markersize = ms)
    lines!(pd_ax_b, [0, 100], [0, 100], color=:black)

    #Energy and Theta
    eg = GridLayout(f[3, 1:2])
    Label(eg[0, 1:2], text = "Energy", fontsize = title_fs)
    theta_ax = Axis(eg[1, 1], title = L"\Theta")
    E_ax = Axis(eg[1, 2], title = L"F_{sol}")

    theta_mark = @lift(Point2f($x, $@lift(thetas[$x])))
    E_mark = @lift(Point2f($x, $@lift(Es[$x])))

    scatter!(theta_ax, xs, thetas, color = :magenta, markersize = ms)
    scatter!(theta_ax, theta_mark, color=:black, markersize = 10)

    scatter!(E_ax, xs, Es, color = :blue,  markersize = ms)
    scatter!(E_ax, E_mark, color=:black, markersize = 10)
    
    # Configuration
    cg = GridLayout(f[1:2, 1:2])
    Label(cg[0, 1:2], text = "Configuration", fontsize = title_fs)
    conf_ax = LScene(cg[1:2, 1:2])
    n_atoms = size(realizations[1])[2]
    n_mol = input["n_mol"]
    points = @lift([Point3f(realizations[$x][1,i], realizations[$x][2,i], realizations[$x][3,i]) for i in 1:n_atoms])
    colors = vcat([[j for _ in 1:1206] for j in 1:n_mol]...)
    conf_ms = 10

    scatter!(conf_ax, points, markersize = conf_ms, color = colors, colormap = :rainbow)

    # for (label, layout) in zip(["A", "B", "C", "D"], [cg, pg, mg, eg])
    #     Label(layout[1, 1, TopLeft()], label,
    #         fontsize = 12,
    #         font = :bold,
    #         padding = (0, 10, 10, 0),
    #         halign = :right)
    # end

    f
end

function visualize_persistence_and_ma(input, output)
    realizations = [[e for e in eachcol(hcat(MorphoMol.Utilities.get_matrix_realization(state, input["template_centers"])...))] for state in output["states"]]
    realizations = [hcat(realizations[i]...) for i in 1:length(realizations)]
    dgms = output["PDGMs"]
    Es = output["Es"]
    pf = MorphoMol.Energies.get_prefactors(input["rs"], input["η"])
    Vs = pf[1] .* output["Vs"]
    As = pf[2] .* output["As"]
    Cs = pf[3] .* output["Cs"]
    Xs = pf[4] .* output["Xs"]
    OLs = input["overlap_slope"] .* output["OLs"]
    mol_type = input["mol_type"]    
    exp_template_centers = MorphoMol.Utilities.TWOTMVSU_EXPERIMENTAL_ASSEMBLY[mol_type]["template_centers"]
    exp_state = MorphoMol.Utilities.TWOTMVSU_EXPERIMENTAL_ASSEMBLY[mol_type]["state"]
    thetas = [MorphoMol.Utilities.average_offset_distance(exp_template_centers, input["template_centers"], exp_state, state) for state in output["states"]]
    xs = [i for i in 1:length(dgms)];
    zero_ps = [MorphoMol.Energies.get_total_persistence(dgm[1], 0.1) for dgm in dgms]
    one_ps = [MorphoMol.Energies.get_total_persistence(dgm[2], -0.1) for dgm in dgms]
    two_ps = [MorphoMol.Energies.get_total_persistence(dgm[3], -0.1) for dgm in dgms]
    tps = [MorphoMol.Energies.get_total_persistence_summed(dgm, [0.1, -0.1, -0.1, 0.0]) for dgm in dgms];

    f = Figure(fontsize = 12)
    cm = :seaborn_muted
    cr = (1, 11)
    ms = 5
    title_fs = 15
    #Slider
    sl_i = Slider(f[5, 1:5], range = 1:length(dgms), startvalue = 1)
    x = sl_i.value

    # Persistence
    cm = :Spectral_4
    cr = (1, 4)

    pg = GridLayout(f[1:2, 3:5])
    Label(pg[0, 1:3], text = "Persistence", fontsize = title_fs)
    zero_p_ax = Axis(pg[2, 1], title = L"0.1 \text{total}(b_0)")
    one_p_ax = Axis(pg[1, 2], title = L"-0.1 \text{total}(b_1)")
    two_p_ax = Axis(pg[2, 2], title = L"-0.1 \text{total}(b_2)")
    tp_ax = Axis(pg[1, 1], title = L"Σ")

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

    pd_ax_a = Axis(pg[1, 3], title = L"b_1", xticks = 0:25:100, yticks = 0:25:100)
    pd_ax_b = Axis(pg[2, 3], title = L"b_2", xticks = 0:25:100, yticks = 0:25:100)
    dgm_points_1 = @lift([Point2f(dgms[$x][2][i,1], dgms[$x][2][i,2]) for i in 1:size(dgms[$x][2])[1]])
    dgm_points_2 = @lift([Point2f(dgms[$x][3][i,1], dgms[$x][3][i,2]) for i in 1:size(dgms[$x][3])[1]])


    scatter!(pd_ax_a, dgm_points_1, color = 2, colormap = cm, colorrange = cr, markersize = ms)
    lines!(pd_ax_a, [0, 100], [0, 100], color=:black)
    scatter!(pd_ax_b, dgm_points_2, color = 3, colormap = cm, colorrange = cr, markersize = ms)
    lines!(pd_ax_b, [0, 100], [0, 100], color=:black)

    #Measures
    mg = GridLayout(f[3, 1:5])
    Label(mg[0, 1:5], text = "Geometric Measures", fontsize = title_fs)
    cm = :Set2_5
    cr = (1, 5)

    Vs_ax = Axis(mg[1, 1], title = L"pV")
    As_ax = Axis(mg[1, 2], title = L"\sigma A")
    Cs_ax = Axis(mg[1, 3], title = L"\kappa C")
    Xs_ax = Axis(mg[1, 4], title = L"\overline{\kappa} X")
    OLs_ax = Axis(mg[1, 5], title = L"OLs")

    Vs_mark = @lift(Point2f($x, $@lift(Vs[$x])))
    As_mark = @lift(Point2f($x, $@lift(As[$x])))
    Cs_mark = @lift(Point2f($x, $@lift(Cs[$x])))
    Xs_mark = @lift(Point2f($x, $@lift(Xs[$x])))
    OLs_mark = @lift(Point2f($x, $@lift(OLs[$x])))

    scatter!(Vs_ax, xs, Vs, color = 1, colormap = cm, colorrange = cr, markersize = ms)
    scatter!(Vs_ax, Vs_mark, color=:black, markersize = 10)

    scatter!(As_ax, xs, As, color = 2, colormap = cm, colorrange = cr, markersize = ms)
    scatter!(As_ax, As_mark, color=:black, markersize = 10)

    scatter!(Cs_ax, xs, Cs, color = 3, colormap = cm, colorrange = cr, markersize = ms)
    scatter!(Cs_ax, Cs_mark, color=:black, markersize = 10)

    scatter!(Xs_ax, xs, Xs, color = 4, colormap = cm, colorrange = cr, markersize = ms)
    scatter!(Xs_ax, Xs_mark, color=:black, markersize = 10)

    scatter!(OLs_ax, xs, OLs, color = 5, colormap = cm, colorrange = cr, markersize = ms)
    scatter!(OLs_ax, OLs_mark, color=:black, markersize = 10)

    #Energy and Theta
    eg = GridLayout(f[4, 1:2])
    Label(eg[0, 1:2], text = "Energy", fontsize = title_fs)
    theta_ax = Axis(eg[1, 1], title = L"\Theta")
    E_ax = Axis(eg[1, 2], title = L"F_{sol}")

    theta_mark = @lift(Point2f($x, $@lift(thetas[$x])))
    E_mark = @lift(Point2f($x, $@lift(Es[$x])))

    scatter!(theta_ax, xs, thetas, color = :magenta, markersize = ms)
    scatter!(theta_ax, theta_mark, color=:black, markersize = 10)

    scatter!(E_ax, xs, Es, color = :blue,  markersize = ms)
    scatter!(E_ax, E_mark, color=:black, markersize = 10)
    
    # Configuration
    cg = GridLayout(f[1:2, 1:2])
    Label(cg[0, 1:2], text = "Configuration", fontsize = title_fs)
    conf_ax = LScene(cg[1:2, 1:2])
    n_atoms = size(realizations[1])[2]
    n_mol = input["n_mol"]
    points = @lift([Point3f(realizations[$x][1,i], realizations[$x][2,i], realizations[$x][3,i]) for i in 1:n_atoms])
    colors = vcat([[j for _ in 1:1206] for j in 1:n_mol]...)
    conf_ms = 10

    scatter!(conf_ax, points, markersize = conf_ms, color = colors, colormap = :rainbow)

    # for (label, layout) in zip(["A", "B", "C", "D"], [cg, pg, mg, eg])
    #     Label(layout[1, 1, TopLeft()], label,
    #         fontsize = 12,
    #         font = :bold,
    #         padding = (0, 10, 10, 0),
    #         halign = :right)
    # end

    f
end