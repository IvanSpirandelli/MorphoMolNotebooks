function visualize_configuration_sequence!(cgl::GridLayout, obs_index, input, output, figure_config, label_text = "Configuration")
    template_centers = input["template_centers"]
    template_radii = input["template_radii"]
    n_atoms_per_mol = length(template_radii)
    n_mol = input["n_mol"]
    configurations = [MorphoMol.Utilities.get_point3f_realization(state, template_centers) for state in output["states"][figure_config["vis_range"]]]

    Label(cgl[0, 1:2], text = label_text, fontsize = figure_config["title_fs"])
    conf_ax = LScene(cgl[1:2, 1:2])
    points = @lift([configurations[$obs_index][i] for i in 1:n_mol*n_atoms_per_mol])
    colors = vcat([[j for _ in 1:n_atoms_per_mol] for j in 1:n_mol]...)
    conf_ms = 10

    scatter!(conf_ax, points, markersize = conf_ms, color = colors, colormap = :rainbow)
    cgl
end

function visualize_interface_sequence!(igl::GridLayout, obs_index, interfaces, show_wireframe, figure_config, label_text = "Interface")
    ifils = interfaces["filtrations"][figure_config["vis_range"]]
    barycentric_subdivisions = interfaces["vertices"][figure_config["vis_range"]]

    min_v = minimum([minimum([e[2] for e in filtration if length(e[1]) == 1]) for filtration in ifils])
    max_v = sqrt(maximum([maximum([e[2] for e in filtration if length(e[1]) == 1]) for filtration in ifils]))

    Label(igl[0, 1:2], text = label_text, fontsize = figure_config["title_fs"])
    i_sc = LScene(igl[1:2, 1:2], show_axis=false, scenekw = (lights = [AmbientLight(RGBf(1.0, 1.0, 1.0))],))

    faces = [[TriangleFace(e[1]) for e in filtration if length(e[1]) == 3] for filtration in ifils]
    meshes = [GeometryBasics.Mesh(bcs, fs) for (bcs, fs) in zip(barycentric_subdivisions, faces)]
    colors = [[e[2] for e in filtration if length(e[1]) == 1] for filtration in ifils]

    l_ms, l_cs = @lift(meshes[$obs_index]), @lift(colors[$obs_index])
    mesh!(i_sc, l_ms, color = l_cs, colorrange = (min_v, max_v), colormap = :viridis)
    if show_wireframe
        wireframe!(i_sc, l_ms, color=:black, linewidth=0.25)
    end
    igl
end

function visualize_ma_measures!(mgl::GridLayout, obs_index, input, output, figure_config, label_text = "Morphometric Measures")
    xs = [i for i in 1:length(output["states"][figure_config["vis_range"]])];
    pf = MorphoMol.Energies.get_prefactors(input["rs"], input["η"])
    Vs = pf[1] .* output["Vs"][figure_config["vis_range"]]
    As = pf[2] .* output["As"][figure_config["vis_range"]]
    Cs = pf[3] .* output["Cs"][figure_config["vis_range"]]
    Xs = pf[4] .* output["Xs"][figure_config["vis_range"]]
    OLs = input["overlap_slope"] .* output["OLs"][figure_config["vis_range"]]
    
    Label(mgl[0, 1:5], text = label_text, fontsize = figure_config["title_fs"])
    cm = :Set2_5
    cr = (1, 5)

    Vs_ax = Axis(mgl[1, 1], title = L"pV")
    As_ax = Axis(mgl[1, 2], title = L"\sigma A")
    Cs_ax = Axis(mgl[1, 3], title = L"\kappa C")
    Xs_ax = Axis(mgl[1, 4], title = L"\overline{\kappa} X")
    OLs_ax = Axis(mgl[1, 5], title = L"OLs")

    Vs_mark = @lift(Point2f($obs_index, $@lift(Vs[$obs_index])))
    As_mark = @lift(Point2f($obs_index, $@lift(As[$obs_index])))
    Cs_mark = @lift(Point2f($obs_index, $@lift(Cs[$obs_index])))
    Xs_mark = @lift(Point2f($obs_index, $@lift(Xs[$obs_index])))
    OLs_mark = @lift(Point2f($obs_index, $@lift(OLs[$obs_index])))

    plot_ms = figure_config["plot_ms"]
    tracker_ms = figure_config["tracker_ms"]

    scatter!(Vs_ax, xs, Vs, color = 1, colormap = cm, colorrange = cr, markersize = plot_ms)
    scatter!(Vs_ax, Vs_mark, color=:black, markersize = tracker_ms)

    scatter!(As_ax, xs, As, color = 2, colormap = cm, colorrange = cr, markersize = plot_ms)
    scatter!(As_ax, As_mark, color=:black, markersize = tracker_ms)

    scatter!(Cs_ax, xs, Cs, color = 3, colormap = cm, colorrange = cr, markersize = plot_ms)
    scatter!(Cs_ax, Cs_mark, color=:black, markersize = tracker_ms)

    scatter!(Xs_ax, xs, Xs, color = 4, colormap = cm, colorrange = cr, markersize = plot_ms)
    scatter!(Xs_ax, Xs_mark, color=:black, markersize = tracker_ms)

    scatter!(OLs_ax, xs, OLs, color = 5, colormap = cm, colorrange = cr, markersize = plot_ms)
    scatter!(OLs_ax, OLs_mark, color=:black, markersize = tracker_ms)
    mgl
end

function visualize_interface_persistence_measures!(ipgl::GridLayout, obs_index, (zero_tps, one_tps, stps), persistence_weights, figure_config, label_text = "Interface Persistence")
    xs = [i for i in 1:length(stps)];
    
    cm = :Spectral_4
    cr = (1, 4)
    
    Label(ipgl[0, 1:3], text = label_text, fontsize = figure_config["title_fs"])
    zero_tp_ax = Axis(ipgl[2, 1], title = "$(persistence_weights[1]) b_0")
    one_tp_ax = Axis(ipgl[1, 2], title = "$(persistence_weights[2]) b_1")

    zero_tp_mark = @lift(Point2f($obs_index, $@lift(zero_tps[$obs_index])))
    one_tp_mark = @lift(Point2f($obs_index, $@lift(one_tps[$obs_index])))

    stp_ax = Axis(ipgl[1, 1], title = L"Σ")
    stp_mark = @lift(Point2f($obs_index, $@lift(stps[$obs_index])))

    plot_ms = figure_config["plot_ms"]
    tracker_ms = figure_config["tracker_ms"]
    
    scatter!(zero_tp_ax, xs, zero_tps, color = 1, colormap = cm, colorrange = cr, markersize = plot_ms)
    scatter!(zero_tp_ax, zero_tp_mark, color=:black, markersize = tracker_ms)

    scatter!(one_tp_ax, xs, one_tps, color = 2, colormap = cm, colorrange = cr, markersize = plot_ms)
    scatter!(one_tp_ax, one_tp_mark, color=:black, markersize = tracker_ms)

    scatter!(stp_ax, xs, stps, color = 4, colormap = cm, colorrange = cr, markersize = plot_ms)
    scatter!(stp_ax, stp_mark, color=:black, markersize = tracker_ms)
end

function visualize_alpha_shape_persistence_measures!(ipgl::GridLayout, obs_index, (zero_tps, one_tps, two_tps, stps), persistence_weights, figure_config, label_text = "Interface Persistence")
    xs = [i for i in 1:length(stps)];
    
    cm = :Spectral_4
    cr = (1, 4)
    
    Label(ipgl[0, 1:3], text = label_text, fontsize = figure_config["title_fs"])
    zero_tp_ax = Axis(ipgl[2, 1], title = "$(persistence_weights[1]) H_0")
    one_tp_ax = Axis(ipgl[1, 2], title = "$(persistence_weights[2]) H_1")
    two_tp_ax = Axis(ipgl[2, 2], title = "$(persistence_weights[3]) H_2")

    zero_tp_mark = @lift(Point2f($obs_index, $@lift(zero_tps[$obs_index])))
    one_tp_mark = @lift(Point2f($obs_index, $@lift(one_tps[$obs_index])))
    two_tp_mark = @lift(Point2f($obs_index, $@lift(two_tps[$obs_index])))

    stp_ax = Axis(ipgl[1, 1], title = L"Σ")
    stp_mark = @lift(Point2f($obs_index, $@lift(stps[$obs_index])))

    plot_ms = figure_config["plot_ms"]
    tracker_ms = figure_config["tracker_ms"]
    
    scatter!(zero_tp_ax, xs, zero_tps, color = 1, colormap = cm, colorrange = cr, markersize = plot_ms)
    scatter!(zero_tp_ax, zero_tp_mark, color=:black, markersize = tracker_ms)

    scatter!(one_tp_ax, xs, one_tps, color = 2, colormap = cm, colorrange = cr, markersize = plot_ms)
    scatter!(one_tp_ax, one_tp_mark, color=:black, markersize = tracker_ms)

    scatter!(two_tp_ax, xs, two_tps, color = 3, colormap = cm, colorrange = cr, markersize = plot_ms)
    scatter!(two_tp_ax, two_tp_mark, color=:black, markersize = tracker_ms)

    scatter!(stp_ax, xs, stps, color = 4, colormap = cm, colorrange = cr, markersize = plot_ms)
    scatter!(stp_ax, stp_mark, color=:black, markersize = tracker_ms)
end

function visualize_energy_and_theta!(etgl::GridLayout, obs_index, input, output, figure_config, label_text = "Energy and Theta")
    Es = output["Es"][figure_config["vis_range"]]
    mol_type = input["mol_type"]    
    exp_template_centers = MorphoMol.Utilities.TWOTMVSU_EXPERIMENTAL_ASSEMBLY[mol_type]["template_centers"]
    exp_state = MorphoMol.Utilities.TWOTMVSU_EXPERIMENTAL_ASSEMBLY[mol_type]["state"]
    thetas = [MorphoMol.Utilities.average_offset_distance(exp_template_centers, input["template_centers"], exp_state, state) for state in output["states"][figure_config["vis_range"]]]
    xs =  [i for i in 1:length(Es)]
    
    #Energy and Theta
    Label(etgl[0, 1:2], text = label_text, fontsize = figure_config["title_fs"])
    E_ax = Axis(etgl[1:2, 1], title = L"F_{sol}")
    theta_ax = Axis(etgl[1:2, 2], title = L"\Theta")

    theta_mark = @lift(Point2f($obs_index, $@lift(thetas[$obs_index])))
    E_mark = @lift(Point2f($obs_index, $@lift(Es[$obs_index])))

    plot_ms = figure_config["plot_ms"]
    tracker_ms = figure_config["tracker_ms"]

    scatter!(theta_ax, xs, thetas, color = :magenta, markersize = plot_ms)
    scatter!(theta_ax, theta_mark, color=:black, markersize = tracker_ms)

    scatter!(E_ax, xs, Es, color = :blue, markersize = plot_ms)
    scatter!(E_ax, E_mark, color=:black, markersize = tracker_ms)
end