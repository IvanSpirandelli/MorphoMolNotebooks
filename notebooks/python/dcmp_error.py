import oineus as oin
import numpy as np
import torch
import diode

def get_points(file):
    points = []
    with open(file) as f:
        for line in f:
            x, y, z = line.split()
            points.append(np.array([float(x), float(y), float(z)]))
    return points

def get_interface_diagram_and_filtration(points, n_atoms_per_mol):
    points = np.asarray(points)
    simplices = diode.fill_alpha_shapes(points)
    fil = oin.Filtration_double([oin.Simplex_double(s[0], s[1]) for s in simplices])

    def is_multi(sigma):
        return len(set(v // n_atoms_per_mol for v in sigma.vertices)) >= 2
    fil = fil.subfiltration(is_multi)

    def recalculate_filtration_value(cell):
        parts = [v // n_atoms_per_mol for v in cell.vertices]
        n_parts = len(set(parts))
        p_agg = [np.array([0.0, 0.0, 0.0]) for _ in range(n_parts)]
        weights = [0 for _ in range(n_parts)]
        for i, p in enumerate(parts):
            p_agg[p] += points[cell.vertices[i]]
            weights[p] += 1
        bcs = [p_agg[i] / weights[i] for i in range(n_parts)]
        filtration_value = sum([np.linalg.norm(bcs[i] - bcs[j]) for i in range(n_parts) for j in range(i + 1, n_parts)]) / n_parts
        cell.value = filtration_value
        return cell

    altered_cells = [recalculate_filtration_value(cell) for cell in fil.cells()]
    fil = oin.Filtration_double(altered_cells)

    dcmp = oin.Decomposition(fil, True)
    params = oin.ReductionParams()
    params.clearing_opt = False
    dcmp.reduce(params)
    dgm = dcmp.diagram(fil, include_inf_points=False)
    return dgm, fil

points = get_points("test.txt")
get_interface_diagram_and_filtration(points, 1206)