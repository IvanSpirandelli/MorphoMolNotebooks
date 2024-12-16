import numpy as np
import chromatic_tda as chro

def compute_six_pack():
    points = np.random.rand(1000,3)
    labels = [i % 2 for i in range(len(points))]
    chro_alpha = chro.ChromaticAlphaComplex(points, labels)
    simplicial_complex = chro_alpha.get_simplicial_complex(
                sub_complex="mono-chromatic",
                full_complex="all",
                relative="mono-chromatic",
    )  # these options make sense for three colors; for two use, e.g., just sub_complex='mono-chromatic'
    six_pack = simplicial_complex.bars_six_pack()

compute_six_pack()