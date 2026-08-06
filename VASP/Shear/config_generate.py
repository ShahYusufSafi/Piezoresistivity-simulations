import numpy as np
from ase.io import read, write
import os, glob

EQ = '../Equilibirium_run/outputs/Yusuf'
atoms0 = read(f'{EQ}/POSCAR')
cell0  = np.array(atoms0.get_cell())
reset_dirs = True
strains = [-0.010, -0.005, -0.002, 0.000, +0.002, +0.005, +0.010]

INCAR_RELAX = """ISTART=0
ICHARG=2
ENCUT=320
ISMEAR=0
SIGMA=0.05
ISYM=0
ISIF=2
IBRION=2
NSW=60
EDIFFG=-0.01
GGA=PE
"""
INCAR_SCF = """ISTART=0
ICHARG=2
ENCUT=320
ISMEAR=0
SIGMA=0.05
ISYM=0
EDIFF=1E-6
LCHARG=.TRUE.
NSW=0
GGA=PE
"""
INCAR_BAND = """ISTART=0
ICHARG=11
ENCUT=320
ISMEAR=0
ISYM=0
LORBIT=11
GGA=PE
"""

def write_line_kpoints(filepath, npts=40):
    segs = [((0.0, 0.5, 0.5), "X_par"),  ((0.0, 0.0, 0.0), "\\Gamma"),
            ((0.0, 0.0, 0.0), "\\Gamma"), ((0.5, 0.0, 0.5), "X_perp")]
    out = ["Delta-valley split path", str(npts), "line", "reciprocal"]
    for i in range(0, len(segs), 2):
        (k1, l1), (k2, l2) = segs[i], segs[i+1]
        out.append(f"{k1[0]:.6f} {k1[1]:.6f} {k1[2]:.6f}  {l1}")
        out.append(f"{k2[0]:.6f} {k2[1]:.6f} {k2[2]:.6f}  {l2}")
        out.append("")
    open(filepath, "w").write("\n".join(out).rstrip() + "\n")

for eps in strains:
    M = np.eye(3) + np.array([[0, eps, eps],
                              [eps, 0, eps],
                              [eps, eps, 0]])          # [111] shear
    a = atoms0.copy()
    a.set_cell(cell0 @ M, scale_atoms=True)

    d = f'eps_{eps:+.4f}'
    os.makedirs(d, exist_ok=True)
    if reset_dirs:
        for f in glob.glob(f'{d}/*'):
            os.remove(f)

    write(f'{d}/POSCAR', a, format='vasp')

    open(f'{d}/INCAR.relax', 'w').write(INCAR_RELAX)   # shear needs ionic relaxation
    open(f'{d}/INCAR.scf',   'w').write(INCAR_SCF)
    open(f'{d}/INCAR.band',  'w').write(INCAR_BAND)

    os.system(f'cp {EQ}/POTCAR {d}/')
    write_line_kpoints(f'{d}/KPOINTS.line', npts=40)
    open(f'{d}/KPOINTS.mesh', 'w').write("auto\n0\nGamma\n8 8 8\n0 0 0\n")