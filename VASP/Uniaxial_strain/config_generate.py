# ============================================================================
#  Uniaxial [100] strain sweep generator.
#
#  CHANNEL: Gamma_12 (tetragonal) + a small Gamma_1 (hydrostatic) part, since
#  diag(1+eps,1,1) is 1/3 hydrostatic + tetragonal. Probes Xi_u via the
#  Delta-valley splitting.
#
#  NO RELAXATION STAGE: for axial [100] strain the two-atom basis is frozen by
#  symmetry (no Kleinman displacement), so affine scale_atoms=True is exact.
#  Only the [111] shear channel needs INCAR.relax. Hence this generator writes
#  INCAR.scf and INCAR.band but NOT INCAR.relax -> jobs.sh auto-skips Stage A.
#
#  DUAL-X PATH: KPOINTS.line visits X_par (along the strain axis) AND X_perp
#  (transverse). The valley splitting lives in the DIFFERENCE between these two;
#  a single cubic X point would hide it. At eps=0 the two must be degenerate.
#
#  STAGED INCARs: SCF (uniform mesh) -> CHGCAR, then non-SCF bands (ICHARG=11)
#  on KPOINTS.line. A cold ICHARG=11 with no CHGCAR gives unconverged, wrong
#  eigenvalues - the bug that produced the bad Xi_u.
# ============================================================================


import numpy as np
from ase.io import read, write
import os, glob

EQ = '../Equilibirium_run/outputs'
atoms0 = read(f'{EQ}/POSCAR')
cell0  = np.array(atoms0.get_cell())
reset_dirs = True
strains = [-0.010, -0.005, -0.002, 0.000, +0.002, +0.005, +0.010]

INCAR_SCF = """ISTART=0
ICHARG=2
ENCUT=320
PREC=Accurate
ISMEAR=0
SIGMA=0.05
ISYM=0
EDIFF=1E-7
LCHARG=.TRUE.
NSW=0
GGA=PE
"""
INCAR_BAND = """ISTART=0
ICHARG=11
ENCUT=320
PREC=Accurate
ISMEAR=0
SIGMA=0.05
ISYM=0
LORBIT=11
GGA=PE
"""



def write_line_kpoints(filepath, npts=40):
    segs = [((0.0, 0.5, 0.5), "X_par"),  ((0.0, 0.0, 0.0), "\\Gamma"),
            ((0.0, 0.0, 0.0), "\\Gamma"), ((0.5, 0.0, 0.5), "X_perp")]
    out = ["X_par - Gamma - X_perp", str(npts), "line", "reciprocal"]
    for i in range(0, len(segs), 2):
        (k1, l1), (k2, l2) = segs[i], segs[i+1]
        out.append(f"{k1[0]:.6f} {k1[1]:.6f} {k1[2]:.6f}  {l1}")
        out.append(f"{k2[0]:.6f} {k2[1]:.6f} {k2[2]:.6f}  {l2}")
        out.append("")
    open(filepath, "w").write("\n".join(out).rstrip() + "\n")

for eps in strains:
    M = np.diag([1+eps, 1.0, 1.0])          # [100] uniaxial
    a = atoms0.copy()
    a.set_cell(cell0 @ M, scale_atoms=True)

    d = f'eps_{eps:+.4f}'
    os.makedirs(d, exist_ok=True)
    if reset_dirs:
        for f in glob.glob(f'{d}/*'):
            os.remove(f)

    write(f'{d}/POSCAR', a, format='vasp')

    open(f'{d}/INCAR.scf',  'w').write(INCAR_SCF)     # no INCAR.relax for uniaxial
    open(f'{d}/INCAR.band', 'w').write(INCAR_BAND)
    os.system(f'cp {EQ}/POTCAR {d}/')
    write_line_kpoints(f'{d}/KPOINTS.line', npts=40)
    open(f'{d}/KPOINTS.mesh', 'w').write("auto\n0\nGamma\n8 8 8\n0 0 0\n")