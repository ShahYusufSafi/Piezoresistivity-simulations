import numpy as np
from ase.io import read, write
import os

# ============================================================================
#  mass_eps0 — effective-mass validation run (UNSTRAINED Si)
#
#  PURPOSE
#  -------
#  Extract the longitudinal (m_l) and transverse (m_t) effective masses of the
#  Si conduction-band Delta-valley from first principles, and validate them
#  against literature (Hamaguchi Ch.4: m_l = 0.98 m_e, m_t = 0.19 m_e). This is
#  a one-off validation of the whole DFT -> band -> curvature-fit workflow BEFORE
#  it is applied to the strained eps_* sweep. It is NOT part of the sweep.
#
#  WHY A DEDICATED RUN (and why it can give BOTH masses)
#  -----------------------------------------------------
#  The effective mass is the curvature of E(k) at the valley minimum, and
#  curvature is DIRECTIONAL: the ellipsoidal valley has one curvature along its
#  axis (-> m_l) and a different one across it (-> m_t). A line-mode run only
#  resolves curvature along the directions it actually samples. So a single line
#  along the valley axis yields m_l alone — it contains no information about the
#  perpendicular curvature.
#
#  This run therefore samples TWO perpendicular legs through the same valley
#  centre, combined into ONE KPOINTS.line:
#      L leg  — along the valley axis (Cartesian y-hat)  ->  m_l
#      T leg  — perpendicular to it   (Cartesian x-hat)  ->  m_t
#  One SCF charge density is computed once on the mesh; both legs are evaluated
#  non-self-consistently on top of it (ICHARG=11), so no second SCF is needed —
#  only the sampled k-points differ between the legs.
#
#  WHY THIS COULD NOT BE DONE "FROM THE START"
#  -------------------------------------------
#  The transverse leg must be centred on the valley and perpendicular to the
#  valley AXIS — but you cannot define "perpendicular to the axis" until a first
#  run has located that axis. The prior eps_+0.0000 run was that reconnaissance:
#  it located the valley (and, as a by-product, gave a first m_l). Its printed
#  valley position (k0_cart below) is the input that lets this run aim the
#  transverse leg correctly. Locate first, then sample perpendicular.
#
#  THE ONE SUBTLE TRAP (Cartesian vs fractional)
#  ---------------------------------------------
#  A step written in FRACTIONAL coordinates (e.g. [0, delta, 0]) is NOT
#  perpendicular in real space, because the FCC reciprocal basis is
#  non-orthogonal. Such a "perpendicular" leg silently re-samples the valley
#  axis and returns m_l TWICE. The fix, used below: define the leg directions in
#  CARTESIAN space (y-hat, x-hat), then map to fractional via B^-1 for the
#  KPOINTS file. The commented tests at the bottom verify (a) both legs share the
#  valley centre as midpoint and (b) cos(angle between legs) ~ 0.
#
#  STAGED INCARs
#  -------------
#  INCAR.scf  : self-consistent on the uniform mesh -> writes CHGCAR.
#  INCAR.band : non-SCF (ICHARG=11) on KPOINTS.line, reads that CHGCAR. A cold
#               ICHARG=11 with no CHGCAR gives unconverged, wrong eigenvalues.
#
#  OUTPUT / CROSS-CHECK
#  --------------------
#  The L leg reproduces m_l independently of eps_+0.0000 — if the two agree, the
#  two runs are sampling equivalent <100> valleys consistently and the leg
#  geometry is trustworthy for m_t. mass_eps0 is then the self-contained source
#  of truth for both masses (one valley, one consistent frame).
# ============================================================================


EQ = '../Equilibirium_run/outputs/Yusuf'
atoms0 = read(f'{EQ}/POSCAR')

d = 'mass_eps0'          # dedicated one-off, NOT part of the strain sweep
os.makedirs(d, exist_ok=True)
write(f'{d}/POSCAR', atoms0, format='vasp')   # unstrained

# ---- INCARs: reuse staged SCF -> non-SCF band recipe ----
open(f'{d}/INCAR.scf', 'w').write(
"ISTART=0\nICHARG=2\nENCUT=320\nISMEAR=0\nSIGMA=0.05\nISYM=0\n"
"EDIFF=1E-6\nLCHARG=.TRUE.\nNSW=0\nGGA=PE\n")
open(f'{d}/INCAR.band', 'w').write(
"ISTART=0\nICHARG=11\nENCUT=320\nISMEAR=0\nSIGMA=0.05\nISYM=0\nLORBIT=11\nGGA=PE\n")
os.system(f'cp {EQ}/POTCAR {d}/')
open(f'{d}/KPOINTS.mesh', 'w').write("auto\n0\nGamma\n8 8 8\n0 0 0\n")

# ---- the mass path: two short legs crossing at the Si Delta-valley ----
# Si CBM (Delta) at fractional ~ (0.42, 0, 0.42) in the FCC primitive recip basis.
# Longitudinal = along Gamma->X (the (h,0,h) line). Transverse = perpendicular.
#k0 = np.array([0.42, 0.00, 0.42])          # valley center (fractional)
#d_long = np.array([0.06, 0.00, 0.06])      # along the axis (small step)
#d_perp = np.array([0.00, 0.06, 0.00])      # perpendicular (ky)

recip = atoms0.cell.reciprocal() * 2*np.pi
B    = np.array(recip)          # frac -> cart
Binv = np.linalg.inv(B)         # cart -> frac

k0_cart = np.array([0.0, 0.9768, 0.0])   # valley center, straight from your print
delta   = 0.08                            # Å⁻¹ half-width

k0     = k0_cart @ Binv

# longitudinal = along the valley axis ŷ ; transverse = along x̂ (⊥ to ŷ)
long_dir = np.array([0.0, 1.0, 0.0])
perp_dir = np.array([1.0, 0.0, 0.0])

L_cart = [k0_cart - delta*long_dir, k0_cart + delta*long_dir]
T_cart = [k0_cart - delta*perp_dir, k0_cart + delta*perp_dir]

L_end = [p @ Binv for p in L_cart]   # -> fractional for KPOINTS
T_end = [p @ Binv for p in T_cart]


def line(a, b, l1, l2):
    return [f"{a[0]:.4f} {a[1]:.4f} {a[2]:.4f}  {l1}",
            f"{b[0]:.4f} {b[1]:.4f} {b[2]:.4f}  {l2}", ""]

npts = 30
out = ["mass legs: longitudinal then transverse", str(npts), "line", "reciprocal"]
out += line(L_end[0], L_end[1], "L-", "L+")   # along ŷ  -> m_l
out += line(T_end[0], T_end[1], "T-", "T+")   # along x̂  -> m_t
open(f'{d}/KPOINTS.line', 'w').write("\n".join(out).rstrip() + "\n")


'''
# tests
print("L midpoint (frac):", np.round((np.array(L_end[0])+L_end[1])/2, 4))
print("T midpoint (frac):", np.round((np.array(T_end[0])+T_end[1])/2, 4))



# these two MUST be equal, and equal to k0_cart @ Binv
L = np.array([L_end[0], L_end[1]])
T = np.array([T_end[0], T_end[1]])
dL = L[1] - L[0]      # L leg direction (fractional)
dT = T[1] - T[0]      # T leg direction (fractional)
# convert both to Cartesian and check they're perpendicular
dL_cart = dL @ B
dT_cart = dT @ B
cos = np.dot(dL_cart, dT_cart) / (np.linalg.norm(dL_cart)*np.linalg.norm(dT_cart))
print("cos(angle between legs) =", round(cos, 4))   # must be ~0
'''