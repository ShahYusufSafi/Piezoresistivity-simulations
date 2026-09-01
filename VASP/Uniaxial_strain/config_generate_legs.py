import numpy as np
from ase.io import read, write
import os, glob


from pymatgen.io.vasp import BSVasprun
from pymatgen.electronic_structure.core import Spin



EQ = './eps_+0.0000'
atoms0 = read(f'{EQ}/POSCAR')

d = 'mass_eps0'          # directory to write files
os.makedirs(d, exist_ok=True)


reset_dirs = False
if reset_dirs:
    for f in glob.glob(f'{d}/*'):
        os.remove(f)
write(f'{d}/POSCAR', atoms0, format='vasp')   # unstrained


os.system(f'cp {EQ}/INCAR.scf {d}')
os.system(f'cp {EQ}/INCAR.band {d}')
os.system(f'cp {EQ}/KPOINTS.mesh {d}')
os.system(f'cp {EQ}/POTCAR {d}')

delta   = 0.039                         #  sampling half-width in Ang^-1: how far we step off k0 to feel the curvature.

recip = atoms0.cell.reciprocal() * 2*np.pi # frac -> cart
B    = np.array(recip)          
Binv = np.linalg.inv(B)         # cart -> frac (we use this map to convert back to fracrtional coordinates)

# We first find the valley minimum
bs = BSVasprun(f'{EQ}/vasprun.band.xml').get_band_structure(
    kpoints_filename=f'{EQ}/kpoints.line'
    )
cbm_ind = bs.get_cbm()["band_index"][Spin.up][0]
E = bs.bands[Spin.up][cbm_ind] - bs.get_cbm()["energy"]

E_min_ind = np.argmin(E)

# now we find that index's kpoint in frac -> cartesian
k_frac = bs.kpoints[E_min_ind].frac_coords
k_cart = np.round(bs.structure.lattice.reciprocal_lattice.get_cartesian_coords(k_frac), 4) # in A^-1


# Let's ensure we are at CBM of Si (0.85 of the way between X and gamma)
a_lat = np.linalg.norm(atoms0.cell[0]) * np.sqrt(2)
k0 = np.linalg.norm(k_cart) / (2 * np.pi/a_lat)
assert abs(k0 - 0.85) < 1e-2, f"Not at Δ minimum: k0={k0:.4f}"

#k0_cart = np.array([0.0, 0.9768, 0.0])   # valley center, from uniaxial.ipynb
# We build 2 orthogonal vectors
par_dir = np.array([1.0, 0.0, 0.0])
perp_dir = np.array([0.0, 1.0, 0.0])

# Now we shift the k0 to the left-right with vector parallel to it, and up-down with vector prependicular to it (because k0 is in (0,1,0) direction) 
L_cart = [k_cart - delta*par_dir, k_cart + delta*par_dir]
T_cart = [k_cart - delta*perp_dir, k_cart + delta*perp_dir]


# now we convert every of those shifts back into fractional coordinate
L_frac = [p @ Binv for p in L_cart]   # -> fractional for KPOINTS
T_frac = [p @ Binv for p in T_cart]

# writing kpoints file
def line(a, b, l1, l2):
    return [f"{a[0]:.4f} {a[1]:.4f} {a[2]:.4f}  {l1}",
            f"{b[0]:.4f} {b[1]:.4f} {b[2]:.4f}  {l2}", ""]

npts = 30
out = ["mass legs: longitudinal then transverse", str(npts), "line", "reciprocal"]
out += line(L_frac[0], L_frac[1], "L-", "L+")   # along ŷ  -> m_l
out += line(T_frac[0], T_frac[1], "T-", "T+")   # along x̂  -> m_t
open(f'{d}/KPOINTS.line', 'w').write("\n".join(out).rstrip() + "\n")


'''
# tests
k0     = k0_cart @ Binv
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