import numpy as np
import sys, os
from pymatgen.io.vasp import BSVasprun
from pymatgen.electronic_structure.core import Spin
from ase.io import read

# we first find the position of k0
#k0 = np.array([0.0, 0.0, 0.0]) # a guess , we will fix it once grid is working

def grid(k0, h = 0.05):

    # basis vectors
    hx = np.array([1.0 ,0.0, 0.0]) * h
    hy = np.array([0.0, 1.0, 0.0]) * h

    # Now we generate few points around k0 with such transformation
    points = []
    for x in range(-1, 2):
        for y in range(-1, 2):                   
            px =  hx * x  
            py =  hy * y  
            p  =  px + py + k0
            points.append(p)

    return np.array(points) 

#pnts= grid(k0)
#print(pnts.shape)

# now I find the actual k0 coordinate from eps_+0.0000
bs = BSVasprun('./eps_+0.0000/vasprun.xml').get_band_structure('./eps_+0.0000/KPOINTS.line', line_mode=True)

cbm1 = bs.get_cbm()
bi1  = cbm1["band_index"][Spin.up][0]
E1   = bs.bands[Spin.up][bi1] - cbm1["energy"]

imin1 = int(np.argmin(E1))



# Original minimum with no fitting
#kfrac_org = bs.kpoints[imin1].frac_coords
#kcart1 = bs.structure.lattice.reciprocal_lattice.get_cartesian_coords(kfrac_org)
#
#grid_pnts = grid(kcart1)
# end of original minimum


# One Newton step 
# instead I use following for matching the longitudinal leg better (but note we have to deal with strain swepps separately)
# Cartesian coords of every line k-point
kcart_line = np.array([
    bs.structure.lattice.reciprocal_lattice.get_cartesian_coords(kp.frac_coords)
    for kp in bs.kpoints
])


# On the X_par->Gamma segment the line is (t,0,0)_cart, so only kx varies.
w  = 3                                  # window half-width (points each side)
sl = slice(max(imin1-w, 0), imin1+w+1)
kx = kcart_line[sl, 0]
Ew = E1[sl]

c        = np.polyfit(kx, Ew, 2)        # E ~ c2 kx^2 + c1 kx + c0
kx_star  = -c[1] / (2*c[0])             # vertex
k0_cart  = np.array([kx_star, 0.0, 0.0])

grid_pnts = grid(k0_cart)               # <- center on the vertex, not argmin
# end Newtonian fit

# Now I convert my cartesian=s back to fractional, since all of my files are in fractional coordinat
#atom0 = read('./eps_+0.0000/POSCAR')
#recip = atom0.cell.reciprocal() * 2*np.pi
#frac = np.linalg.inv(np.array(recip))
#
#k_frac = kcart1 @ frac

# or simply as bellow
k_frac = np.array([bs.structure.lattice.reciprocal_lattice.get_fractional_coords(p) for p in grid_pnts])

k0_frac = bs.structure.lattice.reciprocal_lattice.get_fractional_coords(k0_cart)

# A quick check that fractional coordinated match the actual ones in band structure
assert np.allclose(k_frac.mean(axis=0), k0_frac, atol=1e-15), "wrong conversion"
print("The plane is on CBM and ready to be used as kpoints coordinates")

# now we run a non-SCF VASP, with all previous files except kpoint
DIR = './eps_+0.0000'

new_DIR = f'{DIR}/mass_val'

os.makedirs(new_DIR, exist_ok=True)

# copy files
os.system(f'cp {DIR}/POSCAR {new_DIR}')
os.system(f'cp {DIR}/POTCAR {new_DIR}')
os.system(f'cp {DIR}/INCAR.band {new_DIR}')
os.system(f'cp {DIR}/CHGCAR {new_DIR}')

# create KPOINT
lines = "Explicit k-points for mass patch\n9\nReciprocal\n"
for k in k_frac:
    lines += f"{k[0]:.10f} {k[1]:.10f} {k[2]:.10f} 1\n"
open(f'{new_DIR}/KPOINTS.grid', 'w').write(lines)