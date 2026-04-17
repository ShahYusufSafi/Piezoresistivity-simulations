from ase.build import bulk
from ase.visualize import view
import matplotlib.pyplot as plt
from ase.visualize.plot import plot_atoms
import nglview as ngv


# Create silicon crystal (diamond cubic structure)
si = bulk('Si', 'diamond', a=5.43)

# Repeat to make it bigger (more visible)
si = si.repeat((3, 3, 3))

# Visualize
#view(si)

view = ngv.show_ase(si)
view

#plot_atoms(si)
#plt.savefig("./Paticles Simulation/structure.png")