import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve

def FEM_element_based(sigma_values, L, Nx):
    """Assemble the global stiffness matrix from element contributions.

    The 1D mesh has `Nx` nodes and therefore `Nx - 1` elements. Element `e`
    connects local nodes `e` and `e + 1`; its local matrix is scattered into
    the corresponding global rows/columns.

    Parameters
    ----------
    sigma_values : np.ndarray, shape (Nx - 1,)
        Per-element material coefficient.
    L : float
        Total domain length.
    Nx : int
        Number of nodes (=> Nx - 1 elements).

    Returns
    -------
    K_global : np.ndarray, shape (Nx, Nx)
        Assembled global stiffness matrix.
    h : float
        Uniform element length.
    """
    h = L / (Nx - 1)

    K_global = np.zeros((Nx, Nx))

    for e in range(Nx - 1):
        sigma = sigma_values[e]
        K_e = (sigma / h) * np.array([[1, -1], [-1, 1]])
        K_global[e:e+2, e:e+2] += K_e

    return K_global, h

def BC_Solver(K_Full, L_BC, R_BC):
    """Solve the FEM system K*V = b with Dirichlet boundary conditions.

    Applies the boundary conditions by partitioning the global system:
    the boundary DOFs are fixed (known), and their contributions are moved
    to the RHS before solving for the interior DOFs only.

    Parameters
    ----------
    K_Full : np.ndarray, shape (Nx, Nx)
        Full (unconstrained) global stiffness matrix.
    L_BC : float
        Dirichlet value imposed at the left boundary node (index 0).
    R_BC : float
        Dirichlet value imposed at the right boundary node (index -1).

    Returns
    -------
    V : np.ndarray, shape (Nx,)
        Solution vector with boundary values enforced and interior solved.
    """
    V = np.zeros(Nx)
    V[0]  = L_BC
    V[-1] = R_BC

    # use CSR for efficient sparse solve (stiffness matrix is tridiagonal)
    K = csr_matrix(K_Full[1:-1, 1:-1])

    # RHS: no volumetric source, but fixed boundary values enter via
    # b_i = -K_{i,0}*V[0] for the left and -K_{i,-1}*V[-1] for the right;
    # only the first and last interior nodes couple to the boundary
    b = np.zeros(Nx - 2)
    b[0]  -= K_Full[1,  0 ] * V[0]    # left  boundary contribution
    b[-1] -= K_Full[-2, -1] * V[-1]   # right boundary contribution

    # Solve the reduced sparse system K_interior * V_interior = b
    V_interior = spsolve(K, b)
    V[1:-1] = V_interior               # insert interior solution back
    return V

# Parameters
L = 1.0
Nx = 50
x = np.linspace(0, L, Nx)

# Create a Generator with a fixed seed for reproducibility
rng = np.random.default_rng(seed=42)

#sigma_values = 0.5 * np.ones(Nx) # The conductivity parameter

sigma_values = rng.random((Nx-1,))
# Assemble full matrix
K_full, h = FEM_element_based(sigma_values, L, Nx)

# Incorporate BC, and solve
V = BC_Solver(K_full, 0.0, 1.0)

# ===POSTPROCESSING ====
x_mid = 0.5 * (x[:-1] + x[1:])          # element midpoints
# Derived quantities (BELLOW ARE FOR NODES)
E = -np.diff(V) / h # per elements
J = sigma_values * E # per elements

print("V min/max:", V.min(), V.max())
print("Mean electric field:", np.mean(E))
print("Mean current density:", np.mean(J))

# Exact solution
V_exact = x / L

# Plot
plt.figure(figsize=(8, 4))
plt.plot(x, V, '-o', markersize=2,label='FEM solution', )
plt.plot(x, V_exact, '--', label='Exact solution')
plt.plot(x_mid, J, label='Current')
plt.plot(x_mid, E, label='Electric field')
plt.step(x_mid, E, where='mid', label='E field')
plt.step(x_mid, J.mean() / sigma_values, where='mid', label='J/σ (predicted E)')
# these two should lie exactly on top of each other
plt.xlabel('Position x')
plt.ylabel('Potential V')
plt.title('1D Steady Electrical Conduction')
plt.grid(True)
plt.legend()
plt.savefig("./P1d.png", dpi=200)