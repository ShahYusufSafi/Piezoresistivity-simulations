import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve

def FEM_element_based(sigma_values, L, Nx):
    h = L / (Nx - 1)

    K_global = np.zeros((Nx, Nx))

    for e in range(Nx - 1):
        sigma = sigma_values[e]
        K_e = (sigma / h) * np.array([[1, -1], [-1, 1]])
        K_global[e:e+2, e:e+2] += K_e

    return K_global, h


def BC_Solver(K_Full, L_BC, R_BC):
    # Potential
    V = np.zeros(Nx)
    V[0] =  L_BC
    V[-1] = R_BC

    # Reduced matrix
    K = csr_matrix(K_Full[1:-1, 1:-1])

    # RHS: no source in the PDE, but boundary values contribute
    b = np.zeros(Nx - 2)
    b[0]  -= K_Full[1,0]   * V[0]
    b[-1] -= K_Full[-2,-1] * V[-1]

    # Solve
    V_interior = spsolve(K, b)
    V[1:-1] = V_interior
    return V


# Parameters
L = 1.0
Nx = 50
x = np.linspace(0, L, Nx)

sigma_values = 0.5 * np.ones(Nx) # The conductivity parameter

# Assemble full matrix
K_full, h = FEM_element_based(sigma_values, L, Nx)

# Incorporate BC, and solve
V = BC_Solver(K_full, 0.0, 1.0)

# Derived quantities
E = -np.gradient(V, x)
J = sigma_values * E

print("V min/max:", V.min(), V.max())
print("Mean electric field:", np.mean(E))
print("Mean current density:", np.mean(J))

# Exact solution
V_exact = x / L

# Plot
plt.figure(figsize=(8, 4))
plt.plot(x, V, label='FEM solution')
plt.plot(x, V_exact, '--', label='Exact solution')
plt.plot(x, J, label='Current')
plt.plot(x, E, label='Electric field')
plt.xlabel('Position x')
plt.ylabel('Potential V')
plt.title('1D Steady Electrical Conduction')
plt.grid(True)
plt.legend()
plt.savefig("./FEM/P1d.png", dpi=200)