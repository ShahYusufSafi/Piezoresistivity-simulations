import numpy as np
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import spsolve

def FEM_element_based(sigma, L, Nx):
    h = L / (Nx - 1)

    K_global = np.zeros((Nx, Nx))
    K_e = (sigma / h) * np.array([[1, -1], [-1, 1]])

    for e in range(Nx - 1):
        K_global[e:e+2, e:e+2] += K_e

    return K_global, h

# Parameters
L = 1.0
Nx = 50
sigma = 0.5

x = np.linspace(0, L, Nx)

# Potential
V = np.zeros(Nx)
V[0] = 0.0
V[-1] = 1.0

# Assemble full matrix
K_full, h = FEM_element_based(sigma, L, Nx)

# Reduced matrix
K = csr_matrix(K_full[1:-1, 1:-1])

# RHS: no source in the PDE, but boundary values contribute
b = np.zeros(Nx - 2)
b[0]  -= K_full[1,0]   * V[0]
b[-1] -= K_full[-2,-1] * V[-1]


# Solve
V_interior = spsolve(K, b)
V[1:-1] = V_interior

# Derived quantities
E = -np.gradient(V, x)
J = sigma * E

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


