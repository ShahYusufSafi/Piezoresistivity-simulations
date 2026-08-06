"""
fem_aniso.py

FEM assembly for ANISOTROPIC conductivity.

Why this is needed
------------------
The existing solvers in fem.ipynb take a scalar sigma per element:

    K_e = (sigma_e / h) * [[1,-1],[-1,1]]                     (1D)
    K_e = gamma1*S1 + gamma2*S2 + gamma3*S3                   (2D, isotropic)

That is correct for the HYDROSTATIC (Gamma_1) channel, where strain shifts the
gap and sigma stays a scalar. It cannot represent the UNIAXIAL (Gamma_12)
channel, whose entire physical content is that sigma becomes DIRECTIONAL:
valley repopulation makes sigma_xx != sigma_yy. Feeding a scalar into the
solver discards exactly the effect the channel exists to describe.

This module replaces the element matrices with the tensor form

    K_ij = Area * (grad phi_i)^T . sigma . (grad phi_j)

which reduces exactly to the isotropic case when sigma = sigma0 * I
(verified against the original compute_element_matrix to ~1e-14).

Usage
-----
    from fem_aniso import assemble_2d, solve_dirichlet, element_stiffness_2d

    sig_e = model.sigma_2d(eps_xx_elem, 0.0, sigma0=s0, mode="beam")  # (Ne,2,2)
    K = assemble_2d(points, triangles, sig_e)
    V = solve_dirichlet(K, fixed_idx, fixed_val, n_nodes=len(points))
"""

import numpy as np
from scipy.sparse import lil_matrix, csr_matrix
from scipy.sparse.linalg import spsolve


# =====================================================================
#  1D:  linear elements, sigma_xx only
# =====================================================================
def assemble_1d(sigma_xx, L, Nx):
    """
    1D conduction matrix. Only the xx component of the tensor is meaningful
    in 1D -- pass sigma[...,0,0] from the uniaxial model.

    Returns (K_global (Nx,Nx) dense, h).
    """
    
    
    sigma_xx = np.asarray(sigma_xx, float).ravel()
    if len(sigma_xx) != Nx - 1:
            raise ValueError(f"need {Nx-1} element values, got {len(sigma_xx)}")
    h = L / (Nx - 1)
    K = np.zeros((Nx, Nx))
    for e in range(Nx - 1):
        Ke = (sigma_xx[e] / h) * np.array([[1.0, -1.0], [-1.0, 1.0]])
        K[e:e + 2, e:e + 2] += Ke
    return K, h


# =====================================================================
#  Dirichlet solve
# =====================================================================
def solve_dirichlet(K, fixed_idx, fixed_val, n_nodes=None):
    """
    Solve K V = 0 with Dirichlet nodes fixed.

    K         : (N,N) sparse or dense
    fixed_idx : indices of constrained nodes
    fixed_val : their prescribed potentials
    """
    K = csr_matrix(K)
    N = n_nodes if n_nodes is not None else K.shape[0]

    fixed_idx = np.asarray(fixed_idx, int)
    fixed_val = np.asarray(fixed_val, float)

    V = np.zeros(N)
    V[fixed_idx] = fixed_val

    free = np.setdiff1d(np.arange(N), fixed_idx)
    # move the known columns to the RHS
    b = -np.asarray(K[free][:, fixed_idx] @ fixed_val).ravel()
    V[free] = spsolve(csr_matrix(K[free][:, free]), b)
    return V


# =====================================================================
#  2D:  P1 triangles, anisotropic sigma
# =====================================================================
def make_mesh(L, H, nx, ny):
    x = np.linspace(0.0, L, nx)
    y = np.linspace(0.0, H, ny)
    X, Y = np.meshgrid(x, y)
    points = np.column_stack([X.ravel(), Y.ravel()])

    triangles = []
    for i in range(ny - 1):
        for j in range(nx - 1):
            v0 = i * nx + j
            v1 = (i + 1) * nx + j
            v2 = i * nx + j + 1
            v3 = (i + 1) * nx + j + 1
            triangles.append([v0, v1, v2])
            triangles.append([v1, v3, v2])
    return points, np.array(triangles)




def element_stiffness_2d(p1, p2, p3, sigma):
    """
    P1 triangle stiffness for an anisotropic conductivity tensor.

    p1, p2, p3 : (2,) vertex coordinates
    sigma      : (2,2) conductivity tensor (or scalar for the isotropic case)

    Returns (3,3). Reduces to the standard isotropic stiffness when
    sigma = sigma0 * I.
    """
    sigma = np.asarray(sigma, float)
    if sigma.ndim == 0:
        sigma = sigma * np.eye(2)

    (x1, y1), (x2, y2), (x3, y3) = p1, p2, p3
    detJ = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1)
    if abs(detJ) < 1e-300:
        raise ValueError("degenerate triangle (zero area)")
    area = abs(detJ) / 2.0

    # constant shape-function gradients on a P1 triangle
    b = np.array([y2 - y3, y3 - y1, y1 - y2])
    c = np.array([x3 - x2, x1 - x3, x2 - x1])
    grad = np.column_stack([b, c]) / detJ          # (3,2), row i = grad phi_i

    return area * grad @ sigma @ grad.T


def assemble_2d(points, triangles, sigma_elements):
    """
    Assemble the global conduction matrix.

    points         : (N,2)
    triangles      : (Ne,3) node indices
    sigma_elements : (Ne,2,2) tensor per element, or (Ne,) scalars,
                     or a single (2,2)/scalar applied to every element.
    """
    points = np.asarray(points, float)
    triangles = np.asarray(triangles, int)
    Ne = len(triangles)

    sig = np.asarray(sigma_elements, float)
    if sig.ndim == 0:
        sig = np.broadcast_to(sig * np.eye(2), (Ne, 2, 2))
    elif sig.shape == (2, 2):
        sig = np.broadcast_to(sig, (Ne, 2, 2))
    elif sig.ndim == 1:
        sig = sig[:, None, None] * np.eye(2)
    if sig.shape != (Ne, 2, 2):
        raise ValueError(f"sigma_elements must broadcast to ({Ne},2,2), got {sig.shape}")

    K = lil_matrix((len(points), len(points)))
    for e, tri in enumerate(triangles):
        i1, i2, i3 = tri
        Ke = element_stiffness_2d(points[i1], points[i2], points[i3], sig[e])
        idx = [i1, i2, i3]
        for a, ia in enumerate(idx):
            for b_, ib in enumerate(idx):
                K[ia, ib] += Ke[a, b_]
    return csr_matrix(K)


def element_current_2d(points, triangles, sigma_elements, V):
    """
    Per-element current density J = -sigma . grad V.  Returns (Ne,2).
    Useful as the conservation check: the x-flux through every cross-section
    of a slab must be equal.
    """
    points = np.asarray(points, float)
    triangles = np.asarray(triangles, int)
    Ne = len(triangles)

    sig = np.asarray(sigma_elements, float)
    if sig.ndim == 1:
        sig = sig[:, None, None] * np.eye(2)
    elif sig.shape == (2, 2):
        sig = np.broadcast_to(sig, (Ne, 2, 2))

    J = np.zeros((Ne, 2))
    for e, tri in enumerate(triangles):
        i1, i2, i3 = tri
        (x1, y1), (x2, y2), (x3, y3) = points[i1], points[i2], points[i3]
        detJ = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1)
        b = np.array([y2 - y3, y3 - y1, y1 - y2])
        c = np.array([x3 - x2, x1 - x3, x2 - x1])
        grad = np.column_stack([b, c]) / detJ
        gradV = grad.T @ V[[i1, i2, i3]]
        J[e] = -sig[e] @ gradV
    return J
