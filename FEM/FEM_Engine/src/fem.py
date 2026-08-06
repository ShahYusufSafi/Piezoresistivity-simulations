"""
fem_aniso.py

FEM assembly and post-processing for ANISOTROPIC conductivity.

Why this exists
---------------
A scalar-sigma solver is correct for the HYDROSTATIC (Gamma_1) channel, where
strain shifts the gap and sigma stays a scalar. It cannot represent the UNIAXIAL
(Gamma_12) channel, whose entire physical content is that sigma becomes
DIRECTIONAL: valley repopulation makes sigma_xx != sigma_yy. Every routine here
takes a (2,2) tensor per element, and reduces exactly to the isotropic case when
sigma = sigma0 * I.

Everything is vectorised over elements -- no Python loop over the mesh -- so the
same code stays usable on a refined mesh.

Layout
------
  mesh          : make_mesh
  kernel        : element_geometry            (shared by assembly and fields)
  assembly      : element_stiffness_2d, assemble_2d, assemble_1d
  boundaries    : boundary_nodes, solve_dirichlet, solve_contacts
  post-process  : element_fields_2d, contact_currents, cross_section_flux,
                  device_resistance

Conventions
-----------
  sigma_elements accepts, for a mesh of Ne elements:
      scalar               -> isotropic, same everywhere
      (2,2)                -> anisotropic, same everywhere
      (Ne,)                -> isotropic, per element
      (Ne,2,2)             -> anisotropic, per element
  V is nodal; E and J are element-wise constant (P1).
"""

import numpy as np
from scipy.sparse import coo_matrix, csr_matrix
from scipy.sparse.linalg import spsolve


# =====================================================================
#  Mesh
# =====================================================================
def make_mesh(L, H, nx, ny, diagonal="right"):
    """
    Structured triangulation of the rectangle [0,L] x [0,H].

    diagonal : "right", "left", or "cross" -- how each quad is split.
               "cross" alternates, removing the directional bias that a
               uniform diagonal introduces into anisotropic problems.

    Returns (points (N,2), triangles (Ne,3)).
    """
    X, Y = np.meshgrid(np.linspace(0.0, L, nx), np.linspace(0.0, H, ny))
    points = np.column_stack([X.ravel(), Y.ravel()])

    tris = []
    for i in range(ny - 1):
        for j in range(nx - 1):
            v0 = i * nx + j
            v1 = (i + 1) * nx + j
            v2 = i * nx + j + 1
            v3 = (i + 1) * nx + j + 1
            if diagonal == "right" or (diagonal == "cross" and (i + j) % 2 == 0):
                tris += [[v0, v1, v2], [v1, v3, v2]]
            else:
                tris += [[v0, v1, v3], [v0, v3, v2]]
    return points, np.array(tris, dtype=int)


# =====================================================================
#  Shared kernel: P1 shape-function gradients and areas
# =====================================================================
def element_geometry(points, triangles):
    """
    Constant shape-function gradients, areas and centroids for every P1 triangle.

    Returns
    -------
    grad      : (Ne,3,2)   grad[e,i] = grad(phi_i) on element e
    area      : (Ne,)
    centroids : (Ne,2)

    This is the single place element geometry is computed. Assembly and
    post-processing both call it, so they cannot disagree about the mesh.
    """
    points = np.asarray(points, float)
    triangles = np.asarray(triangles, int)

    p = points[triangles]                                  # (Ne,3,2)
    x1, y1 = p[:, 0, 0], p[:, 0, 1]
    x2, y2 = p[:, 1, 0], p[:, 1, 1]
    x3, y3 = p[:, 2, 0], p[:, 2, 1]

    detJ = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1)   # (Ne,) = 2 * signed area
    if np.any(np.abs(detJ) < 1e-300):
        bad = int(np.argmin(np.abs(detJ)))
        raise ValueError(f"degenerate triangle at index {bad} (zero area)")

    b = np.stack([y2 - y3, y3 - y1, y1 - y2], axis=1)      # (Ne,3)
    c = np.stack([x3 - x2, x1 - x3, x2 - x1], axis=1)      # (Ne,3)
    grad = np.stack([b, c], axis=2) / detJ[:, None, None]  # (Ne,3,2)

    return grad, np.abs(detJ) / 2.0, p.mean(axis=1)


def _as_tensor_field(sigma_elements, Ne):
    """Broadcast any accepted sigma specification to (Ne,2,2)."""
    sig = np.asarray(sigma_elements, float)
    eye = np.eye(2)[None, :, :]
    if sig.ndim == 0:
        return np.broadcast_to(sig * eye, (Ne, 2, 2))
    if sig.shape == (2, 2):
        return np.broadcast_to(sig[None, :, :], (Ne, 2, 2))
    if sig.ndim == 1:
        if len(sig) != Ne:
            raise ValueError(f"sigma has {len(sig)} entries but mesh has {Ne} elements")
        return sig[:, None, None] * eye
    if sig.shape == (Ne, 2, 2):
        return sig
    raise ValueError(
        f"sigma_elements must be a scalar, (2,2), ({Ne},) or ({Ne},2,2); "
        f"got shape {sig.shape}"
    )


# =====================================================================
#  Assembly
# =====================================================================
def element_stiffness_2d(p1, p2, p3, sigma):
    """
    Single-element stiffness,  K_ij = Area * grad(phi_i) . sigma . grad(phi_j).

    Standalone entry point for testing and one-off elements. assemble_2d does
    not call it -- it uses the vectorised path. Reduces exactly to the standard
    isotropic stiffness when sigma = sigma0 * I.
    """
    pts = np.array([p1, p2, p3], float)
    grad, area, _ = element_geometry(pts, np.array([[0, 1, 2]]))
    sig = _as_tensor_field(sigma, 1)
    return float(area[0]) * grad[0] @ sig[0] @ grad[0].T


def assemble_2d(points, triangles, sigma_elements):
    """
    Global conduction matrix, vectorised over elements.

        K_ij = sum_e  Area_e * grad(phi_i) . sigma_e . grad(phi_j)
    """
    triangles = np.asarray(triangles, int)
    Ne, N = len(triangles), len(points)

    grad, area, _ = element_geometry(points, triangles)
    sig = _as_tensor_field(sigma_elements, Ne)

    # all (3,3) local matrices at once
    Ke = area[:, None, None] * np.einsum('eia,eab,ejb->eij', grad, sig, grad)

    rows = np.broadcast_to(triangles[:, :, None], (Ne, 3, 3)).ravel()
    cols = np.broadcast_to(triangles[:, None, :], (Ne, 3, 3)).ravel()
    # coo sums duplicate (row,col) entries on conversion -- that IS the assembly
    return coo_matrix((Ke.ravel(), (rows, cols)), shape=(N, N)).tocsr()


def assemble_1d(sigma_xx, L, Nx):
    """
    1D conduction matrix. Only the xx component is meaningful in 1D --
    pass sigma[...,0,0] if you are holding a tensor field.

    Returns (K (Nx,Nx) dense, h).
    """
    sigma_xx = np.asarray(sigma_xx, float).ravel()
    if len(sigma_xx) != Nx - 1:
        raise ValueError(f"need {Nx-1} element values, got {len(sigma_xx)}")
    h = L / (Nx - 1)
    K = np.zeros((Nx, Nx))
    i = np.arange(Nx - 1)
    np.add.at(K, (i, i),         sigma_xx / h)
    np.add.at(K, (i + 1, i + 1), sigma_xx / h)
    np.add.at(K, (i, i + 1),    -sigma_xx / h)
    np.add.at(K, (i + 1, i),    -sigma_xx / h)
    return K, h


# =====================================================================
#  Boundaries
# =====================================================================
def boundary_nodes(points, axis=0, value=0.0, tol=1e-10):
    """Indices of nodes with points[:, axis] == value, within tol."""
    points = np.asarray(points, float)
    return np.where(np.abs(points[:, axis] - value) < tol)[0]


def solve_dirichlet(K, fixed_idx, fixed_val, n_nodes=None):
    """
    Solve K V = 0 with the listed nodes held at prescribed potentials.
    Standard lifting: the known columns move to the right-hand side.
    """
    K = csr_matrix(K)
    N = n_nodes if n_nodes is not None else K.shape[0]

    fixed_idx = np.asarray(fixed_idx, int)
    fixed_val = np.broadcast_to(np.asarray(fixed_val, float), fixed_idx.shape)

    V = np.zeros(N)
    V[fixed_idx] = fixed_val

    free = np.setdiff1d(np.arange(N), fixed_idx)
    rhs = -np.asarray(K[free][:, fixed_idx] @ fixed_val).ravel()
    V[free] = spsolve(csr_matrix(K[free][:, free]), rhs)
    return V


def solve_contacts(K, points, L, V_left=1.0, V_right=0.0, tol=1e-10):
    """
    Dirichlet contacts at x=0 and x=L, natural (zero-flux) boundaries elsewhere.
    Returns (V, left_nodes, right_nodes).
    """
    left = boundary_nodes(points, axis=0, value=0.0, tol=tol)
    right = boundary_nodes(points, axis=0, value=L, tol=tol)
    if len(left) == 0 or len(right) == 0:
        raise ValueError("no nodes on one of the contacts -- check L and tol")
    idx = np.concatenate([left, right])
    val = np.concatenate([np.full(len(left), V_left), np.full(len(right), V_right)])
    return solve_dirichlet(K, idx, val, n_nodes=len(points)), left, right


# =====================================================================
#  Post-processing
# =====================================================================
def element_fields_2d(points, triangles, V, sigma_elements):
    """
    Element-wise fields for a P1 solution (constant on each element):

        E = -grad V,      J = sigma . E

    Returns (E (Ne,2), J (Ne,2), centroids (Ne,2)).
    """
    triangles = np.asarray(triangles, int)
    Ne = len(triangles)

    grad, _, centroids = element_geometry(points, triangles)
    sig = _as_tensor_field(sigma_elements, Ne)

    E = -np.einsum('eia,ei->ea', grad, np.asarray(V, float)[triangles])
    J = np.einsum('eab,eb->ea', sig, E)
    return E, J, centroids


def contact_currents(K, V, left, right):
    """
    Total current through each contact, from the discrete reaction fluxes
    r = K V summed over the constrained nodes.

    EXACT for the discrete problem and conservative by construction:
    I_left + I_right = 0 to round-off on any mesh, however coarse. Prefer this
    over integrating J across a cross-section, which is only as good as the
    quadrature used.

    Returns (I_left, I_right) -- equal and opposite.
    """
    r = csr_matrix(K) @ np.asarray(V, float)
    return (float(r[np.asarray(left, int)].sum()),
            float(r[np.asarray(right, int)].sum()))


def cross_section_flux(centroids, J, x_positions=None, tol=1e-9):
    """
    Uniformity diagnostic: sum of J_x over the elements at each x-station.

    NOT a current. It sums current DENSITIES without weighting by element
    extent, so its absolute value scales with how many elements sit in a
    cross-section (on a 12-row mesh it comes out 12x the true current). Use
    contact_currents for anything quantitative.

    What it is good for is locating trouble: the RELATIVE spread across
    stations should be at round-off, and if it isn't, this shows you WHERE
    conservation degrades rather than merely that it did.

    Returns (x_positions, flux_sum).
    """
    centroids = np.asarray(centroids, float)
    if x_positions is None:
        x_positions = np.unique(np.round(centroids[:, 0], 9))
    flux = np.array([J[np.abs(centroids[:, 0] - xc) < tol, 0].sum()
                     for xc in x_positions])
    return x_positions, flux


def device_resistance(K, V, left, right, V_applied=None):
    """R = V_applied / |I|, using the exact contact current."""
    I_left, _ = contact_currents(K, V, left, right)
    if V_applied is None:
        Va = np.asarray(V, float)
        V_applied = abs(Va[np.asarray(left, int)][0] - Va[np.asarray(right, int)][0])
    return V_applied / abs(I_left)


# ---------------------------------------------------------------------
#  Backward-compatible alias
# ---------------------------------------------------------------------
def element_current_2d(points, triangles, sigma_elements, V):
    """Deprecated: use element_fields_2d, which also returns E and centroids."""
    _, J, _ = element_fields_2d(points, triangles, V, sigma_elements)
    return J