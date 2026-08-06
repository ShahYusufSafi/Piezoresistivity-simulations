"""
uniaxial_model.py

The Gamma_12 (uniaxial / valley-repopulation) conductivity model, packaged for
the FEM layer.

Why this module exists
----------------------
The Herring-Vogt tensor is ANALYTIC in the strain: given (Xi_u, m_l, m_t, T)
it can be evaluated at any strain the FEM solver produces. Shipping FEM a table
of sigma sampled at a handful of strain points would force it to interpolate,
which throws away exactly the property that makes the handoff correct.

So this module ships the PARAMETERS plus the function, not sampled outputs.

Usage
-----
    from uniaxial_model import UniaxialModel

    model = UniaxialModel.from_csv("../VASP/Uniaxial/uniaxial_scalars.csv")

    sig = model.sigma_2d(eps_xx, eps_yy, mode="plane_stress")   # (...,2,2)
    rho = model.rho_2d(eps_xx, eps_yy, mode="beam")             # (...,2,2)

    # linear (small-strain) coefficients, if the solver wants them:
    model.pi11, model.pi12

    # and the strain beyond which the linear model is no longer trustworthy:
    model.eps_linear_5pct
"""

from dataclasses import dataclass
import numpy as np

KB = 8.617333262e-5  # eV/K


@dataclass
class UniaxialModel:
    """Valley-repopulation conductivity model with all first-principles inputs."""

    Xi_u: float          # uniaxial deformation potential [eV]
    m_l: float           # longitudinal effective mass [m_e]
    m_t: float           # transverse effective mass [m_e]
    C11: float           # elastic stiffness [GPa]
    C12: float           # elastic stiffness [GPa]
    T: float = 300.0     # temperature [K]

    # optional, carried through for reference / reporting
    pi11: float = None           # [Pa^-1]
    pi12: float = None           # [Pa^-1]
    eps_linear_5pct: float = None  # strain where linear model hits 5% error

    # ---------------- derived ----------------
    @property
    def beta(self):
        """Xi_u / kT -- the Boltzmann argument per unit strain [1/strain]."""
        return self.Xi_u / (KB * self.T)

    @property
    def inv_mc(self):
        """Isotropic conductivity mass 1/m_c = (1/3)(1/m_l + 2/m_t)."""
        return (1.0 / 3.0) * (1.0 / self.m_l + 2.0 / self.m_t)

    @property
    def nu(self):
        """Poisson ratio C12/(C11+C12)."""
        return self.C12 / (self.C11 + self.C12)

    @property
    def dS_Pa(self):
        """(S11 - S12) = 1/(C11 - C12), in Pa^-1."""
        return 1.0 / ((self.C11 - self.C12) * 1e9)

    # ---------------- the model ----------------
    def sigma_tensor(self, e_xx, e_yy=0.0, e_zz=0.0, sigma0=1.0):
        """
        Herring-Vogt conductivity tensor under diagonal strain.
        Analytic -> valid at any strain. Returns (...,3,3), normalised so
        sigma(0) = sigma0 * I.
        """
        e = np.stack(np.broadcast_arrays(e_xx, e_yy, e_zz), axis=-1).astype(float)
        dE = self.beta * e
        dE = dE - dE.min(axis=-1, keepdims=True)   # stability; cancels in the ratio
        w = np.exp(-dE)
        n = w / (2.0 * w.sum(axis=-1, keepdims=True))

        d = np.arange(3)
        inv_m = np.full(e.shape + (3,), 1.0 / self.m_t)
        inv_m[..., d, d] = 1.0 / self.m_l

        sig = 2.0 * np.einsum('...i,...ia->...a', n, inv_m)
        sig = sig / self.inv_mc * sigma0

        out = np.zeros(e.shape[:-1] + (3, 3))
        out[..., d, d] = sig
        return out

    def sigma_2d(self, e_xx, e_yy=0.0, sigma0=1.0, mode="plane_stress"):
        """
        In-plane 2x2 conductivity for a 2D FEM mesh.

        mode selects the out-of-plane closure:
          plane_stress : e_zz = -(C12/C11)(e_xx + e_yy)   thin film, free surface
          plane_strain : e_zz = 0                          thick / constrained
          beam         : e_yy = e_zz = -nu * e_xx          free uniaxial bar
        """
        e_xx = np.asarray(e_xx, float)
        e_yy = np.asarray(e_yy, float)
        if mode == "plane_stress":
            e_zz = -(self.C12 / self.C11) * (e_xx + e_yy)
        elif mode == "plane_strain":
            e_zz = np.zeros_like(e_xx)
        elif mode == "beam":
            e_yy, e_zz = -self.nu * e_xx, -self.nu * e_xx
        else:
            raise ValueError("mode must be plane_stress | plane_strain | beam")
        return self.sigma_tensor(e_xx, e_yy, e_zz, sigma0)[..., :2, :2]

    def rho_2d(self, e_xx, e_yy=0.0, sigma0=1.0, **kw):
        """Resistivity tensor (diagonal here, so just reciprocals)."""
        s = self.sigma_2d(e_xx, e_yy, sigma0, **kw)
        r = np.zeros_like(s)
        r[..., 0, 0] = 1.0 / s[..., 0, 0]
        r[..., 1, 1] = 1.0 / s[..., 1, 1]
        return r

    # ---------------- I/O ----------------
    @classmethod
    def from_csv(cls, path):
        """Build the model from the long-format scalar results file."""
        from results_io import load_result

        def opt(name):
            try:
                return load_result(name, path=path)
            except KeyError:
                return None

        return cls(
            Xi_u=load_result("Xi_u", path=path),
            m_l=load_result("m_l", path=path),
            m_t=load_result("m_t", path=path),
            C11=load_result("C11", path=path),
            C12=load_result("C12", path=path),
            T=load_result("T", path=path),
            pi11=opt("pi11"),
            pi12=opt("pi12"),
            eps_linear_5pct=opt("eps_linear_5pct"),
        )

    def summary(self):
        return (
            f"UniaxialModel(Xi_u={self.Xi_u:.3f} eV, m_l={self.m_l:.3f}, "
            f"m_t={self.m_t:.3f}, C11={self.C11:.1f} GPa, C12={self.C12:.1f} GPa, "
            f"T={self.T:.0f} K)\n"
            f"  beta = Xi_u/kT      = {self.beta:.1f} per unit strain\n"
            f"  1/m_c               = {self.inv_mc:.6f}\n"
            f"  Poisson ratio nu    = {self.nu:.4f}\n"
            f"  (S11-S12)           = {self.dS_Pa:.4e} Pa^-1"
        )
