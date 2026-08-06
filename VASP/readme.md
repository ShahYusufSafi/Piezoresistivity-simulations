# DFT Strain Sweeps with VASP

This directory computes how silicon's electronic structure responds to strain,
using density-functional theory (DFT). The output is a set of **deformation
potentials** — the slopes that tell us how the band edges move and split as the
crystal is deformed. Those slopes are what turn the conductivity in the FEM
model from a fitted guess into a physics-derived quantity $\sigma(\varepsilon)$.

The logic of the whole directory is one sentence: **strain changes where the
atoms sit, that changes the potential the electrons feel, that changes the band
energies, and the change in the gap is the deformation potential.** Everything
below is the machinery that makes each link in that chain trustworthy.

---

## 1. What VASP actually computes

DFT finds the ground-state electron density $n(\mathbf{r})$, written as a sum
over the occupied Kohn–Sham (KS) orbitals $\psi_i$:

$$n(\mathbf{r}) = \sum_{i=1}^{N} |\psi_i(\mathbf{r})|^2 .$$

The orbitals come from a single-particle Schrödinger-like equation with an
effective potential:

$$\hat{H}_\text{KS} = -\frac{\hbar^2}{2m}\nabla^2 + V_\text{eff}(\mathbf{r}),
\qquad
V_\text{eff} = V_\text{ext} + V_H + V_\text{xc}.$$

The three pieces:

- **External potential $V_\text{ext}(\mathbf{r})$** — the attraction the valence
  electrons feel from the positive ion cores. In VASP the cores are PAW
  pseudopotentials (nucleus + frozen core electrons), not bare nuclei.
- **Hartree term $V_H(\mathbf{r})$** — the classical electrostatic repulsion of
  an electron against the whole electron cloud:
  $$V_H(\mathbf{r}) = e^2 \int \frac{n(\mathbf{r}')}{|\mathbf{r}-\mathbf{r}'|}\,d^3r'.$$
- **Exchange–correlation $V_\text{xc}(\mathbf{r})$** — all the many-body quantum
  effects beyond Hartree, approximated by a functional (we use PBE, a GGA):
  $$V_\text{xc}(\mathbf{r}) = \frac{\delta E_\text{xc}[n]}{\delta n(\mathbf{r})}.$$

### Where the crystal structure actually lives

This is the part that matters for strain, and it is easy to get wrong.
$V_\text{ext}$ is **not** stored in the POTCAR. POTCAR only gives the *shape* of
one ion's pseudopotential, per chemical species. The actual potential is a sum
of those shapes, one copy centred on every ion in the (infinitely repeated)
crystal:

$$V_\text{ext}(\mathbf{r}) = \sum_{I}\sum_{\mathbf{R}}
  v^{\text{PP}}_{s(I)}\!\left(\mathbf{r}-\boldsymbol{\tau}_I-\mathbf{R}\right).$$

Read against the input files:

- $v^{\text{PP}}_{s}$ — the per-species pseudopotential shape → **POTCAR**.
- $\boldsymbol{\tau}_I$ — the atomic positions → **POSCAR**.
- $\mathbf{R}=n_1\mathbf{a}_1+n_2\mathbf{a}_2+n_3\mathbf{a}_3$ — the lattice
  tiling → **POSCAR** (the three lattice vectors).

So POTCAR is the stencil and POSCAR decides where the stencil is stamped and how
space is tiled. **Strain only ever enters the KS problem through $V_\text{ext}$**,
because strain only moves $\boldsymbol{\tau}_I$ and $\mathbf{a}_i$. $V_H$ and
$V_\text{xc}$ are not where strain lives — they only respond afterwards, because
the density relaxes to the new $V_\text{ext}$.

### Self-consistent field (SCF) cycle

VASP cannot build $V_\text{eff}$ without knowing $n$, and cannot get $n$ without
$V_\text{eff}$, so it iterates:

1. guess $n(\mathbf{r})$,
2. build $V_\text{eff}$,
3. solve $\hat{H}_\text{KS}\psi_i = \varepsilon_i\psi_i$,
4. rebuild $n$ from the new $\psi_i$,
5. repeat until the energy and density stop changing.

The converged eigenvalues $\varepsilon_{n\mathbf{k}}$ are the band structure.

---

## 2. From bands to conductivity

The band energies give the density of states,

$$g(E) = \sum_{n,\mathbf{k}} \delta(E - \varepsilon_{n\mathbf{k}}),$$

and the carrier density follows from filling those states up to the Fermi level,
$n \propto \int g(E)\,f(E)\,dE$ with $f$ the Fermi–Dirac distribution. When
strain shifts the band edges, $g(E)$ and hence the carrier density change, and
that is the microscopic origin of $\sigma(\varepsilon)$:

$$\Delta E_g(\varepsilon) \approx (a_c + a_v)\,\varepsilon_h
\;\Longrightarrow\;
\sigma(\varepsilon)\ \text{from first principles},$$

replacing the phenomenological $\sigma_0(1 + \pi\varepsilon + \alpha\varepsilon^2)$.
The coefficients $a_c+a_v$ and their shear counterparts are exactly the
deformation potentials this directory extracts.

---

## 3. The three strain channels

A strain is applied as a deformation gradient $F = \mathbb{1}+\varepsilon$ acting
on the lattice vectors (`cell = cell0 @ M` in the generator). Cubic symmetry
splits any strain into three independent pieces, and each one isolates a
different deformation potential:

| Channel | Matrix `M` | What it does to the box | Probes |
|---|---|---|---|
| Hydrostatic ($\Gamma_1$) | `diag(1+ε, 1+ε, 1+ε)` | uniform scaling, volume changes | $a = a_c + a_v$ |
| Uniaxial ($\Gamma_{12}$) | `diag(1+ε, 1, 1)` | one axis stretched, becomes tetragonal | $\Xi_u$ (valley splitting) |
| Shear ($\Gamma_{25'}$) | `I + ε·[[0,1,1],[1,0,1],[1,1,0]]` | cell **angles** tilt, volume fixed | $d$ (valence splitting) |

Why they look different in the band structure:

- **Hydrostatic** keeps every direction equivalent, so no degeneracy lifts — the
  gap just shifts a little. This is the expected "boring" result and is itself a
  validation (it confirms the strain carries no symmetry-breaking part).
- **Uniaxial** makes the strained axis different from the other two, splitting
  the six conduction $\Delta$ valleys into a 2-fold set (along the axis) and a
  4-fold set (transverse).
- **Shear** tilts the bonds unequally and splits the triply degenerate valence
  top at $\Gamma$. This is the clearest degeneracy lift to see, because it sits
  right at $\Gamma$.

---

## 4. The calculation workflow — and why each piece exists

The original setup copied one `INCAR` and one `KPOINTS` into every directory.
That is correct only for the simplest case. Here is what changed and the reason
behind each change.

### 4.1 Three stages per directory: relax → SCF → bands

A band-structure plot is a **non-self-consistent** run (`ICHARG=11`): it reads a
*fixed* charge density and just diagonalises the Hamiltonian along a k-path. So
that density has to already exist and be converged. If you run `ICHARG=11` cold,
VASP uses an unconverged starting density and the eigenvalues are computed
against the wrong potential — wrong by tens of meV, which is exactly the size of
the band shifts we are trying to measure. So each directory runs in stages:

- **`INCAR.scf`** (`ICHARG=2`) — self-consistent run on a uniform mesh; converges
  the density and writes `CHGCAR`.
- **`INCAR.band`** (`ICHARG=11`) — non-self-consistent bands on the k-path,
  reading that `CHGCAR`.
- **`INCAR.relax`** — *shear directories only* (see 4.4).

The SCF stage is repeated in **every** directory, never shared: each strained
cell has its own potential and its own converged density, so a strained run must
not read another cell's `CHGCAR`.

### 4.2 `ISYM = 0` in every strained run

By default VASP exploits crystal symmetry. With a small strain it can decide the
cell is still "close enough" to cubic (within `SYMPREC`) and symmetrise it back,
which re-imposes the very degeneracy the strain is supposed to lift — giving a
splitting of exactly zero. Turning symmetry off (`ISYM=0`) keeps the
inequivalent points distinct so the splitting survives.

### 4.3 KPOINTS: a mesh for SCF, a line for bands

These are two genuinely different sampling jobs, so they need two files.

- **`KPOINTS.mesh`** — an `8×8×8` Γ-centred grid. The SCF stage needs to sample
  the *whole* Brillouin zone to build an accurate charge density (the density is
  an integral over all occupied states across the zone). A uniform grid does
  that.
- **`KPOINTS.line`** — a line through chosen high-symmetry points. The band stage
  only needs the eigenvalues *along a path*, to draw the dispersion. A line is
  not enough to converge a density (that is why it must read the mesh-derived
  `CHGCAR`), but it is exactly what a band plot wants.

### 4.4 The dual-X path — why these specific coordinates

This is the subtle one. The valley splitting under uniaxial strain only shows up
if the path visits **both** of the now-inequivalent X points. A standard cubic
path visits only one X, so it hides the effect entirely.

Under `diag(1+ε,1,1)` we stretch the Cartesian **x** axis. That makes the
$\Delta$ valley along x (call it $X_\parallel$) shift differently from the
valleys along y and z ($X_\perp$). To see the difference, the path goes
$X_\parallel \to \Gamma \to X_\perp$, with these fractional reciprocal
coordinates:

```
X_par  = (0,   1/2, 1/2)     -> along the strained x axis
Gamma  = (0,   0,   0  )
X_perp = (1/2, 0,   1/2)     -> transverse (y axis)
```

Why exactly those triples? Fractional reciprocal coordinates are coefficients of
the reciprocal lattice vectors $\mathbf{b}_i$. For the FCC primitive cell the
reciprocal lattice is BCC, with
$\mathbf{b}_1=\tfrac{2\pi}{a}(-1,1,1)$,
$\mathbf{b}_2=\tfrac{2\pi}{a}(1,-1,1)$,
$\mathbf{b}_3=\tfrac{2\pi}{a}(1,1,-1)$. Then

$$(0,\tfrac12,\tfrac12)\!:\ \tfrac12\mathbf{b}_2+\tfrac12\mathbf{b}_3
   = \tfrac{2\pi}{a}(1,0,0)\quad\text{(Cartesian x — the strain axis)},$$
$$(\tfrac12,0,\tfrac12)\!:\ \tfrac12\mathbf{b}_1+\tfrac12\mathbf{b}_3
   = \tfrac{2\pi}{a}(0,1,0)\quad\text{(Cartesian y — transverse)}.$$

So those two triples are precisely the X points along the strained axis and
perpendicular to it. And because they are written in *fractional* (basis-tied)
coordinates, they keep pointing at the strained-x and transverse faces even
after the cell is deformed — they track the strain automatically, with no need
to recompute them per strain.

**Built-in sanity check:** at $\varepsilon=0$ the conduction minima on the
$X_\parallel$ and $X_\perp$ sides must be at the *same* energy (cubic, still
degenerate). If they differ at zero strain, the path or the run is wrong before
any deformation potential is trusted.

### 4.5 Internal (Kleinman) relaxation — shear only

For hydrostatic and uniaxial strain, symmetry pins the two-atom basis at
$(\tfrac14,\tfrac14,\tfrac14)$, so simply scaling the atoms with the cell
(`scale_atoms=True`) is exact — no relaxation stage. Under $[111]$ shear that
symmetry is gone: the two sublattices shift relative to each other along $[111]$
(the Kleinman internal displacement). This must be relaxed inside the fixed
sheared cell, or the valence splitting — and therefore $d$ — comes out too
small. Hence the shear directories carry an extra `INCAR.relax` (`ISIF=2`,
`IBRION=2`), run first, with `CONTCAR → POSCAR` carried into the SCF stage.

### 4.6 `jobscript` and clean re-runs

The remote jobscript runs relax → SCF → bands in order, swapping in the matching
INCAR and KPOINTS at each stage. A guard, `if [ -f INCAR.relax ]`, runs the
relaxation only when that file is present — so the *same* script serves all
three channels: shear triggers the relax stage, hydrostatic and uniaxial skip
it. Each remote directory is wiped and recreated before upload (`rm -rf … &&
mkdir -p …`) so a stale `POSCAR` or `CHGCAR` from a previous channel cannot be
re-run by accident. `vasprun.xml` is overwritten by each stage, so the retrieved
file is always the band run.

---

## 5. Extracting the deformation potentials

The potentials are slopes of band features versus strain:

$$a = a_c + a_v = \frac{dE_g}{d\,\mathrm{Tr}(\varepsilon)}
\quad\text{(hydrostatic gap shift)},$$

$$\Xi_u = \frac{d}{d\varepsilon}\big[E_c(X_\parallel) - E_c(X_\perp)\big]
\quad\text{(uniaxial valley splitting)},$$

$$d = \frac{\Delta E_{\Gamma_{25'}}}{3\sqrt{3}\,\varepsilon_{xy}}
\quad\text{(shear valence splitting; symmetric }M\Rightarrow \varepsilon_{xy}=\epsilon).$$

| Deformation potential | This work | Literature (Si) |
|---|---|---|
| $a = a_c + a_v$ | $+1.83$ eV | $\sim +1.5$ |
| $\Xi_u$ | _in progress_ | $9.0 \pm 2$ |
| $\lvert d\rvert$ | $5.46$ eV | $\sim 5.3$ ($d<0$) |
| $dE_g/dP = -a/B$ | $-18.7$ meV/GPa | $\sim -15$ |

The pressure coefficient $dE_g/dP$ is a convention-independent cross-check: it
should be negative (the indirect gap closes under pressure) and near the
measured value, which it is.

---

## 6. VASP settings

| Parameter | Value | Notes |
|---|---|---|
| Code | VASP (PAW) | not GPAW |
| `ENCUT` | 320 eV | plane-wave cutoff |
| `KPOINTS` (SCF) | `8×8×8` Γ-centred | uniform mesh for the density |
| `KPOINTS` (bands) | line mode | path for the dispersion |
| `XC` | PBE (GGA) | not LDA |
| `ISMEAR` / `SIGMA` | 0 / 0.05 | Gaussian smearing |
| `ISYM` | 0 | keep strain-broken symmetry |

All runs in a sweep share the same POTCAR and INCAR settings; only the lattice
vectors (and, for shear, the relaxed internal positions) differ between
directories. Any change in the band structure is therefore attributable to
strain alone — which is the whole point.