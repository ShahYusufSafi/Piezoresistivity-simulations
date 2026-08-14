# VASP and phonopy file reference

Documentation of every input and output file used in the piezoresistivity pipeline: what each one is, what it controls, why the settings are what they are, and which failure it produces when it is wrong.

---

## 1. How VASP organises work

VASP has no command-line arguments and no configuration outside the working directory. It reads exactly four files from the directory it is launched in — `POSCAR`, `POTCAR`, `INCAR`, `KPOINTS` — and writes everything else into that same directory. This has one consequence that shapes the whole pipeline: **a calculation is a directory**. There is no way to run two calculations in the same folder without one overwriting the other's output.

That is why the strain sweep is one subdirectory per strain (`eps_-0.0100/`, `eps_+0.0000/`, …), each holding its own complete set of four inputs. It is also why multi-stage runs inside one directory — relax, then SCF, then bands — must copy their outputs to distinct names between stages, because the second `vasp` invocation silently destroys the first one's `OUTCAR`.

The four inputs divide cleanly by what they describe. `POSCAR` is the geometry — where the atoms are, what shape the cell is. `POTCAR` is the chemistry — how the ion cores interact with valence electrons. `KPOINTS` is the sampling — which points of the Brillouin zone the integrals run over. `INCAR` is everything else: what kind of calculation, how accurately, and what to write out.

**Strain lives in `POSCAR`, and only in `POSCAR`.** This is worth stating explicitly because it is the single most common conceptual error at the start. Deforming the crystal means multiplying the lattice vectors by a deformation matrix and writing the result as a new `POSCAR`. Nothing in `POTCAR`, `INCAR`, or `KPOINTS` knows or needs to know that the cell is strained.

---

## 2. Input files

### 2.1 POSCAR — the structure

The geometry file. Line 1 is a free-text comment, line 2 a global scale factor, lines 3–5 the three lattice vectors, then species names, species counts, the coordinate mode (`Direct` or `Cartesian`), and the atomic positions.

For diamond silicon the cell is two atoms at fractional `(0,0,0)` and `(¼,¼,¼)` on an FCC lattice with `a = 5.4303 Å`. Note that the VASP wiki's example Si POSCAR is a *fictitious single-atom FCC silicon at a = 3.9 Å* — it is a tutorial artefact, not real silicon, and using it produces a plausible-looking but physically wrong band structure.

Strain is applied by transforming the cell matrix. In ASE:

```python
M = np.diag([1+eps, 1.0, 1.0])        # [100] uniaxial
a.set_cell(cell0 @ M, scale_atoms=True)
```

`scale_atoms=True` matters: it carries the fractional positions along with the deformation so the atoms move affinely with the lattice rather than staying put in absolute space. For shear channels the matrix is off-diagonal instead, and the resulting cell no longer has cubic symmetry — which is what forces `ISYM=0` downstream (§2.4).

Two derived filenames appear throughout the pipeline. `POSCAR.relaxed` is the ionically-relaxed primitive cell, saved after the relaxation stage so later stages have an unambiguous reference; `POSCAR-001`, `POSCAR-002`, … are phonopy's displaced supercells (§5.2). These are all POSCAR-format files — the suffix is bookkeeping, not a different format.

### 2.2 POTCAR — the pseudopotential

Contains the PAW projector data for each species: how the frozen core is replaced by a smooth pseudopotential, and how valence pseudo-wavefunctions map back onto all-electron ones. This project uses PAW-PBE silicon throughout.

Three things about `POTCAR` are practically important. The species blocks must be concatenated in **the same order as the species appear in `POSCAR`** — VASP matches them positionally and will happily compute nonsense if they are swapped. The file also carries a recommended `ENCUT` (`ENMAX`), which is the source of the default when `ENCUT` is not set in `INCAR`; for Si this is around 245 eV, so the production value of 320 eV sits comfortably above it. And the number of valence electrons per species (`ZVAL`) determines `NELECT` — 4 per Si atom, hence `NELECT = 8` for the two-atom cell, which is a useful one-line check that the right POTCAR was picked up.

### 2.3 KPOINTS — Brillouin-zone sampling

Two modes are in use, and the distinction is not about accuracy but about *purpose*.

**Mesh mode** samples the BZ on a regular grid for integration. Every quantity that is an integral over occupied states — total energy, charge density, stress — needs this.

```
Regular 8 x 8 x 8 mesh centered at Gamma
0
Monkhorst-Pack
 8 8 8
 0 0 0
```

The 8×8×8 grid reduces to 29 irreducible k-points under cubic symmetry, which is why a converged unstrained OUTCAR shows 29 eigenvalue blocks. Under shear with `ISYM=0` there is no reduction and the count rises accordingly — the run gets slower, which is expected and not a bug.

**Line mode** samples along a path between high-symmetry points for plotting dispersion. It is *not* suitable for integration — the sampling is deliberately non-uniform — and is always run non-self-consistently on a charge density produced by a prior mesh-mode run.

```
X_par - Gamma - X_perp
40
line
reciprocal
0.000000 0.500000 0.500000  X_par
0.000000 0.000000 0.000000  \Gamma

0.000000 0.000000 0.000000  \Gamma
0.500000 0.000000 0.500000  X_perp
```

Segments are given as endpoint pairs, blank-line separated, with 40 points per segment. The `X_par`/`X_perp` path is the diagnostic one for the uniaxial channel — it cuts through one parallel and one perpendicular valley in a single plot so the Ξ_u splitting is directly visible. Earlier work used a 160-point L–Γ–X–K–Γ path for general band structure figures.

**Why line mode is authoritative for the gap.** An 8×8×8 mesh has no k-point at silicon's true Δ-minimum (~0.85 along Γ–X); the nearest sampled point on that mesh is `(0.375, 0.375, 0)`, and the CBM read off the mesh is therefore systematically too high. For deformation potentials the *shift* of a consistently-sampled point is what matters, so the mesh is adequate for slopes — but any absolute gap quoted in the thesis must come from the line-mode run.

The four named variants in a phonon directory are `KPOINTS.mesh` (8×8×8, SCF stage), `KPOINTS.line` (band stage), `KPOINTS.prim` (8×8×8 on the 2-atom primitive cell) and `KPOINTS.sc` (4×4×4 on the 16-atom supercell). The last pair are the ones that get confused: the supercell is eight times larger in real space, so its reciprocal cell is eight times smaller and needs a proportionally coarser mesh. **Copying `KPOINTS.prim` onto a supercell force run makes it 8× denser than intended** — the run does not fail, it just takes vastly longer than it should. That symptom (a "seconds" job running for minutes) is diagnostic.

### 2.4 INCAR — calculation control

Everything else. Parameters group naturally by what they govern.

**Basis set and accuracy**

| Tag | Value | What it does |
|---|---|---|
| `ENCUT` | 320 eV | Plane-wave kinetic-energy cutoff — the basis-set size. **Must be identical across every run in a sweep**, or energies are not comparable. |
| `PREC` | `Accurate` | Sets FFT grid density and internal tolerances. `Normal` is adequate for eigenvalues; stress needs `Accurate`. |
| `EDIFF` | 1E-6 … 1E-8 | SCF convergence threshold in eV. See below. |
| `LREAL` | `.FALSE.` | Projection in reciprocal space. Correct for small cells; the real-space alternative trades accuracy for speed on large ones. |
| `GGA` | `PE` | PBE exchange–correlation. The only physical approximation in the DFT step. |

`EDIFF` is not one value across the project because different quantities converge at different rates. Eigenvalues are the easiest, and 1E-4 was enough for the first band-edge extractions. **Stress converges considerably more slowly than energy**, so anything feeding the elastic constants uses 1E-7. Phonon force runs use 1E-8, because force constants are second derivatives and inherit the noise of the underlying forces twice over.

**Smearing**

| Tag | Value | What it does |
|---|---|---|
| `ISMEAR` | 0 | Gaussian smearing. Correct for semiconductors; the tetrahedron and Methfessel–Paxton schemes are for metals. |
| `SIGMA` | 0.05 eV | Smearing width. Broadens the occupation step for numerical stability. |

**Ionic degrees of freedom** — this is where the strain either survives or is destroyed.

| Tag | Value | What it does |
|---|---|---|
| `NSW` | 0 | Number of ionic steps. **0 freezes the atoms**, preserving the applied strain. |
| `IBRION` | −1 | Ion-update algorithm. −1 = no movement, paired with `NSW=0`. 2 = conjugate gradient relaxation. 6 = finite-difference Hessian/elastic. 8 = DFPT. |
| `ISIF` | 2 | Stress/relaxation control. **2 computes stress but holds cell shape fixed.** 3 relaxes the cell too. |

`ISIF=2` is the single most important tag in the strained runs. With `ISIF=3` the cell would relax and the applied strain would simply disappear — the calculation would converge cleanly and report the unstrained answer. Nothing in the output flags this; the deformation potential just comes out near zero.

**Symmetry**

| Tag | Value | What it does |
|---|---|---|
| `ISYM` | 0 | Switch off symmetry detection and symmetrisation. |

Required for any sheared cell. VASP detects the point group from `POSCAR` and uses it to reduce the k-mesh and symmetrise the charge density. A sheared cell has lower symmetry than cubic, but numerical noise in the cell vectors can make VASP's detector round back up to a higher group — at which point it symmetrises away the very distortion under study. `ISYM=0` costs time and removes the failure mode.

**Restart and charge**

| Tag | Value | What it does |
|---|---|---|
| `ISTART` | 0 | Start from scratch, ignore `WAVECAR`. |
| `ICHARG` | 2 | Initial charge from superposed atomic densities (self-consistent run). |
| `ICHARG` | 11 | Read `CHGCAR` and hold it fixed — **non-self-consistent**, for band-path runs. |

`ICHARG=11` with no `CHGCAR` present is a silent catastrophe: VASP starts from an empty or atomic density, never converges it, and writes eigenvalues that look structurally reasonable but are wrong. This produced a bad Ξ_u before it was caught. The band stage must always follow a completed mesh SCF stage in the same directory.

**Output control**

| Tag | Value | What it does |
|---|---|---|
| `LCHARG` | `.TRUE.` | Write `CHGCAR` — required if a band stage follows. |
| `LWAVE` | `.TRUE.`/`.FALSE.` | Write `WAVECAR`. Large; only needed for restarts. |
| `LVHAR` | `.TRUE.` | Write the Hartree potential to `LOCPOT`, for cross-run alignment. |
| `LORBIT` | 11 | Site- and orbital-projected DOS/band character. |

**Elastic constants** use a distinct combination: `LEPSILON = .TRUE.`, `IBRION = 6`, `ISIF = 3`, `NFREE = 4`, with `ENCUT` raised to 500 and `EDIFF` to 1E-8. This is a single calculation at one geometry, not a sweep, and it produces the full elastic tensor directly — the source of C₁₁ = 161.9 GPa and C₁₂ = 65.2 GPa. It belongs in its own `elastic/` directory, kept separate from the strain sweep, because the `LEPSILON` machinery slows every step and its output is the stress tensor rather than clean band-edge data.

**The INCAR variants in use**, all sharing `ENCUT=320`, `ISMEAR=0`, `SIGMA=0.05`, `ISYM=0`, `GGA=PE`:

```
INCAR.scf     ISTART=0  ICHARG=2   PREC=Accurate  EDIFF=1E-7  NSW=0   LCHARG=.TRUE.
INCAR.band    ISTART=0  ICHARG=11  LORBIT=11
INCAR.relax   IBRION=2  NSW=40     ISIF=2         EDIFF=1E-7
INCAR.force   IBRION=-1 NSW=1      EDIFF=1E-8     LREAL=.FALSE.
```

`INCAR.force` deserves a warning. It once carried `IBRION=8` from an abandoned DFPT attempt, and running that against each displaced supercell means computing *a full DFPT Hessian of an already-displaced cell* — physically meaningless, hours per displacement, and it still writes forces. `FORCE_SETS` then assembles from the wrong kind of run with no error anywhere. Finite-displacement force runs are single static evaluations: `IBRION=-1`, `NSW=1`, always.

---

## 3. The staged run pattern

A single strain directory typically runs three VASP invocations in sequence, driven by a jobscript. The cluster has no scheduler, so this is plain bash executed over SSH.

```bash
source /etc/profile.d/modules.sh 2>/dev/null || true
module load intel-oneapi-mkl

# Stage A - ionic relaxation (shear channels only)
if [ -f INCAR.relax ]; then
    cp INCAR.relax INCAR; cp KPOINTS.mesh KPOINTS
    vasp > log.relax
    cp CONTCAR POSCAR
fi

# Stage B - SCF on uniform mesh -> CHGCAR.  STRESS TENSOR LIVES HERE.
cp INCAR.scf INCAR; cp KPOINTS.mesh KPOINTS
vasp > log.scf
grep -q "reached required accuracy" log.scf || echo "WARNING: SCF not converged"
cp OUTCAR OUTCAR.scf; cp vasprun.xml vasprun.scf.xml

# Stage C - non-SCF bands on the k-path, reads CHGCAR
cp INCAR.band INCAR; cp KPOINTS.line KPOINTS
vasp > log.band
cp OUTCAR OUTCAR.band; cp vasprun.xml vasprun.band.xml
```

The stage structure follows from the physics. Stage A moves ions to their force-free positions *within* the imposed strain — necessary for shear, where the sublattices displace relative to one another (the Kleinman internal displacement), and unnecessary for hydrostatic and uniaxial strain, where symmetry fixes the internal coordinates. Stage B is the self-consistent run that produces both the converged charge density and the stress tensor. Stage C rides on that density along the k-path.

**The copy lines after each stage are not optional.** Stage C's `vasp` call overwrites `OUTCAR` and `vasprun.xml`. Without the preservation step, the stress tensor computed correctly in Stage B is destroyed before it can be read — which is exactly what caused C₁₂(yy) and C₁₂(zz) to disagree by 30% until the stages were separated. The stress was never wrong; it was clobbered.

`bash --login` is used when invoking remotely so the module environment is actually sourced, and `vasp` must be called *without* a trailing `&` — backgrounding it makes the script exit before the calculation finishes.

---

## 4. Output files

### 4.1 OUTCAR — the human-readable record

Everything VASP knows, in the order it learned it. The file is long and repetitive across ionic steps; for a static run only the last block matters.

**Convergence.** The line to look for is:

```
------------------------ aborting loop because EDIFF is reached ----------------
```

Its absence means the SCF did not converge, and nothing below should be trusted.

**Structure echo.** Near the top, VASP reports what it thinks it was given:

```
volume of cell :      40.0258
position of ions in fractional coordinates:
   0.000  0.000  0.000
   0.250  0.250  0.250
NELECT = 8.0000
The static configuration has the point symmetry T_d
```

Two atoms, diamond positions, eight valence electrons, tetrahedral point symmetry. This is the cheapest possible check that the intended structure was actually read — and the point-symmetry line is the one that catches a strain that failed to apply.

**Fermi level.**

```
E-fermi :   6.1354     XC(G=0):  -9.4103     alpha+bet :-11.9845
```

E_F is a reference level *within this run only*. It is not a valid cross-run alignment reference, for reasons in §4.6.

**Eigenvalues.** One block per irreducible k-point:

```
 k-point     1 :       0.0000    0.0000    0.0000
  band No.  band energies     occupation
      1      -6.1292      2.00000
      2       5.8325      1.99998
      3       5.8325      1.99998
      4       5.8325      1.99998
      5       8.3996      0.00000
      6       8.3996      0.00000
      7       8.3996      0.00000
      8       9.2095      0.00000
```

At Γ, bands 2–4 sit at 5.8325 eV and are occupied — the triply degenerate Γ₂₅′ valence-band maximum. (The spin–orbit split-off would separate from these, but there is no SOC in this calculation.) Bands 5–7 at 8.3996 eV are the conduction states at Γ, which are not the global minimum.

Scanning all 29 k-points gives VBM as the highest occupied energy anywhere and CBM as the lowest empty energy anywhere:

```
VBM = 5.8325 eV  at Γ
CBM = 6.4490 eV  at k-point 18, (0.375, 0.375, 0.000)
gap = 0.6165 eV  INDIRECT
```

VBM at Γ and CBM away from it — the indirect gap, with k-point 18 in the Δ region confirming the valley character. On the caveat about this being a mesh artefact rather than the true minimum, see §2.3.

**Total energy.**

```
free  energy   TOTEN  =       -10.83681854 eV
energy(sigma->0) =      -10.83681850
```

Deformation potentials use band edges, not total energy, but TOTEN versus volume gives the bulk modulus as an independent cross-check.

**Stress tensor.**

```
  FORCE on cell =-STRESS in cart. coord.  units (eV):
  Direction    XX          YY          ZZ          XY          YZ          ZX
  Total       0.47537     0.47537     0.47537    -0.00000    -0.00000     0.00000
  in kB      19.02847    19.02847    19.02847    -0.00000    -0.00000     0.00000
  external pressure =       19.03 kB
```

Diagonal and isotropic, off-diagonals zero — cubic symmetry intact at ε = 0. **The kB values are printed with pressure sign convention and must be negated** to become stress. The 19 kB residual is deliberate; see §6.

**Average electrostatic potential.**

```
 average (electrostatic) potential at core
       1 -83.2420       2 -83.2420
```

Both Si atoms see the same value, as symmetry requires. This is the cross-run alignment reference — the single most important line in the file for deformation potential work (§4.6).

### 4.2 CONTCAR — the final structure

Same format as `POSCAR`, holding the geometry at the end of the run. For a static run it is a copy of the input. After a relaxation it holds the relaxed cell, and the idiom `cp CONTCAR POSCAR` promotes it to the input for the next stage.

Two cautions. With `ISIF=2` only the ions move, so `CONTCAR` retains the imposed cell shape — which is the intended behaviour for strained runs. And `CONTCAR` written during a still-running calculation is a partial structure; copying it mid-run gives a half-relaxed cell with no warning.

### 4.3 vasprun.xml — the machine-readable record

Contains essentially everything in `OUTCAR` plus the k-point coordinates and band structure in parseable XML. This is what pymatgen reads (`BSVasprun`, `BSPlotter`) and what phonopy consumes in both phonon lanes.

Prefer it over regex-parsing `OUTCAR` wherever a library exists to read it. `OUTCAR` parsing is fine for scalars — Fermi level, volume, pressure — but fragile for eigenvalue blocks, whose formatting shifts between VASP versions and spin settings.

### 4.4 EIGENVAL, DOSCAR, IBZKPT, OSZICAR, XDATCAR

`EIGENVAL` is eigenvalues per k-point in a fixed-width plain-text format — the same data as in `OUTCAR` and `vasprun.xml`, simply easier to parse than the former. `DOSCAR` holds the density of states, site-projected when `LORBIT=11`. `IBZKPT` lists the irreducible k-points VASP actually used after symmetry reduction, which is the direct way to confirm that `ISYM=0` did what was intended — with symmetry off, the count equals the full mesh. `OSZICAR` is the one-line-per-SCF-step convergence log, useful for watching a run live. `XDATCAR` is the ionic trajectory, relevant only for relaxations.

### 4.5 CHGCAR, CHG, WAVECAR — the restart chain

`CHGCAR` is the converged charge density on the FFT grid, written when `LCHARG=.TRUE.` and read back by `ICHARG=11`. This is the hinge of the two-stage SCF→band pattern: the expensive self-consistency happens once on the mesh, and the band path reuses the result. `CHG` is a coarser running copy. `WAVECAR` holds the wavefunctions and is by far the largest file — useful for restarts, deletable otherwise.

### 4.6 LOCPOT and the alignment problem

Written when `LVHAR=.TRUE.`, containing the Hartree potential on the real-space grid.

The reason this exists at all is the deepest convention issue in the pipeline. **Kohn–Sham eigenvalues have no absolute zero.** Each run fixes its energy origin by its own average electrostatic potential, so a band-edge energy from a strained run and one from an unstrained run are not directly comparable — they are measured from different, unrelated zeros. Comparing them without alignment gives a deformation potential contaminated by an arbitrary offset.

The fix is to subtract a common reference from both. **That reference is V_avg — the average electrostatic potential at the cores — not E_Fermi.** The Fermi level is not a valid reference because it sits inside the gap of an intrinsic semiconductor, where its precise position is determined by smearing and temperature rather than by anything physical, and it therefore moves between runs for reasons unrelated to band shifts.

So:

```
VBM_aligned = VBM - V_avg
CBM_aligned = CBM - V_avg
```

and deformation potentials are the slopes of the *aligned* quantities with strain. Using `E_Fermi` here is a silent error that produces plausible but wrong numbers.

---

## 5. Phonon files (phonopy)

Phonons enter the project as the input to strain-dependent electron–phonon coupling — the route to π₄₄, which is identically zero in the rigid-valley picture. The dispersion itself is not the goal; ω(q,ε) is EPW input.

### 5.1 The two lanes

There are two ways to get force constants, and they use disjoint file sets. Mixing them is the most common failure in this part of the pipeline.

**Finite displacement (FDM)** — the validated lane. Build a supercell, displace one atom at a time, compute forces in each displaced configuration with a plain static run, and assemble the force-constant matrix from the force response. Many cheap calculations. Force constants arrive as `FORCE_SETS`.

**DFPT** (`IBRION=8`) — one expensive calculation that computes the Hessian by perturbation theory. Force constants arrive inside `vasprun.xml` and must be extracted with `phonopy --fc`, which writes `FORCE_CONSTANTS`, and the band step then needs `FORCE_CONSTANTS = READ` in `band.conf`.

DFPT was abandoned for this sweep on infrastructure grounds — a 64-atom DFPT run took 59 minutes on a saturated single-core shared login node, which does not scale to seven strains. FDM is the production lane. **The lane choice must be consistent across all three of `INCAR.force`, the collection command, and `band.conf`**; a half-switched lane produces either a hard error or, worse, a clean-looking wrong answer.

### 5.2 Displacement generation

```bash
phonopy -d --dim 2 2 2 -c POSCAR.relaxed
```

Reads the relaxed primitive cell and writes three things. `SPOSCAR` is the undisplaced 2×2×2 supercell — 16 atoms from the 2-atom primitive cell. `POSCAR-001`, `POSCAR-002`, … are the displaced supercells, one per symmetry-inequivalent displacement. `phonopy_disp.yaml` records which displacement each one carries, and is what lets phonopy reassemble the force constants later.

The `-c` flag matters. Without it phonopy reads whatever is in `POSCAR`, and at this point in the script `POSCAR` is often a leftover from a previous stage.

### 5.3 Force runs and collection

Each displaced supercell gets a static VASP run, its `vasprun.xml` preserved under a matching name:

```bash
for f in POSCAR-0*; do
    n=${f#POSCAR-}
    cp "$f" POSCAR; cp INCAR.force INCAR; cp KPOINTS.sc KPOINTS
    vasp > log.force-$n 2>&1
    cp vasprun.xml vasprun-$n.xml
done
phonopy -f vasprun-*.xml          # -> FORCE_SETS
```

`FORCE_SETS` is a plain-text file listing, for each displacement, the displacement vector and the resulting force on every atom in the supercell. Its header carries the atom count, which is what phonopy checks against the supercell it builds at the band step.

### 5.4 band.conf — the band-structure configuration

```
ATOM_NAME = Si
DIM = 2 2 2
BAND = 0.0 0.0 0.0  0.5 0.0 0.5  0.375 0.375 0.75  0.0 0.0 0.0  0.5 0.5 0.5
BAND_LABELS = \Gamma X K \Gamma L
BAND_POINTS = 101
```

`DIM` must match the `--dim` used to generate the displacements, or the atom counts will not reconcile. `BAND` is the q-path as a flat list of fractional coordinates. `BAND_LABELS` uses single-backslash LaTeX; if Python's string escaping doubles it to `\\Gamma` in the file on disk, phonopy's label parser can abort before writing any output — `cat band.conf` to check the actual bytes.

For the FDM lane this file must contain **no** `FORCE_CONSTANTS = READ` line. Without it, phonopy finds and reads `FORCE_SETS` automatically. With it, phonopy looks for a `FORCE_CONSTANTS` file that the FDM lane never wrote, and fails.

### 5.5 Assembly and output

```bash
cp POSCAR.relaxed POSCAR                      # see the warning below
phonopy --fc-symmetry band.conf
phonopy-bandplot --gnuplot band.yaml > band.dat
```

`--fc-symmetry` enforces the acoustic sum rule, which is what drives the three acoustic branches to exactly zero at Γ.

`band.yaml` is the full result: q-points, distances along the path, and frequencies per branch. `band.dat` is the two-column gnuplot flattening (distance, frequency) — convenient for a quick imaginary-mode scan, but `band.yaml` is what to plot from.

**The reference-cell trap.** At the band step phonopy reads `POSCAR` as the reference *unit* cell and expands it by `DIM`. If `POSCAR` is still a displaced supercell left over from the force loop, phonopy builds 2×2×2 of a 16-atom cell — 128 atoms — against a `FORCE_SETS` containing 16, and fails with:

```
RuntimeError: Number of forces is not consistent with supercell setting.
```

The upstream tell appears earlier in the same output: `Spacegroup: Cm (8)`. That is the broken symmetry of a cell with an atom pushed off its site. **The check that the reference cell is right is `Spacegroup: Fd-3m (227)`.** Fix by restoring `cp POSCAR.relaxed POSCAR` before the band call, or by passing `-c POSCAR.relaxed` explicitly.

### 5.6 Validation criteria

An unstrained run is trusted when three acoustic branches go to zero at Γ, the optical branches sit near 15.5 THz, no genuine imaginary frequencies appear (a threshold of −0.05 THz filters Γ numerical noise), and the dispersion overlays the independently validated curve from `Vibrational_props.ipynb`. Against Yu–Cardona Fig. 3.1 the expected discrepancy is a 2–3% PBE underestimate of the optical frequencies.

Under shear, the physics check is the Γ optical triplet splitting: 15.3804, 15.3804, 15.4412 THz against a degenerate triplet unstrained.

### 5.7 Strain magnitude for the shear channel

Frequencies are even in ε, so ±ε dispersions should overlay. Measured mirror-asymmetry is 2.1% at ±0.2% strain and 10.7% at ±1% — a 5× reduction tracking a 5× strain reduction, which identifies the asymmetry as genuine odd higher-order contamination scaling linearly with ε rather than as numerical noise. **±0.2% is therefore the clean derivative regime for the shear channel**, in contrast to the ±1% used for the linear-response channels.

---

## 6. The reference geometry: residual pressure is deliberate

Calculations use the **experimental** lattice constant a = 5.4303 Å rather than the PBE-relaxed value (≈ 5.466 Å), to match the conditions of the reference measurement (Smith 1954, room-temperature Si). This leaves a residual hydrostatic pressure of ~1.9 GPa in every cell. The SCF is fully converged; the pressure is a property of the chosen geometry, not a convergence failure.

The consistency check closes exactly: P/B = 1.905/97.4 = 2.0% volume compression → Δa/a = 0.65%, recovering a_PBE = 5.466 Å. (B = (C₁₁+2C₁₂)/3 = 97.4 GPa from this work, against 97.9 GPa experimental.)

**Impact.** It cancels exactly in the valley splitting, since hydrostatic terms are common to all six valleys — Ξ_u is unaffected. It contributes a few percent to C₁₁ and C₁₂ individually, partially cancelling in the difference C₁₁ − C₁₂, which is the only combination entering π.

**Why this departs from usual practice.** Standard DFT convention is to relax to the theoretical equilibrium so that all quantities are internally consistent and the residual pressure vanishes. That is the right choice when comparing against other calculations. It is the wrong choice here: PBE overestimates a by ~0.7%, so relaxing would evaluate the strain response at a volume the real crystal never occupies, while the benchmark is experimental.

**To relax anyway**, if reproduction at PBE equilibrium is wanted:

```
ISIF   = 3
IBRION = 2
NSW    = 40
EDIFF  = 1E-7
EDIFFG = -1E-3
PREC   = Accurate
ENCUT  = 500
```

Run twice, restarting from `CONTCAR`, since the plane-wave basis is fixed at the initial cell and one pass leaves residual Pulay stress. Converged when `external pressure` < ~0.5 kB. All downstream quantities must then be re-extracted.

`ENCUT = 500` here matters more than it looks. The production runs at 320 eV are fine because volume is nearly constant across the sweep and Pulay stress cancels in the slope. A volume relaxation changes the cell enough that 320 eV would give a spuriously converged geometry. For that run the higher cutoff is not optional.

---

## 7. Failure modes and their signatures

Collected because each of these produces plausible output rather than an error.

| Signature | Cause |
|---|---|
| Deformation potential ≈ 0 | `ISIF=3` on a strained run — the strain relaxed away |
| Strain response symmetrised away | `ISYM` not set to 0 on a sheared cell |
| Wrong Ξ_u, bands look fine | `ICHARG=11` with no `CHGCAR` — unconverged eigenvalues |
| C₁₂(yy) ≠ C₁₂(zz) by ~30% | Stage C overwrote Stage B's `OUTCAR` before the stress was read |
| Deformation potential offset by a constant | Aligned to `E_Fermi` instead of `V_avg` |
| Absolute gap too large | Read from the SCF mesh instead of the line-mode run |
| Sign error in stress | kB block not negated (printed as pressure) |
| `Number of forces is not consistent` | `POSCAR` at the band step is a displaced supercell |
| `Spacegroup: Cm (8)` | Same cause — reference cell is displaced, should be `Fd-3m (227)` |
| Force run takes minutes not seconds | `KPOINTS.prim` copied onto a supercell run — 8× too dense |
| Displacement runs take hours each | `IBRION=8` left in `INCAR.force` — DFPT on every displacement |
| Elastic constants a few % off | Pulay stress at fixed `ENCUT` under volume change |
| π extracted from a fit is biased | Nonlinearity above ~0.3% strain; use the analytic derivative at ε=0 |

---

## 8. Grep cheat-sheet

```bash
grep "reached" OUTCAR                              # convergence check
grep "E-fermi" OUTCAR                              # Fermi level
grep "free  energy   TOTEN" OUTCAR | tail -1       # total energy
grep "external pressure" OUTCAR | tail -1          # stress
grep "volume of cell" OUTCAR | tail -1             # volume
grep -A30 "k-point     1 :" OUTCAR                 # first k-point eigenvalues
grep -A3 "average (electrostatic)" OUTCAR          # alignment reference
grep "NELECT" OUTCAR                               # valence electron count
grep "point symmetry" OUTCAR                       # detected point group
```

---

## 9. Extraction

One function pulling everything needed from a single run.

```python
import re
from pathlib import Path
import numpy as np

def extract_run(outcar_path):
    text = Path(outcar_path).read_text()

    E_f   = float(re.search(r"E-fermi\s*:\s*([-\d.]+)", text).group(1))
    V     = float(re.findall(r"volume of cell :\s*([\d.]+)", text)[-1])
    TOTEN = float(re.findall(r"free  energy   TOTEN\s*=\s*([-\d.]+)", text)[-1])
    P     = float(re.findall(r"external pressure =\s*([-\d.]+)", text)[-1])

    # average electrostatic potential at cores -> cross-run alignment reference
    pot_block = re.search(
        r"average \(electrostatic\) potential at core.*?\n.*?\n((?:\s+\d+\s+[-\d.]+)+)",
        text, re.DOTALL)
    V_avg = None
    if pot_block:
        nums = re.findall(r"-?\d+\.\d+", pot_block.group(1))
        V_avg = np.mean([float(x) for x in nums])

    # eigenvalues, last block only
    kpat = re.compile(
        r"k-point\s+\d+\s*:\s*([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s*\n"
        r"\s*band No\..*?\n((?:\s*\d+\s+[-\d.]+\s+[-\d.]+\s*\n)+)")
    matches = list(kpat.finditer(text))
    nk = int(re.search(r"NKPTS\s*=\s*(\d+)", text).group(1))

    eig, occ = [], []
    for m in matches[-nk:]:
        e, o = [], []
        for line in m.group(4).strip().split("\n"):
            p = line.split()
            e.append(float(p[1])); o.append(float(p[2]))
        eig.append(e); occ.append(o)
    eig = np.array(eig); occ = np.array(occ)

    VBM = eig[occ > 0.5].max()
    CBM = eig[occ <= 0.5].min()

    return {
        "E_fermi": E_f, "volume": V, "TOTEN": TOTEN, "pressure": P,
        "V_avg": V_avg, "VBM": VBM, "CBM": CBM, "gap": CBM - VBM,
        "VBM_aligned": VBM - V_avg if V_avg else None,
        "CBM_aligned": CBM - V_avg if V_avg else None,
    }
```

Expected output on the ε = 0 reference run:

```
E_fermi         = 6.1354
volume          = 40.0258
TOTEN           = -10.83681854
pressure        = 19.03
V_avg           = -83.2420
VBM             = 5.8325
CBM             = 6.4490
gap             = 0.6165
VBM_aligned     = 89.0745
CBM_aligned     = 89.6910
```

Sweeping and taking slopes:

```python
strains = [-0.010, -0.005, -0.002, 0.000, 0.002, 0.005, 0.010]
rows = [dict(extract_run(f"eps_{e:+.4f}/OUTCAR.scf"), eps=e) for e in strains]

eps_arr = np.array([r["eps"]         for r in rows])
VBM_al  = np.array([r["VBM_aligned"] for r in rows])
CBM_al  = np.array([r["CBM_aligned"] for r in rows])
gap     = np.array([r["gap"]         for r in rows])

a_c   = np.polyfit(eps_arr, CBM_al, 1)[0]
a_v   = np.polyfit(eps_arr, VBM_al, 1)[0]
a_gap = np.polyfit(eps_arr, gap,    1)[0]
print(f"a_c = {a_c:.3f} eV,  a_v = {a_v:.3f} eV,  d(gap)/dε = {a_gap:.3f} eV")
```

`d(gap)/dε` should equal `a_c − a_v` — a free internal consistency check.

Note that these slopes come from `polyfit` over the full sweep, which is appropriate for the *deformation potentials*, where the band edges are linear in ε over this range. It is **not** the right method for extracting the π coefficients themselves: those require the analytic derivative at ε = 0, because fit-based extraction is severely biased by nonlinearity above ~0.3% strain.

---

## 10. Directory layout

```
project/
├── elastic/                       # one LEPSILON calculation
│   ├── POSCAR  POTCAR  KPOINTS
│   ├── INCAR                      # LEPSILON=.TRUE., IBRION=6, ISIF=3
│   └── OUTCAR                     # -> C11, C12, C44
│
├── deformation_potentials/        # strain sweep, one dir per strain
│   └── eps_+0.0000/
│       ├── POSCAR  POTCAR
│       ├── INCAR.scf  INCAR.band
│       ├── KPOINTS.mesh  KPOINTS.line
│       ├── OUTCAR.scf   vasprun.scf.xml     # stress, alignment
│       └── OUTCAR.band  vasprun.band.xml    # band edges, masses
│
└── disp_curve_sweeps/             # phonons, one dir per strain
    └── eps_+0.0000/
        ├── POSCAR.relaxed                   # reference cell
        ├── INCAR.relax  INCAR.force
        ├── KPOINTS.prim  KPOINTS.sc
        ├── SPOSCAR  POSCAR-001 ...          # phonopy -d
        ├── phonopy_disp.yaml
        ├── vasprun-001.xml ...              # per-displacement forces
        ├── FORCE_SETS                       # phonopy -f
        ├── band.conf
        └── band.yaml  band.dat              # phonopy --fc-symmetry
```
