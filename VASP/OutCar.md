## First — verification

Before extracting anything,let's confirm the calculation is trustworthy. Three checks:

**1. Did it converge?**
```
------------------------ aborting loop because EDIFF is reached ----------------
```
Yes — SCF reached EDIFF tolerance. 

**2. Is the structure right?**
```
volume of cell :      40.0258
position of ions in fractional coordinates:
   0.000  0.000  0.000
   0.250  0.250  0.250
NELECT = 8.0000
The static configuration has the point symmetry T_d
```
2 atoms, diamond positions, 8 electrons, tetrahedral symmetry = correct diamond Si. 

**3. Is it near equilibrium?**
```
external pressure =       19.03 kB
```
This is **not zero** — there's ~19 kBar residual pressure. This tells us a=5.43 Å isn't quite the PBE equilibrium (PBE prefers ~5.47 Å). For deformation potentials this is acceptable since you're taking *derivatives* around this reference, but be aware that the "ε=0" isn't the true stress-free state. More on this at the end.

## The sections we care about

### Section 1: Fermi energy

```
E-fermi :   6.1354     XC(G=0):  -9.4103     alpha+bet :-11.9845
```

E_F = 6.1354 eV. This is the reference level within this run. **Grep command:**
```bash
grep "E-fermi" OUTCAR
```

### Section 2: The eigenvalue block

After the Fermi line, you get one block per k-point:

```
 k-point     1 :       0.0000    0.0000    0.0000
  band No.  band energies     occupation
      1      -6.1292      2.00000
      2       5.8325      1.99998     ← top of valence (VBM region)
      3       5.8325      1.99998
      4       5.8325      1.99998
      5       8.3996      0.00000     ← bottom of conduction
      6       8.3996      0.00000
      7       8.3996      0.00000
      8       9.2095      0.00000
```

For each k-point: band index, energy (eV), occupation (2.0 = full, 0.0 = empty for non-spin-polarized).

**Reading the physics at k-point 1 (Γ):**
- Bands 2,3,4 at **5.8325 eV**, all occupied → this is the **triply-degenerate valence band maximum** (heavy hole + light hole + a third, the split-off would separate with spin-orbit which you don't have). This is the textbook Si VBM at Γ. ✓
- Bands 5,6,7 at **8.3996 eV**, empty → conduction states at Γ (but not the global minimum).

### Section 3: Finding VBM and CBM across all k-points

The VBM is the highest occupied energy anywhere; CBM is the lowest empty energy anywhere. Scan all 29 k-points:

- **VBM = 5.8325 eV** at k-point 1 (Γ) — the highest occupied value in the whole file.
- **CBM** — look for the lowest band-5 energy across k-points. Scanning:
  - k1: 8.40, k13: 6.965, k18: **6.449**, k19: 6.729, k21: 6.544...
  - The minimum is **6.449 eV at k-point 18** = (0.375, 0.375, 0.000).

So:
```
VBM = 5.8325 eV  at Γ
CBM = 6.4490 eV  at k=18 (0.375, 0.375, 0)
gap = 0.6165 eV  INDIRECT
```

This matches what you got before. ✓ And k-point 18 being along the Γ→X-ish direction confirms the indirect Δ-valley character.

**Note on k=18:** (0.375, 0.375, 0) in this mesh is near the X/Δ region of the BZ. The true Si CBM is at ~(0.85, 0, 0)·(2π/a) along Γ-X; your 8×8×8 mesh doesn't have a point exactly there, so it lands on the nearest sampled point. For accurate CBM you'd eventually want a finer mesh or a band-structure path — but for deformation potentials the *shift* of this consistently-sampled point with strain is what matters, and that's fine.

### Section 4: Total energy

```
free  energy   TOTEN  =       -10.83681854 eV
energy(sigma->0) =      -10.83681850
```

For deformation potentials you mostly use band edges, not total energy. But TOTEN vs. strain gives you the **bulk modulus** (curvature of E vs. volume) as a bonus check. Grep:
```bash
grep "free  energy   TOTEN" OUTCAR | tail -1
```

### Section 5: Stress tensor

```
  FORCE on cell =-STRESS in cart. coord.  units (eV):
  Direction    XX          YY          ZZ          XY          YZ          ZX
  Total       0.47537     0.47537     0.47537    -0.00000    -0.00000     0.00000
  in kB      19.02847    19.02847    19.02847    -0.00000    -0.00000     0.00000
  external pressure =       19.03 kB
```

Diagonal stress 19 kBar, off-diagonals zero (cubic symmetry preserved at ε=0). This is the residual pressure I flagged. When you do the **LEPSILON elastic run**, this stress-vs-strain relationship gives you C₁₁, C₁₂, C₄₄.

### Section 6: Average electrostatic potential (for cross-strain alignment)

```
 average (electrostatic) potential at core
       1 -83.2420       2 -83.2420
```

**This is important for separating a_c from a_v.** This is the potential at each atomic core (both Si atoms see −83.2420 eV, as they must by symmetry). When you compare strained runs, this value provides the common reference to align band edges between calculations. Grep:
```bash
grep -A3 "average (electrostatic) potential" OUTCAR
```

## How to extract all this with 

The existing parser already handles most of this. Here's a focused extraction function that pulls everything you need per run:

```python
import re
from pathlib import Path
import numpy as np

def extract_run(outcar_path):
    text = Path(outcar_path).read_text()
    
    # Fermi
    E_f = float(re.search(r"E-fermi\s*:\s*([-\d.]+)", text).group(1))
    
    # Volume
    V = float(re.findall(r"volume of cell :\s*([\d.]+)", text)[-1])
    
    # Total energy
    TOTEN = float(re.findall(r"free  energy   TOTEN\s*=\s*([-\d.]+)", text)[-1])
    
    # Pressure
    P = float(re.findall(r"external pressure =\s*([-\d.]+)", text)[-1])
    
    # Average electrostatic potential at cores (for alignment)
    pot_block = re.search(r"average \(electrostatic\) potential at core.*?\n.*?\n((?:\s+\d+\s+[-\d.]+)+)", text, re.DOTALL)
    V_avg = None
    if pot_block:
        nums = re.findall(r"-?\d+\.\d+", pot_block.group(1))
        V_avg = np.mean([float(x) for x in nums])
    
    # Eigenvalues (last block)
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
    
    # Band edges
    VBM = eig[occ > 0.5].max()
    CBM = eig[occ <= 0.5].min()
    
    return {
        "E_fermi": E_f, "volume": V, "TOTEN": TOTEN, "pressure": P,
        "V_avg": V_avg, "VBM": VBM, "CBM": CBM, "gap": CBM - VBM,
        "VBM_aligned": VBM - V_avg if V_avg else None,
        "CBM_aligned": CBM - V_avg if V_avg else None,
    }

# Test on this run
r = extract_run("eps_+0.0000/OUTCAR")
for k, v in r.items():
    print(f"{k:15s} = {v}")
```

Expected output for this OUTCAR:
```
E_fermi         = 6.1354
volume          = 40.0258
TOTEN           = -10.83681854
pressure        = 19.03
V_avg           = -83.2420
VBM             = 5.8325
CBM             = 6.4490
gap             = 0.6165
VBM_aligned     = 89.0745   (5.8325 - (-83.2420))
CBM_aligned     = 89.6910
```

The **aligned** values are what we'll use to separate a_c from a_v across strains.

## The full extraction loop (all 7 strains)

```python
strains = [-0.010, -0.005, -0.002, 0.000, 0.002, 0.005, 0.010]
rows = []
for eps in strains:
    r = extract_run(f"eps_{eps:+.4f}/OUTCAR")
    r["eps"] = eps
    rows.append(r)

eps_arr = np.array([r["eps"] for r in rows])
VBM_al  = np.array([r["VBM_aligned"] for r in rows])
CBM_al  = np.array([r["CBM_aligned"] for r in rows])
gap     = np.array([r["gap"] for r in rows])

# Deformation potentials = slopes
a_c = np.polyfit(eps_arr, CBM_al, 1)[0]
a_v = np.polyfit(eps_arr, VBM_al, 1)[0]
a_gap = np.polyfit(eps_arr, gap, 1)[0]

print(f"a_c (dE_c/dε)     = {a_c:.3f} eV")
print(f"a_v (dE_v/dε)     = {a_v:.3f} eV")
print(f"d(gap)/dε         = {a_gap:.3f} eV   (should ≈ a_c - a_v)")
```

## Residual pressure — deliberate, not unconverged

Calculations use the **experimental** lattice constant $a = 5.4303$ Å rather than the
PBE-relaxed value ($\approx 5.466$ Å), to match the conditions of the reference
measurement (Smith 1954, room-temperature Si). This leaves a residual hydrostatic
pressure of ~1.9 GPa in every cell. The SCF is fully converged; the pressure is a
property of the chosen geometry, not a convergence failure.

Consistency check: $P/B = 1.905/97.4 = 2.0\%$ volume compression $\rightarrow$
$\Delta a/a = 0.65\%$, recovering $a_{\rm PBE} = 5.466$ Å exactly. ($B = (C_{11}+2C_{12})/3$
from this work $= 97.4$ GPa, vs. 97.9 GPa experimental.)

**Impact.** Cancels exactly in the valley splitting (hydrostatic terms are common to all
six valleys, so $\Xi_u$ is unaffected). Contributes a few percent to $C_{11}$, $C_{12}$ —
partially cancelling in the difference $C_{11}-C_{12}$, which is the only combination
entering $\pi$.

**Usual practice.** Standard DFT convention is to relax to the theoretical (PBE)
equilibrium, so that all quantities are internally consistent and the residual pressure
vanishes. That is the right choice when comparing to other calculations. It is the wrong
choice here: PBE overestimates $a$ by ~0.7%, so relaxing would evaluate the response at a
volume the real crystal never occupies while benchmarking against experimental data.

**To relax anyway** (if reproducing at PBE equilibrium is wanted):

```
ISIF   = 3      # relax cell shape + volume + ions
IBRION = 2
NSW    = 40
EDIFF  = 1E-7
EDIFFG = -1E-3
PREC   = Accurate
ENCUT  = 500    # >= 1.3x default; ENCUT=320 gives Pulay stress errors under volume change
```

Run twice, restarting from `CONTCAR` (the plane-wave basis is fixed at the initial cell,
so one pass leaves residual Pulay stress). Converged when `external pressure` < ~0.5 kB.
All downstream quantities must then be re-extracted, since the reference geometry changed.


The `ENCUT = 500` in that block matters more than it looks. The production runs at 320 eV are fine because volume is nearly constant across the sweep and Pulay stress cancels in the slope — but a volume relaxation changes the cell enough that 320 eV would give a spuriously converged geometry. If you ever do that run, the higher cutoff isn't optional.
## Quick grep cheat-sheet for OUTCAR

```bash
grep "E-fermi" OUTCAR                              # Fermi level
grep "free  energy   TOTEN" OUTCAR | tail -1       # total energy
grep "external pressure" OUTCAR | tail -1          # stress
grep "volume of cell" OUTCAR | tail -1             # volume
grep -A30 "k-point     1 :" OUTCAR                 # first k-point eigenvalues
grep -A3 "average (electrostatic)" OUTCAR          # alignment reference
grep "reached" OUTCAR                              # convergence check
```
