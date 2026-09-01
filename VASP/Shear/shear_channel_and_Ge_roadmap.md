# Closing the Shear Channel, and the Si–Ge Paper

*Research briefing and roadmap. Prepared 2 August 2026.*

---

## 0. The headline

**Your π₄₄ = 0 theorem is correct, but the conclusion you drew from it is not.**

You concluded that the nonzero experimental π₄₄ in n-Si "lives entirely in strain-dependent
scattering, invisible to RTA by construction," and therefore that closing the shear channel
requires EPW and beyond-RTA physics.

That is almost certainly wrong. The dominant mechanism for π₄₄ in n-Si is a **band-structure**
effect, not a scattering effect, and it is fully accessible inside constant-τ RTA. It is missing
from your pipeline for a specific and fixable reason: the Herring–Vogt Hamiltonian you are using
keeps only the *diagonal* part of the strain coupling, and shear enters the conduction band
*off-diagonally*.

A back-of-envelope estimate using the corrected physics gives

> **π₄₄(n-Si) ≈ −13.4 × 10⁻¹¹ Pa⁻¹**, against Smith's measured **−13.6 × 10⁻¹¹ Pa⁻¹**.

The agreement to ~2% is partly luck (one input parameter is uncertain at the 10–15% level), but
the *order of magnitude and the sign are not luck*. This means the shear channel is closable in
weeks with your existing VASP setup, not months with EPW.

Section 1 explains the physics from scratch. Section 2 gives the derivation and the numbers you
can check against. Section 3 is the roadmap for closing the channel. Section 4 is Ge — where the
same calculation, run with the *original* rigid-valley model, reproduces Smith's n-Ge π₄₄ to 5%.
Section 5 is the resulting paper. Section 6 is annotated literature.

---

## 1. The physics, from intuition

### 1.1 What your theorem actually proves

The Herring–Vogt Hamiltonian is

$$H_{\mathrm{HV}}^{(v)} = \Xi_d\,\mathrm{Tr}(\varepsilon) + \Xi_u\,\hat{m}^{(v)\mathsf{T}}\varepsilon\,\hat{m}^{(v)}$$

where $\hat{m}^{(v)}$ is the unit vector along the axis of valley $v$. This describes a **rigid
valley**: the whole ellipsoid moves up or down in energy by a scalar amount, keeping its shape,
its curvature, and its position in $k$-space fixed.

For Si's six ⟨100⟩ valleys and a pure shear strain (all diagonal components zero), the bracket
$\hat{m}^{\mathsf{T}}\varepsilon\hat{m}$ picks out only diagonal components of $\varepsilon$,
which are zero. All six valleys shift by nothing. No splitting, no repopulation, no signal.
π₄₄ = 0, exactly.

**This is a theorem about the rigid-valley approximation, not about relaxation-time approximation.**
Those are two independent assumptions and you have been conflating them. Constant-τ RTA is an
assumption about *scattering*. The rigid valley is an assumption about *band structure*. Your
proof only uses the second.

Hamaguchi states the physical content directly in §6.3.9: uniaxial stress along [111] acts
equivalently on all six Si valleys, so the degeneracy is not removed. Your derivation is cleaner
and more general than the textbook statement, but a referee will recognise it, so it cannot be
sold as the novelty.

### 1.2 Why the rigid valley fails under shear — the nonsymmorphic story

Here is the part that connects to something you already understand.

**Vocabulary.** A *symmorphic* space group can be built entirely from point-group operations
(rotations, reflections, inversions) about a single fixed origin, plus lattice translations. A
*nonsymmorphic* space group additionally requires operations that combine a point operation with
a translation by a *fraction* of a lattice vector. A **glide plane** is one such operation: reflect
through a plane, then slide by a fractional translation. Diamond structure (Fd-3m) is
nonsymmorphic, and its glide operations are exactly what relate your two-atom basis (SL1 and SL2
at τ = (¼,¼,¼)a) to each other.

You already noted in your own working that the two-atom basis is what produces the *X*-point
degeneracies. That is precisely the hook.

In Si, the two lowest conduction bands — labelled **Δ₁** and **Δ₂′** in BSW notation — are
**degenerate exactly at the X point**, and this degeneracy is *enforced by the diamond glide
symmetry*. It is not accidental. The conduction band minimum sits at
$k_0 = 0.15 \times (2\pi/a)$ *away* from X (i.e. at 85% of the Γ→X line), which is close enough
that the presence of the second band a fraction of an eV away controls the shape of the valley.

Now: **shear strain ε_xy destroys the glide plane at z = a₀/8**, while leaving the glide planes at
x = a₀/8 and y = a₀/8 intact. The symmetry protecting the X-point degeneracy is gone, so the two
bands split. The size of the split is $2D\varepsilon_{xy}$, where $D$ is a new deformation
potential.

**Vocabulary.** $D$ is the **interband (off-diagonal) shear deformation potential**. It is a
different object from Ξ_u and Ξ_d. Ξ_u and Ξ_d are *diagonal* matrix elements — they tell you how
much a band moves. $D$ is an *off-diagonal* matrix element — it tells you how strongly strain
*mixes two bands together*. Rigid-valley Herring–Vogt has no slot for it, which is exactly why
your pipeline cannot see shear.

### 1.3 What the two-band model gives you

Sverdlov, Ungersboeck, Kosina & Selberherr (ESSDERC 2007, and IEEE TED 54, 2183, 2007) wrote down
the minimal two-band k·p Hamiltonian around X and solved it. Near the [001] valley pair:

$$E(\mathbf{k}) = \frac{\hbar^2 k_z^2}{2m_l} + \frac{\hbar^2(k_x^2+k_y^2)}{2m_t} + \delta E_C
- \sqrt{\left(\frac{\hbar k_z p}{m_0}\right)^2 + \left(D\varepsilon_{xy} - \frac{\hbar^2 k_xk_y}{M}\right)^2}$$

where $\delta E_C$ is the usual Herring–Vogt diagonal shift, $p$ is a momentum matrix element,
and $M$ is a coupling mass parameter.

Define the **dimensionless shear** $\eta \equiv 2D\varepsilon_{xy}/\Delta$, where $\Delta$ is the
Δ₁–Δ₂′ energy gap evaluated at the valley minimum. Then shear does three things, and they scale
differently:

| Effect | Scaling in η | Contributes to *linear* π₄₄? |
|---|---|---|
| Valley minimum moves toward X: $k_{\min} = -k_0\sqrt{1-\eta^2}$ | quadratic | no |
| Valley energy shift: $\Delta E = -(\Delta/4)\eta^2$ | quadratic | **no** |
| Transverse mass splits into two branches | **linear** | **yes** |
| Longitudinal mass: $m_l(\eta)/m_l = (1-\eta^2)^{-1}$ | quadratic | no |
| Nonparabolicity α increases | quadratic | no |

The transverse mass branches are

$$\frac{m_{t1}(\eta)}{m_t} = \left(1 - \eta\frac{m_t}{M}\right)^{-1}, \qquad
\frac{m_{t2}(\eta)}{m_t} = \left(1 + \eta\frac{m_t}{M}\right)^{-1}$$

for the directions across and along [110] respectively.

**This is the whole story.** The energy shift under shear *is* second order — so your theorem
survives at linear order, and you were right that repopulation cannot produce π₄₄. But the
**mass anisotropy is first order**, and that is what produces π₄₄.

### 1.4 The intuitive picture

Draw the constant-energy surface of the [001] valley. Unstrained, it is an ellipsoid: long along
z (mass $m_l$), and *circular* in the x–y cross-section (mass $m_t$ in both directions).

Apply ε_xy. The circular cross-section becomes an **ellipse with its principal axes along [110]
and [1̄10]** — i.e. rotated 45° from the crystal axes. The valley has not moved in energy and has
not repopulated. It has been *squashed sideways*.

An ellipse whose principal axes are rotated 45° from your coordinate frame has an **off-diagonal
inverse-mass tensor component** $(1/m)_{xy}$. That off-diagonal inverse mass is exactly what
produces an off-diagonal conductivity σ_xy, which is exactly what π₄₄ measures.

So: **π₄₄ in n-Si is a valley-shape effect, not a valley-population effect.** Repopulation is
forbidden by symmetry; shape distortion is not.

---

## 2. Verifiable numbers

### 2.1 Closed-form result

Working through the algebra (derivation sketched below, **you must verify it independently** —
see §2.4):

$$\boxed{\;\pi_{44} \simeq -\frac{D\,m_c}{3\,\Delta\,M\,C_{44}}\;}$$

Note what is **absent**: there is no $1/k_BT$. This is a mass-modulation effect, so to leading
order it is **temperature-independent**. Contrast with the repopulation channels (π₁₁, π₁₂), which
carry $\Xi_u/k_BT$ and therefore fall roughly as 1/T. This is a sharp, falsifiable prediction and
there is a classic paper to test it against (Morin, Geballe & Herring, *Phys. Rev.* **105**, 525,
1957 — temperature dependence of piezoresistance in Si and Ge).

### 2.2 Input parameters and where each comes from

| Symbol | Meaning | Value | Source | Can you get it yourself? |
|---|---|---|---|---|
| $k_0$ | CBM offset from X | $0.15 \times 2\pi/a$ | standard | yes — read off your band structure |
| $\Delta$ | Δ₁–Δ₂′ gap at CBM | ≈ 0.48 eV | $=2\hbar^2k_0^2/m_l$ | **yes — $\Delta = 4\times[E_c(X) - E_{\rm CBM}]$** |
| $D$ | interband shear def. pot. | ≈ 14 eV | Hensel/Sverdlov — **verify** | **yes — slope of X-point splitting vs ε_xy** |
| $M$ | coupling mass | ≈ 0.24 $m_0$ | $M \approx m_t/(1-m_t/m_0)$ | yes — from mass anisotropy at finite shear |
| $m_c$ | conductivity mass | 0.264 $m_0$ | $1/m_c=\frac13(1/m_l+2/m_t)$ | yes — you have $m_l$, $m_t$ |
| $C_{44}$ | shear elastic constant | 79.6 GPa | experiment | **yes — you already run IBRION=6** |

Using your own $m_l = 0.959$, $m_t = 0.194$:

- $\Delta = 2\hbar^2k_0^2/m_l$ with $\hbar^2/m_0 = 7.62$ eV·Å², $k_0 = 0.1736$ Å⁻¹ → **Δ ≈ 0.479 eV**
- Equivalently, X should sit **≈ 0.12 eV above your CBM** — check this on your existing band structure right now, it costs nothing
- $M = 0.194/(1-0.194) = $ **0.241 $m_0$**
- $1/m_c = \frac13(1.043 + 10.31) = 3.784$ → **$m_c$ = 0.264 $m_0$**
- $\eta = 2D\varepsilon_{xy}/\Delta \approx 58\,\varepsilon_{xy}$

Then:

$$\pi_{44} = -\frac{14}{0.479}\times\frac{0.264}{0.241}\times\frac{1}{3\times79.6\times10^9}
= -\frac{29.2 \times 1.098}{2.388\times10^{11}} = -1.34\times10^{-10}\ \mathrm{Pa}^{-1}$$

**π₄₄ = −13.4 × 10⁻¹¹ Pa⁻¹.  Smith: −13.6 × 10⁻¹¹ Pa⁻¹.**

### 2.3 An independent consistency check on Δ and D

This one is worth doing because it validates two parameters at once without any transport
calculation.

For uniaxial stress $X$ along [110], the shear component is $\varepsilon_{xy} = X/(4C_{44})$.
Setting $\eta = 1$ (the point where the two valleys merge at X and the model breaks down) requires
$\varepsilon_{xy} = 1/58 = 1.7\%$, hence

$$X = 4 \times 79.6\ \mathrm{GPa} \times 0.017 = 5.4\ \mathrm{GPa}$$

Sverdlov et al.'s Fig. 6 inset plots η against [110] stress and reaches η = 1 at ≈ 5–6 GPa.
**Independent confirmation that Δ ≈ 0.48 eV and D ≈ 14 eV are mutually consistent.**

### 2.4 Derivation sketch — verify this yourself

Do not take the boxed formula on trust. Here is the chain; check every step, especially the
Voigt factors.

1. Linearise the mass branches:
   $1/m_{t1} = 1/m_t - \eta/M$, $\;1/m_{t2} = 1/m_t + \eta/M$.
2. Principal axes are $\hat u = (1,1,0)/\sqrt2$ and $\hat v = (1,-1,0)/\sqrt2$.
   Build $(1/m)_{ij} = (1/m_{t2})\hat u_i\hat u_j + (1/m_{t1})\hat v_i\hat v_j$ in the xy block.
   Since $\hat u_x\hat u_y = +\tfrac12$ and $\hat v_x\hat v_y = -\tfrac12$:
   $$(1/m)_{xy} = \eta/M, \qquad (1/m)_{xx} = (1/m)_{yy} = 1/m_t \;\text{(unchanged at first order)}$$
   The second equality is an important check: it guarantees shear makes **no** first-order
   contribution to π₁₁ or π₁₂, as symmetry demands.
3. Only the [001] valley pair responds to ε_xy. The [100] pair couples to ε_yz, the [010] pair to
   ε_zx — both zero here. So only $n/3$ of the carriers contribute, and their population is
   unchanged at first order (energy shift is O(η²)).
4. $\sigma_{xy} = \frac{n}{3}e^2\tau\,\eta/M$; $\;\sigma_0 = ne^2\tau/m_c$;
   $\;\sigma_{xy}/\sigma_0 = \frac{\eta}{3}\frac{m_c}{M}$.
5. $\Delta\rho_{xy}/\rho_0 = -\sigma_{xy}/\sigma_0$ to first order.
6. **Voigt:** $\varepsilon_6 = 2\varepsilon_{xy}$, $\;\sigma_6 = C_{44}\varepsilon_6 = 2C_{44}\varepsilon_{xy}$,
   $\;\Delta\rho_6 = \pi_{44}\sigma_6$. Substituting $\eta = 2D\varepsilon_{xy}/\Delta$ gives the box.

⚠ **The factor-of-two trap.** $\varepsilon_4 = 2\varepsilon_{yz}$, $S_{44} = 1/C_{44}$, and
$\Delta\rho_4 = \pi_{44}\sigma_4$. Cross any one of these conventions and π₄₄ comes out exactly
2× or ½× — which looks like a physics discrepancy and is not. Write an assertion into the
notebook that checks the shear conversion against a hand-computed case before trusting any output.

⚠ **Sign of the branch assignment.** Which of $m_{t1}$/$m_{t2}$ is "along" versus "across" [110]
sets the sign of π₄₄. Determine this from your own DFT rather than from my labelling.

---

## 3. Roadmap: closing the shear channel

Everything here uses your existing VASP + notebook stack. No EPW.

### Step 1 — Look at your existing band structure again (½ day)

You already have the three-panel figure. Two checks, zero new computation:

- **Is X ≈ 0.12 eV above your CBM?** If yes, Δ ≈ 0.48 eV is confirmed from your own data.
- **In your shear panel, is the conduction band degeneracy split at X?** Your note says
  "Δ-valley minima remain degenerate" — that is correct and expected, because the *energy* shift is
  second order. But you were looking at the wrong place. Look at **X**, and look at **curvature**,
  not at the minimum energy.

You also flagged that you had mislabelled the ~2.5 eV conduction splitting as Γ₂₅′ when it is
Γ₁₅. Fix that at the same time; it is the same figure.

### Step 2 — Extract D from DFT (2–3 days)

Apply pure ε_xy at several magnitudes (suggest ±0.25%, ±0.5%, ±1.0%). At each, read the splitting
of the two lowest conduction bands **at the X point**. Plot splitting vs ε_xy. It should be linear
through the origin with

$$\text{splitting} = 2D\varepsilon_{xy}$$

Slope/2 gives you **D, from your own calculation, parameter-free**. Compare against ~14 eV.

This is the single highest-value calculation available to you right now. It is cheap, it is
decisive, and it converts a literature parameter into an own-work number.

### Step 3 — Extract the strained mass tensor (1 week)

At each shear magnitude, fit the curvature of the lowest conduction band in the $k_x$–$k_y$ plane
at $k_z = k_{\min}$. You want the full 2×2 inverse-mass tensor, not just the trace. Confirm:

- principal axes rotate to [110] / [1̄10]
- $(1/m)_{xy}$ is **linear** in ε_xy with slope $1/(M) \times (2D/\Delta)$
- $(1/m)_{xx}$ is flat at first order
- $k_{\min}$ moves toward X as $\sqrt{1-\eta^2}$ — quadratic, so nearly flat at small strain

Each of these is a falsifiable prediction of the two-band model against your DFT. If they hold,
you have validated the model with your own data, which is a publishable figure on its own.

### Step 4 — Wire into the transport module (3–4 days)

Your conductivity sum currently assumes a rigid inverse-mass tensor per valley. Generalise it to
accept a **strain-dependent** inverse-mass tensor $(1/m)^{(v)}(\varepsilon)$ with off-diagonal
components. This is a modest refactor and it makes the module correct for arbitrary strain, not
just shear.

Then compute π₄₄ by the same analytic-derivative-at-zero method you used for π₁₁ and π₁₂ (do not
fit finite-strain curves — the nonlinearity argument that bit you before still applies).

### Step 5 — Temperature sweep (2 days)

Compute π₁₁, π₁₂, π₄₄ over 200–400 K. Prediction: π₁₁ and π₁₂ scale roughly as 1/T; π₄₄ is
roughly flat. Compare against Morin, Geballe & Herring (1957). If this holds it is the single most
convincing figure in the paper, because it is a *qualitative* discrimination between two
mechanisms rather than a number that happens to match.

### Step 6 — Honest scope statement

What remains genuinely beyond RTA: the *residual* discrepancy after mass modulation is included,
plus the ~15–20% shortfall you already see in π₁₁/π₁₂. Strain-dependent intervalley scattering is
real — Poncé's group shows intervalley scattering drops by about a third as the valleys split —
it is just not the *leading* term for π₄₄. Say so plainly and the EPW discussion becomes a
well-posed outlook rather than an unfulfilled promise.

---

## 4. Germanium

### 4.1 Why Ge is the perfect complement — and it's cleaner than Si

Here is the thing that makes the paper: **in Ge, π₄₄ comes from the mechanism Si cannot use, and
Si's π₁₁/π₁₂ mechanism is the one Ge cannot use.** They are exact mirror images.

Ge's conduction minima are four L valleys along ⟨111⟩, at the Brillouin zone boundary.

**Under [100] uniaxial strain.** For every L valley, $\hat m^{\mathsf T}\varepsilon\hat m =
(\varepsilon_1 + 2\varepsilon_2)/3$ — identical for all four. No splitting, no repopulation, so
$\pi_{11} - \pi_{12}$ must vanish in the valley channel. Smith measured
π₁₁ = −4.7 and π₁₂ = −5.0, i.e. π₁₁ − π₁₂ = **+0.3**, zero inside experimental scatter. The
common part (both ≈ −5) is the hydrostatic/carrier-density channel.

**Under ε_xy shear.** $\hat m^{\mathsf T}\varepsilon\hat m = 2\varepsilon_{xy}m_xm_y = \pm\tfrac23\varepsilon_{xy}$.
Two valleys up, two down, splitting $\tfrac43\Xi_u\varepsilon_{xy}$. **Nonzero, first order,
rigid-valley.** Ordinary Herring–Vogt repopulation, no two-band correction needed.

### 4.2 Ge π₄₄ estimate — check this

Same machinery as your existing Si π₁₁ calculation, just with ⟨111⟩ valleys:

$$\pi_{44}^{\rm Ge} \simeq -\frac{2}{9}\left(\frac{1}{m_t}-\frac{1}{m_l}\right)m_c
\frac{\Xi_u}{2k_BT\,C_{44}}$$

With Ge values $m_l = 1.59$, $m_t = 0.082$ (→ $m_c = 0.120\,m_0$, $1/m_t - 1/m_l = 11.57$),
$\Xi_u^L = 16.3$ eV, $C_{44} = 67.1$ GPa, $k_BT = 0.02585$ eV:

**π₄₄ ≈ −145 × 10⁻¹¹ Pa⁻¹.  Smith: −137.9 × 10⁻¹¹ Pa⁻¹.  ≈ 5% agreement.**

From constant-τ RTA and literature deformation potentials. Nothing fitted.

### 4.3 The 2×2 that is the paper

| | **Si** (⟨100⟩, 6 valleys) | **Ge** (⟨111⟩, 4 valleys) |
|---|---|---|
| **[100] uniaxial** | splits 2:4 → strong repopulation | degenerate → **no repopulation** |
| Smith π₁₁ − π₁₂ | **−155.6** | **+0.3** |
| **xy shear** | degenerate → **no repopulation** | splits 2:2 → strong repopulation |
| Smith π₄₄ | **−13.6** (residual, mass modulation) | **−137.9** |
| shear/uniaxial ratio | 0.087 | ~460 |
| mechanism of π₄₄ | Δ₁–Δ₂′ mass modulation, T-independent | valley repopulation, ∝ 1/T |

All units 10⁻¹¹ Pa⁻¹. One criterion — whether $\hat m^{\mathsf T}\varepsilon\hat m$ is degenerate
across the valley star — predicts every cell, and the two materials are complements across a
factor of ~5000 in the shear/uniaxial ratio.

### 4.4 The Ge go/no-go: functional choice

**This is the one thing that can kill the Ge extension, and you should test it in week one.**

Plain PBE does not describe Ge's conduction band. QuantumATK's benchmark study is blunt about it:
PBE and PBEsol close the direct gap at Γ and give a qualitatively wrong, semi-metallic Ge band
structure, and should not be used for Ge band gaps. TB09-MGGA (mBJ), pps-PBE, and HSE preserve the
correct ordering. Experimentally the Γ minimum is only **0.14 eV** above L, so the ordering is
fragile even in a good calculation.

Worse: a very recent paper (arXiv:2512.08857, Dec 2025) shows HSE results for Ge are *scattered
even between groups using the same code, pseudopotentials, mixing and screening* — and that HSE
often fails to get the Γ–L and Γ–Γ gaps right simultaneously. They propose a DFT+α correction on
the 4s-like orbitals. Read this before choosing.

**Recommended path:**
1. Relax with PBE (lattice constant is fine).
2. Single-point band structure with **mBJ (TB09)** — cheap, known to work for Ge, available in VASP
   as METAGGA = MBJ.
3. **Spin-orbit coupling is not optional.** Ge's Δ₀ ≈ 0.29 eV versus Si's 0.044 eV.
4. **Confirm L is below Γ.** If it is not, stop and re-plan.

**The escape hatch worth testing explicitly.** Ξ_u is an energy *derivative* with respect to
strain, and derivatives are often far better behaved in PBE than absolute band positions. Compute
$\Xi_u^L$ in both PBE and mBJ/HSE at two or three strain points. If they agree to a few percent,
you can run the full sweep in PBE and cite the spot-check. That single paragraph converts your
biggest vulnerability into a demonstrated result — and a referee who knows Ge will look for
exactly it.

### 4.5 A free simplification

Ge's L points sit at (½,½,½) in reciprocal lattice units — half a reciprocal lattice vector, hence
a **time-reversal-invariant momentum**. The band gradient vanishes there by symmetry for arbitrary
strain, so the L minima stay pinned at fixed fractional coordinates no matter how you deform the
cell.

You do not have the problem you had with Si, where the Δ minima drift along the line under strain
and must be re-located. Say this in the paper — it makes the Ge shear calculation cleaner than the
Si one, and it signals you know what you are doing.

Also worth noting: L valleys have **no Δ₁–Δ₂′ analogue**. There is no second conduction band
degenerate with them at L requiring a two-band treatment. Ge's π₄₄ is therefore *simpler* than
Si's, not harder — which is a nice inversion of the expected narrative.

### 4.6 Ge reference values

| Quantity | Value | Source |
|---|---|---|
| $\Xi_u^L$ | 16.2 ± 0.4 eV (80 K) | Balslev, *Phys. Rev.* **143**, 636 (1966) |
| $\Xi_u^L$ | 16.3 (80 K), 19.3 (4 K) | Hamaguchi Table 6.1 |
| $\Xi_u^L$ | 16.8 eV | Grundmann Table 6.8 |
| $\Xi_d^L$ | −4.43 eV | Grundmann Table 6.8 |
| $\Xi_d^L$ | −12.3 (4 K), −4.4 (Fischetti) | Hamaguchi Table 6.1 |
| $m_l$, $m_t$ | 1.59, 0.082 | standard |
| $C_{11}$, $C_{12}$, $C_{44}$ | 128.5, 48.3, 67.1 GPa | standard — compute your own |
| Smith n-Ge π₁₁, π₁₂, π₄₄ | −4.7, −5.0, −137.9 | Grundmann Table 8.3 (after Smith 1954) |
| Smith p-Ge π₁₁, π₁₂, π₄₄ | −10.6, +5.0, +98.6 | same |

⚠ The Ξ_u spread (16.2 to 19.3 eV depending on temperature and source) is much wider than Si's.
Quote which one you used and show the sensitivity. Also note the **compliance conversion must be
rebuilt** — Ge's $S_{11}-S_{12}$ conversion factors differ substantially from Si's near-unity
values.

---

## 5. What the paper becomes

### 5.1 Structure

1. **Si as the control, briefly.** State π₁₁, π₁₂. Cite Smith (1954) for experiment and Roisin
   et al. (2024) for the state-of-the-art first-principles calculation. Be plain that their
   iterative BTE with full electron–phonon coupling is the more accurate calculation — 3% from
   experiment versus your ~17%. Do not argue otherwise.
2. **The symmetry criterion.** One bracket, $\hat m^{\mathsf T}\varepsilon\hat m$, and the 2×2
   table. This is Figure 1.
3. **Si π₄₄: where the rigid valley fails.** The nonsymmorphic X-point degeneracy, the off-diagonal
   coupling D, mass modulation, the closed form, and the temperature-independence prediction.
4. **Ge: the mirror.** Rigid-valley repopulation recovers π₄₄ to ~5%, and π₁₁ − π₁₂ ≈ 0 falls out
   for free.
5. **Device layer.** See below.
6. **Scope and limits.** What is left for beyond-RTA.

### 5.2 The claim to make

Do not claim to beat EPW. Claim this instead:

> A single symmetry criterion, applied to the valley star, predicts which piezoresistance
> coefficients are accessible to rigid-valley transport and which are not. Where the criterion says
> "accessible," constant-τ RTA recovers the measured coefficient to within 5–20% at negligible cost.
> Where it says "inaccessible," the coefficient is not zero in nature but arises from a *different*
> and identifiable mechanism, which we isolate and quantify. Si and Ge are complementary test cases
> across a factor of ~5000 in shear-to-uniaxial response ratio.

A method with an adjustable τ can fit any of these numbers. A method that predicts zero and gets
zero, predicts large and gets large, and predicts *different temperature scalings for different
coefficients of the same tensor* is testing something. That is the asset.

### 5.3 The device layer — use torsion

The loading state that makes the Si–Ge contrast visible is **torsion**. A torsion bar is in pure
shear throughout its cross-section. Within rigid-valley transport, a Si gauge gives identically
zero response while a Ge gauge gives a large one — same code, same approximation, same mesh.

Then the mass-modulation correction restores a small, temperature-independent Si signal. So torsion
also becomes an experimental discriminator for the mechanism, not just for the material.

This complements rather than competes with your Si cantilever result: bending discriminates linear
versus exponential constitutive laws; torsion discriminates the two materials and the two
mechanisms.

### 5.4 Where this stands relative to prior art

- Ungersboeck et al. (2007) did the shear band structure and full-band Monte Carlo mobility, using
  **empirical pseudopotentials and k·p with fitted parameters**. They touched piezoresistance
  ("the limits of the linear piezoresistance model can be determined") but did not present a
  parameter-free π₄₄.
- Roisin/Poncé et al. (2024) did π₁₁ and π₁₂ from full first principles, and **explicitly state
  that π₄₄ is not accessible with uniaxial [100] strain in their setup**.
- Nobody, as far as this search goes, has done the parameter-free π₄₄ from DFT, the closed form,
  the temperature-scaling discrimination, or the Si–Ge symmetry pairing.

That is the gap. It is narrow but real, and it sits precisely where the leading group flagged
their own limit.

### 5.5 Venue

Not PRB for a Si-only result. With the Si–Ge pairing plus the device layer:
**J. Appl. Phys.**, **Journal of Computational Electronics**, **Solid-State Electronics**, or
**Computer Physics Communications** if you lead with the code. **SISPAD** is a fast, well-matched
conference and it is the Selberherr/Kosina community — which is also the Smirnov lineage you are
already benchmarking against.

---

## 6. Annotated literature

### Read first — these change what you do next

**Sverdlov, Ungersboeck, Kosina & Selberherr**, "Effects of Shear Strain on the Conduction Band in
Silicon: An Efficient Two-Band k·p Theory," *Proc. ESSDERC 2007*, 386–389.
→ https://www.iue.tuwien.ac.at/pdf/ib_2008/hashed_links/5d4PPhr.mQJunWmrCY_us.pdf
The four-page paper containing every formula in §1.3. Free PDF. **Read this today.**

**Ungersboeck, Dhar, Karlowatz, Sverdlov, Kosina & Selberherr**, "The Effect of General Strain on
the Band Structure and Electron Mobility of Silicon," *IEEE Trans. Electron Devices* **54**, 2183
(2007).
→ https://www.iue.tuwien.ac.at/pdf/ib_2008/JB2007_Ungersboeck_1.pdf
The journal version. Analytical expressions for both mass change and valley splitting, validated
against nonlocal EPM with spin-orbit. Includes full-band Monte Carlo and an explicit discussion of
where the linear piezoresistance model breaks down.

**Sverdlov, Karlowatz, Kosina & Selberherr**, "Two-band k·p model for the conduction band in
silicon."
→ https://www.iue.tuwien.ac.at/pdf/ib_2008/CP2007_Sverdlov_4.pdf
Source of $M \approx m_t/(1-m_t/m_0)$, and notes a weak stress dependence $D(\eta) = D + \beta\eta^2$
with β = 0.7 eV needed to improve agreement near η ≈ 1.

**Hensel, Hasegawa & Nakayama**, "Cyclotron resonance in uniaxially stressed silicon. II. Nature of
the covalent bond," *Phys. Rev.* **138**, A225 (1965).
The original source for the off-diagonal strain coupling and for $M$. Everything above traces back
here.

### The scoop, and your best conversation-opener

**Roisin, Brunin, Rignanese, Flandre, Raskin & Poncé**, "Phonon-limited mobility for electrons and
holes in highly-strained silicon," *npj Comput. Mater.* **10**, 242 (2024). Open access.
→ https://www.nature.com/articles/s41524-024-01425-0
Full IBTE with EPW-level e-ph coupling, plus their own four-point bending experiment. Electron
values: π₁₁ = −989, π₁₂ = +502 TPa⁻¹ (IBTE); SERTA −876/+475; MRTA −841/+514. Smith: −1022/+534.
**Your constant-τ result (−850/+430) sits essentially on top of their SERTA/MRTA** — worth a
sentence in the paper. Code and data on Materials Cloud. Note the Author Correction (Feb 2025).

### Temperature dependence — your discriminating test

**Morin, Geballe & Herring**, "Temperature dependence of the piezoresistance of high-purity silicon
and germanium," *Phys. Rev.* **105**, 525 (1957).
The benchmark for §3 Step 5. If Si's π₄₄ is flat in T while π₁₁ falls as 1/T, this is where you
check it.

### Ge functional problem

**QuantumATK Ge benchmark study** — https://docs.quantumatk.com/tutorials/germanium/germanium.html
Direct side-by-side of PBE, PBEsol, TB09-MGGA, pps-PBE, HSE on Ge band ordering, lattice constant,
and elastic constants. Concise and decisive.

**"An accurate alternative to hybrid functionals for germanium: DFT+α,"** arXiv:2512.08857 (Dec 2025).
→ https://arxiv.org/html/2512.08857
Shows HSE results for Ge are scattered even under identical settings and often fail to get Γ–L and
Γ–Γ right simultaneously. Read before committing to a functional.

### Foundational

- **Smith**, *Phys. Rev.* **94**, 42 (1954) — the measurements, both materials.
- **Herring & Vogt** (1956) — the many-valley transport theory.
- **Balslev**, *Phys. Rev.* **143**, 636 (1966) — deformation potentials for Si and Ge from the
  indirect absorption edge under shear. Si 295 K: Ξ_u = 9.2 ± 0.3 eV. Ge 80 K: Ξ_u = 16.2 ± 0.4 eV.
- **Kanda**, *IEEE TED* **29**, 64 (1982) — the P(N,T) doping/temperature factor.
- **Fischetti & Laux**, *J. Appl. Phys.* **80**, 2234 (1996) — deformation potentials and mobility
  in strained Si, Ge, SiGe.
- **Hamaguchi**, *Basic Semiconductor Physics* 4th ed., §6.3.9 and Table 6.1. ⚠ Ξ_u = 10.0 eV in
  §6.3.9 versus 9.0 ± 2 in Table 6.1 — check before citing. m_l/m_t are in **Ch. 2 p. 74**.
- **Grundmann**, *Physics of Semiconductors* — Table 8.3 (π for Si, Ge, GaAs), Table 6.8
  (deformation potentials Si and Ge, both Δ and L).

---

## 7. Risks and how a referee will push

| Risk | Response |
|---|---|
| "The two-band shear physics is Ungersboeck 2007" | True. Your contribution is the *parameter-free DFT* extraction of D, the closed form for π₄₄, the temperature-scaling discrimination, and the Ge pairing. Cite them prominently and early. |
| "π₁₁/π₁₂ was done better by Poncé" | Agree in the text. Frame constant-τ as the instrument that *isolates* the repopulation channel, not as a cheaper attempt at the same target. |
| "Your π₄₄ agreement to 2% is fortuitous" | It is. Show the sensitivity to D (π₄₄ ∝ D linearly), and lead with the temperature-scaling prediction, which is qualitative and therefore robust. |
| "Ge with PBE is meaningless" | Pre-empt with the mBJ/HSE spot-check on Ξ_u described in §4.4. |
| "The FEM layer is bolted on" | Torsion makes it load-bearing: it is the geometry where the symmetry statement becomes an engineering prediction. |

---

## 8. Sequence

**Now (shear channel, your stated priority)**
1. Re-inspect existing band structures — X-point position and splitting. Fix Γ₁₅ label. *(½ day)*
2. Extract D from shear sweep. *(2–3 days)* ← highest value per unit effort
3. Extract strained mass tensor; validate the two-band predictions against DFT. *(1 week)*
4. Generalise transport module to off-diagonal inverse-mass tensors. *(3–4 days)*
5. Temperature sweep, compare to Morin–Geballe–Herring. *(2 days)*

**Then (Ge)**
6. Functional go/no-go: mBJ + SOC, confirm L below Γ. *(2–3 days — binary outcome)*
7. Ξ_u^L PBE-vs-mBJ spot-check. *(2 days)*
8. Deformation potentials and elastic constants, both loading modes. *(1 week)*
9. Transport, compare all three coefficients to Smith. *(1 week)*

**Then (device)**
10. Torsion FEM, Si versus Ge, with and without mass modulation. *(1 week)*

**Throughout**
- Raise authorship and framing with Franchini *before* drafting.
- Even unfinished, the D extraction plus the 2×2 table is a concrete thing to put in front of Poncé
  ahead of the September–October FNRS promoter window.

---

*One closing note. The most valuable thing in this document is not the π₄₄ number — it is that*
*"rigid valley" and "constant-τ RTA" are independent approximations. Keep them separate in the*
*thesis text. Several of your defence notes currently conflate them, and a committee member who*
*knows the Sverdlov work will find that immediately.*
