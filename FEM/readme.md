# FEM Formulation of Electrical Conduction (Static Case)

## Problem Statement

We consider the steady-state electrical conduction problem in 1D:

$$\nabla \cdot (\sigma \nabla V) = 0$$

In 1D this reduces to:

$$\frac{d}{dx}\left(\sigma \frac{dV}{dx}\right) = 0$$


## Boundary Conditions

We impose Dirichlet boundary conditions:

$$V(0) = V_0, \quad V(L) = 0$$

	•	These represent fixed voltages (electrodes).
	•	Known values → removed from system (modify matrices accordingly).



## Material Model
$\sigma$ = electrical conductivity
	
- First case:
$$\sigma = \text{constant}$$

Later extension:
$$\sigma = \sigma(\varepsilon(u))
\quad \text{(strain-dependent conductivity)}$$



## Weak (Variational) Formulation

- Start from:

$$\frac{d}{dx}\left(\sigma \frac{dV}{dx}\right) = 0$$

Multiply by test function v(x) and integrate:

$$\int_0^L \frac{d}{dx}\left(\sigma \frac{dV}{dx}\right) v \, dx = 0$$



- Integration by Parts

$$\int_0^L \sigma \frac{dV}{dx} \frac{dv}{dx} \, dx
= \int_{\partial \Omega} v \, (\sigma \nabla V \cdot n)$$

- For Dirichlet BCs (fixed voltage), boundary term vanishes:

$$\int_0^L \sigma \frac{dV}{dx} \frac{dv}{dx} \, dx = 0$$


- FEM Approximation (Galerkin).Approximate the solution:

$$V(x) \approx \sum_{j=1}^{N} V_j \phi_j(x)$$

Test functions:
$$v(x) = \phi_i(x)$$


## Substitution into Weak Form

$$\int_0^L \sigma \left( \sum_{j=1}^{N} V_j \frac{d\phi_j}{dx} \right) \frac{d\phi_i}{dx} dx = 0$$

Rearranging:

$$\sum_{j=1}^{N} V_j \int_0^L \sigma \frac{d\phi_j}{dx} \frac{d\phi_i}{dx} dx = 0$$


- Matrix Formulation

- Stiffness Matrix (Electrical)

$$K_{ij} = \int_0^L \sigma \frac{d\phi_i}{dx} \frac{d\phi_j}{dx} dx$$



## Linear System

$$K V = 0$$

After applying boundary conditions:

$$K_{\text{reduced}} V_{\text{unknown}} = f$$


## Interpretation


-  Recover Physical Quantities

After solving for V:

Electric Field

$$E = -\frac{dV}{dx}$$

Current Density

$$J = \sigma E = -\sigma \frac{dV}{dx}$$

⸻

- Extension: Coupling with Mechanics

Once displacement u(x) is known:

Strain:

$$\varepsilon(u) = \frac{du}{dx}$$

Conductivity becomes:

$$\sigma = \sigma(\varepsilon(u))$$


-  Updated PDE

$$\frac{d}{dx}\left(\sigma(\varepsilon(u)) \frac{dV}{dx}\right) = 0$$

⸻

### Important

Now:

$$\frac{d}{dx}(\sigma \frac{dV}{dx})
= \sigma \frac{d^2 V}{dx^2} + \frac{d\sigma}{dx} \frac{dV}{dx}$$

This introduces true coupling physics.

⸻

🧩 Final System Structure

Step 1: Solve Mechanics

$$\nabla \cdot \sigma_{\text{mech}}(u) = 0$$

Step 2: Compute strain

$$\varepsilon(u)$$

Step 3: Solve Electrical

$$\nabla \cdot (\sigma(\varepsilon(u)) \nabla V) = 0$$



The electrical problem becomes a material-dependent PDE, where deformation modifies conductivity, which in turn reshapes the electric field.

⸻

✅ Final Takeaway
	•	FEM reduces PDE → linear system
	•	Electrical conduction = Laplace-type problem
	•	Coupling enters through:
\sigma = \sigma(\varepsilon(u))
	•	This creates a multiphysics system (mechanics + electricity)

