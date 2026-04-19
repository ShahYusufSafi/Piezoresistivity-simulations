# Piezoresistivity Simulations

This project explores the piezoresistive effect, where mechanical deformation leads to changes in electrical conductivity and resistance in materials.

The main idea is to understand and model the relation between applied force, strain, and electrical response (voltage/current), starting from simple physical models and extending toward more realistic descriptions.

## Current Status

- Implemented a 1D finite element model for steady electrical conduction with constant conductivity  
- Verified numerical solution against the analytical solution (linear potential profile, constant field)  
- Established a baseline for further extensions  

## Next Steps

- Introduce strain-dependent conductivity: $\sigma = \sigma(\varepsilon(u))$  
- Couple mechanical deformation with electrical conduction  
- Extend to higher dimensions and more realistic geometries  
- Explore physically informed models for $\sigma(\varepsilon)$ based on transport theory  

## Structure

- `FEM/` – finite element implementation and test cases  
- `Particle Simulation/` – An initial geometric view of crystal structure

## Goal

The long-term goal is to build a computational framework that connects mechanical deformation, electronic transport, and continuum modeling in piezoresistive materials.