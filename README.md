# Moments-Based Regularized LBM with Curved Boundary Conditions

This repository contains a CUDA implementation of a **Moments-Based Regularized Lattice Boltzmann Method (MRLBM)** featuring a novel **curved boundary treatment** for accurate and efficient simulation of flows past complex geometries on Cartesian grids.

## 🧩 Project Overview

We develop a robust curved boundary condition for the 2D MRLBM, combining the computational efficiency of uniform Cartesian grids with accurate curved-boundary representation through **moment regularization** and **extrapolation**.

### Key Features:
- **Memory-Efficient MRLBM**: Stores macroscopic moments instead of full distribution functions, reducing memory usage by **33% (2D)** compared to conventional LBM.
- **Novel Curved Boundary Treatment**: Implements two methods:
  - **Method A**: Cartesian coordinate approach.
  - **Method B (Recommended)**: Rotated coordinate system aligned with boundary normals for enhanced stability.
- **Extrapolation-Based Technique**: Uses quadratic extrapolation to impose wall conditions accurately on stair-case approximated surfaces.
- **Validated Benchmark**: Flow past a circular cylinder at Re = 10-200.

🧩 Governing Model
- Lattice: D2Q9 model
- Collision Operator: Regularized single-relaxation-time (BGK-type) in moment space
- **Boundary Conditions:**
    - Inlet: Constant velocity
    - Outlet: Neumann (zero-gradient) condition
    - Cylinder: No-slip curved boundary condition
    - Walls: Bounce-back or symmetry

## 📁 Repository Structure



```text
.
├── src/                 # Source code
│   ├── core/           # Core LBM classes (collision, streaming)
│   ├── boundary/       # Boundary condition implementations
│   ├── moments/        # Moments-based regularization
│   └── utils/          # Helper functions (geometry, IO)
├── examples/           # Example simulations
│   └── cylinder_flow/ # Primary validation case
├── docs/              # Documentation
└── results/           # Output data and visualization scripts

