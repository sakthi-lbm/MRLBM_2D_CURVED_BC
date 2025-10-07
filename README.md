# Moments-Based Regularized LBM with Curved Boundary Conditions

This repository contains a CUDA implementation of a **Moments-Based Regularized Lattice Boltzmann Method (MRLBM)** featuring a novel **curved boundary treatment** for accurate and efficient simulation of flows past complex geometries on Cartesian grids.

## 🧩 Project Overview

We develop a robust curved boundary condition for the 2D MRLBM, combining the computational efficiency of uniform Cartesian grids with accurate curved-boundary representation through **moment regularization** and **extrapolation**.

### ✳️ Key Features:
- **Memory-Efficient MRLBM**: Stores macroscopic moments instead of full distribution functions, reducing memory usage by **33% (2D)** compared to conventional LBM.
- **Novel Curved Boundary Treatment**: Implements two methods:
  - **Method A**: Cartesian coordinate approach.
  - **Method B (Recommended)**: Rotated coordinate system aligned with boundary normals for enhanced stability.
- **Extrapolation-Based Technique**: Uses quadratic extrapolation to impose wall conditions accurately on stair-case approximated surfaces.
- **Validated Benchmark**: Flow past a circular cylinder at Re = 10-200.

### ✳️ Governing Model
- **Lattice:**  D2Q9 model
- **Collision Operator :** Regularized single-relaxation-time (BGK-type) in moment space
- **Boundary Conditions :**
    - *Inlet :* Constant velocity
    - *Outlet :* Neumann (zero-gradient) condition
    - *Cylinder :* No-slip curved boundary condition
    - *Walls :* No-slip or periodic

## 📁 Repository Structure
```text
MR_LBM/
├── post/                             # Post-processing scripts or analysis tools
│
├── src/                              # Source files for MRLBM solver
│   ├── cases/                        # Case setup files (boundary/geometry definitions)
│   ├── colrec/                       # Color/recording utilities (if any)
│   ├── includeFiles/                 # Shared include files
│   │
│   ├── arrayIndex.h                  # Indexing macros and utilities
│   ├── boundaryCondition.cuh         # Boundary condition implementations (GPU)
│   ├── compile.sh                    # Build script
│   ├── cylinder_mrlbm.cu             # Flow past a cylinder main case
│   ├── definitions.h                 # Global constant and type definitions
│   ├── errorDef.h                    # Error handling macros
│   ├── globalFunctions.h             # Global utility functions
│   ├── globalStructs.h               # Global data structures
│   ├── lbmInitialization.cu          # Initialization routines (GPU)
│   ├── lbmInitialization.cuh         # Header for initialization kernels
│   ├── main.cu                       # Main driver file
│   ├── main.cuh                      # Header for main solver
│   ├── mlbm.cu                       # MRLBM solver implementation
│   ├── mlbm.cuh                      # Header for MRLBM kernels
│   ├── nodeTypeMap.h                 # Node-type mapping (fluid/solid/boundary)
│   ├── saveData.cu                   # Data output routines
│   ├── saveData.cuh                  # Header for output routines
│   ├── treat_data.cu                 # Post-treatment or averaging kernels
│   ├── treat_data.cuh                # Header for post-processing kernels
│   └── var.h                         # Variable declarations and parameters
│
├── x64/                              # Build directory (compiled objects/binaries)
│
├── .gitignore                        # Git ignore rules
└── sim_D2Q9_sm86                     # Compiled executable (D2Q9 GPU architecture)


