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

### 🧮 Model Description
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
├── post/                                 # Post-processing and visualization scripts
│
├── src/                                  # Source files for the MRLBM solver
│   ├── cases/                            # Simulation cases
│   │   └── cylinder/                     # Flow past a circular cylinder
│   │       ├── aux_functions.cuh         # Auxiliary GPU device functions
│   │       ├── boundaries.cuh            # Curved boundary condition implementations
│   │       ├── constants.h               # Physical and numerical constants
│   │       ├── interpolation_utilities.cuh # Interpolation and geometry utilities
│   │       ├── numerical_solutions.cuh   # Numerical solver functions
│   │       └── outputs_and_model.h       # Output and model configuration
│   │
│   ├── colrec/                           # Collision and reconstruction modules
│   │   ├── 2nd_order/
│   │   │   └── collision_and_reconstruction.cuh # Second-order regularization model
│   │   ├── 3rd_order/                    # Third-order regularization model (future)
│   │   └── 4th_order/                    # Fourth-order regularization model (future)
│   │
│   ├── includeFiles/                     # Common include files
│   │   ├── interface_handling.cuh        # Interface and memory handling utilities
│   │   ├── interface.h                   # Global interface definitions
│   │   ├── popLoad.inc                   # Population loading routines
│   │   └── popSave.inc                   # Population saving routines
│   │
│   ├── arrayIndex.h                      # Macros for indexing and grid mapping
│   ├── boundaryCondition.cuh             # Generic boundary condition functions
│   ├── compile.sh                        # Build script for compilation with NVCC
│   ├── cylinder_mrlbm.cu                 # Cylinder flow main case
│   ├── definitions.h                     # Global constants and definitions
│   ├── errorDef.h                        # Error definitions and macros
│   ├── globalFunctions.h                 # General-purpose GPU/CPU functions
│   ├── globalStructs.h                   # Structs for grid, field, and LBM data
│   ├── lbmInitialization.cu              # Initialization routines
│   ├── lbmInitialization.cuh             # Headers for initialization kernels
│   ├── main.cu                           # Main driver file
│   ├── main.cuh                          # Header for main solver
│   ├── mlbm.cu                           # MRLBM collision/streaming implementation
│   ├── mlbm.cuh                          # Header for MRLBM kernels
│   ├── nodeTypeMap.h                     # Node-type (fluid/solid/boundary) mapping
│   ├── saveData.cu                       # Output and file writing routines
│   ├── saveData.cuh                      # Header for output routines
│   ├── treat_data.cu                     # Post-treatment/averaging kernels
│   ├── treat_data.cuh                    # Header for post-treatment routines
│   └── var.h                             # Variable declarations and parameters
│
├── x64/                                  # Compiled objects and executables
│
├── .gitignore                            # Git ignore file
└── sim_D2Q9_sm86                         # Compiled binary for D2Q9 (sm_86 GPU)
```
## 🛠️ Compilation & Installation
git clone https://github.com/sakthi-lbm/MRLBM_2D_CURVED_BC.git
cd MR_LBM/src
bash compile.sh
./../sim_D2Q9_sm86


## 📊 Output Files & Visualization
```text
The simulation generates:
    - CYLINDER/SIM_ID/grid.x: Computational grid geometry
    - CYLINDER/SIM_ID/data_XXX.f: Time-dependent flow fields (bindary format)
    - CYLINDER/SIM_ID/master.p3d: Paraview master file to load all the output files including grid.x
    
CYLINDER/SIM_ID/: also Contains:
    - forces_SIM_ID.dat: Time history of drag and lift forces
    - pressure_SIM_ID.dat: Time history of surface pressure at different theta
    - rho_inlet_SIM_INLET.dat: Inlet average density (free-stream density)
```

## 📜 How to Cite
If you use or modify this solver in your research, please cite the following article:

*Sakthivel, M., B.Y. dos Anjos, L.A. Hegele Jr. A moments-based regularized lattice Boltzmann method (MRLBM):
Incompressible curved boundary implementation and validation in two-dimensions. Physics of Fluids (2025).*

```text
bibtex:

@misc{Sakthivel2025MRLBM,
  author = {Sakthivel, M., B.Y. dos Anjos, L.A. Hegele Jr},
  title = {Moment-based Regularized Lattice Boltzmann Solver (MRLBM) for Flow Past a Cylinder},
  year = {2025},
}
```



## 🧩 Future Work
- Extension to 3D (D3Q19, D3Q27) lattice models
- Implementation of higher-order regularization (4th-order and 6th order)
- Implementation of *cut-cell* approach for the curved boundary




