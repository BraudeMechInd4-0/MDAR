# Bi-level Evolutionary Framework for Multi-Debris Active Removal

This repository contains the implementation code for the paper:

**"Bi-level Evolutionary Framework for Designing Multi-Debris Active Removal Missions"**  
By Elad Denenberg and Adham Salih

*Full citation will be added upon publication.*

---

## Overview

This framework implements a bi-level evolutionary optimization approach for Multi-Debris Active Removal (MDAR) mission design. The framework supports:

- **High-fidelity force models**: Zonal harmonics up to J6 and atmospheric drag
- **Bi-level optimization**: Upper level (debris sequence + timing) and lower level (trajectory parameters)
- **Multiple low-level optimizers**: Random search, Genetic Algorithm (GA), and CMA-ES
- **Direct transfer architecture**: Proof-of-concept on Iridium-33 debris cloud

---

## Installation

### 1. Clone the Repository

```bash
git clone <repository-url>
cd <repository-name>
```

### 2. Initialize Submodules

This project requires several submodules for parsing, orbital mechanics utilities, and Lambert solvers:

```bash
git submodule init
git submodule update
```

The required submodules are:
- `ParseGP/` - TLE parsing utilities
- `Satellites/` - Satellite data storage
- `vallado/` - Vallado's astrodynamics utilities (includes Lambert solvers)

### 3. MATLAB Requirements

- MATLAB R2020b or later
- Parallel Computing Toolbox (optional, for faster computation)

---

## Usage

### Step 1: Generate Debris Propagation Models

Before running optimization, you must generate Chebyshev propagation models for the debris cloud:

```matlab
createModelsScript
```

**What this does:**
- Parses TLE data from `Satellites/Iridium.xml` using ParseGP utilities
- Generates initial states (positions, velocities, masses, cross-sections) for debris
- Integrates debris trajectories using J6+drag dynamics over 2-year mission horizon
- Fits Chebyshev polynomials for fast propagation
- Saves models to `Satellites/Iridium33Model.mat`

**Requirements:**
- TLE data file: `Satellites/Iridium.xml` (Iridium-33 debris cloud)
- Reference epoch and duration are configurable in the script

**Note:** This step may take several minutes to hours depending on debris count and mission duration.

### Step 2: Run Bi-level Optimization

Execute the main optimization framework:

```matlab
HiLevel
```

**What this does:**
- Loads pre-computed models from `Satellites/Iridium33Model.mat`
- Runs bi-level evolutionary optimization with three low-level solver types:
  - Random search
  - Genetic Algorithm (GA)
  - CMA-ES (Covariance Matrix Adaptation Evolution Strategy)
- Performs 3 independent runs per solver type
- Saves results to files: `Run1Rand.mat`, `Run1GA.mat`, `Run1CMA.mat`, etc.

**Key parameters** (configurable in `HiLevel.m`):
- `PopSizeHi = 30` - High-level population size
- `MaxGenHi = 50` - Number of generations
- `NtoRemove = 10` - Number of debris to remove
- `T1Max = 2*60*60*24*30` - Maximum wait time between maneuvers (≈60 days)
- `wM = 1` - Weight for removed mass
- `wDV = 8` - Weight for total ΔV

---

## Repository Structure

```
.
├── HiLevel.m                          # Main optimization script
├── createModelsScript.m               # Model generation script
│
├── High-Level GA Components/
│   ├── initPopHi.m                    # Population initialization
│   ├── EvaluateHi.m                   # High-level evaluation
│   ├── tournamentSelection.m          # Selection operator
│   ├── CrossOver.m                    # Crossover coordinator
│   ├── Crossover_Ordered_Operator.m   # Ordered crossover (OX)
│   ├── Mutation.m                     # Mutation coordinator
│   ├── PolyMutation.m                 # Polynomial mutation
│   ├── ExchangeMutation.m             # Exchange mutation
│   ├── EliteProcedure.m               # Elitism handling
│   └── SBX.m                          # Simulated binary crossover
│
├── Low-Level Optimizers/
│   ├── LowLevelOptimizationRandom.m   # Random search
│   ├── LowLevelOptimizationGA.m       # GA optimizer
│   └── LowLevelOptimizationCMA.m      # CMA-ES optimizer
│
├── Orbital Mechanics/
│   ├── orbit_eq_J6_drag.m             # J6+drag dynamics
│   ├── orbit_eq_J2_drag.m             # J2+drag dynamics
│   ├── prtlambertbf.m                 # Perturbed Lambert solver
│   ├── hitearth.m                     # Earth collision detection
│   ├── odeMPCI.m                      # Modified Picard-Chebyshev integrator
│   ├── propCheb.m                     # Chebyshev propagation
│   ├── kep_elements.m                 # Keplerian element conversion
│   └── drag_accel.m                   # Atmospheric drag model
│
└── Submodules/
    ├── ParseGP/                       # TLE parsing
    ├── Satellites/                    # Debris data (.mat files)
    └── vallado/                       # Astrodynamics utilities (Lambert, etc.)
```

---

## Data Files

### Required Input Data

**`Satellites/Iridium.xml`** - TLE data for Iridium-33 debris cloud
- Standard XML format from space-track.org or similar sources
- Parsed automatically by `generateList()` function in ParseGP utilities

### Generated Output Files

- **`Satellites/Iridium33Model.mat`** - Pre-computed Chebyshev propagation models
  - Generated by `createModelsScript.m`
  - Contains: `model` structure with Chebyshev coefficients and time segments
- **`Run<N><Type>.mat`** - Optimization results for each run and solver type
  - Contains: `BestSol`, `wDV`, `wM`, `LowLevelType`, `GenParam`

---

## Key Features

### Force Model
- **Zonal harmonics**: J2 through J6 (Earth oblateness)
- **Atmospheric drag**: Exponential density model (Vallado 2022, Table A.6)
- **Extensible**: Framework supports additional perturbations (SRP, third-body)

### Optimization Framework
- **Bi-level structure**: Decouples sequence planning from trajectory optimization
- **Modular design**: Easy to swap low-level optimizers
- **Realistic dynamics**: No analytical approximations during optimization

### Mission Architecture
- **Direct transfer**: Impulsive maneuvers between debris
- **Bi-impulse transfers**: Departure + arrival velocity matching
- **Flexible timing**: Optimizes wait times between maneuvers

---

## Computational Notes

- **Parallel evaluation**: Low-level optimizers use `parfor` (requires Parallel Computing Toolbox)
- **Runtime**: Full 3×3 runs (3 solvers × 3 repetitions) may take several hours
- **Memory**: Pre-computing models requires sufficient RAM (~1-2 GB for 98 debris over 2 years)

---

## Citation

If you use this code in your research, please cite:

```bibtex
@article{denenberg2025mdar,
  title={Bi-level Evolutionary Framework for Designing Multi-Debris Active Removal Missions},
  author={Denenberg, Elad and Salih, Adham},
  journal={Acta Astronautica},
  year={2026},
  note={In press}
}
```

*(Citation details will be updated upon publication)*

---

## License

[Add your license here]

---

## Contact

For questions or issues, please contact:
- Elad Denenberg: eladd@braude.ac.il
- Adham Salih: adhamsalih@braude.ac.il

---

## Acknowledgments

This research uses:
- Modified Picard-Chebyshev Integration (Woollands & Junkins, 2019; Bai & Junkins, 2011)
- Perturbed Lambert solver (Thompson et al., 2018)
- Vallado's astrodynamics utilities (Vallado, 2022)