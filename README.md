# Bi-level Evolutionary Framework for Multi-Debris Active Removal

This repository contains the implementation code for two papers on bi-level optimization for
Multi-Debris Active Removal (MDAR) mission design — a **single-objective** framework and its
**multi-objective** extension. Both share one code base; the entry point you run selects which.

**Paper 1 — single-objective (published):**
Elad Denenberg, Adham Salih, **"Bi-level Evolutionary Framework for Designing Multi-Debris Active
Removal Missions"**, Acta Astronautica, Volume 246, 2026, Pages 247-257, ISSN 0094-5765,
https://doi.org/10.1016/j.actaastro.2026.04.002.

**Paper 2 — multi-objective (submitted):**
Adham Salih, Elad Denenberg, **"Multi-objective Bi-level Optimization Framework for Designing
Multi-Debris Active Removal Missions"**, submitted to Acta Astronautica.

---

## Overview

A bi-level evolutionary optimization approach for MDAR mission design. An **upper level** searches the
debris visiting sequence and maneuver timing; a **lower level** resolves each transfer's trajectory
(Δv) with realistic dynamics. The framework supports:

- **Single-objective (paper 1):** weighted-sum fitness (removed mass vs. total Δv), with three
  interchangeable lower-level optimizers — Random search, Genetic Algorithm (GA), and CMA-ES.
- **Multi-objective (paper 2):** an NSGA-II upper level returning the Pareto front of the two
  objectives — **maximize removed mass** and **minimize total Δv**.
- **High-fidelity force models:** zonal harmonics through J6 and atmospheric drag; the multi-objective
  runs additionally include **solar radiation pressure (SRP) and lunar third-body** perturbations.
- **Perturbed-Lambert lower level:** Thompson's perturbed-Lambert method (`prtlambertT`) integrated
  with a Modified Picard-Chebyshev propagator.
- **Scenarios:** the Iridium-33 debris cloud and a synthetic cloud.

---

## MATLAB Requirements

- MATLAB R2020b or later
- Parallel Computing Toolbox (optional, for faster computation — the lower level uses `parfor`)

---

## Installation

### 1. Clone the repository
```bash
git clone https://github.com/BraudeMechInd4-0/MDAR.git
cd MDAR
```

### 2. Initialize submodules
The project uses two submodules (see `.gitmodules`):
- `ParseGP/` — TLE / GP catalog parsing utilities (small)
- `fundamentals-of-astrodynamics/` — Vallado's astrodynamics library (provides the Lambert seed
  `lambertb` and related utilities)

The Vallado repository is **very large** — it ships C++/Python/Fortran sources, STK data, and test
datasets we don't need. Use a **sparse checkout** to pull only its MATLAB sources
(`software/matlab/`). Requires **Git 2.25+** (on Windows, run these in **Git Bash**):

```bash
# ParseGP is small — initialize it normally
git submodule update --init ParseGP

# Vallado — register without checking out, then sparse-checkout only the MATLAB sources
git submodule update --init --no-checkout fundamentals-of-astrodynamics
cd fundamentals-of-astrodynamics
git sparse-checkout init --cone
git sparse-checkout set software/matlab
git checkout
cd ..
```

This leaves only `fundamentals-of-astrodynamics/software/matlab/` on disk — which the scripts add to
the path via `addpath(genpath("fundamentals-of-astrodynamics"))`. (The `Satellites/` folder ships
with the repo and holds the Iridium-33 TLE data; it is not a submodule.)

---

## Usage

### Step 1 — Generate debris propagation models
Before any optimization, build the fast Chebyshev propagation models for the debris cloud:
```matlab
createModelsScript
```
Parses the catalog, integrates each object under the selected force model, fits Chebyshev
polynomials, and saves the models under `models/`. This step can take minutes to hours depending on
cloud size and mission horizon.

### Step 2 — Run an optimizer (pick one)

**Multi-objective (paper 2):**
```matlab
MDAR_NSGAII_main
```
NSGA-II upper level, two objectives (removed mass, total Δv). Runs the Iridium-33 and synthetic-cloud
scenarios under the J6+Drag+SRP+Moon force model, repeated for statistical pooling into a reference
Pareto front. Saves per-run results (non-dominated set + per-generation logs).

**Single-objective (paper 1):**
```matlab
HiLevel
```
Weighted-sum fitness with a selectable lower-level optimizer (`LowLevelType = "Rand" | "GA" | "CMA"`).

**Weighted-sum sweep baseline (paper 2 comparison):**
```matlab
HiLevelVaringM
```
Sweeps the objective weights to trace a front via scalarization.

**Key upper-level parameters** (editable at the top of each main): `PopSizeHi`, `MaxGenHi`,
`NtoRemove`, `T1Max` (max inter-maneuver wait), and — single-objective only — the weights `wM`, `wDV`.

---

## Reproducing paper 1 (original single-objective results)

Paper 1 used the **J6+drag** force model and the original Battin/Thompson Lambert settings. The shipped
code standardizes on the newer physics and the consolidated `prtlambertT` solver, so to reproduce the
paper-1 numbers:

1. **Force model** — in the lower-level optimizers use `orbit_eq_J6_drag` instead of
   `orbit_eq_J6_drag_SRP_moon`. `LowLevelOptimizationGA` and `LowLevelOptimizationCMA` already use
   `orbit_eq_J6_drag`; change the handle in `LowLevelOptimizationRandom` to match.
2. **Solver settings** — set `options.maxIter = 20;` and `options.delta = 8;` in the lower-level
   evaluators (`options.Tolr = 1e-5` is already the default). With these, `prtlambertT` reproduces the
   behavior of the original `prtlambertbf` solver (same Thompson algorithm; the newer defaults are 100
   iterations / `delta` 16).

---

## Repository structure

```
.
├── MDAR_NSGAII_main.m        # Paper 2: multi-objective (NSGA-II) main
├── HiLevel.m                 # Paper 1: single-objective main (Rand/GA/CMA)
├── HiLevelVaringM.m          # Weighted-sum sweep baseline
├── createModelsScript.m      # Debris propagation-model generation
│
├── Upper-level operators/
│   ├── initPopHi.m           # Population initialization
│   ├── EvaluateHi.m          # Single-objective evaluation
│   ├── EvaluateHiVarM.m      # Multi-objective (varying-mass) evaluation
│   ├── tournamentSelection.m, SelectionByRank.m        # selection (SO / NSGA-II)
│   ├── CrossOver.m, Crossover_Ordered_Operator.m, SBX.m # crossover
│   ├── Mutation.m, PolyMutation.m, ExchangeMutation.m   # mutation
│   ├── EliteProcedure.m, EliteFullSorting.m             # elitism (SO / NSGA-II)
│   ├── NDSort.m, CrowdingDistance.m, CalcRankAndDistance.m  # NSGA-II ranking (PlatEMO)
│
├── Lower-level optimizers/
│   ├── LowLevelOptimizationRandom.m, LowLevelOptimizationGA.m,
│   │   LowLevelOptimizationCMA.m, LowLevelOptimizationOneVar.m
│   └── EvaluateModel.m       # per-leg cost evaluation
│
├── Orbital mechanics/
│   ├── orbit_eq_J6_drag_SRP_moon.m  # J6 + drag + SRP + Moon (paper 2)
│   ├── orbit_eq_J6_drag.m, orbit_eq_J2_drag.m  # J6/J2 + drag
│   ├── drag_accel.m, moon_v.m, sun_v.m         # perturbation models / ephemerides
│   ├── prtlambertT.m, prtlambertTRB.m, prtlambertMPS.m, prtlambertMPSRB.m  # Lambert solvers
│   ├── odeMPCI.m, odeMPCIrev.m  # Modified Picard-Chebyshev integrators (forward / backward)
│   ├── propCheb.m, kep_elements.m, hitearth.m
│
└── Submodules/
    ├── ParseGP/                        # TLE / GP parsing
    └── fundamentals-of-astrodynamics/  # Vallado astrodynamics utilities
```

---

## Data files

- **`Satellites/Iridium.xml`** — Iridium-33 TLE data (parsed by ParseGP).
- **`models/*.mat`** — generated Chebyshev propagation models (produced by `createModelsScript.m`;
  gitignored — regenerate locally, do not commit).
- **Run outputs** (`*Run*.mat`, `*VM_MO_*.mat`) — per-run results (best/non-dominated solutions plus
  logs); gitignored.

---

## Citations

If you use this code, please cite the relevant work:

**Framework (paper 1, single-objective):**
Elad Denenberg, Adham Salih, "Bi-level Evolutionary Framework for Designing Multi-Debris Active
Removal Missions," Acta Astronautica, 246 (2026) 247-257. https://doi.org/10.1016/j.actaastro.2026.04.002

**Framework (paper 2, multi-objective):**
Adham Salih, Elad Denenberg, "Multi-objective Bi-level Optimization Framework for Designing
Multi-Debris Active Removal Missions," submitted to Acta Astronautica.

**Perturbed-Lambert solvers** (`prtlambertT`, `prtlambertTRB`, `prtlambertMPS`, `prtlambertMPSRB`,
`odeMPCIrev`):
Elad Denenberg, "Improved Convergence and Efficiency of Perturbed Lambert Solvers via Backward
Chebyshev-Picard Integration," in review, Advances in Space Research.

**NSGA-II ranking utilities** (`NDSort`, `CrowdingDistance`):
Ye Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, "PlatEMO: A MATLAB platform for evolutionary
multi-objective optimization [educational forum]," IEEE Computational Intelligence Magazine, 2017,
12(4): 73-87.

---

## Acknowledgments / third-party components

- **PlatEMO** (Tian et al., 2017) — non-dominated sorting and crowding distance (`NDSort`,
  `CrowdingDistance`).
- **Vallado's astrodynamics utilities** (Vallado, 2022) — via the `fundamentals-of-astrodynamics`
  submodule (Lambert seed, etc.).
- **Modified Picard-Chebyshev Integration** (Woollands & Junkins, 2019; Bai & Junkins, 2011).
- The multi-objective force model (`orbit_eq_J6_drag_SRP_moon`, `moon_v`, `sun_v`) and the
  perturbed-Lambert solvers are the work of E. Denenberg.

---

## License

Distributed under the **GNU Affero General Public License v3.0** — see [`LICENSE`](LICENSE).

---

## Contact

- Elad Denenberg — eladd@braude.ac.il
- Adham Salih — adhamsalih@braude.ac.il
