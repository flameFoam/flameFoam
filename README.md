# flameFoam

Hydrogen-air premixed turbulent combustion **solver** for **OpenFOAM-9**.

- Progress-variable approach (`c`)
- TFC and ETFC models for RANS (Zimont, Bradley, Bray correlations)
- FSD model for LES (Charlette or Pitsch-Duchamp wrinkling factor)
- Laminar burning velocity: constant, Malet (lean mixtures), or ANN correlations
- Quenching and wall-quenching models (RANS)
- Multi-region capable (fluid + solid)

flameFoam is a standalone solver (not a library). It is executed directly with
the `flameFoam` command after compilation.

flameFoam is not approved or endorsed by the OpenFOAM Foundation or OpenCFD.

## Requirements

- [OpenFOAM-9](https://openfoam.org/release/9/)
- Linux environment with `wmake`

## Project structure

```
flameFoam/
├── flameFoam.C                 Solver main (time loop, multi-region control)
├── EEqn.H                      Energy transport equation
├── UEqn.H                      Momentum transport equation
├── pEqn.H                      Pressure transport equation
├── cEqn.H                      Progress variable (c) transport equation
├── createFields.H              Top-level field creation include
├── createMeshes.H              Top-level mesh creation include
├── createMeshesPostProcess.H   Post-processing mesh helpers
├── initCombustion.H            Combustion model setup, parameter reading, field
|                               initialization
├── combustion/
│   ├── RANS.H                  TFC / ETFC reaction-rate models + Turbulent
|   |                           burning velocity correlations + quenching and
|   |                           wall-quenching models
│   ├── LES.H                   FSD model + Wrinkling factor correlations
│   └── ANN.H                   Artificial neural network laminar burning
|                               velocity estimation
├── fluid/
│   ├── compressibleCourantNo.C/.H          Courant-number calculation
│   ├── compressibleMultiRegionCourantNo.H  Multi-region Courant helper
│   ├── createFFields.H                     Fluid-region field creation
│   ├── createFluidMeshes.H                 Fluid mesh creation
│   ├── setRegionFFields.H                  Per-region fluid field pointers
│   └── solveFluid.H                        Fluid-region solve sequence
├── solid/
│   ├── createSolidFields.H                 Solid-region field creation
│   ├── createSolidMeshes.H                 Solid-region mesh creation
│   ├── readSolidTimeControls.H             Solid-region time-step controls
│   ├── setRegionSolidFields.H              Per-region solid field pointers
│   ├── solidRegionDiffNo.C/.H              Diffusion number calculation
│   ├── solidRegionDiffusionNo.H            Diffusion-number helper
│   └── solveSolid.H                        Solid-region solve sequence
├── include/
│   ├── setInitialMultiRegionDeltaT.H       Initial multi-region time-step
│   └── setMultiRegionDeltaT.H              Multi-region time-step adjustment
├── tutorial/                   Example case (flame-acceleration tube)
│   ├── 0.orig/Fluid/           Initial fields (T, U, p, c, ...)
│   ├── constant/Fluid/         Modeling options, parameters
│   ├── system/                 Solution control, OpenFOAM utility dicts
│   ├── Allmesh                 Mesh generation
│   ├── Allrun-serial           Runs simulation in one core
│   ├── Allrun-parallel         Runs simulation on multiple cores
│   └── Allclean                Cleans simulation results from the case
├── Make/
│   ├── files                   Source list and compiled executable name
│   └── options                 Include paths and linked libraries
├── CHANGELOG.md                Version history
├── LICENSE                     GPL-3.0 licence definition
└── README.md                   This file
```

## Build

From the solver root (this directory):

```bash
wmake
```

The executable is installed as `$FOAM_USER_APPBIN/flameFoam`.

## Execute program

In the case directory, simply run:

```bash
flameFoam
```

(or use the provided `Allrun` scripts in the tutorial).

## Tutorial

The tutorial case is a flame acceleration tube with annular obstacles, filled
with premixed 13% hydrogen-air mixture, demonstrating ETFC + Zimont + Malet
simulation.

Create the mesh:

```bash
cd tutorial
./Allmesh
```

Run simulation in serial

```bash
./Allrun-serial
```

or parallel:

```bash
./Allrun-parallel
```

Expected behaviour: progress variable `c` advances, temperature/pressure rise,
sample p and T data is written to postProcessing/.

Delete the simulation results:

```bash
cd tutorial
./Allclean
```

## Contributors

- Lead developer: Mantas Povilaitis, 2019 – current, mantas.povilaitis@lei.lt
- Co-developer: [Julius Venckus](https://github.com/jolonas), 2022 – current
- DNN model development and implementation: Andrius Ambrutis, 2022 – current
- Initial co-developer: Justina Jaseliūnaitė, 2019 – 2022

## Publications

## Publications
- General paper on the initial version:
  - flameFoam: An open source CFD solver for turbulent premixed combustion
  https://www.sciencedirect.com/science/article/pii/S0029549321003137
- Experiment simulations:
  - The role of CFD combustion modelling in hydrogen safety management — IX:
  Validation of ETFC model implementation in flameFoam for large-scale
  hydrogen-air-steam deflagration
  https://www.sciencedirect.com/science/article/abs/pii/S0029549324008379
  - Simulation of Hydrogen-Air-Diluents Mixture Combustion in an Acceleration
  Tube with FlameFoam Solver https://www.mdpi.com/1996-1073/14/17/5504
  - RANS- and TFC-Based Simulation of Turbulent Combustion in a Small-Scale
  Venting Chamber https://www.mdpi.com/1996-1073/14/18/5710
- Presentation of the deep neural network developed for the estimation of
laminar burning velocity:
  - Development of a CFD-Suitable Deep Neural Network Model for Laminar Burning
  Velocity https://www.mdpi.com/2076-3417/12/15/7460
```
