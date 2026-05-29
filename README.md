# flameFoam-0.12.3

Hydrogen-air premixed turbulent combustion solver for OpenFOAM v9:
- Progress variable approach
- TFC and ETFC models for RANS with Zimont, Bradley and Bray correlations
- Laminar burning velocity can be set by user (constant value) or estimated using Malet correlation (for lean mixtures only) or ANN model
- FSD model for LES with Charlette correlation

## Compilation
The solver is compiled using **wmake** command. [OpenFOAM-9](https://openfoam.org/release/9/) needs to be installed.

## Activation
The solver is executed using **flameFoam** command.

## Example case/tutorial
Tutorial case is available in the `tutorial` folder. The case demonstrates modeling flame acceleration in a tube with annular obstacles.

## Contributors
- Lead developer: Mantas Povilaitis, 2019 - current, mantas.povilaitis@lei.lt
- Co-developer: [Julius Venckus](https://github.com/jolonas), 2022 - current
- DNN model development and implementation: Andrius Ambrutis, 2022 - current
- Co-development of detonation model: [Stephen Adesina](https://github.com/Adstefnum), 2025
- Initial implementation of LBVPower class: [Mariia Nikolaieva](https://github.com/MariiaNikolaieva), 2025
- Porting to OpenFOAM-11: [Ilaryon Saladkou](https://github.com/IlaryonSaladkou), 2024
- Initial co-developer: Justina Jaseliūnaitė, 2019 - 2022

## Publications
- General paper on the initial version:
  - flameFoam: An open source CFD solver for turbulent premixed combustion https://www.sciencedirect.com/science/article/pii/S0029549321003137
- Experiment simulations:
  - The role of CFD combustion modelling in hydrogen safety management — IX: Validation of ETFC model implementation in flameFoam for large-scale hydrogen-air-steam deflagration https://www.sciencedirect.com/science/article/abs/pii/S0029549324008379
  - Simulation of Hydrogen-Air-Diluents Mixture Combustion in an Acceleration Tube with FlameFoam Solver https://www.mdpi.com/1996-1073/14/17/5504
  - RANS- and TFC-Based Simulation of Turbulent Combustion in a Small-Scale Venting Chamber https://www.mdpi.com/1996-1073/14/18/5710
- Presentation of the deep neural network developed for the estimation of laminar burning velocity:
  - Development of a CFD-Suitable Deep Neural Network Model for Laminar Burning Velocity https://www.mdpi.com/2076-3417/12/15/7460
