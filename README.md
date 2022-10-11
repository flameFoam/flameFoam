# flameFoam
Hydrogen-air premixed turbulent combustion solver for OpenFOAM
- Progress variable approach
- TFC and ETFC models for RANS with Zimont, Bradley and Bray correlations
- Laminar burning velocity can be set by user (constant value) or estimated using Malet correlation (for lean mixtures only) or custom DNN model (for dry mixtures only)
- FSD model for LES with Charlette correlation

## Contributors
- Lead developer: Mantas Povilaitis, 2019 - current, mantas.povilaitis@lei.lt
- Co-developer: Julius Venckus, 2022 - current
- DNN model development and implementation: Andrius Ambrutis, 2022 - current
- Initial co-developer: Justina Jaseliūnaitė, 2019 - 2022

## Publications
- General paper on the initial version (Elsevier paywall):
  - flameFoam: An open source CFD solver for turbulent premixed combustion https://www.sciencedirect.com/science/article/pii/S0029549321003137
- Presentation of the deep neural network developed for the estimation of laminar burning velocity (Open access):
  - Development of a CFD-Suitable Deep Neural Network Model for Laminar Burning Velocity https://www.mdpi.com/2076-3417/12/15/7460
- Experiment simulations (Open access):
  - Simulation of Hydrogen-Air-Diluents Mixture Combustion in an Acceleration Tube with FlameFoam Solver https://www.mdpi.com/1996-1073/14/17/5504
  - RANS- and TFC-Based Simulation of Turbulent Combustion in a Small-Scale Venting Chamber https://www.mdpi.com/1996-1073/14/18/5710
