Hydrogen-air premixed turbulent combustion solver for OpenFOAM

- Progress variable approach
- TFC and ETFC models for RANS with Zimont, Bradley and Bray correlations
- Laminar burning velocity can be set by user (constant value) or estimated using Malet correlation (for lean mixtures only) or custom DNN model (for dry mixtures only)
- Charlette model for LES

General paper on in initial version:
flameFoam: An open source CFD solver for turbulent premixed combustion
https://www.sciencedirect.com/science/article/pii/S0029549321003137

Experiment simulations:
Simulation of Hydrogen-Air-Diluents Mixture Combustion in an Acceleration Tube with FlameFoam Solver
https://www.mdpi.com/1996-1073/14/17/5504

RANS- and TFC-Based Simulation of Turbulent Combustion in a Small-Scale Venting Chamber
https://www.mdpi.com/1996-1073/14/18/5710
