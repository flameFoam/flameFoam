# flameFoam

Hydrogen-air premixed turbulent combustion library for **OpenFOAM-13**.

- Progress-variable approach (`c` / `b`)
- TFC and ETFC models for RANS (Zimont, Bradley, Bray)
- FSD model for LES (Charlette wrinkling factor)
- Laminar burning velocity: Malet, power-law (`LBVPower`), or ANN (DNN)

flameFoam is a *combustion model library*, not a standalone solver. It is loaded into `foamRun` with `solver multicomponentFluid`.

flameFoam is not approved or endorsed by the OpenFOAM Foundation or OpenCFD.

## Requirements

- [OpenFOAM-13](https://openfoam.org/release/13/)
- Linux environment with `wmake` (WSL is fine)

## Project structure

```
flameFoam/
  flameFoam.C / flameFoam.H     Combustion model entry point
  reactionRateModels/           TFC, ETFC, FSD + LBV / TBV / wrinkling
  ThermophysicalTransportModels/  non-unity Lewis diffusivity models
  tutorial/                     Tutorial / test case
  Make/files, Make/options      wmake inputs (not the linux*Gcc* object dir)
  version.H                     Library version (13.0.0)
```

**Migration note:** the transport-model directory was renamed from `ThermopysicalTransportModels` to `ThermophysicalTransportModels`. Update any local scripts that used the old path.

## Build

From the library root (this directory):

```bash
wmake
```

The shared library is installed as `$FOAM_USER_LIBBIN/flameFoam-13.0.0.so`.

## Activate in a case

`system/controlDict`:

```
application     foamRun;
solver          multicomponentFluid;

libs
(
    "flameFoam-13.0.0.so"
);
```

`constant/combustionProperties`:

```
combustionModel flameFoam;

flameFoamCoeffs
{
    debug           false;   // also sent to all sub-models
    debugFields     false;

    X_H2O           0.0;
    X_H2_0          0.13;
    ...

    reactionRate
    {
        model   ETFC;

        laminarBurningVelocity
        {
            model   Malet;
        }

        turbulentBurningVelocity
        {
            model   Zimont;
            Zimont
            {
                ZimontA     0.52;
            }
        }
    }
}
```

To use the ANN laminar burning velocity, set `model ANN`. Training-scale constants and network weights stay in `ANN.C` and must not be changed.

Critical startup lines are written both to the solver log (`Info`) and to `flameFoam.<region>.combustionInfo`.

ETFC reads `Sct` from the already-registered `thermophysicalTransport` object (RAS or LES sub-dictionary).

## Tutorial / test case

```bash
cd tutorial
./Allrun
```

This runs `Allmesh`, then `Allrun-serial` (`foamRun`), then checks the log for:

- `flameFoam combustion model selected`
- `Reaction rate model: ETFC`

Expected log file: `tutorial/log.foamRun`.

Clean:

```bash
cd tutorial
./Allclean
```

Parallel (after mesh):

```bash
cd tutorial
./Allrun-parallel
```

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
