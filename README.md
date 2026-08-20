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
  reactionRateModels/           TFC, ETFC, FSD, DDT + LBV / TBV / wrinkling / auto-ignition
  ThermophysicalTransportModels/  non-unity Lewis diffusivity models
  tutorial/                     Smoke-test case (OpenFOAM tutorial layout)
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

ANN optional coefficients (defaults match the published trained network). Layer sizes follow the compiled weight matrices; the weights themselves stay in `ANN.C`.

```
laminarBurningVelocity
{
    model   ANN;
    ANN
    {
        pRef        3970000;
        ERRef       7.16;
        TRef        864;
        X_H2_dryAir 0.705;
        X_O2_dryAir 0.295;
    }
}
```

Critical startup lines are written both to the solver log (`Info`) and to `flameFoam.<region>.combustionInfo`. Model-selection text is collected on `reactionRate` (the old `infoPass` helper class was removed).

ETFC / aITransport read `Sct` from the already-registered `thermophysicalTransport` object (RAS or LES sub-dictionary). They do not construct a second `IOdictionary`.

## DDT / auto-ignition table

```
reactionRate
{
    model   DDT;
    tIgn    1.5e-4;

    wrinklingFactor
    {
        model   Charlette;
    }

    autoIgnition
    {
        model   aITransport;
        aITransport
        {
            ADTDir  constant/ADT;
        }
    }
}
```

Put one `*.ADT` file per pressure in `constant/ADT/` (filename = pressure, e.g. `100000.ADT`; 4 header lines then `T  tau` columns). Lookup uses bilinear interpolation in (p, T), not integer rounding.

## Tutorial / smoke test

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

- flameFoam: An open source CFD solver for turbulent premixed combustion https://www.sciencedirect.com/science/article/pii/S0029549321003137
- Validation of ETFC model implementation in flameFoam https://www.sciencedirect.com/science/article/abs/pii/S0029549324008379
- Simulation of Hydrogen-Air-Diluents Mixture Combustion in an Acceleration Tube https://www.mdpi.com/1996-1073/14/17/5504
- RANS- and TFC-Based Simulation of Turbulent Combustion in a Small-Scale Venting Chamber https://www.mdpi.com/1996-1073/14/18/5710
- Development of a CFD-Suitable Deep Neural Network Model for Laminar Burning Velocity https://www.mdpi.com/2076-3417/12/15/7460
