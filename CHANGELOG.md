# Changelog

## [12.3.0] - 2025-05-02
### Added
- Initial version of generalized LBV power law correlation LBVPower with user-defined powers

## [12.2.0] - 2024-11-07
### Changed
- Added new switches to control optional correction of unburnt mixture temperature and debug field writing
- Implemented support for calculation restart from non-zero times
- Various clean-ups and fixes, including for regressions after porting

## [12.1.0] - 2024-10-15
### Changed
- Implemented proper diffusion coefficient for the ETFC model
- Improved check logic for too low turbulent burning velocities in the TFC model

## [12.0.0] - 2024-10-07
### Changed
- Ported to OpenFOAM 12 - minor changes

## [11.0.0] - 2024-03-26
### Changed
- Ported to OpenFOAM 11 by a complete rewrite of the solver into a collection of classes now to be used with standard OpenFOAM modules, same capabilities maintained (except ETFC and ANN LBV models, porting pending)
- Version numbering scheme changed to indicate target OpenFOAM version

## [0.12.2] - 2022-06-14
### Added
- Artificial neuron network based model of laminar burning velocity
- LES combustion models
- ETFC combustion model
### Changed
- Ported to OF9
- Refactored combustion modeling
- Minor changes, fixes and clean-ups

## [0.8] - 2021-01-14
### Changed
- Initial public release
