# CELES Changelog

<!--
  format inspired by [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
  
  please use the following types where appropriate
    Added       for new features.
    Changed     for changes in existing functionality.
    Deprecated  for soon-to-be removed features.
    Removed     for now removed features.
    Fixed       for any bug fixes.
-->


## [Testing]

## [Unreleased]
### Changed
- subclass celes classes from `matlab.System`
- implement `setProperties` methods for name-value style initialization
- implement `validatePropertiesImpl` methods for validation
- limit use of `Dependent` properties to avoid redundant calculations
- implement `setupImpl` where appropriate for one-time calculations to improve performance
- compute and set maximal particle distance inside particles class
- prefer implicit expansion over `bsxfun` (requires MATLAB >= R2016b)
- fancier, faster plotting functions
- use convex hull in compute_maximal_particle_distance
- provide robust, fallback method for compute_maximal_particle_distance

### Fixed
- don't try to compute initial power for plane waves
- fix bug in scattered field plots
- fix bug in computeTotalFieldPWP

## [2.1] - 2017-10-25
### Changed
- define only trigonometric versions of legendre and spherical functions
- precalculate coefficients in spherical functions @tkfryett
- case-insensitive flags and types

### Removed
- `disperse` flag (radii must be specified)

### Fixed
- performance issue with polydisperse samples @fragkrag

## [2.0] - 2017-09-22
### Added
- support for polydisperse sphere size and refractive index @fragkrag

### Changed
- patch MATLAB's built-in GMRES to monitor progress
- unified log messages from iterative solvers
- gather arrays from GPU to save memory when using GMRES
- precalculate maximal particle distance

### Removed
- monitor flag for iterative solvers

### Fixed
- fix wrong sign for downward propagation of PWP Gaussian beam


## 1.0 - 2017-02-24

[Testing]: https://github.com/disordered-photonics/celes/compare/HEAD...develop
[Unreleased]: https://github.com/disordered-photonics/celes/compare/v2.1...HEAD
[2.1]: https://github.com/disordered-photonics/celes/compare/v2.0...v2.1
[2.0]: https://github.com/disordered-photonics/celes/compare/v1.0...v2.0
