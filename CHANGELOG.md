# CELES Changelog

<!--
  format inspired by [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
  
  please use the following types where appropriate
    Added     	for new features.
    Changed   	for changes in existing functionality.
    Deprecated  for soon-to-be removed features.
    Removed   	for now removed features.
    Fixed     	for any bug fixes.
-->


## [Unreleased]
### Changed
- precalculate coefficients in spherical_functions_angular @tkfryett
- case-insensitive flags and types
- remove distinction between mono- and poly-disperse (radii must be specified)


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

[Unreleased]: https://github.com/disordered-photonics/celes/compare/v2.0...HEAD
[2.0]: https://github.com/disordered-photonics/celes/compare/v1.0...v2.0
