# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### To add
- CT based source
- CT detector tally

## [1.1.2] - YYYY-MM-DD

### Added
- DICOM based spatial source
- DICOM dose distribution tally

### Fix
- DICOM read bugs on selected data representation

## [1.1.1] - 2020-05-18

### Added
- Kerma track length estimator tally
- Load balance description at documentation
- Kerma track length estimator description at documentation
- This changelog file

### Changed
- Spherical dose tally now accepts azimuthal and polar bins
- Cylindrical dose tally now accepts azimuthal bins
- Updated documentation for spherical and cylindrical tallies

### Fix
- Fixed an error that duplicate particles at phase space file sources when splitting is set to 1
- Fixed an error affecting external sources on voxelized geometries

## [1.1.0] - 2020-05-13

### Added
- Load balance system for multi-threading and MPI executions
- Missage errors for parse library
- Main program prints now a text missage instead of a number identifier for configuration parsing errors

### Changed
- Tallies, sources and the main program have been changed to count histories as integers instead of doubles
- Sources incorporates an optional Load Balance system based on [RUPER-LB](https://github.com/PenRed/RUPER-LB)

### Fixed
- Minor bugs at parse library
- Explicit multi-threading compatibility check on MPI initialization at main program

### Removed
- Compiled folder

## [1.0.0] - 2020-03-02

### Added
- First stable version of PenRed code system
- Documentation at "doc" folder
- Complete examples with sample results at *examples* folder
- README.md file with basic description
- Docker container to execute PenRed simulations
