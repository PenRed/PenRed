# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### To add
- CT based source
- CT detector tally


## [1.2.3] - 2021-04-21

### Added
- Enable PSF translations and rotations
- PSF translation and rotation example (quadrics/5-accelerator-3)
- Utility to convert the IAEA PSF format to PenRed PSF format
- Utility to convert the PenRed PSF format to IAEA PSF format

## [1.2.2b] - 2021-04-08

### Fix
- Bug: Errenoeus calculation of the ECUTR variable when WCC is greater than minimum EABS
- Bug: Erroneous history skipping at dump resume when no specific samplers have been created

### Added
- Soft energy deposition corrections in mutren utility

## [1.2.2] - 2021-03-28

### Fix
- Corrected incorrect access of the IED and IEU variables.

### Added
- Mutren utility based on the mutren code in the original PENELOPE package

### Changed
- Changed IPHF and IPHL variables of the element data base class to store the indexes in a C like format (starting from 0) instead of FORTRAN like (starting from 1)
- The mu_en data for the kerma track length tally must be specified in the input file in eV instead of keV

## [1.2.1] - 2021-03-11

### Added
- File spectrum energy sampler	

## [1.2.0] - 2021-03-08

### Fix
- Main program prints Bremsstrahlung splitting number (IBRSPL) instead of x-ray splitting number (IXRSPL)
- Generic sources configuration error for verbose levels lesser than 2
- Compilation warnings on load balance server utility

### Added
- Variance reduction (VR) modules to implement custom variance reduction techniques
- Splitting VR module
- Russian Roulette VR module
- Radial splitting VR module (experimental and undocumented)

### Changed
- X-Ray splitting converted to VR module
- Increased the maximum number of intervals in the "energy intervals" energetic source

## [1.1.3] - 2020-11-27

### Fix
- DICOM based spatial source bugs
- DICOM dose distribution tally bugs
- DICOM parser bugs

### Added
- Extended load balance capabilities via TCP/HTTP/HTTPS communication

## [1.1.2] - 2020-06-05

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
