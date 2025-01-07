# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.12.1] -

### Added
- CMake instructions to download and build all necessari dependencies for DICOM support
- Option to force the download and compilation of the DICOM package instead of using the system-installed version
- Dirent header file for Windows
- Set default simulation configuration file name to 'config.in' when no file is provided

### Changed
- The embedded database generation has been refactored to save memory usage during compilation

### Fixed
- Bug: DICOM geometry configuration return success even if enclosure parameters are not provided
- Bug: Forced termination of simulation is not functioning correctly in asynchronous simulation mode
- Compilation: Python binding compilation fails when DCMTK is enabled and fetched
- The PSF sampler triggers a segmentation fault when it encounters a particle with energy exceeding the defined maximum

## [1.12.0b] - 2024-11-06

### Changed
- Updated the documentation for the variance reduction techniques section
- Logs that were previously flushed only at the end of the simulation are now flushed immediately when printed
- VR logs now display body aliases instead of body indexes


## [1.12.0] - 2024-10-18

### Added
- Global log system
- Better support to perform penred simulations from python
- Detection singles/pulses tally
- Utility to convert penred data format files to YAML. These YAML files can be used by the penred python wrapper
- Script to automatically compile and install the pyPenred python wrapper

### Changed
- The dockerfile to create containers has been updated
- The code has been refactored to limit stack usage to under 1MB, preventing issues on Windows systems
- The execution of register functions is now forced by code and not by linker directives (like --whole-archive)

### Fixed
- Bug: Only tallies which implement the 'begin_sim' function gets particle stacks pointers
- Bug: Erroneous composition log print during material creation

## [1.11.0b] - 2024-07-11

### Added

- ASCII input format for voxelized geometries

### Fixed

- Materials creation via the configuration file fails for materials with a single element

## [1.11.0] - 2024-07-01

### Added

- Auxiliary *measurement* and *results* classes to the penred::measurements namespace to automate the tallying of magnitudes with arbitrary dimensions, as well as the generation and reading of formatted files

- Auxiliary *aliasing* class to the penred::sampling namespace to automate the sampling of distributions of arbitrary dimension using Walker's aliasing. Distributions can be obtained from the *measurement* and *results* classes

- *DETECTION_SPATIAL_DISTRIB* tally to generate 1D, 2D, 3D, or 4D spatial and energy distributions of impinging particles in a specific detector

- *1D_MEASURE*, *2D_MEASURE*, and *3D_MEASURE* spatial samplers to sample particle positions according to 1D, 2D, and 3D distributions, respectively. These distributions can be generated with the *DETECTION_SPATIAL_DISTRIB*

- Configuration example for utilizing the new tallies and sampler. It is located in *examples/quadrics/9-measure-source*

- x-ray based utilities to streamline simulations of this kind of devices

- Profile functions and utility to create n-dimensional profiles from *measurement* and *results* classes.

- Databases of compositions, along with utilities to list them

- Methods to create materials from compositions databases

- Finish simulation function to stop simulations in the *simulation* library functions

- *simulator* class to wrapp the main program functioanlities. For now, MPI and load balance support is not implemented in this class.

- New tally to register energy spectrums of emerging particles in spherical bins

- Additional constraints for DICOM segmentation

- First version of a Python wrapper using the *simulator* class and pybind11

### Changed

- The compilation of the TCP library is now optional, allowing the code to be compiled without the ASIO library

- Erroneous material ID assignation on automatic generated material section by *pen_context* configuration.

- Integrated data bases uses now split string literals to allow to fit specific compiler limitations on string literals length.

- Particle stacks are now allocated dynamically to avoid stackoverflows in some systems.

- Removed number of bins limit in cylindrical and spherical dose distribution tallies

### Fix

- In some cases, *COMBO* geometry skips low priority geometries when a higher geometry contains void bodies.

- Erroneous verbose levels in geometry prints.

- Kerma track length tally does not write the automatically generated absortption coefficients until the end of the simulation.

## [1.10.0] - 2024-04-15

### Added

- By default, the material database is embedded in the required executables, allowing to construct materials automatically at runtime. This compilation option can be disabled to reduce the executable size
- New library to streamline the material creation by code
- Helper classes to handle configuration definitions, reading and documentation generation in the internal data library
- Utility named *penContextInfo* which provides an automatically generated documentation of PENELOPE context configuration
- Library for automatically handling context simulations with automatic particle instantiation

### Changed

- C++ required standard has been upgrated to C++14
- Old and depecrated eabs_e+, eabs_e- and eabs_gamma options to specify absorption energies in configuration files have been removed
- Contexts require now to implement two new pure virtual functions, *config* and *configWithGeo*
- PENELOPE context configuration has been moved into *config* and *configWithGeo*, streamlining main programs construction
- The PENELOPE context does not require an element database as a constructor argument
- "TRACK" tally is now less verbose and can be used to plot simulated tracks

## [1.9.4] - 2024-02-23

### Added
- x-ray utilities
- DICOM class reads and process CT specific tags
- New geometry type "filters"
- New PSF sampler which stores the whole PSF in memory
- Specific samplers can now share data after configuration, as was done for tallies
- Materials expose their composition and mean Z
- Simulation library to streamline the creation and handling of simulations
- Particles can now update their state directly from the internal stack
- Stream support to parse internal data sections/configurations
- Triangular surface geoemtries can now define vertex groups and transform them via the configuration
- Triangular surface geoemtries can now be loaded from a internal string instead of a external file
- Blender plugin now exports vertex groups
- Geometries vector of *combo geometry* instances is now accessible
- Main stack of particle instances is now accessible
- Capability to request simulation status information during execution
- Defined virtual methods to access internal geometries for combined geometry types (like *COMBO*)
- Seeds rotation in BRACHY source can be enabled/disable
- Errors on context initialization are now printed via the standard output by default

### Changed
- PSF tally is now able to select which particles must be recorded
- Mathematical template classes have been moved to a independent header file
- A specfic Macro to declare specific samplers have been created (DECLARE_SPECIFIC_SAMPLER)
- Replaced context virtual method "getIF" by "setForcing"
- BRACHY source accepts now DICOM geometries or a combined geometry type with a DICOM inside it

### Fix
- Bug: Blender addon crash on quadric geometry export
- Bug: Combo geometry overwerrites the global bodies material number with detector number
- Bug: Maximum energy on specific samplers is overwritten when the parameter is not used in the specific configuration

## [1.9.3b] - 2023-09-15

### Added
- Geometry bodies name is now accessible through the geometry wrapper class
- Default global material energy absorption energies can be assigned for each particle
- Configuration parameter to set the default object octree partition to all bodies for mesh based geoemtries

### Changed
- Energy body deposition tally now prints body names
- Increased maximum number of materials to 200

### Fix
- Bug: Triangular mesh geometries can't handle void worlds correctly
- Bug: Triangular mesh geometries can't handle particles sampled in void regions within the geometry

## [1.9.3] - 2023-08-06

### Added
- Database ID or source filename to base material class
- Number of threads as configure function parameter for specific sources
- Library to calculate muen data

### Changed
- PSF based sources no longer need to specify the number of partitions in configuration
- Kerma track length tally will create the mutren data and file if it does not exist
- Material energy deposition tally no longer need to specify the number of materials

### Fix
- Bug: Kerma track length tally does not check if the material is used nor created
- Bug: PSS tally does not implement the shared configuration function
## [1.9.2] - 2023-06-11

### Added
- Geo to wavefront mesh utility for visualization
- Capability to chain dumps

### Fixed
- Bug: Tallies composed by subtallies (DICOMkerma and PSS) do not dump the subtally information
- Energy body deposition compilation warnings

## [1.9.1] - 2023-06-04

### Added
- Mesh test capability to find bodies intersections
- Options in the Blender plugin to export only active or non hided objects

### Changed
- Energy deposition in body tally: No longer needs to specify the number of bodies
- Energy deposition in body tally: Begins body number from 0 instead of 1, to be consistent with geometry packages

### Fixed
- Bug: DSMAX parameter in mesh geometries only accepts body names with 5 or less characters
- Bug: Some mesh geometries lose triangles when the object is subdivided in regions

## [1.9.0] - 2023-05-23

### Added
- Now a particle range can be specified instead of absorption energy
- Default parameters to material parameters (Eabs, C1, C2, etc.)
- A tally postconfiguration step to enable sharing information between threads. This avoids repeating huge configuration calculus and allow their parallelization
- Capability to finish all workers on load balance server

### Changed
- The default maximum number of materials have been increased to 100
- The configuration parameter names to specify absorption energies have been changed (the old ones have not been removed yet for retrocompatibility)
- Spatial dose distribution tally now calculates the voxel mass once and more efficiently

### Fixed
- Bug: Mesh geometry file read fails sometimes
- Bug: In the load balance server, the speed of new workers with no reports is considered as 0 on the balance step

## [1.8.1b] - 2023-05-12

### Added
- Load balance server instruction to force worker finish

### Fix
- DICOM geometry config not converting contour name to lower case
- DICOM contour ranges assign
- Results dump was flagged as incomplete simulation

## [1.8.1] - 2023-04-17

### Added
- Access to particle stacks as constants objects from tallies
- PSS generic tally to registrer energy deposition from primary, scatter and multi-scatter contributions

## [1.8.0] - 2023-03-19

### Added
- DICOM: Capability to specify intensity and density ranges for specific segmentation in each contour
- Cylindrical shell spatial sampling
- Geometry "COMBO" type: Allow to combine different geometries in a single one regardless of its type

### Fix
- Geometry view library compilation with DICOM support

## [1.7.1b] - 2023-02-15

### Fix
- Bug: Symbols not exported on viewer C binding DLL for Windows

## [1.7.1] - 2023-02-14

### Added
- Default compilation options for clang compiler

### Changed
- Material Range utility is now compiled by default when utilities compilation is enabled

### Fix
- Compilation error on clang compilers
- Several warning messages on clang compiler

## [1.7.0] - 2023-02-11

### Added
- Triangular mesh geometry simulation support
- Blender plugin to construct Quadric and Triangular Mesh based geometries

## [1.6.0] - 2022-07-14

### Added
- Image exporter library for tally results. Actually, only MHD and gnuplot formats are implemented
- Support for CT source to use generic samplers instead of a PSF
- Library with C bindings to render geometries in 2D and 3D

### Changed
- Voxel based geometries now use a box enclosure instead of a spherical enclosure
- spatial, direction, energy and time sampler pointers are not passed as arguments of the specific sampler configuration anymore
- Spatial dose now do not print depth dose results by default, but can be enabled.
- DICOM contour detection voxels algorithm

### Fix
- Bug: Errors on tracking transport for voxelized geometries when the particle source is outside the enclosure
- Bug: Specific sources do not update their generic samplers pointers (spatial, direction, energy and time) before the configuration call

## [1.5.2] - 2022-06-14

### Added
- Enclosures to mesh based geometries to consider backscattering effects at boundaries.

### Changed
- Mesh based geometries IBODY 0 is reserved for the enclosure
- Body identifiers of mesh based geometries are now equal to the material index, i.e. (IBODY = MAT) instead of (IBODY = MAT + 1)

## [1.5.1] - 2022-05-03

### Added
- Support to redirect output files to user defined directories
- Support to disable context report creation
- Complete DICOM example consisting of a internal-rt treatment. Can be found in *examples/internal_rt*.

## [1.5.0] - 2022-03-08

### Added
- Capability to simulate brachytherapy treatments via a specific source named "BRACHY"
- Tally to score kerma on DICOM geometries and obtain Dose Volume Cummulative Histogram (DVH) for each provided contour
- Capability to specify source material on PENNUC specific source
- DICOM parser now creates also individual masks for each contour

### Fix
- Bug: DICOM parser erroneous seed information reading
- Bug: DICOM structures, such as contours and seeds, were displaced half voxel in each axis

## [1.4.0] - 2021-12-04

### Added
- Capability to add dump files
- Capability to specify a different number of threads for each MPI rank
- Cabability to dump final results of each thread
- Automatic number of threads selection
- PENNUC specific source

### Fix
- Bug: Erroneous history count when dumping a resumed simulation
- Bug: When a specific sampler produce a secondary particle (dhist=0) in a void region and if the particle reaches a non void region, its kinetic energy is scored as deposited energy. However, this bug does not affect the particle tracking and other components, only energy counters of the first reached material.

### Changed
- Tally funcion "beginHist" has been replaced by the function "sampledPart". This one is called when a source samples a particle.

## [1.3.0b] - 2021-10-15

### Fix
- Bug: Contour voxel plane assign was inconsistent when the contour Z coordinate is located on the edge between two voxels 

## [1.3.0] - 2021-10-06

### Added
- CT based source
- CT detector tally
- Geometry information is now accessible for all sampler types via the virtual method *updateGeometry*
- Image sampler can, optionally, ajust the source position automatically if it is used with DICOM based geometries
- Function to access DICOM seed data
- The main program will check the returned value of the geometry configuration function, in addition to the configure status variable

### Fix
- Bug: DICOM geometry does not update the configure status variable during the configuration function call
- Bug: Kerma track length tally fails on configuration when multiple materials are used

### Changed
- Documentation has been split and expanded

## [1.2.4b] - 2021-07-02

### Changed
- Optimisation of the quadric gometry library, specially for complex geometries

## [1.2.4] - 2021-06-21

### Added
- Negative interpretation of the IF factor as the PENELOPE's penmain program does, i.e. the IF is calculated to produce an average number of interactions equal to abs(factor)

### Changed
- Changes in the CMake files to be albe to compile the program using the MSVS

### Removed
- The creation of error files during the compilation has been removed, as it not require more changes.

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
