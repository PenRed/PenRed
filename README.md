
<p align="center">
  <a href="https://github.com/PenRed/PenRed">
    <img src="logo/PenRed_color.png" width="400" alt="PenRed logo">
  </a>
</p>

[![CMake](https://github.com/PenRed/PenRed/actions/workflows/cmake.yml/badge.svg)](https://github.com/PenRed/PenRed/actions/workflows/cmake.yml) 

Parallel Engine for Radiation Energy Deposition


## Introduction

The PenRed package provides a powerful parallel engine for radiation transport Monte Carlo simulations. Comprehensive information about installation, usage, and details of each component can be found in the provided documentation (*PenRed_user_manual.pdf*), located in the *doc* folder.

For further reading, users can refer to a descriptive [academic article](https://arxiv.org/abs/2003.00796) available on ArXiv, which has been accepted and published in the Computer Physics Communications journal. The article can be accessed at https://doi.org/10.1016/j.cpc.2021.108065. If you use PenRed in research that leads to publications, please cite this article.


## Installation

PenRed can be compiled with any compiler supporting the C++14 standard and [CMake](https://cmake.org/download/) version 3.4 or greater. Once a suitable C++ compiler is available, to compile PenRed, first download the source code from this repository. This can be done manually from the repository page or using *git* with the following command line instruction:

```
git clone https://github.com/PenRed/PenRed.git
```

With the code downloaded, navigate to the *src* folder to compile it. This folder includes a CMake file and two scripts to compile it automatically, depending on the operating system:

  + **Linux/Unix**: Run the bash script *compile.sh* with

    > bash compile.sh

  + **Windows**: Run the batch script *compile.bat*

The compilation scripts will create a *build* folder where all the temporary files required by the compilation will be stored. Finally, a *compilation* folder will be created in the *src* directory with all generated executables. Among all the executables, the PenRed main program, the one used to run the simulation examples, is located at:

```
src/compiled/mains/pen_main
```

Editing the corresponding script, allows the user to enable/disable the following main optional features:

  + **DICOMs**: If enabled, PenRed's capabilities to read and simulate DICOM images will be active. This option requires the DICOM toolkit (dcmtk) library.
  + **Multi-threading**: This option enables multi-threading capabilities. PenRed implements multi-threading via the standard thread library specified in the C++ standard. Thus, no extra library is required to enable this option.
  + **MPI**: This option enables MPI simulations. It requires a library with an implementation of the MPI standard, such as *OpenMPI* or *MPICH*.
  + **Load Balance**: Enables load balancing between threads and MPI processes. This option requires at least multi-threading capabilities, as it uses threads to handle MPI communications.

These options, except for multi-threading capabilities which are included in the C++ standard, require some external libraries. These libraries are typically available in most Linux and Unix package repositories, or they must be downloaded separately for Windows. Additionally, PenRed offers several compilation options for specific tools and bindings. These options can be configured using CMake utilities, as explained in the documentation.

### Basic usage

To execute the main program, the user needs a configuration file and, likely, the required database files, such as material and geometry files. Their paths should be specified in the configuration file, leaving the configuration file path as the only program argument. All details regarding the simulation configuration can be found in the documentation section *Framework Usage*. Additionally, the *examples* folder contains several configuration file examples with the corresponding material and geometry files ready to be executed. These ones are described in the documentation *Examples* section. To execute the program, the user must use the command line instruction:

```
./pen_main path/to/configuration/file
```

Therefore, to run the provided examples, copy the main program executable in the corresponding example folder and run the instruction:

```
./pen\_main config.in
```

If MPI capabilities have been enabled during compilation, the code should be executed as any MPI program. For example:

```
mpirun -np Nprocesses ./pen\_main path/to/configuration/file
```

where *Nprocesses* specify the number of MPI processes to use. Of course, the user can use any other options of the *mpirun* command, such us specify the hosts where execute the code via the *hostfile* option.

## Docker containers

The *containers* directory provides basic Dockerfiles to generate containers with an entry point ready to run PenRed simulations.
