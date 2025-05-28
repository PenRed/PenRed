
<p align="center">
  <a href="https://github.com/PenRed/PenRed">
    <img src="logo/PenRed_color.png" width="400" alt="PenRed logo">
  </a>
</p>

[![CMake](https://github.com/PenRed/PenRed/actions/workflows/cmake.yml/badge.svg)](https://github.com/PenRed/PenRed/actions/workflows/cmake.yml) 

Parallel Engine for Radiation Energy Deposition


## Introduction

PenRed is a powerful parallel engine for radiation transport Monte Carlo simulations. Comprehensive information about installation, usage, and component details can be found in the provided documentation (*PenRed_user_manual.pdf*), located in the *doc* folder.

For further reading, users can refer to a descriptive [academic article](https://arxiv.org/abs/2003.00796) available on ArXiv, which has been accepted and published in the Computer Physics Communications journal. The article can be accessed at [https://doi.org/10.1016/j.cpc.2021.108065](https://doi.org/10.1016/j.cpc.2021.108065). If you use PenRed in research that leads to publications, please cite this article.

## pyPenred

A Python module named [pyPenred](https://pypi.org/project/pyPenred/) has been developed to run PenRed simulations through **Python**. The easiest way to install the pyPenred package is via pip:

```bash
pip install pyPenred
```

On some systems, you may also need to install the *pyYAML* package manually:

```bash
pip install pyYAML
```

The complete [API documentation](https://penred.github.io/PenRed/pyPenred) can be found online and included in this package at *src/bindings/python/pyPenred/docs/*

### Compiling pyPenred

To compile and install the pyPenred package manually (to enable native optimizations or include custom-developed modules), use the appropriate compilation script from the src folder:

Linux/macOS:

```bash
./installPyPenred.sh
```

Windows:

```bash
installPyPenred.bat
```

### Usage

A brief explanation of pyPenred features and usage can be found in the [pyPenred readme](https://pypi.org/project/pyPenred/), which is also included in this repository at:

```
src/bindings/python/pyPenred/README.md
```

For a detailed description, see the *Python Wrapper (pyPenred)* section in the PenRed documentation.

## Blender Integration

A [Blender](https://www.blender.org/) plugin has been developed to integrate geometry construction, simulation configuration, and execution within the Blender environment, providing a user-friendly graphical interface. The plugin can be installed via the zip file:

```
src/utilities/Blender/4.2.3/penred.zip
```

The plugin [documentation is available online](https://penred.github.io/PenRed/Blender), and both the plugin and documentation source code can be found in the same folder as the zipped plugin.

## Code Compilation

PenRed can be compiled with any compiler supporting the C++14 standard and [CMake](https://cmake.org/download/) version 3.11 or later. The steps to compile the code are described following.

### Download the Source Code

Download the source code either manually from the repository page or using git:

```bash
git clone https://github.com/PenRed/PenRed.git
```

### Compilation

Navigate to the src folder and use the appropriate compilation script for your operating system:

  + **Linux/Unix**: Run the bash script *compile.sh* with

    ```bash
    bash compile.sh
    ```

  + **Windows**: Run the batch script *compile.bat*

The scripts will create:

  1. A build folder for temporary compilation files

  2. A compilation folder containing all generated executables

The main PenRed executable (used to run simulation examples) is located at:  

```
src/compiled/mains/pen_main
```

### Compilation Options

Edit the compilation scripts to enable/disable multiple features. Some of them are listed below:

  + **DICOMs**: Enables DICOM image reading/simulation capabilities (requires DCMTK library). If the DCMTK library is not installed, it will be downloaded and compiled automatically. This step can take several minutes.
  + **Multi-threading**: Enables multi-threading via C++ standard thread library (no additional dependencies).
  + **MPI**: Enables MPI simulations (requires MPI implementation library like OpenMPI or MPICH).
  + **Load Balance**: Enables load balancing between threads/MPI processes (requires multi-threading)

Additional compilation options for specific tools and bindings can be configured using CMake utilities (see documentation).  

### Basic usage

To execute the main program, the user needs to specify the configuration file path:

```bash
./pen_main path/to/configuration/file
```

The configuration file should specify paths to required database files (materials, geometries, etc.). See the *Framework Usage* section in the documentation for details.

The examples folder contains ready-to-run configuration files with corresponding material and geometry files (described in the *Examples* documentation section). To run an example:

1. Copy the executable to the example folder
2. Run:

  ```bash
  ./pen\_main config.in
  ```
  
  or simply double click the executable. This will assume a configuration file named *config.in* exists in the directory.

For MPI-enabled builds, execute as any MPI program:  

```
mpirun -np Nprocesses ./pen\_main path/to/configuration/file
```

where *Nprocesses* specifies the number of MPI processes. Standard mpirun options (like hostfile) are supported.

## Docker containers

The *containers* directory provides basic Dockerfiles to generate containers with an entry point ready to run PenRed simulations.
