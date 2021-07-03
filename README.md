# PenRed
Parallel Engine for Radiation Energy Deposition


## Introduction

PenRed package provides a parallel engine for radiation transport Monte-Carlo simulations. The information about installation, usage and description of source files can be found at the provided documentation. This one is distributed within the PenRed package at

> doc/PenRed_doc.pdf

In addition, users can find a descriptive [academic article](https://arxiv.org/abs/2003.00796) at ArXiv which has been accepted and published in the Computer Physics Communications journal (https://doi.org/10.1016/j.cpc.2021.108065). If PenRed is used for research conducting to publications, please, cite the published article.



## Installation

PenRed compilation has been tested in several Linux distributions with gcc 5 to 10 versions, clang and icc, and in Windows 10 environment using the MSVS 2019. However, the support on windows is not complete tested yet. For example, DICOM geometries have not been tested, but multithreading and MPI capabilities are tested and expected to work.

### Linux code compilation

To install PenRed via the code compilation, download PenRed sources from this GitHub repository,

```
git clone https://github.com/PenRed/PenRed.git
```

The code must be compiled in the *src* folder. To simplify the installation, PenRed includes a CMake file and a bash script (*compile.sh*) to compile the code automatically. In this script, you can enable/disable the following optional features,

1. DICOMs: If it is enabled, PenRed capabilities to read and simulate DICOM images will be active. This option requires the library dicom toolkit (dcmtk).
2. Multi-threading: This option enables multi-threading capabilities. PenRed implements multi-threading via the standard thread library specified in the C++11 standard. Thus, it is not required any extra library to enable this option.
3. MPI: This option enables MPI parallelism. It requires a library with an implementation of the MPI standard, such as openmpi or mpich.
4. Load balance: Enables load balancing system between threads and MPI processes. This option requires, at least, multi-threading capabilities, because the MPI balance procedure uses threads to handle MPI communications.

It is possible also to enable many other compilation options, which are described in the documentation. Notice that all previous dependencies are optional, thus PenRed can be executed with no external libraries besides the standard C++ library. Moreover, the dependencies libraries can be found at most linux package repositories. For example, to compile PenRed with DICOM support in Fedora, you can use the *dnf* command to install the *dcmtk* library,

```
sudo dnf install dcmtk dcmtk-devel
```

So, to compile PenRed, you can use the provided script,

```
bash compile.sh
```

or do it yourself using the ccmake tool,

```
cd /path/to/PenRed/repository/src
mkdir build
cd build
ccmake ../
make install
```

With *ccmake* you can configure the optional PenRed features with a more friendly interface. But, of course, you can use directly *cmake* defining the appropriate flags like,

```
cmake -DWITH_DICOM="ON" -DWITH_MULTI_THREADING="ON" -DWITH_MPI="OFF" -DWITH_LB="OFF" -DDEVELOPMENT_WARNINGS="OFF" ../
```

Notice that the *cmake* version 3 is the minimum required to perform the compilation.

### Windows code compilation

To compile the PenRed code using MSVS on windows, first, select "Clone a repository" from the Visual Studio start window.

![newRepo](./img/MSVS/NewCloneRepo.png)

Then, set the PenRed repository url and push on the clone button.

![clone](./img/MSVS/ClonePenRed.png)

The download will start automatically. Once the Cmake configuration ends, to avoid compiling with a debug profile, add a new configuration with the configuration manager.

![newRepo](./img/MSVS/manageConfigurations.png)

In the configuration panel to the left, click on the button with the green *plus* sign to add a new configuration and select the release depending on your system. In the following image, we selected a release for 64-bit system.

![newRepo](./img/MSVS/addConfiguration.png)

Once you have selected the new configuration, and push on "Save and generate CMake cache to load variables" to be able to change the CMake variables for this configuration and compile it.

![newRepo](./img/MSVS/selectConfiguration.png)

Now, we can change all the Cmake configuration variables to enable or disable MPI support and other features.

![newRepo](./img/MSVS/setCompileConfig.png)

Finally, build and install PenRed,

![newRepo](./img/MSVS/buildAndInstall.png)

If the compilation finishes successfully, a new folder named *compiled* will be created inside the *src* folder containing the PenRed's main program and the executables for all other enabled utilities.

![newRepo](./img/MSVS/compiled.png)

### Usage

Once the code has been compiled, the user can found the executable of our provided main program ready to simulate at,

> src/compiled/mains/pen_main

To execute the program, the user needs a configuration file and, probably, some required data base files, such as material and geometry files. To introduce how PenRed should be used, the package includes several examples ready to be executed in the *examples* directory. These examples include all the materials and geometry files as well as the corresponding configuration file. Notice that, in some examples, the user must specify the number of histories to simulate and the dump time interval in seconds. In addition, all the available paramterers to be used in the configuration file are summarised in the provided documentation.

To execute the program, only the path to the configuration file is required, so the user can use the following instruction in a Unix environment,

```
./pen_main path/to/configuration/file
```

or with the inverted slashes (\\) if is executed in Windows. On the other hand, if the MPI capabilities have been enabled at the compilation, the code should be executed as any MPI program, for example,

```
mpirun -np Nprocesses ./pen_main path/to/configuration/file
```

where *Nprocesses* specifies the number of MPI processes to use. Of course, the user may use any other options of the *mpirun* command, such as specifying the hosts where the code will be executed via the *hostfile* option.

## Docker containers

As the code grows, several compilation options are included in the *CMakeLists* files. Therefore, providing a prebuilded container for each possible parameter combination is no longer feasible.

Nevertheless, to obtain a containerised version of PenRed, we provide three Dockerfiles in the folder *containers* to generate them. The first one, named *complete* includes all the packages to compile and run PenRed in the Fedora environment. For size optimised versions, we provide also two alpine based Dockerfiles, named *Alpine* and *Alpine-dicom*, which include only the necessary runtime libraries and the PenRed executable for executions without and with DICOM capabilities. However, notice that *Alpine* uses the *musl* library instead of the *glib*, which causes the multithreading PenRed capabilities to don't work properly. 

Notice that each Dockerfile includes a set of predefined compilation options that should be modified by the user according to their needs. For example, if the container will be executed on different architectures, the user should disable the native architecture optimisation specifying the flag

```
-DWITH_NATIVE="OFF"
```

A complete list of all available compilation options can be found in the documentation.
