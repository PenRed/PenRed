# PenRed
Parallel Engine for Radiation Energy Deposition


## Introduction

PenRed package provides a parallel engine for radiation transport Monte-Carlo simulations. The information about installation, usage and description of source files can be found at the provided documentation. This one is distributed within the PenRed package at

> doc/PenRed_doc.pdf

In addition, users can find a descriptive [academic article](https://arxiv.org/abs/2003.00796) at ArXiv which has been sent to be published.



## Installation

To install PenRed, the user has two options. Compile the code by itself or install it using a package manager.

### Code compilation

To install PenRed directly by the code compilation, download PenRed sources from this GitHub repository,

```
git clone https://github.com/PenRed/PenRed.git
```

The code must be compiled in the *src* folder. Thus, move into this folder. To simplify the installation, PenRed includes a CMake file and a bash script (*compile.sh*) to compile the code automatically. In this script, you can enable/disable the following optional features,

1. DICOMs: If it is enabled, PenRed capabilities to read and simulate DICOM images will be active. This option requires the library dicom toolkit (dcmtk).
2. Multi-threading: This option enables multi-threading capabilities. PenRed implements multi-threading via the standard thread library specified in the C++11 standard. Thus, it is not required any extra library to enable this option.
3. MPI: This option enables MPI simulations. It requires a library with an implementation of the MPI standard, such as openmpi or mpich.

Notice that all previous dependencies are optional. The corresponding libraries can be found at most linux package repositories. For example, to compile PenRed with DICOM support in Fedora, you can use the *dnf* command to install the dependencies,

```
sudo dnf install dcmtk dcmtk-devel
```

After that simple step, you can launch the compilation using the provided script,
```
bash compile.sh
```

or doing it yourself,

```
cd /path/to/PenRed/repository/src
mkdir build
cd build
ccmake ../
make install
```

With *ccmake* you can configure the optional PenRed features with a more friendly interface. But, of course, you can use directly *cmake* defining the appropriate flags like,

```
cmake -DWITH_DICOM="ON" -DWITH_MULTI_THREADING="ON" -DWITH_MPI="OFF" -DDEVELOPMENT_WARNINGS="OFF" ../
```

Once the code has been compiled, the user can found the executable of our provided main program ready to simulate at,

> src/compiled/mains/pen_main

To execute the program, the user needs a configuration file and, probably, the required data base files, such as material and geometry files. Their path should be specified at the configuration, remaining the configuration file path as the only program argument. So, to execute the program, the user must use 

```
./pen_main path/to/configuration/file
```

Otherwise, if the MPI capabilities have been enabled at the compilation, the code should be executed as any MPI program, for example,

```
mpirun -np Nprocesses ./pen_main path/to/configuration/file
```

where *Nprocesses* specify the number of MPI processes to use. Of course, the user can use any other options of the *mpirun* command, such us specify the hosts where execute the code via the *hostfile* option.

### Use package managers

Actually, PenRed can be installed via a *rpm* package using the *yum* and *dnf* package managers. A *deb* package is under development and will be provided soon. 

All these packages can be found at the repository 

> https://github.com/PenRed/packages

You can download a single package instead of the whole repository, using the curl comand. For example, to 
download some of the available *rpm* files, use,

```
curl https://raw.githubusercontent.com/PenRed/packages/master/rpm/PenRed-XXXX.x86_64.rpm --output PenRed.rpm
```

where XXXX depends on the desired file to download. Finally, install the package using

```
sudo dnf install PenRed.rpm
```

or

```
sudo yum install PenRed.rpm
```

Once the installation has been finished, PenRed can be executed directly using,

```
penred path/to/configuration/file
```

Also, some utilities will be installed too. They are explained at the package repository.
