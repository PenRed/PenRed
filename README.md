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
4. Load balance: Enables load balancing system between threads and MPI processes. This option requires, as least, multi-threading capabilities, because uses threads to handle MPI communications.

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
cmake -DWITH_DICOM="ON" -DWITH_MULTI_THREADING="ON" -DWITH_MPI="OFF" -DWITH_LB="OFF" -DDEVELOPMENT_WARNINGS="OFF" ../
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

## Docker containers

Another method to run PenRed is using a docker container. We provide several docker images with different characteristics. The first one **vigial/penred-compiled** contains the whole github package and all the required dependencies to compile the source code, including the dcmtk library for DICOM suport. To download the image and execute the simulation, use the following command in the simulation folder,

```
sudo docker run -v $PWD:/home/penred:z vigial/penred-compiled path/to/configuration/file
```

For size optimized containers, we provide two Alpine based images with only the compiled executable and the runtime libraries. The first one has not DICOM support enabled and has a compressed size lesser than 5 MB. Can be downloaded using the instruction,

```
sudo docker push vigial/penred-alpine
```

The second one has DICOM support enabled and its compressed size is lesser than 20 MB. Can be downloaded using

```
sudo docker pull vigial/penred-alpine-dicom
```

To execute a simulation with that images, use the corresponding instruction on the simulation folder,

```
sudo docker run -v $PWD:/home:z vigial/penred-alpine path/to/configuration/file
sudo docker run -v $PWD:/home:z vigial/penred-alpine-dicom path/to/configuration/file
```
