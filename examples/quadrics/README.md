
# Quadric geomety based examples

This folder contain several examples which uses quadric geometries. Each folder contains 
an example with the corresponding configuration file `config.in`, geometry and required material
files. However, configuration files must be modified to be executed with PenRed. The user must
specify the number of histories to simulate, the number of threads to use and, finally, the time
interval between dumps. To do that, replace the following texts with the corresponding values,

* __NHISTS__ : Replace by the number of histories, for example, `1.0e9`
* __NTHREADS__ : Replace by the number of simulation threads, for example `2`
* __DUMPTIME__ : Replace by the time between dumps (in seconds), for example, `300`

## Results

To folder `results` contain a results sample for all examples.
