
# Internal rt

This example consists on an internal RT treatment adapted from the Gate examples to be able to compare the results (CITA). Indeed, the results of this example executed with the Gate toolkit are provided in the *Gate* folder. The remaining folder structure is described following:

 * *data*: Contains all the required input data for the simulation, including a tarfile with the Mhd images produced by Gate, that can be converted automatically ASCII format via the scripts provided in *scripts* folder. Notice that the *WalkerAliasing* folder contains a preprocessed data of the walkers algorithm. This is not required to run the simulation, as PenRed will compute that file. However, once computed, it can be stored automatically to be used in simulations with the same distributed source and save computation time. Therefore, the data in this folder is provided to speed-up the computation. In addition, the folders contain:
 
    + *Dicom*: A tarfile with the patient CT and SPECT DICOM images, which are directly processed by PenRed
    
    + *Materials*: A tarfile with all the required materials. To be able to compare with the Gate code, according to their metodology, each material has been defined multiple times as independent materials with different densities instead of handling automatically voxels with the same material but different densities.
    
    + *PennucSpectrum*: Contains the Pennuc required data to produce the source spectrum. Alternatively, a configuration file with a approximated spectrum is provided.
    
 * *infiles*: Contains two configuration files, one to use the Pennuc source and another one to use an aproximated spectrum.
 
 * *Output* and *outputs*: Auxiliary folders to be filled with the simulation results.
 
 * *scripts* Postprocessing scripts to plot and compare the results of both codes, PenRed and Gate.

# How to run

1. Compile PenRed with DICOM support (with the DCMTK library). The compilation process is described in the documentation. Once compiled, copy the executable to this folder.

2. Extract/decompress the dicom, material and Walker algorithm tar and xz files:

        cd data/Dicom/; tar -xvf dicoms.tar.xz; cd ../..
        cd data/Materials/; tar -xvf materials.tar.xz; cd ../..
        cd data/Mhd/; tar -xvf Mhd.tar.xz; cd ../..
        xz -dk data/WalkerAliasing/ActivityInitialization.dat.xz 

3. Extract also the gate tar in *Gate* folder to compare the results

        cd Gate/; tar -xvf gate.tar.xz; cd ..      

4. Execute *pen_main* executable with one of the configuration files of the folder "infiles". The default number of threads to use in the configuration file is 1, but can be changed modifing the line

        simulation/threads 1

    Then, run pen_main
      
        ./pen_main infiles/Internal-rt-Y90LineSpectrum.in &> outputs/Internal-rt-Y90LineSpectrum.log &

    or
       
        ./pen_main infiles/Internal-rt-Y90PennucSpectrum.in &> outputs/Internal-rt-Y90PennucSpectrum.log & 

    The tally and output files generated will be saved in the output folder. 
      
5. Use the bash script *Check-Gate_PenRed.sh* to verify the PenRed results against Gate. To execute the comparison with verbose output use:

        ./scripts/Activity-Dose_DistributionChecker.sh -c all -q all -r Gate -t output.gcc9.3.1-pennuc -R Gate -T PenRed
     instead, for a less verbose output use:

        ./scripts/Activity-Dose_DistributionChecker.sh -c all -q all -r Gate -t output.gcc9.3.1-pennuc -R Gate -T PenRed -s
        
     For further help details, enter 
     
        ./scripts/Activity-Dose_DistributionChecker.sh [or optional -c help]
