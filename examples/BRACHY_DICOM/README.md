
# BRACHYTHERAPY DICOM EXAMPLE

This example consists of a complete Brachytherapy treatment. For this purpose, DICOM files of a computational patient are used. These ones mimick a clinical scenario of interstitial High Dose Rate breast case. The example is subdivided in two stages. First, the geometry of  MBDCA Ir-192 source is build and a Phase Space File (PSF) at the exit of the encapsulation is stored. Then, the PSF is used as a source for each seed location which is read automatically from the DICOM plan.

## Phase Space File generation
The simultation of the `PSF` folder must be executed first to obtain the Phase Space File at the exit of the encapsulated seed. This one contains the following files:

* aapm-Ir-192-mbdca_psf.geo: geometry file of the generic MBDCA Ir-192 source. The geometry information is available at https://www.uv.es/radiofisica/
* Ir.spc: Ir-192 spectrum
* config.in: configuration file 
* material files (.mat) 


## Treatment
The `treatment` folder includes the needed files to execute a Breast Test Case treatment. DICOM files and all the information about this case can be found at http://rpc.mdanderson.org/rpc/BrachySeeds/BrachySeeds/index.html  Test Case Data → Breast test case


This folder contains the files/folder listed below:

* breastCaseHDR.zip: DICOM files of the computational patient
* mutrenData folder: with the mu-tren information for each material. Have been obtained using the *mutren* utility included in the penRed package.
* config.in: configuration file
* material files (*mat)

Notice that, to be able to execute the treatment, the path where the `psf-merged.dat` file of the first simultaion was created, must be added to the config.in file.
