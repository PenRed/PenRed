
## World Tallies

World tallies are not associated with any object or spatial region of the geometry system. This tallies can be created and removed from the panel **world > Tallies Properties**:

<img src="../../simulation-configuration/images/worldTalliesProperties.png" alt="World Tallies" width="300" style="display: block; margin: 0 auto"/>

The available tallies of this type are listed below and described in this section:

- [DICOM Kerma](#dicom-kerma)
- [Emerging Particles](#emerging-particles)
- [Tracks](#tracks)

---

### DICOM Kerma

This tally has been implemented to obtain the Kerma Spatial Distribution and Dose Volume Cumulative Histogram (DVH) for each contoured structure of a DICOM image. Therefore, this tally requires a DICOM based geometry to be used. To obtain the absorbed dose in Gy for each volume, this tally calculates the kerma, using the kerma track length estimator. The configuration tally parameters are described below.

Within the **Energy** box, two parameters must be specified, **Min** and **Max**, which define the energy range where the kerma is tallied. Any gamma particle whose energy falls outside that range will not contribute to the measurement. If its energy is larger, once it decreases due to interactions and falls within the range, the subsequent paths will be accounted for.

In addition, this tally requires a list of \(\mu_{en}\) values for the simulated energy range and for each material in the DICOM region. These values are automatically calculated by penRed and saved in files for use in subsequent simulations. The user can provide a prefix to specify where these files should be saved using the **Data Prefix** parameter. For example:

- To save \( \mu_{en} \) files in a folder named `muen`, set the **Data Prefix** value to `muen/` on Unix systems or `muen\` on Windows.

Finally, the **DVH** box defines the required parameters to convert Monte Carlo kerma values to real dose and tally the DVH curves.

- The **MC to Dose Factor** is used to convert the resulting simulated kerma in *eV / g history* to dose in *Gy*. To obtain this value for brachytherapy sources, check the *DICOM kerma tally* section of the penRed documentation, specifically the **MCkerma2Dose** parameter.

- The **Max Dose** parameter defines the absorbed dose limit, in Gy, for the DVH graph.

- The **Bins** parameter specifies the number of dose bins used in the DVH representation.

---

### Emerging Particles

This tally measures the energy distribution of particles that leave the geometry within a specific energy range, expressed in particles/(eV·history). The tally considers two cases:

1. **Down-bound emerging particles**: Particles whose Z-axis direction is less than or equal to zero.
2. **Up-bound emerging particles**: Particles whose Z-axis direction is greater than zero.

Additionally, the tally measures the number of particles that leave the geometry according to their direction, dividing the spectrum into polar and azimuthal intervals in spherical coordinates. In this case, measurements are expressed in particles/(sr·history).

The configuration parameters for this tally are as follows:

- **`Energy`**  
  Energy parameters specify the energy range and bins to define the spectrum.
    - **`Min`**  
      Minimum energy to be recorded, in eV.
    - **`Max`**  
      Maximum energy to be recorded, in eV.
    - **`Bins`**  
      Number of energy bins.

- **`Spatial`**  
  Defines the spherical mesh to measure the emerging particles according to their direction.
    - **`Polar bins`**  
      Number of polar bins.
    - **`Azimuthal bins`**  
      Number of azimuthal bins.

---
      
### Tracks

The **Tracks** tally stores the simulated particle tracks in separate files by particle type. Note that a single history can produce long tracks, especially for electrons, and the generated files can consume a large amount of disk space if the number of tracked histories is large. It is recommended not to track more than 50 histories, depending on the energy, simulated particles, and geometry. Only one tally of this type can be defined. The tally parameters are:

- **`Enabled`**  
  Enable or disable the tracks tally
- **`Histories`**  
  Specify the number of histories to be tallied per thread
  
Take into account that only particles reaching the geometry system will be tracked. That is, if a particle is sampled in a void region and never reaches a non-void body, it will not appear in the corresponding track file.

In the generated files, particle tracks are separated by two blank lines, and each line contains the complete state of a particle at the specified point.

