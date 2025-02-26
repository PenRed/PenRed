
## World Tallies

World tallies are not associated with any object or spatial region of the geometry system. This tallies can be created and removed from the panel **world > Tallies Properties**:

<img src="/simulation-configuration/images/worldTalliesProperties.png" alt="World Tallies" width="300" style="display: block; margin: 0 auto"/>

The available tallies of this type are listed below and described in this section:

- [Emerging Particles](#emerging-particles)
- [Tracks](#tracks)

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
