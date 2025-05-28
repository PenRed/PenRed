## Detector Tallies

Like [Object tallies](object-tallies.md), detector tallies are associated with Blender objects, from which the position and region size (if needed) are obtained according to the corresponding object's bounding box. Additionally, multiple tallies of the same type can be defined within the same Blender object, for example, to use different scoring mesh bins. However, this type of tally must necessarily be associated with a material body that has the **detector** flag enabled in the [Body properties](body-parameters.md). Therefore, the tally is enabled for all bodies sharing the same detector index.

For tally types that affect a limited spatial region (instead of the entire detector volume), note that the **defined region includes only the body bounding box of the object where the tally is defined**.

Like object-based tallies, these tallies can be created and removed from the panel **object > Tally Properties**:

<img src="../../simulation-configuration/images/detectorTallyProperties.png" alt="Detector Tallies" width="300" style="display: block; margin: 0 auto"/>

Regarding the naming of tallies, the final name is created according to the following pattern:

{Object Name}\_{index}\_{Tally Name}

here:
- **Object Name** is the name of the Blender object.
- **index** corresponds to the tally index within the specific type.
- **Tally Name** is the name specified by the user for this tally.

The following object tallies can be defined within the plugin:

- [Imact Detector](#impact-detector)
- [Spatial Distribution](#spatial-distribution)
- [PSF](#psf)
- [Angular Detector](#angular-detector)

---

### Impact Detector

This tally is capable of measuring different magnitudes: fluence, particle energy, energy deposition, and particle age spectra. Each magnitude is described below:

- **`Fluence Spectrum`**  
  Measures the spectral fluence integrated over the detector volume, expressed in cm/eV. The output file provides the fluence for each particle type as well as the total fluence within the specified detector.

- **`Energy Spectrum`**  
  Reports the energy spectrum of particles entering the specified detector volume. Particles generated inside the detector, such as secondary particles, are not included in this tally. The units are given in 1/(eV·history). The output file contains the energy spectrum for each particle type.

- **`Age Spectrum`**  
  Records the age distribution of particles impacting the detector, with units in 1/(seconds·history). The output file contains the probability distribution across different time intervals for all simulated particles. Notice the particle age will be calculated only if it is enabled in [Source Parameters](particle-sources.md).

- **`Energy Deposition Spectrum`**  
  Reports the energy deposition spectrum measured in the specified detector. The units are given in 1/(eV·history). The output file contains the probability distribution across energy intervals, considering all simulated particles.

Each magnitude can be enabled or disabled using the corresponding toggle in the tally panel. Additionally, a linear or logarithmic scale can be used for each case. To configure the recording histograms, the following parameters must be set:

- **`Energy`**  
  Energy parameters specify the energy range and bins for all energy-based magnitudes. Additionally, energy range limits also the particles to be recorded for age spectrum.
    - **`Min`**  
      Minimum energy to be recorded, in eV.
    - **`Max`**  
      Maximum energy to be recorded, in eV.
    - **`Bins`**  
      Number of energy bins.

- **`Age`**  
  Age parameters specify the time range and bins for the **Age Spectrum**. These are only required if the **Age Spectrum** is enabled.
    - **`Min`**  
      Minimum particle age to be recorded, in seconds.
    - **`Max`**  
      Maximum particle age to be recorded, in seconds.
    - **`Bins`**  
      Number of age bins.

---
    
### Spatial Distribution

This tally records the energy spectrum of particles reaching a specified detector within a spatial 3D grid. Like other impact detectors, particles created inside the detector will not be recorded. The available configuration parameters are as follows:

- **`Energy`**  
  Energy parameters specify the energy range and bins to score the spectrum.
    - **`Min`**  
      Minimum energy to be recorded, in eV.
    - **`Max`**  
      Maximum energy to be recorded, in eV.
    - **`Bins`**  
      Number of energy bins.

- **`Spatial`**  
  Spatial parameters specify the number of bins, in each axis, for the 3D spatial grid.
    - **`X Bins`**  
      Number of bins in the X-axis.
    - **`Y Bins`**  
      Number of bins in the Y-axis.
    - **`Z Bins`**  
      Number of bins in the Z-axis.

- **`Particle Type`**  
  Specifies the particle type to be recorded. Other types will be ignored by the tally.

- **`Coordinates`**  
  If enabled, the coordinates (X, Y, and Z) will be printed alongside the results, in cm.

- **`Print Bins`**  
  If enabled, the bin indexes in the X, Y, and Z axes will be printed alongside the results.

---
  
### PSF

This tally generates a particle phase space file (PSF) storing all particles impacting the specified detector. Additionally, the impacting particles can be filtered by energy range and particle type. The configuration parameters are described below:

- **`Min`**  
  Minimum particle energy to be recorded, in eV.

- **`Max`**  
  Maximum particle energy to be recorded, in eV.

- **`Gammas`**  
  Enable/disable gamma recording.

- **`Electrons`**  
  Enable/disable electron recording.

- **`Positrons`**  
  Enable/disable positron recording.
  
---

### Angular Detector

This tally reports the angular energy spectrum in a specified detector. In this tally, the polar and azimuthal angles of impacting particles are defined relative to the Z-axis (0, 0, 1) and the particle's direction. Therefore, these angles are **not related to the detector's orientation**. The recorded particles are filtered by energy and angular intervals. The energy spectra of particles are tallied in units of 1/(eV·sr·particle). Below are the parameters used to configure this tally:

- **`Energy`**  
  Energy parameters specify the energy range and bins to score the spectrum.
    - **`Min`**  
      Minimum particle energy to be recorded, in eV.
    - **`Max`**  
      Maximum particle energy to be recorded, in eV.
    - **`Bins`**  
      Number of energy bins.

- **`Polar Angle`**  
  Defines the polar angle interval to be recorded.
    - **`Min`**  
      Minimum polar angle, in degrees.
    - **`Max`**  
      Maximum polar angle, in degrees.

- **`Azimuthal Angle`**  
    Defines the azimuthal angle interval to be recorded.
    - **`Min`**  
      Minimum azimuthal angle, in degrees.
    - **`Max`**  
      Maximum azimuthal angle, in degrees.
