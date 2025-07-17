## Materials

Material properties, such as cross-sections, density, and composition, are defined in material files. These files can be generated at runtime by specifying the composition in the configuration, as discussed later. Materials are identified by a **unique name and index**, which are used to assign materials to each body.

A material can be added or removed in the **World > Material Properties** panel, as shown in the following image:

<img src="../../images/materials.png" alt="Material Panels" width="300" style="display: block; margin: 0 auto"/>

In the image above, two materials **water** and **Bone** have been defined, corresponding to material indexes `1` and `2`, respectively. The material index is assigned according to the order in which materials appear in the panel, and the index `0` is reserved for **void** regions.

For a single simulation, up to **200 different materials** can be defined. A new material can be added using the *Add Material* button, and the last one can be removed using the *Remove Material* button.

Each material must define simulation-specific parameters, which are described in the following sections. 

### Cutoffs

Cutoff parameters specify when the simulation of a specific particle type ends. For each particle type, the following cutoff parameters must be set:

- **`Type`**  
  The cutoff type to be applied for this particle type in the current material. The options are:
  - **`Energy`**: An energy limit, in eV, is specified. A particle with energy below this limit will be locally absorbed.
  - **`Range`**: A particle with a range below the specified value will be absorbed locally. The meaning of this range depends on the particle type:
    - **`Gammas`**: The range is defined as the particle's mean free path (inverse attenuation coefficient), in cm.
    - **`Electrons and Positrons`**: The range is interpreted as the Continuous Slowing Down Approximation (CSDA), in cm.
        
### Composition

Material composition can be defined using either of two methods: selecting a predefined composition from a database or manually specifying the composition. The desired method can be chosen by selecting either **Database** or **Composition** from the Definition dropdown menu.

<img src="../../simulation-configuration/images/materialCompositionType.png" alt="Material Composition Definition" width="350" style="display: block; margin: 0 auto"/>

#### Manual Definition

To define a material manually, specify:

- The atomic number (*Z*) and corresponding weight fraction for each element in the compound.
- The material density in \( \text{g/cm}^3 \)

Users can add elements using the *Add Element* button or remove the last added element with the *Remove Element* button. Note that weight fractions do not require normalization.

#### Material Databases

Alternatively, users can select a predefined material from either:

- The internal PENELOPE material database
- An external database. Currently the following are included:
    + ICRP Adult Female (*ICRP_AF*)
    + ICRP Adult Male (*ICRP_AM*)

<img src="../../simulation-configuration/images/materialDBSelection.png" alt="Material DB Selection" width="350" style="display: block; margin: 0 auto"/>

After selecting the database, the available materials appear in a list below the selection menu. A search function allows users to filter materials by name.

<img src="../../simulation-configuration/images/materialsDBSearch.png" alt="Material DB Selection" width="350" style="display: block; margin: 0 auto"/>

### Advanced Parameters

Electron and positron transport use the Class II scheme to group interactions and speed up the simulation, which can be tunned setting its configuration parameters. If the advanced parameter settings are disabled, these values will be set automatically by **penRed**, ensuring a good compromise between accuracy and simulation speed in most cases. However, these parameters can be specified by the user if needed:

- **`C1`**  
  The average angular deflection, \( C1 \approx 1 - \langle \cos \theta \rangle \), produced by multiple elastic scattering along a path length equal to the mean free path between consecutive hard elastic events. The maximum allowed value is \( 0.2 \), but a value of \( 0.05 \) is usually adequate.

- **`C2`**  
  The maximum average fractional energy loss between consecutive hard elastic events. Like \( C1 \), the maximum allowed value is \( 0.2 \), and a value of about \( 0.05 \) is usually adequate.

- **`WCC`**  
  The cutoff energy loss, in eV, for hard inelastic collisions.

- **`WCR`**  
  The cutoff energy loss, in eV, for hard bremsstrahlung emission.
