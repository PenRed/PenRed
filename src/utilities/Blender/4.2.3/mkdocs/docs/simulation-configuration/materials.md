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

To define a material, its composition must be specified via the atomic number (*Z*) and the corresponding weight fraction for each element present in the compound. Additionally, the material density must be set in \( \text{g/cm}^3 \). The user can include as many elements as needed to define the compound using the *Add Element* button or remove the last one with the *Remove Element* button.

The weight fractions of elements composing the material do not need to be normalized.

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
