## Body Parameters

This section outlines the parameters used to define body properties in the configuration file. Note that body parameters related to geometry definition are covered in the [Geometry Construction](../geometry-construction/index.md) section. These parameters can be accessed under **Object > Body Properties**:

<img src="../../simulation-configuration/images/bodyParameters.png" alt="Body Properties" width="300" style="display: block; margin: 0 auto"/>

The behavior of each body in the geometry can be configured with the following parameters:

- **`Material Index`**  
  Specifies the material index assigned to the body. The material index corresponds to the order in which materials are defined in the world properties (see the [Materials](materials.md) section). Note that the first defined material corresponds to material index `1`, as index `0` is reserved for void regions.

- **`Detector`**  
  When enabled, the body is treated as a detector part. To associate the body with a specific detector, the **`index`** parameter must be set accordingly. Certain tallies are exclusive to detector bodies, requiring this flag to be enabled (see [Detector Tallies](detector-tallies.md)).
  
---
  
### Advanced Parameters

Advanced parameters allow users to define body-specific Class II behavior and apply variance reduction techniques:

- **`Class II Limit Distance`**  
  Defines the maximum distance between hard interactions for particles using the Class II transport scheme (applicable to electrons and positrons). This parameter is useful for very thin bodies (relative to the particle range), as the scheme requires sufficient interactions within the body. As a guideline, ensuring at least 10 interactions is typically sufficient. In most cases, this parameter can be left disabled.

- **`X-Ray Splitting`**  
  When enabled, generated X-rays are split into multiple particles according to the specified `Factor`. The weight of each resulting particle is reduced by a factor of $1/\text{Factor}$.

- **`Bremsstrahlung Splitting`**  
  When enabled, generated bremsstrahlung photons are split into multiple particles according to the specified `Factor`. The weight of each resulting particle is reduced by a factor of $1/\text{Factor}$.

- **`Interaction Forcing (IF)`**  
  Allows users to increase the probability of specific interactions for a given particle type. Multiple *IF* modules can be defined, with each module affecting a single interaction type for a specific particle. The available parameters for defining interaction forcing modules are:

    - **`Particle`**  
    Specifies the particle type for which interactions will be forced.

    - **`Interaction`**  
    Specifies the interaction type to be forced.

    - **`Factor Type`**  
    Determines how the forcing factor is applied. The available options are:
        - **`Multiply`**  
        The selected interaction becomes `Factor` times more probable.
        - **`Average`**  
        The `Factor` represents the average number of interactions for a particle with the maximum simulation energy until it comes to rest (for electrons and positrons) or traverses a mean free path (for gammas).

    - **`Weight Window`**  
    Restricts the application of interaction forcing to particles whose weights fall within the specified range.
