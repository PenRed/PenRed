## Body Parameters

This section covers the parameters used to define body properties in the configuration file. Note that body parameters related to geometry definition are covered in the [Geometry Construction](../geometry-construction/index.md) section. These parameters can be found in the **Object > Body Properties**:

<img src="../../simulation-configuration/images/bodyParameters.png" alt="Body Properties" width="300" style="display: block; margin: 0 auto"/>

The behavior of each body in the geometry can be configured with the following parameters:

- **`Material Index`**  
  Sets the material index filling this body. The material index corresponds to the order in which materials have been defined in the world properties (see section [Materials](materials.md)). Note that the first defined material corresponds to material index `1`, as index `0` is reserved for void regions.

- **`Detector`**  
  If enabled, this material body will be considered a detector. To specify the detector the body belongs to, the **`index`** parameter must be set accordingly. Some tallies can only be used on detector bodies, which require this flag to be enabled (see [Detector Tallies](detector-tallies.md)).
  
- **`Class II Advanced Parameters`**  
  The **`Limit Distance`** parameter within this block allows the user to specify a maximum distance between hard interactions for particles using the Class II transport scheme (electrons and positrons). This parameter should be tuned for very thin bodies (compared to the particle range), as the scheme requires particles to interact enough times within the crossed body. As a reference, ensuring at least 10 interactions should be sufficient. In most cases, this parameter can be disabled.
