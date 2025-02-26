## Quadric Surfaces

This plugin is **not** capable of transforming a generic mesh surface into quadric geometries. Instead, it provides a set of predefined objects for constructing the geometry. If more flexibility is needed, the user can use the [Triangular Meshes](triangular-meshes.md) geometry format. These objects are listed in a submenu named **Quadric** in the **Add** menu, as shown in the following image:

<img src="/images/blenderPluginQuadAdd.png" alt="Quadric Objects Menu" width="600" style="display: block; margin: 0 auto"/>

Note that **only the objects listed in the *Quadric* submenu can be exported as quadric surfaces**. Other types of objects will be ignored by the quadrics exporter. To determine if an object is quadric-compatible, the user can check if the **Body Properties** panel has a properties block named **Quadric Properties**:

<img src="/images/cylQuad.png" alt="Quadric Properties Panel" width="400" style="display: block; margin: 0 auto"/>

Regardless of the quadric object type, the following parameters must be specified to define the object in the geometry:

- **`Material Index`**  
  Sets the material index filling this body. Note that the first defined material corresponds to material index `1`, as index `0` is reserved for void regions.

- **`Module`**  
  If enabled, the geometry package will assume this object is a module. This means none of its surfaces are crossed by other bodies or module surfaces. Therefore, all descendant bodies and modules must be completely included inside the module volume. If used, this option speeds up the particle tracking calculations, but the user must ensure the requirements are fulfilled.

The only quadric objects requiring additional parameters are the **CONE** and the **CONE_SHELL**. For both, the top radius and bottom radius must be specified.
