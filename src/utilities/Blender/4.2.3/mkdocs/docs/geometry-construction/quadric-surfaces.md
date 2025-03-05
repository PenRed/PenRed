## Quadric Surfaces

This plugin is **not** capable of transforming a generic mesh surface into quadric geometries. Instead, it provides a set of predefined objects for constructing the geometry. If more flexibility is needed, the user can use the [Triangular Meshes](triangular-meshes.md) geometry format. These objects are listed in a submenu named **Quadric** in the **Add** menu, as shown in the following image:

<img src="../../images/blenderPluginQuadAdd.png" alt="Quadric Objects Menu" width="400" style="display: block; margin: 0 auto"/>

Note that **only the objects listed in the *Quadric* submenu can be exported as quadric surfaces**. Other types of objects will be ignored by the quadrics exporter. To determine if an object is quadric-compatible, the user can check if the **Body Properties** panel has a properties block named **Quadric Properties**:

<img src="../../images/cylQuad.png" alt="Quadric Properties Panel" width="400" style="display: block; margin: 0 auto"/>

Regardless of the quadric object type, the following parameters must be specified to define the object in the geometry:

- **`Material Index`**  
  Sets the material index filling this body. Note that the first defined material corresponds to material index `1`, as index `0` is reserved for void regions.

- **`Module`**  
  If enabled, the geometry package will assume this object is a module. This means none of its surfaces are crossed by other bodies or module surfaces. Therefore, all descendant bodies and modules must be completely included inside the module volume. If used, this option speeds up the particle tracking calculations, but the user must ensure the requirements are fulfilled.

The only quadric objects requiring additional parameters are **CONE** bodies. For this type, the top radius and bottom radius must be specified.

### Cutting Planes

Optionally, any quadric object can be cut using planes to modify its shape. These *cuts* are defined as boolean operations, which can be added or removed from the object using the **Quadric Properties** panel. This ensures the structure is created as intended or removed entirely.

When a cutting plane is created, both a boolean modifier and the corresponding plane are generated. You can select these planes by clicking their names in the **Cutting Planes** panel:

<img src="../../images/selectCuttingPlane.png" alt="Select Cutting Plane" width="300" style="display: block; margin: 0 auto"/>

Once selected, the cutting plane becomes visible in the 3D viewport and can be moved or rotated to achieve the desired shape. Note that cutting planes are only visible when selected:

<img src="../../images/cuttingPlaneExample.png" alt="Cutting Plane Example" width="500" style="display: block; margin: 0 auto"/>

To remove a specific cutting plane, click the trash icon next to its name in the panel.
