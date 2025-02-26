## Triangular Meshes

Unlike the [Quadric Surfaces](quadric-surfaces.md) exporter, the exporter for triangular mesh surfaces can convert any mesh constructed with Blender to the PenRed format. Therefore, there is no need to use the predefined quadric objects. Additionally, Blender can be used to import meshes constructed with other programs and export them to the PenRed format.

The only parameter required for each body is the **`Material Index`**, which can be found in the **Body Properties** panel within the **Object Properties**. 

<img src="/images/materialIndexMesh.png" alt="Object transforms" width="300" style="display: block; margin: 0 auto"/>

Additionally, if vertex groups are defined for any object, the plugin will also export them according to the expected PenRed format inside the geometry file. These vertex groups can then be transformed during the simulation configuration at runtime (See the PenRed documentation for more details).

<img src="/images/vertexGroupsMesh.png" alt="Object transforms" width="300" style="display: block; margin: 0 auto"/>


For this type of geometry, mesh surfaces **cannot be crossed by other objects or the same object's mesh surfaces**. This restriction is enforced to achieve better performance in particle tracking. To determine which object is inside another, the plugin uses Blender's parent system. Finally, note that triangular mesh geometries require a *world* object that includes all other defined objects.
