## Utilities

This plugin provides some utilities to streamline the geometry construction. 



### Object Positioning

Focusing on relative object posicioning, specially in [Quadric Surfaces](quadric-surfaces.md) based geometries, stacking objects is a common operation needed to construct complex devices. To simplify this process, a set of utilities have been added to object transforms:

<img src="../../images/utilitiesTransform.png" alt="Position Utilities" width="500" style="display: block; margin: 0 auto"/>

This operations are applied along the *Z* axis, regardless the object rotations, and take into account the objects bounding box to determine the final object position. The operations are described following:

- **Put Down**: Move the selected objects under the active object bounding box.

- **Put Top**: Move the selected objects on the active object bounding box.

- **Put Inside Top**: Move the selected objects inside the active object, snapping them to the top face of the active object bounding box.

- **Put Inside Down**: Move the selected objects inside the active object, snapping them to the bottom face of the active object bounding box.

- **Put Inside Centered**: Move the selected objects inside the active object, making their bounding box center coincide with the active object center.

### Revolution Quadrics

On construction...
