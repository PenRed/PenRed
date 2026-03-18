# Dynamic Movement Example: The "Demilune" Scan

This example demonstrates the time-dependent geometry capabilities in PenRed by simulating a simple scanning pattern. A small radiation source moves relative to a fixed detector, tracing a path that combines a 180-degree rotation with a linear translation, forming a shape reminiscent of a **demilune** or "D".

The primary goal is to illustrate different implementation approaches for achieving the same relative motion, highlighting the flexibility of the system and the concept of transformation propagation in a geometry tree.

## Understanding Transformations in PenRed

Before diving into the examples, it's important to understand how PenRed applies time-dependent transformations to geometry objects:

1. **Rotation interpolation**: When an object has rotation keyframes, PenRed interpolates between them using **SLERP** (Spherical Linear Interpolation). This ensures smooth and physically meaningful rotations regardless of the keyframe spacing.

2. **Translation interpolation**: The position of an object between keyframes can be interpolated in two ways:
   - **Spline interpolation**: Provides smooth, curved trajectories.
   - **Linear interpolation**: Available as an option for simpler, straight-line motion.

This distinction becomes critical when designing complex motions, as the choice of interpolation method affects the final path, especially when combining rotation and translation.

## The Scenario

*   **Components**: Two simple geometric bodies: a small cylindrical **source** and a larger cuboid **detector**.
*   **Desired Motion**: The source should trace a path relative to the center of the detector, consisting of two phases:
    1.  A 180-degree rotation around the detector's center.
    2.  A linear translation along the X-axis, starting from the end position of the rotation back to the beginning.
*   **Resulting Path**: This combined motion creates a "D" shape (or demilune).

We have implemented this relative motion using three distinct methods, as detailed below.

## The Three Approaches

Each approach is contained in its own subdirectory (`Case1`, `Case2`, `Case3`) along with two gnuplot scripts for plotting, allowing for direct comparison of the setup and results. Additionally, the corresponding Blender files can be found within the *examples/triangualr_surfaces/blender/animations* folder.

### 1. `Case1` - Moving the Detector

*   **Concept**: The source remains fixed in the global coordinate system. The detector is moved (rotated and translated) in such a way that, from the source's perspective, the relative motion is the desired demilune path.
*   **Implementation**: Time-dependent transformations are applied directly to the detector's geometry object.
*   **Key Insight**: Motion is relative. This approach can be simpler if the source's position is complex to define or if multiple sources need to move relative to a single object.

### 2. `Case2` - Using an Auxiliary Void Object (Recommended)

*   **Concept**: An "auxiliary" object is placed in the geometry tree as a parent to the source. All transformations (rotation, translation) are applied to this auxiliary object. The source inherits this motion because transformations propagate to all descendants.
*   **Implementation**:
    1.  Create a dummy `void` object (`sourceHead`) positioned at the detector's center.
    2.  Place the source as a child of this void object, offset to its starting position.
    3.  Apply the desired sequence of rotations and translations to the auxiliary void object.
*   **Key Insight**: This is often the most **convenient and intuitive method**, especially when the desired motion is naturally described relative to a point in space (like the detector's center). It avoids complex, nested transformations on the moving object itself. Notice that a local transform can also be applied to the source, combining parent and child transformations.

### 3. `Case3` - Moving the Source Directly

*   **Concept**: The source is moved directly in the global coordinate system, following the absolute coordinates of the demilune path.
*   **Implementation**: Time-dependent transformations are applied directly to the source geometry object, defining its translation and rotation at each time step.
*   **Key Insight**: This is the most straightforward approach conceptually. However, as implemented here, it introduces differences. This is because the motion between defined time points is interpolated. For a path combining rotation and translation, a direct spline interpolation of the source's position may not perfectly replicate the combined rotational and linear motion achieved by the other methods, leading to variations in the trajectory. To illustrate this effect, this case uses the same number of time frames as the previous cases with the default Blender spline configuration.

## Results and Key Takeaways

*   **Approaches 1 and 2** produce the expected results. This confirms that applying motion to an object or to its parent (with appropriate offsets) are equivalent strategies, thanks to the propagation of transformations.
*   **Approach 3** produces a slightly different path due to the difference between interpolating rotation with the SLERP method and interpolating spatial positions with splines. This can be mitigated introducing more keyframes and adjusting the spline definition.
    > **When defining complex motions involving rotation, it is often more accurate and simpler to use an auxiliary object placed at the center of rotation.**

This example serves as a practical guide for users to choose the right method for their own dynamic simulations in PenRed. The auxiliary object method (`Case2`) is generally the recommended starting point for movements that involve rotation around a point.
