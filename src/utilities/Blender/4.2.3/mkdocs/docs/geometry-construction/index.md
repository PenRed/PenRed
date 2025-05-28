# Geometry Construction

This section covers the created tools, in addition to those provided by Blender, to construct geometries to simulate on. When constructing a geometry in Blender for export to PenRed, it is important to note that the default Blender units are already considered as centimeters. Therefore, no further transformations are necessary. The specific considerations for each geometry type will be discussed in the following sections. However, those sections does not contain an introduction to Blender usage, as it is already well-documented, and users can find many tutorials on the internet.

Therefore, to use this plugin effectively, we assume that the user has a basic understanding of how to use Blender. The key features to be familiar with include navigating in the 3D viewport, creating objects and adjusting their properties, establishing relationships between objects, editing meshes, and installing addons.

Notice that, regardless the geometry type, **the plugin utilizes Blender's parentage mechanics to determine which objects are included within others**. Therefore, the user must construct the object tree according to the geometry requirements.

The available geometry types to build with the plugin are:

- [Quadric Surfaces](quadric-surfaces.md)
- [Triangular Meshes](triangular-meshes.md)

