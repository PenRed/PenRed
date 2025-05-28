## Tallies

Tallies are the components responsible for extracting information from the simulation. PenRed provides many tallies to obtain different magnitudes. These tallies can be classified into three groups within the Blender plugin:

- [World Tallies](world-tallies.md)
- [Object Tallies](object-tallies.md)
- [Detector Tallies](detector-tallies.md)

All tallies have a **Name** parameter to identify them. Depending on the type, some prefixes will be added to the user specified name.

In addition to the tallies defined by the user, all generated configurations includes the tallies to score the energy deposition in each material and body in the system.

Some tallies score magnitudes in a spacial ragion. Those regions are displayed in blue within the 3D viewport:

<img src="../../simulation-configuration/images/talliesWireframe.png" alt="Tallies Wireframe" width="400" style="display: block; margin: 0 auto"/>
