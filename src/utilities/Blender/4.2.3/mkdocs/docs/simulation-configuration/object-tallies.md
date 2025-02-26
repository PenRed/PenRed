## Object Tallies

Object tallies are defined in a spatial region within the geometry system. To define this region, these tallies are associated with a Blender object, from which the position and region size are obtained. For mesh-based Blender objects, both the center of the tallied region and its size are derived from the object's bounding box (a box surrounding all vertices of the Blender object). If the object is a non-mesh object (e.g., a light, camera, or empty), the center of the tallied region corresponds to the object's **location**, and the region size is set to the **scale** property of the object. Therefore, object tallies can be defined in any position within the geometry, and it is not necessary to create a material body to define them. Note that multiple tallies of the same type can be defined within the same Blender object, for example, to use different scoring mesh bins.

These tallies can be created and removed from the panel **object > Tally Properties**:

<img src="../../simulation-configuration/images/objectTallyProperties.png" alt="Object Tallies" width="300" style="display: block; margin: 0 auto"/>

Regarding the naming of tallies, the final name is created according to the following pattern:

{Object Name}\_{index}\_{Tally Name}

where:
- **Object Name** is the name of the Blender object.
- **index** corresponds to the tally index within the specific type.
- **Tally Name** is the name specified by the user for this tally.

The following object tallies can be defined within the plugin:

- [Cylindrical Dose Distribution](#cylindrical-dose-distribution)
- [Spatial Dose Distribution](#spatial-dose-distribution)
- [Spherical Dose Distribution](#spherical-dose-distribution)
- [Kerma](#kerma)

---

### Cylindrical Dose Distribution

This tally scores the absorbed dose in a cylindrical mesh oriented along the Z-axis. The cylindrical region is defined to include the entire object. For mesh objects, the entire bounding box is included within the cylinder, with its radius calculated as:

$$
r = \sqrt{\left(\frac{dx}{2}\right)^2 + \left(\frac{dy}{2}\right)^2}
$$

where:
- \( dx \) and \( dy \) are the bounding box dimensions along the X and Y axes, respectively.

The configuration parameters for this tally type are:

- **`Radial Bins`**  
  The number of radial bins in the cylindrical region.
- **`Z Bins`**  
  The number of height bins along the Z-axis.
- **`Angular Bins`**  
  The number of polar bins dividing the cylindrical region.

---
  
### Spatial Dose Distribution

This tally scores the absorbed dose in the voxels of a Cartesian mesh that fills the object's bounding box.

The configuration parameters for this tally type are:

- **`X Bins`**  
  The number of bins along the X-axis.
- **`Y Bins`**  
  The number of bins along the Y-axis.
- **`Z Bins`**  
  The number of bins along the Z-axis.

---

### Spherical Dose Distribution

This tally scores the absorbed dose in a spherical mesh. The sphere's radius is calculated as:

$$
r = \sqrt{\left(\frac{dx}{2}\right)^2 + \left(\frac{dy}{2}\right)^2 + \left(\frac{dz}{2}\right)^2}
$$

where:
- \( dx \), \( dy \), and \( dz \) are the bounding box dimensions along the X, Y, and Z axes, respectively.

The configuration parameters for this tally type are:

- **`Radial Bins`**  
  The number of radial bins in the spherical region.
- **`Polar Bins`**  
  The number of polar bins in the spherical region.
- **`Azimuth Bins`**  
  The number of azimuthal bins in the spherical region.

---

### Kerma

This tally is based on the equivalence of particle fluence and the total photon path length per unit volume. The estimator can be tallied using three types of meshes: voxel (Cartesian-based), cylindrical, and spherical meshes. The parameters required to define each mesh type are equivalent to those described for the [Spatial](#spatial-dose-distribution), [Cylindrical](#cylindrical-dose-distribution), and [Spherical](#spherical-dose-distribution) dose distribution tallies.

In addition, this tally requires a list of \( \mu_{en} \) values for the simulated energy range and for each material in the tallied region. These values are automatically calculated by penRed and saved in files for use in subsequent simulations. The user can provide a prefix to specify where these files should be saved using the **Data Prefix** parameter. For example:

- To save \( \mu_{en} \) files in a folder named `muen`, set the **Data Prefix** value to `muen/` on Unix systems or `muen\` on Windows.
