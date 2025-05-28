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

This tally scores the absorbed dose in a cylindrical mesh oriented along the Z-axis. The cylindrical region radius depends on the object's bounding box and configuration parameters.

The configuration parameters for this tally type are:

- **`Inside Bounding Box`**  
  Controls how the radius of the cylindrical volume is selected:
    - **Fit inside the bounding box**: Enable the **`Inside Bounding Box`** option to constrain the cylinder within the object's bounding box. The radius is set as:
    
        $$
        r = min\left(dx, dy\right)/2
        $$

        where \( dx \), and \( dy \) are the bounding box dimensions along the X, and Y axes, respectively.
        
    - **Enclose the bounding box**: Disable the **`Inside Bounding Box`** option to ensure the cylinder fully encloses the bounding box. In this case, the cylinder radius is calculated as:
    
        $$
        r = \sqrt{\left(\frac{dx}{2}\right)^2 + \left(\frac{dy}{2}\right)^2}
        $$

        where \( dx \) and \( dy \) are the bounding box dimensions along the X and Y axes, respectively.

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

This tally scores the absorbed dose in a spherical mesh. The sphere's radius depends on the object's bounding box and configuration parameters.

The configuration parameters for this tally type are:

- **`Inside Bounding Box`**  
  Controls how the radius of the spherical volume is selected:
    - **Fit inside the bounding box**: Enable the **`Inside Bounding Box`** option to constrain the sphere within the object's bounding box. The radius is set as:
    
        $$
        r = min\left(dx, dy, dz\right)/2
        $$

        where \( dx \), \( dy \), and \( dz \) are the bounding box dimensions along the X, Y, and Z axes, respectively.
    
    - **Enclose the bounding box**: Disable the **`Inside Bounding Box`** option to ensure the sphere fully encloses the bounding box. The radius is calculated as:
    
        $$
        r = \sqrt{\left(\frac{dx}{2}\right)^2 + \left(\frac{dy}{2}\right)^2 + \left(\frac{dz}{2}\right)^2}
        $$

        where: \( dx \), \( dy \), and \( dz \) are the bounding box dimensions along the X, Y, and Z axes, respectively.    

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

---

### CT

This tally is designed to work with a **CT-like source** and must be attached to a Blender *Empty* object with an enabled CT source. It collects particles reaching a virtual detector positioned behind the CT circle.  

For each projection in the CT circle, a virtual detector (shown in blue in the figure below) is placed in front of the corresponding source position:

  <img src="../../simulation-configuration/images/ctTally.png" alt="CT Tally" width="400" style="display: block; margin: 0 auto"/>

Each detector **only records particles** with ages falling within the projection time (defined in the source). This ensures that a detector exclusively captures particles originating from the source directly in front of it. The result is a **sinogram** containing measurements from all projections.  

#### **Configuration Parameters**  

##### **Energy Parameters**  
Defines the energy range, in eV, of particles to be recorded.  
- **`Minimum`**  
  Lower energy threshold for particle detection.  
- **`Maximum`**  
  Upper energy threshold for particle detection.

##### **Detector Parameters**  
Configures the detector's geometric properties.  
- **`Pixels`**  
  Number of subdivisions (pixels) in the detector.  
- **`Depth`**  
  Radial depth of the detector, in cm.  
- **`Aperture`**  
  Angular width of the detector, in degrees, subtended from the isocenter.  

##### **Detection Options**  
Controls which events are detected.  
- **`Detect Scatter`**  
  If disabled, only **non-scattered particles** are recorded.  
- **`Particle`**  
  Specifies the particle type to detect (other types are ignored).  
