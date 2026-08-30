# DICOM

The DICOM geometry module in penRed converts DICOM medical image series into voxelized geometry for radiation transport simulations. Utilizing this geometry type requires pyPenred to be installed in Blender's Python environment. Check the installation details in the [pyPenred install](pyPenred.md) section.

## Segmentation

First, you must specify how materials and densities are assigned to each voxel. This process can be performed using different approaches, which can also be combined.

### Default values

First, a default density (in $\text{g/cm}^3$) and material index are assigned to every voxel. These default values will be overwritten if other segmentation methods are defined.

<img src="../images/dicomDefaultMatDens.png" alt="Default material and density parameters" width="300" style="display: block; margin: 0 auto"/>
dicomDefaultMatDens.png

### Intensity ranges

Multiple intensity ranges can be defined to assign a density and material index based on the voxel intensity value. For each range, a name, density value, material index, and lower/upper intensity bounds must be provided, as shown below:

<img src="../images/dicomIntensityRange.png" alt="Intensity ranges" width="300" style="display: block; margin: 0 auto"/>

Every voxel whose value falls within the interval $[\text{lower}, \text{upper})$ will be assigned the material and density specified for that range. The range name is merely for identification and does not influence the segmentation logic.

### Calibration

A calibration polynomial can be provided to convert voxel Hounsfield Units (HU) to density ($\text{g/cm}^3$) in CT-based images. It consists of an array defining a polynomial calibration curve. If calibration is not specified for CT images, raw data will be used, and densities must be assigned using alternative techniques described in this section. The transform applied to voxel intensities (after transformation by the slope and intercept parameters in the DICOM file) is:

$$
\rho = \sum_{i=0}^n a_ix^i
$$

where $\rho$ is the final density, $n$ is the polynomial grade and $x$ is the voxel intensity value.

### Density ranges

This option is analogous to *Intensity ranges*, but uses density values instead of raw voxel intensity. Therefore, this method requires a *calibration* step to calculate densities before ranges are applied. Each range requires a name, assigned material index, and lower/upper density bounds.

<img src="../images/dicomDensRange.png" alt="Density ranges" width="300" style="display: block; margin: 0 auto"/>

### Contours

Finally, DICOM contours stored in RTSTRUCT files can be used for segmentation. The specified contour name must match the name stored in the DICOM file exactly. Within each contour, the previous segmentation methods (default material and density, intensity ranges, and density ranges) can be applied locally. To manage overlapping or nested contours, each contour must be assigned a priority; contours with lower priority values will be overwritten by contours with higher priority values.

<img src="../images/dicomContours.png" alt="Contours" width="300" style="display: block; margin: 0 auto"/>

Note that specifying default material and density values is optional for individual contours. This allows you to selectively overwrite only the materials and/or densities of voxels that fall within specified ranges inside the contour.

### Method Priority

Because multiple methods can assign material indices and densities to voxels, a method hierarchy determines precedence. Preferential methods overwrite lower-priority ones according to the following graph:

<img src="../images/dicomMatDensAssign.jpg" alt="Method priority" width="150" style="display: block; margin: 0 auto"/>

Note that this same hierarchy is applied locally when processing a contour's internal default values and ranges.

## Enclosure

An enclosure volume with homogeneous material and density surrounds the voxelized DICOM geometry to account for border effects, such as backscatter.

<img src="../images/voxelEnclosure.png" alt="Enclosure" width="200" style="display: block; margin: 0 auto"/>

To define the enclosure, specify the margin (in cm) and enclosure material:

<img src="../images/DICOMenclosure.png" alt="Enclosure parameters" width="300" style="display: block; margin: 0 auto"/>

Note that enclosure density is automatically assigned according to the material's nominal density and cannot be configured manually within the DICOM parameters.

## Extract data

Enabling the *Print ASCII files* option generates a set of text files containing exported DICOM data. The exported files provide detailed information about voxel materials, densities, contours, and seeds:

- **dicomASCII.rep**: Contains the voxelized geometry data, including the material ID and density assigned to each voxel.
- **dicomContourMask.dat**: Contains a mask where each voxel is assigned the index of the contour it belongs to. Voxels not within any contour are assigned a value of -1.
- **dicomContours.dat**: Contains the 3D coordinates of the points defining each contour, with all penRed-specific transformations applied.
- **roi_contour.dat**: One file is created for each DICOM contour (Region of Interest). Each file contains a binary mask, assigning a value of 1 to voxels inside the specific contour and 0 to voxels outside of it. Notice that this masks don’t take into account other contours.

## Loading

Once all parameters are configured, click *Load DICOM* to select the directory containing your DICOM files. Processing may take several minutes depending on hardware performance.

After processing, the DICOM structure is converted to a set of meshes representing segmented materials. Note that these generated meshes are **not** exported as standard geometry objects during penRed export and will be ignored during export. They serve primarily for visualization, checking segmentation accuracy, and defining particle sources, tallies, or combining optional mesh/quadric geometries.

<img src="../images/DICOMloaded.png" alt="Enclosure parameters" width="500" style="display: block; margin: 0 auto"/>

## Exporting

With the DICOM loaded and simulation settings configured, export the DICOM setup via PenRed export using either the **DICOM** or **DICOM+GEO** export types:

- **DICOM**: Exports only the DICOM configuration, ignoring mesh or quadric geometry.
- **DICOM+GEO**: Configures a **COMBO** geometry combining DICOM data with **MESH** or **QUADRIC** geometries. Non-void materials in secondary geometries will overwrite intersecting DICOM voxels.

<img src="../images/DICOMExport.png" alt="Enclosure parameters" width="200" style="display: block; margin: 0 auto"/>
