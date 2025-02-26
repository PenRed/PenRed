# PenRed Blender Plugin Documentation

Welcome to the documentation for the **PenRed Blender Plugin**. This plugin integrates PenRed's simulation capabilities directly into Blender, allowing users to construct geometries, configure simulations, and export everything to the PenRed format—all within the Blender environment.

This documentation covers the plugin's functionalities in detail. However, some PenRed-specific components, such as geometry file formats, are not described extensively here. For a comprehensive understanding of these components, please refer to the official [PenRed documentation](https://github.com/PenRed/PenRed/tree/master/doc) available in the GitHub repository.

---

## Sections

### **Geometry Construction**
Create simulation geometries using Blender's powerful modeling tools. The plugin supports:

- **Quadric Surfaces**: Predefined quadric objects for constructing geometries.

- **Triangular Meshes**: Convert any Blender mesh into a PenRed-compatible format.

For more details, see the [Geometry Construction](geometry-construction/index.md) section.

---

### **Simulation Configuration**
Configure the entire simulation within Blender, including:
- **Particle Sources**: Define particle sources with customizable energy, direction, and time sampling.

- **Tallies**: Set up detectors and scoring regions to extract simulation data.

  - **Object Tallies**: Measure quantities within specific spatial regions.
  
  - **Detector Tallies**: Record particle interactions in detector bodies.
  
  - **World Tallies**: Capture global simulation data.
  
- **Materials**: Define material properties and assign them to geometry bodies.

For more details, see the [Simulation Configuration](simulation-configuration/index.md) section.

---

### **Export**
Export geometries and simulation configurations to the PenRed format. The plugin supports:

- **Quadric Surfaces**: Export predefined quadric objects.

- **Triangular Meshes**: Export any Blender mesh as a triangular mesh surface.

- **Configuration Files**: Save simulation settings for use in PenRed.

For more details, see the [Export](export.md) section.

---

## Getting Started

1. **Install the Plugin**:

    - Go to Blender's **Preferences** > **Add-ons**.
    
    - Install the plugin located in the **src>utilities>Blender** folder in the PenRed package and enable it.

2. **Explore the Plugin**:

    - Use the Blender and **Geometry Construction** tools to build your simulation geometry.
    
    - Configure your simulation using the **Simulation Configuration** options.
    
    - Export your work to the PenRed format using the **Export** functionality.

3. **Refer to the Documentation**:

    - For detailed instructions, explore the sections listed above.
    
    - For advanced PenRed features, consult the official [PenRed documentation](https://github.com/PenRed/PenRed/tree/master/doc).

---

## Additional Resources
- **GitHub Repository**: [PenRed/PenRed](https://github.com/PenRed/PenRed)
- **Blender Documentation**: [Blender Manual](https://docs.blender.org/manual/en/latest/)

---
