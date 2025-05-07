# Simulate

In addition to constructing the geometry and configuring the simulation, you can execute the simulation directly within the Blender environment. 

## Prerequisites

To run simulations with penRed, the **[pyPenred](https://pypi.org/project/pyPenred/)** package must be installed in Blender's Python environment. Note that this **may differ from your system's default Python installation**.

### Finding Blender's Python Environment

You can identify Blender's Python executable by running the following code in Blender's *Scripting* workspace:

```python
import sys
print(sys.executable)
```

<img src="../images/getPythonVersion.png" alt="Blender Python version" width="500" style="display: block; margin: 0 auto"/>

### Automatic Installation

If **pyPenred** is not installed when initiating a simulation:

1. The plugin will attempt to automatically download and install it
2. If installation fails (e.g., due to OS permissions), you'll need to manually install it following the PenRed package documentation

## Running a Simulation

With **pyPenred** properly installed:

1. Locate the *Simulate* button in the *Simulation Properties* panel (within the *World* tab)

<img src="../images/simulationButton.png" alt="Simulation Button" width="300" style="display: block; margin: 0 auto"/>

2. Before simulation begins:
    - An export window will appear for selecting the working directory
    - All geometry and configuration files will be exported to this location

3. During simulation:
    - Progress is displayed via a progress bar in the 3D viewport
    - The button text changes to *Cancel Simulation* for interruption
   
<img src="../images/simulationProgress.png" alt="Simulation Progress" width="600" style="display: block; margin: 0 auto"/>
