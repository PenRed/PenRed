# Simulate

In addition to constructing the geometry and configuring the simulation, you can execute the simulation directly within the Blender environment. 

## Prerequisites

To run simulations with penRed, the **[pyPenred](https://pypi.org/project/pyPenred/)** package must be installed in Blender's Python environment. Check the installation in the [pyPenred install](pyPenred.md) section.

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
