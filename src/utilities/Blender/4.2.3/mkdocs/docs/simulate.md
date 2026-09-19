# Simulate

In addition to constructing the geometry and configuring the simulation, you can execute the simulation directly within the Blender environment. 

## Prerequisites

To run simulations with penRed, the **[pyPenred](https://pypi.org/project/pyPenred/)** package must be installed in Blender's Python environment. It is already bundled with this addon, but check the [pyPenred install](pyPenred.md) section if you need to install a custom version.

## Running a Simulation

With **pyPenred** properly installed:

1. Locate the *Simulate* button in the *Simulation Properties* panel within the side panel (Press *N* to open it)

<img src="../images/simulationButton.png" alt="Simulation Button" width="300" style="display: block; margin: 0 auto"/>

2. Before simulation begins:
    - An export window will appear, asking you to select the working directory
    - All geometry and configuration files will be exported to this location
    
    **Note**: The working directory will be changed to the selected one, since the simulation uses relative paths to write and read files.

3. During simulation:
    - Progress is displayed via a progress bar in the same side panel
    - The *Cancel* button can be used to interrupt the simulation
   
<img src="../images/simulationProgress.png" alt="Simulation Progress" width="600" style="display: block; margin: 0 auto"/>

4. Simulation finished:
    - Results containing 1D or 2D distributions are plotted automatically and saved within Blender as images
    - To visualize the plots, switch to the *Image Editor*:
    
<img src="../images/openResults.png" alt="To Results" width="600" style="display: block; margin: 0 auto"/>

Within the image editor, every plotted result is listed under the available images. The name of each one contains the tally name along with a brief description suffix:
    
<img src="../images/resultsList.png" alt="To Results" width="600" style="display: block; margin: 0 auto"/>

Notice that these plots are generated automatically. For advanced plotting, the results files can be found in the simulation folder specified during configuration.
