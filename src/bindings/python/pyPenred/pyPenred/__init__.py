# pyPenred/__init__.py
import importlib
import yaml
import time
import os
import sys
import numpy as np
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import locale
from pyPenred import data
from pyPenred import geometry
from pyPenred import psf
from pyPenred import simulation

simulation.create = simulation.simulator
geometry.create = geometry.geometry

def _plotRes1D(self, filename = "", with_errors=True, nsigma=2):
    if not isinstance(self, data.results1D):
        raise TypeError(f"Expected a data.results1D object, got {type(self)}")

    values, sigma, xinfo, vheader, _ = self.data()
    xPoints = np.linspace(xinfo[0], xinfo[1], values.shape[0])

    # Create a Figure
    fig = Figure(figsize=(8, 5), dpi=300)
    ax = fig.add_subplot(111)

    if with_errors:
        ax.errorbar(xPoints, values, yerr=nsigma*sigma,
                    capsize=3, fmt="r--x", ecolor="black")
    else:
        ax.plot(xPoints, values)

    if self.title():
        ax.set_title(self.title())
    if xinfo[2]:
        ax.set_xlabel(xinfo[2])
    if vheader:
        ax.set_ylabel(vheader)

    # Check if the image must be saved to a file
    if filename:
        fig.savefig(filename, dpi=300)
        
    # Restore locale after building the figure
    locale.setlocale(locale.LC_TIME, 'C')
    
    return fig

def _plotRes2D(self, filename="", with_errors=True, nsigma=2, max_rel_error=0.20):
    """
    Plot a 2D result.

    Parameters
    ----------
    with_errors : bool
        If True, add a second panel with the relative error.
    nsigma : float
        Number of standar deviations used to calculate the relative error.
    max_rel_error : float
        Maximum absolute relative error shown on the color scale.
        Default 0.20 (i.e. ±20%). If the data maximum is smaller,
        the range is shrunk automatically to preserve contrast.
    """
    if not isinstance(self, data.results2D):
        raise TypeError(f"Expected a data.results2D object, got {type(self)}")

    values, sigma, yinfo, xinfo, vheader, _ = self.data()
    extent = [xinfo[0], xinfo[1], yinfo[0], yinfo[1]]

    if with_errors:
        fig = Figure(figsize=(14, 6), dpi=300)
        ax_val, ax_err = fig.subplots(1, 2)
    else:
        fig = Figure(figsize=(8, 6), dpi=300)
        ax_val = fig.add_subplot(111)
        ax_err = None

    # --- Value panel ---
    im_val = ax_val.imshow(
        values,
        origin="lower",
        extent=extent,
        cmap="viridis",
        aspect="equal",
    )
    if ax_err is not None:
        ax_val.set_title("Values")
    if xinfo[2]:
        ax_val.set_xlabel(xinfo[2])
    if yinfo[2]:
        ax_val.set_ylabel(yinfo[2])
    if vheader:
        fig.colorbar(im_val, ax=ax_val, label=vheader)

    # --- Relative error panel ---
    if ax_err is not None:
        with np.errstate(divide="ignore", invalid="ignore"):
            rel_err = np.where(
                values != 0,
                nsigma * sigma / values,
                np.nan,
            )

        # Cells with zero value -> show as zero error, not NaN
        rel_err = np.nan_to_num(rel_err, nan=0.0)

        # Determine the color range

        # - vmax
        data_max = np.max(rel_err) if rel_err.size else 0.0
        if data_max == 0.0:
            vmax = max_rel_error
        elif data_max < max_rel_error:
            # Shrink the range to the data, with 1% headroom
            vmax = data_max * 1.01
        else:
            vmax = max_rel_error

        # - vmin
        data_min = np.min(rel_err) if rel_err.size else 0.0
        if data_min == 0.0:
            vmin = 0.0
        if data_min > -max_rel_error:
            vmin = data_min * 0.99
        else:
            vmin = -max_rel_error
            
        im_err = ax_err.imshow(
            100.0*rel_err,
            origin="lower",
            extent=extent,
            cmap="magma",
            vmin=vmin*100.0,
            vmax=vmax*100.0,
            aspect="equal",
        )
        ax_err.set_title(f"Relative Uncertainty ({nsigma} $\sigma$)")
        if xinfo[2]:
            ax_err.set_xlabel(xinfo[2])
        if yinfo[2]:
            ax_err.set_ylabel(yinfo[2])
        fig.colorbar(im_err, ax=ax_err, label="Relative Error (%)")

    if self.title():
        fig.suptitle(f"{self.title()}")
        fig.tight_layout(rect=[0, 0, 1, 0.95])
    else:
        fig.tight_layout()

    # Check if the image must be saved to a file
    if filename:
        fig.savefig(filename, dpi=300)
    
    # Restore locale after building the figure
    locale.setlocale(locale.LC_TIME, 'C')
    
    return fig
        
data.results1D.plot = _plotRes1D
data.results2D.plot = _plotRes2D

def readConfigFile(filename):
    '''
    Read a configuration file in both, YAML or penRed
    internal format, and returns the associated dictionary.
    '''
    try:
        f = open(configFile)
        d = yaml.load(f, Loader=yaml.SafeLoader)
        return d
    except:
        yamlstr = simulation.configFile2YAML(filename)
        return yaml.safe_load(yamlstr)

def writeConfigFile(filename, d):
    '''
    Write a configuration file in penRed internal format,
    from the data stored in the dictionary.
    '''
    configText = simulation.dict2SectionString(d)
    with open(filename, "w") as f:
        f.write(configText)
    
def runFromFile(configFile = "config.in",
                statusTime = 20,
                configLog = "config.log",
                simulationLog = "simulation.log"):
    '''
    Configures and runs a simulation from the specified file.
    Reports the simulation status periodically.
    '''

    #Set logs
    simulation.setConfigurationLog(str(configLog))
    simulation.setSimulationLog(str(simulationLog))

    #Create simulation object
    simu = simulation.create()

    #Try to get the configuration from yaml or penred internal format
    try:
        f = open(configFile)
        d = yaml.load(f, Loader=yaml.SafeLoader)
        simu.configure(d)
    except:
        simu.configFromFile(configFile)
        
    print("Configuration set\n")

    #Start the simulation asynchronously
    simu.simulate(True)

    #Simulation started, check status every 30 seconds
    print("Simulation started\n")
    while simu.isSimulating():
        try:
            time.sleep(statusTime)
        except:
            time.sleep(20)
        status = simu.stringifyStatus()
        for e in status:
            print(e)

    print("Simulation finished")

    return simu


def plotResults(simObj, tallyName : str,
                prefix : str = "",
                titles : list = [],
                fontSize : int = 20,
                labelPad : int = 20,
                width : int = 12,
                height : int = 9,
                markerSize : int = 6,
                capSize : int = 5,
                aspect : str = "auto",
                vrange : list = [],
                erange : list = [],
                grid : bool = False,
                show : bool = False):
    '''
    Read and plot the results from a simulation object with a finished simulation.
    This function requires the matplotlib module, which will be installed in the
    system if it is not.
    '''
    
    if not isinstance(simObj, simulation.simulator):
        raise TypeError("Expected a Simulation object. Create it with 'simulation.create()'")    
    if not tallyName:
        raise ValueError("Empty tally name")

    # Import gnuplot
    try:
        import matplotlib.pyplot as plt
    except:
        # Not found, try to install it
        import subprocess
        subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'matplotlib'])
        import matplotlib.pyplot as plt

    # Create the plot folder
    plotDir = os.path.join(prefix, "plots")
    if not os.path.exists(plotDir):
        os.makedirs(plotDir)
        
    # Get simulation results with dimensions information
    results = simObj.getResults(tallyName, True)

    # Iterate over results
    something2Plot = False
    resultIndex = 0
    for result in results:
        #Check if this result has something to plot
        if len(result[0]) > 1:
            # Plot this results
            something2Plot = True

            # Get the shape
            shape = result[0].shape

            # Get the dimensions
            nDim = len(shape)

            # Get title, if provided
            defaultTitle = f"Result {resultIndex}"
            title = defaultTitle
            if len(titles) > resultIndex:
                title = str(titles[resultIndex])

            # Configure the plot title and output
            if title:
                filename = os.path.join(plotDir,title.replace(" ","_"))
            else:
                filename = os.path.join(plotDir,defaultTitle.replace(" ","_"))                

            if nDim == 1:

                x = np.linspace(result[2][0], result[2][1], shape[0], endpoint=False)

                # One dimension plot
                plt.figure(figsize=(width, height))

                plt.title(title, fontsize=fontSize)
                plt.xlabel(result[2][2], fontsize=fontSize, labelpad=labelPad)
                plt.ylabel(result[3], fontsize=fontSize, labelpad=labelPad)
                if grid:
                    plt.grid(True, linestyle='--', alpha=0.7)
                else:
                    plt.grid(False)

                plt.errorbar(x, result[0], yerr=result[1],
                             fmt='o', color='red', ecolor='red',
                             label='', markersize=markerSize,
                             capsize=capSize)

                
                plt.savefig(filename)
                if show:
                    plt.show()
                plt.close()

            else:

                # Get X and Y extents for 2D ploting
                xDim = result[2+nDim-1]                
                yDim = result[2+nDim-2]

                *outerDims, h, w = shape

                for iPlane in np.ndindex(*outerDims):

                    # + Print values
                    plt.figure(figsize=(width, height))

                    if len(vrange) > 1:
                        img = plt.imshow(result[0][iPlane],
                                         extent=[xDim[0], xDim[1], yDim[0], yDim[1]],
                                         aspect=aspect,
                                         origin="lower",
                                         cmap="viridis",
                                         vmin=vrange[0],
                                         vmax=vrange[1])
                    else:
                        img = plt.imshow(result[0][iPlane],
                                         extent=[xDim[0], xDim[1], yDim[0], yDim[1]],
                                         aspect=aspect,
                                         origin="lower",
                                         cmap="viridis")

                    plt.xlim(xDim[0], xDim[1])
                    plt.ylim(yDim[0], yDim[1])

                    # Set colorbar label
                    cbar = plt.colorbar(img)
                    cbar.set_label(result[2+nDim], fontsize = fontSize)

                    if iPlane:
                        plt.title(title + f" {iPlane}", fontsize=fontSize)
                    else:
                        plt.title(title, fontsize=fontSize)
                
                    plt.xlabel(xDim[2], fontsize=fontSize, labelpad=labelPad)
                    plt.ylabel(yDim[2], fontsize=fontSize, labelpad=labelPad)
                
                    if iPlane:
                        plt.savefig(filename + f"_{iPlane}")
                    else:
                        plt.savefig(filename)
                        
                    if show:
                        plt.show()
                    plt.close()                

                    # + Plot relative errors
                    plt.figure(figsize=(width, height))
                    relErr = np.divide(result[1][iPlane],
                                       result[0][iPlane],
                                       where=result[0][iPlane] != 0)
                    relErr[result[0][iPlane] == 0.0] = 0.0

                    if len(erange) > 1:
                        img = plt.imshow(100.0*relErr,
                                         extent=[xDim[0], xDim[1], yDim[0], yDim[1]],
                                         aspect=aspect,
                                         origin="lower",
                                         cmap="viridis",
                                         vmin=erange[0],
                                         vmax=erange[1])
                    else:
                        img = plt.imshow(100.0*relErr,
                                         extent=[xDim[0], xDim[1], yDim[0], yDim[1]],
                                         aspect=aspect,
                                         origin="lower",
                                         cmap="viridis")

                    plt.xlim(xDim[0], xDim[1])
                    plt.ylim(yDim[0], yDim[1])

                    # Set colorbar label
                    cbar = plt.colorbar(img)
                    cbar.set_label("Error (%)", fontsize = fontSize)

                    if iPlane:
                        plt.title(title + f" {iPlane}", fontsize=fontSize)
                    else:
                        plt.title(title, fontsize=fontSize)
                
                    plt.xlabel(xDim[2], fontsize=fontSize, labelpad=labelPad)
                    plt.ylabel(yDim[2], fontsize=fontSize, labelpad=labelPad)
                
                    if iPlane:
                        plt.savefig(filename + f"_error_{iPlane}")
                    else:
                        plt.savefig(filename + f"_error")
                        
                    if show:
                        plt.show()
                    plt.close()                
                    
        resultIndex += 1
