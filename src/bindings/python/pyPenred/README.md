# pyPenred - Python Bindings for PenRed

<p align="center">
  <a href="https://github.com/PenRed/PenRed">
    <img src="logo.png" width="200" alt="PenRed logo">
  </a>
</p>
https://penred.github.io/PenRed/pyPenred/
Python interface for the [PenRed](https://github.com/PenRed/PenRed) Monte Carlo radiation transport simulator, enabling programmatic control and automation.

## Key Features

Currently pyPenred offers:

- **Direct Simulation Control**  
Run PenRed simulations directly from Python scripts

- **Configuration management**  
Convert between PenRed native configuration format, YAML and Python dicts

- **Workflow Automation**  
Chain multiple simulations with parameter sweeps

- **External Software Integration**  
Embed PenRed in custom applications and pipelines

## Installation

### From PyPI (recommended)

pyPenred can be installed directly from pip as a regular package:

```bash
pip install pyPenred
```
### From Source

If the user wants to modify the source code or compile it with native optimizations, another option consists on compiling the binding locally using the *installPyPenred* scripts in the PenRed repository [src](https://github.com/PenRed/PenRed/tree/master/src) folder.

The requirements to compile the code are:
- Python 3.8+ development headers
- CMake 3.12+
- C++14 compatible compiler

If the requirements are fulfilled, just run the corresponding script.

Linux/macOS:

```bash
./installPyPenred.sh
```

Windows:

```bash
installPyPenred.bat
```

## Configuration Format Conversion

pyPenred implements a bi-directional conversion between PenRed's native configuration format and Python dictionaries. Parameter paths with slashes map to nested dictionary structures:

| Native Format | Python Dictionary Equivalent |
|---------------|------------------------------|
| `/geometry/phantom/material` | `config['geometry']['phantom']['material']` |

## Basic Usage

### Running a Simulation

To run a simulation using a configuration file located in the current folder, named *config.in*, just use:

```python
import pyPenred
pyPenred.runFromFile("config.in")  # Supports both native and YAML formats
```

### Modifying Configurations

To modify an existing configuration file (*config.in*) either in regular PenRed or YAML formats, the following code can be used:

```python
config = pyPenred.readConfigFile("config.in") # Load the file as a dictionary
config['path']['to']['modify'] = 42  # Modify parameters
pyPenred.writeConfigFile("modified.in", config)  # Export to native format
```

### Complete Example

A more complete example can be found in the *runFromFile* function code. This one configures the simulation and runs it asynchronously reporing the simulation state in regular intervals. Take into account that this function is defined within the pyPenred module. Therefore, the *pyPenred* module prefix must be used in regular scripts:

```python
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
```

## Documentation

📚 [pyPenred API](https://penred.github.io/PenRed/pyPenred/)

📚 [PenRed Documentation](https://github.com/PenRed/PenRed/blob/master/doc/PenRed_usage_guide.pdf)  

## Support

🐞 [Report Issues](https://github.com/PenRed/PenRed/issues)  
💬 [Discussion Forum](https://github.com/PenRed/PenRed/discussions)

## API

To build the authomatic API documentation, you need the sphinx package:

```python
pip install sphinx furo sphinx-copybutton sphinx-autodoc-typehints sphinxcontrib-napoleon myst-parser
```

once installed, build the documentation running the following instructions inside the *docs* folder

```bash
make html
```

Finally, the documentation will be available opening the file *docs/build/html/index.html* with a web browser.
