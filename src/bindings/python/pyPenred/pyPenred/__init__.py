# pyPenred/__init__.py
import importlib
import yaml
import time
from pyPenred import simulation

simulation.create = simulation.simulator

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
