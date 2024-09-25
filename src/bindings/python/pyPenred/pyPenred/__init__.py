# pyPenred/__init__.py
import importlib
import yaml
import time

simulation = importlib.import_module(".simulation", __name__)

simulation.create = simulation.simulator

def runFromFile(configFile = "config.in",
                statusTime = 120,
                configLog = "config.log",
                simulationLog = "simulation.log"):

    #Set logs
    simulation.setConfigurationLog(str(configLog))
    simulation.setSimulationLog(str(simulationLog))

    #Create simulation object
    simu = simulation.create()

    #Try to get the configuration from yaml or penred internal format
    try:
        f = open(configFile)
        d = yaml.load(f, Loader=yaml.SafeLoader)
        confSetErr = simu.configure(d)
    except:
        confSetErr = simu.configFromFile(configFile)

    if(confSetErr != 0):
        print("Invalid configuration file format. See config log\n")
        return 1
        
    print("Configuration set\n")

    #Start the simulation asyncronously
    err = simu.simulate(True)
    if(err != 0):
        print("Error on simulation configuration. See config.log\n")
        return 2

    #Simulation started, check status every 30 seconds
    print("Simulation started\n")
    while simu.isSimulating():
        try:
            time.sleep(statusTime)
        except:
            time.sleep(120)
        status = simu.stringifyStatus()
        for e in status:
            print(e)

    print("Simulation finished")
