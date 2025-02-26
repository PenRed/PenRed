def createDumps(writeFile, wirteInterval, readFile):
    if writeFile or readFile:
        f.write(f"# Dumps configuration\n")

    if writeFile:
        f.write(f"simulation/dump2write \"{writeFile}\"\n")
        f.write(f"simulation/dump-interval \"{wirteInterval}\"\n")
        f.write("\n")
    if readFile:
        f.write(f"simulation/dump2read \"{readFile}\"\n")        
        f.write("\n")

def createResults(finalDump, asciiResults, partialResults):

    f.write(f"# Results configuration\n")

    if finalDump:
        f.write(f"simulation/finalDump true\n")
    else:
        f.write(f"simulation/finalDump false\n")

    if asciiResults:
        f.write(f"simulation/ascii-results true\n")
        if partialResults:
            f.write(f"simulation/partial-results true\n")            
    else:
        f.write(f"simulation/ascii-results false\n")

    f.write("\n")

def createThreads(nThreads, seedPair):
    f.write(f"# Thread configuration\n")

    f.write(f"simulation/threads {nThreads}\n")
    f.write(f"simulation/seedPair {seedPair}\n")
    f.write("\n")

def createSimTime(limit, amount):
    if limit:
        f.write(f"# Simulation time configuration\n")
        f.write(f"simulation/max-time {amount}\n")
        f.write("\n")
    
def createGeo(filename, geoType, detectors, dsmax):
    f.write(f"# Simulation geometry configuration\n")
    
    if geoType == "QUAD":
        f.write("geometry/type \"PEN_QUADRIC\"")
        f.write(f"geometry/input-file \"{filename}\"")
        f.write("geometry/processed-geo-file \"report.geo\"")
        
    elif geoType == "MESH":
        f.write("geometry/type \"MESH_BODY\"")
        f.write(f"geometry/input-file \"{filename}\"")        

    #Detectors
    if len(detectos > 0):
        f.write("\n# - Detectors\n")
    for det in detectos:
        f.write(f"geometry/kdet/{det[0]} {det[1]}\n")

    #DSMAX
    if len(dsmax > 0):
        f.write("\n# - DSMAX\n")
    for ds in dsmax:
        f.write(f"geometry/dsmax/{ds[0]} {ds[1]}\n")

    f.write("\n")
    
