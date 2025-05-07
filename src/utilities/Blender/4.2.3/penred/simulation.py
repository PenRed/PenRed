#
#
#    Copyright (C) 2022-2025 Universitat de València - UV
#    Copyright (C) 2022-2025 Universitat Politècnica de València - UPV
#    Copyright (C) 2024-2025 Vicent Giménez Alventosa
#
#    This file is part of PenRed: Parallel Engine for Radiation Energy Deposition.
#
#    PenRed is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Affero General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    PenRed is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Affero General Public License for more details.
#
#    You should have received a copy of the GNU Affero General Public License
#    along with PenRed.  If not, see <https://www.gnu.org/licenses/>. 
#
#    contact emails:
#
#        vicent.gimenez.alventosa@gmail.com (Vicent Giménez Alventosa)
#        sanolgi@upvnet.upv.es              (Sandra Oliver Gil)
#        vicente.gimenez@uv.es              (Vicent Giménez Gómez)
#

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
    
