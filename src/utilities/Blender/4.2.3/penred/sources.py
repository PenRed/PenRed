
def createPSF(f, name, nhist, psfFile, psfMaxE,
              split, window, shift, orientation):
    f.write(f"# PSF Source configuration for '{name}'\n")

    prefix = f"sources/generic/{name}"
    
    # Write number of histories
    f.write(f"{prefix}/nhist {nHists:.5e}\n")
    
    prefix = f"{prefix}/specific"
    f.write(f"{prefix}/type \"PSF\"\n")
    f.write(f"{prefix}/filename \"{psfFile}\"\n")
    f.write(f"{prefix}/Emax \"{psfMaxE}\"\n")

    f.write(f"{prefix}/nsplit {split}\n")
    f.write(f"{prefix}/wght-window [{weight[0]:.5e}, {weight[1]:.5e}]\n\n")

    if shift:
        f.write(f"{prefix}/translation/dx {shift[0]}\n")
        f.write(f"{prefix}/translation/dy {shift[1]}\n")
        f.write(f"{prefix}/translation/dz {shift[2]}\n")

    if orientation:
        f.write(f"{prefix}/rotation/omega {orientation[0]}\n")
        f.write(f"{prefix}/rotation/theta {orientation[1]}\n")
        f.write(f"{prefix}/rotation/phi {orientation[2]}\n")
    

def createGeneric(f, name, nhist, particleType,
                  energyType, energy, spcFile,
                  aperture, direction,
                  position, size):
    f.write(f"\n# Generic Source configuration for '{name}'\n")

    prefix = f"sources/generic/{name}"
    
    # Write number of histories
    f.write(f"{prefix}/nhist {nhist:.5e}\n")

    # Write particle type
    if particleType == "PART_GAMMA":
        f.write(f"{prefix}/kpar \"gamma\"\n")
    elif particleType == "PART_ELECTRON":
        f.write(f"{prefix}/kpar \"electron\"\n")
    elif particleType == "PART_POSITRON":
        f.write(f"{prefix}/kpar \"positron\"\n")
    else:
        f.write(f"{prefix}/kpar \"unknown\"\n")

    # Energy
    f.write(f"\n# - Energy parameters\n")
    if energyType == "ENERGY_MONO":
        f.write(f"{prefix}/energy/type \"MONOENERGETIC\"\n")
        f.write(f"{prefix}/energy/energy {energy:.5e}\n")
    else:
        f.write(f"{prefix}/energy/type \"FILE_SPECTRUM\"\n")
        f.write(f"{prefix}/energy/energy \"{spcFile}\"\n")

    # Direction
    f.write(f"\n# - Direction parameters\n")
    
    f.write(f"{prefix}/direction/type \"SOLID_ANGLE\"\n\n")
    f.write(f"{prefix}/direction/u {direction[0]}\n")
    f.write(f"{prefix}/direction/v {direction[1]}\n")
    f.write(f"{prefix}/direction/w {direction[2]}\n\n")

    f.write(f"{prefix}/direction/theta0 {-aperture/2.0}\n")
    f.write(f"{prefix}/direction/theta1 { aperture/2.0}\n\n")
    
    f.write(f"{prefix}/direction/phi0 0.0\n")
    f.write(f"{prefix}/direction/dphi 360.0\n\n")

    # Spatial
    f.write(f"\n# - Spatial parameters\n")
    
    f.write(f"{prefix}/spatial/type \"BOX\"\n\n")
    f.write(f"{prefix}/spatial/position/x {position[0]}\n")
    f.write(f"{prefix}/spatial/position/y {position[1]}\n")
    f.write(f"{prefix}/spatial/position/z {position[2]}\n\n")

    f.write(f"{prefix}/spatial/size/dx {size[0]}\n")
    f.write(f"{prefix}/spatial/size/dy {size[1]}\n")
    f.write(f"{prefix}/spatial/size/dz {size[2]}\n\n")

def createTime(f, name, timeType, halfLife, timeWindow):

    f.write(f"\n# - Time parameters\n")

    prefix = f"sources/generic/{name}"
    f.write(f"{prefix}/record-time true\n\n")

    if timeType == "TIME_DECAY":
        f.write(f"{prefix}/time/type \"DECAY\"\n\n")
        f.write(f"{prefix}/time/halfLife {halfLife}\n")
        f.write(f"{prefix}/time/time/time0 {timeWindow[0]}\n")
        f.write(f"{prefix}/time/time/time1 {timeWindow[1]}\n\n")
    
    
