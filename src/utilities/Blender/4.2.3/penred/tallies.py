
def createTallyCylDoseDistrib(f, name, output, rmax, nr, zmin, zmax, nz, nphi, toRound):
    f.write(f"# Cylindrical dose tally configuration for '{name}'\n")
    f.write(f"tallies/{name}/type \"CYLINDRICAL_DOSE_DISTRIB\"\n")
    f.write(f"tallies/{name}/rmin 0.0\n")
    f.write(f"tallies/{name}/rmax {round(rmax,toRound)}\n")
    f.write(f"tallies/{name}/nbinsr {nr}\n")
    f.write(f"tallies/{name}/zmin {round(zmin,toRound)}\n")
    f.write(f"tallies/{name}/zmax {round(zmax,toRound)}\n")
    f.write(f"tallies/{name}/nbinsz {nz}\n")
    f.write(f"tallies/{name}/nbinsPhi {nphi}\n")
    if output:
        f.write(f"tallies/{name}/outputdir \"{output}\"\n")
    f.write("\n")

def createTallyImpactDetector(f, name, output, det,
                              emin, emax, ebins,
                              fluence, flog,
                              spectrum, slog,
                              edep, edlog,
                              age, agelog,
                              agebins, agemin, agemax):
    
    f.write(f"# Impact detector tally configuration for '{name}'\n")
    f.write(f"tallies/{name}/type \"IMPACT_DET\"\n")
    f.write(f"tallies/{name}/detector {det}\n")
    
    if fluence:
        f.write(f"tallies/{name}/fluence true\n")
        if flog:
            f.write(f"tallies/{name}/linearScale-fln true\n")
    else:
        f.write(f"tallies/{name}/fluence false\n")
        
    if spectrum:
        f.write(f"tallies/{name}/spectrum true\n")
        if slog:
            f.write(f"tallies/{name}/linearScale-spc true\n")
    else:
        f.write(f"tallies/{name}/spectrum false\n")
        
    if edep:
        f.write(f"tallies/{name}/energy-dep true\n")
        if edlog:
            f.write(f"tallies/{name}/linearScale-edep true\n")
    else:
        f.write(f"tallies/{name}/energy-dep false\n")
        
    if age:
        f.write(f"tallies/{name}/age true\n")
        if agelog:
            f.write(f"tallies/{name}/linearScale-age true\n")
        f.write(f"tallies/{name}/nbin-age {agebins}\n")
        f.write(f"tallies/{name}/age-min {agemin:.4}\n")
        f.write(f"tallies/{name}/age-max {agemax:.4}\n")
    else:
        f.write(f"tallies/{name}/age false\n")

    if fluence or spectrum or edep:
        f.write(f"tallies/{name}/emin {emin:.3e}\n")
        f.write(f"tallies/{name}/emax {emax:.3e}\n")
        f.write(f"tallies/{name}/nbin-energy {ebins}\n")
    if output:
        f.write(f"tallies/{name}/outputdir \"{output}\"\n")
    f.write("\n")

def createTallySpatialDoseDistrib(f, name, output,
                                  nx, ny, nz,
                                  xmin, xmax,
                                  ymin, ymax,
                                  zmin, zmax,
                                  toRound):
    f.write(f"# Spatial dose distribution tally configuration for '{name}'\n")
    f.write(f"tallies/{name}/type \"SPATIAL_DOSE_DISTRIB\"\n")
    f.write(f"tallies/{name}/xmin {round(xmin, toRound)}\n")
    f.write(f"tallies/{name}/xmax {round(xmax, toRound)}\n")
    f.write(f"tallies/{name}/nx {nx}\n")
    f.write(f"tallies/{name}/ymin {round(ymin, toRound)}\n")
    f.write(f"tallies/{name}/ymax {round(ymax, toRound)}\n")
    f.write(f"tallies/{name}/ny {ny}\n")
    f.write(f"tallies/{name}/zmin {round(zmin, toRound)}\n")
    f.write(f"tallies/{name}/zmax {round(zmax, toRound)}\n")
    f.write(f"tallies/{name}/nz {nz}\n")
    if output:
        f.write(f"tallies/{name}/outputdir \"{output}\"\n")
    f.write("\n")

def createTallySphericalDoseDistrib(f, name, output, rmax, nr, ntheta, nphi, toRound):
    f.write(f"# Spherical dose distribution tally configuration for '{name}'\n")
    f.write(f"tallies/{name}/type \"SPHERICAL_DOSE_DISTRIB\"\n")
    f.write(f"tallies/{name}/rmin 0.0\n")
    f.write(f"tallies/{name}/rmax {round(rmax, toRound)}\n")
    f.write(f"tallies/{name}/nr {nr}\n")
    f.write(f"tallies/{name}/ntheta {ntheta}\n")
    f.write(f"tallies/{name}/nphi {nphi}\n")
    if output:
        f.write(f"tallies/{name}/outputdir \"{output}\"\n")
    f.write("\n")

def createTallyPSF(f, name, output,
                   det, emin, emax, gamma, electron, positron):
    f.write(f"# Phase Space File (PSF) tally configuration for '{name}'\n")
    f.write(f"tallies/{name}/type \"PSF\"\n")
    f.write(f"tallies/{name}/detector {det}\n")
    f.write(f"tallies/{name}/emin {emin:.3e}\n")
    f.write(f"tallies/{name}/emax {emax:.3e}\n")

    if gamma:
        f.write(f"tallies/{name}/particles/gamma true\n")
    else:
        f.write(f"tallies/{name}/particles/gamma false\n")

    if electron:
        f.write(f"tallies/{name}/particles/electron true\n")
    else:
        f.write(f"tallies/{name}/particles/electron false\n")

    if positron:
        f.write(f"tallies/{name}/particles/positron true\n")
    else:
        f.write(f"tallies/{name}/particles/positron false\n")
    if output:
        f.write(f"tallies/{name}/outputdir \"{output}\"\n")
    f.write("\n")

def createTallyKerma(f, name, output,
                     emin, emax,
                     dataPrefix, materialList, meshType,
                     n1, n2, n3,
                     min1, min2, min3,
                     max1, max2, max3,
                     toRound):
    f.write(f"# Kerma tally configuration for '{name}'\n")
    f.write(f"tallies/{name}/type \"KERMA_TRACK_LENGTH\"\n")
    f.write(f"tallies/{name}/emin {emin:.3e}\n")
    f.write(f"tallies/{name}/emax {emax:.3e}\n")

    if meshType == 0: # Cartesian
        f.write(f"tallies/{name}/cartesian/nx {n1}\n")
        f.write(f"tallies/{name}/cartesian/xmin {round(min1,toRound)}\n")
        f.write(f"tallies/{name}/cartesian/xmax {round(max1,toRound)}\n")
        
        f.write(f"tallies/{name}/cartesian/ny {n2}\n")
        f.write(f"tallies/{name}/cartesian/ymin {round(min2,toRound)}\n")
        f.write(f"tallies/{name}/cartesian/ymax {round(max2,toRound)}\n")

        f.write(f"tallies/{name}/cartesian/nz {n3}\n")
        f.write(f"tallies/{name}/cartesian/zmin {round(min3,toRound)}\n")
        f.write(f"tallies/{name}/cartesian/zmax {round(max3,toRound)}\n")
        
    elif meshType == 1: # Cylindrical
        f.write(f"tallies/{name}/cylindrical/nr {n1}\n")
        f.write(f"tallies/{name}/cylindrical/rmin 0.0\n")
        f.write(f"tallies/{name}/cylindrical/rmax {round(max1,toRound)}\n")
        
        f.write(f"tallies/{name}/cylindrical/nphi {n2}\n")

        f.write(f"tallies/{name}/cylindrical/nz {n3}\n")
        f.write(f"tallies/{name}/cylindrical/zmin {round(min3,toRound)}\n")
        f.write(f"tallies/{name}/cylindrical/zmax {round(max3,toRound)}\n")
        
    elif meshType == 2: # Spherical
        f.write(f"tallies/{name}/spherical/nr {n1}\n")
        f.write(f"tallies/{name}/spherical/rmin 0.0\n")
        f.write(f"tallies/{name}/spherical/rmax {round(max1,toRound)}\n")
        
        f.write(f"tallies/{name}/spherical/ntheta {n2}\n")
        f.write(f"tallies/{name}/spherical/nphi {n3}\n")
        
    for mat in materialList:
        f.write(f"tallies/{name}/dataFiles/{mat[0]} \"{dataPrefix}{mat[1]}.muen\"\n")

    if output:
        f.write(f"tallies/{name}/outputdir \"{output}\"\n")
    f.write("\n")

def createTallySpatialDistrib(f, name, output, det,
                              nx, ny, nz,
                              xmin, ymin, zmin,
                              xmax, ymax, zmax,
                              ebins, emin, emax,
                              particle,
                              coordinates,
                              bins, toRound):
    f.write(f"# Spatial distribution tally configuration for '{name}'\n")
    f.write(f"tallies/{name}/type \"DETECTION_SPATIAL_DISTRIB\"\n")
    f.write(f"tallies/{name}/detector {det}\n")
    
    f.write(f"tallies/{name}/spatial/nx {nx}\n")
    f.write(f"tallies/{name}/spatial/xmin {round(xmin, toRound)}\n")
    f.write(f"tallies/{name}/spatial/xmax {round(xmax, toRound)}\n")
    f.write(f"tallies/{name}/spatial/ny {ny}\n")
    f.write(f"tallies/{name}/spatial/ymin {round(ymin, toRound)}\n")
    f.write(f"tallies/{name}/spatial/ymax {round(ymax, toRound)}\n")
    f.write(f"tallies/{name}/spatial/nz {nz}\n")
    f.write(f"tallies/{name}/spatial/zmin {round(zmin, toRound)}\n")
    f.write(f"tallies/{name}/spatial/zmax {round(zmax, toRound)}\n")

    f.write(f"tallies/{name}/energy/nbins {ebins}\n")
    f.write(f"tallies/{name}/energy/emin {emin:.3e}\n")
    f.write(f"tallies/{name}/energy/emax {emax:.3e}\n")
    
    f.write(f"tallies/{name}/particle \"{particle}\"\n")

    if coordinates:
        f.write(f"tallies/{name}/printCord true\n")
    else:
        f.write(f"tallies/{name}/printCord false\n")

    if bins:
        f.write(f"tallies/{name}/printBins true\n")
    else:
        f.write(f"tallies/{name}/printBins false\n")

    if output:
        f.write(f"tallies/{name}/outputdir \"{output}\"\n")
    f.write("\n")

def createTallyAngularDet(f, name, output, det, emin, emax, ebins,
                          theta1, theta2, phi1, phi2, logSale):
    f.write(f"# Angular detector tally configuration for '{name}'\n")
    f.write(f"tallies/{name}/type \"ANGULAR_DET\"\n")
    f.write(f"tallies/{name}/emin {emin:.3e}\n")
    f.write(f"tallies/{name}/emax {emax:.3e}\n")
    f.write(f"tallies/{name}/nBinsE {ebins}\n")
    f.write(f"tallies/{name}/theta1 {theta1:.3f}\n")
    f.write(f"tallies/{name}/theta2 {theta2:.3f}\n")
    f.write(f"tallies/{name}/phi1 {phi1:.3f}\n")
    f.write(f"tallies/{name}/phi2 {phi2:.3f}\n")
    if logSale:
        f.write(f"tallies/{name}/linearScale false\n")
    else:
        f.write(f"tallies/{name}/linearScale true\n")        
        
    if output:
        f.write(f"tallies/{name}/outputdir \"{output}\"\n")
    f.write("\n")

def createTallyEmergingParticles(f, name, output, emin, emax, ebins, thetaBins, phiBins):
    f.write(f"# Emerging particles distribution tally configuration for '{name}'\n")
    f.write(f"tallies/{name}/type \"EMERGING_PART_DISTRIB\"\n")
    f.write(f"tallies/{name}/emin {emin:.3e}\n")
    f.write(f"tallies/{name}/emax {emax:.3e}\n")
    f.write(f"tallies/{name}/nBinsE {ebins}\n")
    f.write(f"tallies/{name}/nBinsTheta {thetaBins}\n")
    f.write(f"tallies/{name}/nBinsPhi {phiBins}\n")
    if output:
        f.write(f"tallies/{name}/outputdir \"{output}\"\n")
    f.write("\n")
    
def createTallyTrack(f, name, output, nHists):
    f.write(f"# Track tally configuration for '{name}'\n")
    f.write(f"tallies/{name}/type \"TRACK\"\n")
    f.write(f"tallies/{name}/nhists {nHists}\n\n")
    if output:
        f.write(f"tallies/{name}/outputdir \"{output}\"\n")
    f.write("\n")
