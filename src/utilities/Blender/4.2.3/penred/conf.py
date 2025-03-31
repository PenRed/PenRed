import bpy

from . import sources, utils, tallies, materials, vr
from math import sqrt

def createSim(context, f):
    
    f.write("## Simulation Configuration ##\n\n")
    
    # World
    world = context.scene.world
    
    if hasattr(world, "penred_settings"):

        simSett = world.penred_settings.simulation

        f.write(f"# Simulation configuration\n")
        if simSett.enableDumps:
            f.write(f"simulation/dump-interval {simSett.dumpInterval:.2f}\n")
            f.write(f"simulation/dump2write \"{simSett.dumpWriteFile}\"\n")

        if simSett.readDumps:
            f.write(f"simulation/dump2read \"{simSett.dumpReadFile}\"\n")

        if simSett.finalDump:
            f.write(f"simulation/finalDump true\n")
        else:
            f.write(f"simulation/finalDump false\n")
            
        if simSett.asciiResults:
            f.write(f"simulation/ascii-results true\n")
            if simSett.partialResults:
                f.write(f"simulation/partial-results true\n")
        else:
            f.write(f"simulation/ascii-results false\n")

        if simSett.threadSelType == "AUTO":
            f.write("simulation/threads 0\n")
        else:
            f.write(f"simulation/threads {simSett.nThreads}\n")

        f.write(f"simulation/seedPair {simSett.seedPair}\n")

        if simSett.limitSimTime:
            f.write(f"simulation/max-time {simSett.maxSimTime:.2f}\n")
        
        f.write("\n")

def createMaterials(context, f):

    f.write("## Materials ##\n\n")

    # World tallies
    world = context.scene.world
    
    if hasattr(world, "penred_settings"):

        for i, item in enumerate(world.penred_settings.materials):

            # Get cutoffs
            if item.gammaCutoffType == "CUTOFF_EABS":
                gCutoff = item.gammaEabs
            else:
                gCutoff = item.gammaRange

            if item.electronCutoffType == "CUTOFF_EABS":
                eCutoff = item.electronEabs
            else:
                eCutoff = item.electronRange

            if item.positronCutoffType == "CUTOFF_EABS":
                pCutoff = item.positronEabs
            else:
                pCutoff = item.positronRange

            # Create composition list
            compo = []
            for e in item.composition:
                compo.append([e.z, e.wFraction])

            materials.createMaterial(f, item.name, i+1,
                                     item.gammaCutoffType, gCutoff,
                                     item.electronCutoffType, eCutoff,
                                     item.positronCutoffType, pCutoff,
                                     item.density, compo,
                                     item.enableAdvanced,
                                     item.C1, item.C2, item.WCC, item.WCR)


def createSources(context, f, toRound):

    f.write("## Sources ##\n\n")
    
    for obj in context.scene.objects:
        if hasattr(obj, "penred_settings"):
            source = obj.penred_settings.source
            if source.enabled:

                # Get object properties
                x,y,z,dx,dy,dz,sx,sy,sz,omega,theta,phi,name,bsize = utils.getObjInfo(obj)
                
                if source.particleType == "PART_PSF":

                    # Create the PSF source
                    sources.createPSF(f, name, source.nHists,
                                      source.sourcePSF, source.psfMaxE,
                                      source.split, source.psfWindow,
                                      (x,y,z), (omega, theta, phi), toRound)
                else:
                    if source.spatialType == "SPATIAL_CYL":
                        if source.spatialBBFit:
                            radius = min((bsize[0], bsize[1]))/2.0
                        else:
                            radius = sqrt(bsize[0]*bsize[0] + bsize[1]*bsize[1])/2.0
                                         
                        spatialSize = (0.0, radius, 0.0, bsize[2])
                    else:
                        spatialSize = bsize

                    if source.enableSourceMat:
                        sourceMat = source.sourceMat
                    else:
                        sourceMat = -1
                    sources.createGeneric(f, name, source.nHists,
                                          source.particleType,
                                          source.energyType,
                                          source.energy, source.spcFile,
                                          source.aperture, source.direction,
                                          source.spatialType, (x,y,z), spatialSize,
                                          sourceMat, toRound)

                # Check time sampling
                if source.timeRecord:
                    sources.createTime(f, name, source.timeType,
                                       source.decayHalf, source.timeWindow)

def createTallies(context, f, toRound):

    f.write("## Tallies ##\n\n")
    
    # World tallies
    world = context.scene.world
    
    if hasattr(world, "penred_settings"):
        
        # Get output folder
        outputPrefix = world.penred_settings.simulation.outputPrefix

        # Emerging particle tallies
        for i, item in enumerate(world.penred_settings.talliesEmergingParticle):
            tallyName = f"{i}_{item.name}"
            tallies.createTallyEmergingParticles(f, tallyName, outputPrefix,
                                                 item.emin, item.emax, item.ebins,
                                                 item.nTheta, item.nPhi)

        # Track tally
        tracks = world.penred_settings.tracksTally
        if tracks.enable:
            tallyName = "track"
            tallies.createTallyTrack(f, tallyName, outputPrefix, tracks.nHists)

    else:
        outputPrefix = ""

    # Object tallies
    for obj in context.scene.objects:
        if hasattr(obj, "penred_settings"):

            # Get object properties
            x,y,z,dx,dy,dz,sx,sy,sz,omega,theta,phi,name,bsize = utils.getObjInfo(obj)

            xmin = x-bsize[0]/2.0
            xmax = x+bsize[0]/2.0

            ymin = y-bsize[1]/2.0
            ymax = y+bsize[1]/2.0

            zmin = z-bsize[2]/2.0
            zmax = z+bsize[2]/2.0
            
            rIn = min((bsize[0],bsize[1],bsize[2]))/2.0
            rcylOut = sqrt(bsize[0]*bsize[0] + bsize[1]*bsize[1])/2.0
            rsphOut = sqrt(bsize[0]*bsize[0] + bsize[1]*bsize[1] + bsize[2]*bsize[2])/2.0
            
            # Cylindrical dose tallies
            for i, item in enumerate(obj.penred_settings.talliesCylDose):
                tallyName = f"{obj.name}_{i}_{item.name}"
                if item.spatialBBFit:
                    radius = rIn
                else:
                    radius = rcylOut
                
                tallies.createTallyCylDoseDistrib(f, tallyName, outputPrefix,
                                                  radius, item.nr,
                                                  zmin, zmax,
                                                  item.nz, item.nPhi,
                                                  toRound)

            # Spatial Dose Distribution
            for i, item in enumerate(obj.penred_settings.talliesSpatialDoseDistrib):
                tallyName = f"{obj.name}_{i}_{item.name}"
                tallies.createTallySpatialDoseDistrib(f, tallyName, outputPrefix,
                                                      item.nx, item.ny, item.nz,
                                                      xmin, xmax,
                                                      ymin, ymax,
                                                      zmin, zmax,
                                                      toRound)

            # Spherical dose distribution
            for i, item in enumerate(obj.penred_settings.talliesSphericalDoseDistrib):
                tallyName = f"{obj.name}_{i}_{item.name}"
                if item.spatialBBFit:
                    radius = rIn
                else:
                    radius = rsphOut
                
                tallies.createTallySphericalDoseDistrib(f, tallyName, outputPrefix,
                                                        radius,
                                                        item.nr, item.ntheta, item.nphi,
                                                        toRound)

            # Kerma
            if len(obj.penred_settings.talliesKerma) > 0:
                # Get world
                world = context.scene.world
                
                # Get material list
                matsIndex = utils.overlappingMats(obj)
                matPairs = []
                for m in matsIndex:
                    if len(world.penred_settings.materials) >= m:
                        matName = world.penred_settings.materials[m-1].name
                    else:
                        matName = f"mat{m}"
                    matPairs.append([m, matName])
            
            for i, item in enumerate(obj.penred_settings.talliesKerma):
                tallyName = f"{obj.name}_{i}_{item.name}"


                if item.meshType == "MESH_CART":

                    meshType = 0
                    
                    n1 = item.nx
                    min1 = xmin
                    max1 = xmax

                    n2 = item.ny
                    min2 = ymin
                    max2 = ymax

                    n3 = item.nz
                    min3 = zmin
                    max3 = zmax
                    
                elif item.meshType == "MESH_CYL":
                    
                    meshType = 1

                    if item.spatialBBFit:
                        radius = rIn
                    else:
                        radius = rcylOut
                    
                    n1 = item.nr
                    min1 = 0.0
                    max1 = radius
                    
                    n2 = item.nphi
                    min2 = 0.0
                    max2 = 0.0
                    
                    n3 = item.nz
                    min3 = zmin
                    max3 = zmax
                    
                elif item.meshType == "MESH_SPH":
                    
                    meshType = 2
                    
                    if item.spatialBBFit:
                        radius = rIn
                    else:
                        radius = rsphOut
                    
                    n1 = item.nr
                    min1 = 0.0
                    max1 = radius
                    
                    n2 = item.ntheta
                    min2 = 0.0
                    max2 = 0.0
                    
                    n3 = item.nphi
                    min3 = 0.0
                    max3 = 0.0
                    
                tallies.createTallyKerma(f, tallyName, outputPrefix,
                                         item.emin, item.emax,
                                         item.dataPath, matPairs,
                                         meshType,
                                         n1, n2, n3,
                                         min1, min2, min3,
                                         max1, max2, max3,
                                         toRound)

            # Detector based tallies

            if not obj.penred_settings.isDetector:
                continue

            # Get detector index
            det = obj.penred_settings.detector
                
            # Impact detector tallies
            for i, item in enumerate(obj.penred_settings.talliesImpactDetector):
                tallyName = f"{obj.name}_{i}_{item.name}"
                tallies.createTallyImpactDetector(f, tallyName, outputPrefix, det,
                                                  item.emin, item.emax, item.ebins,
                                                  item.fluence, item.fluenceLogScale,
                                                  item.spectrum, item.spectrumLogScale,
                                                  item.edep, item.edepLogScale,
                                                  item.age, item.ageLogScale,
                                                  item.ageBins, item.ageMin, item.ageMax)

            # Spatial distribution tally
            for i, item in enumerate(obj.penred_settings.talliesSpatialDistrib):
                tallyName = f"{obj.name}_{i}_{item.name}"
                if item.particleType == "PART_GAMMA":
                    particle = "gamma"
                elif item.particleType == "PART_ELECTRON":
                    particle = "electron"
                elif item.particleType == "PART_POSITRON":
                    particle = "positron"
                else:
                    particle = "unknown"
                    
                tallies.createTallySpatialDistrib(f, tallyName, outputPrefix, det,
                                                  item.nx, item.ny, item.nz,
                                                  xmin, ymin, zmin,
                                                  xmax, ymax, zmax,
                                                  item.ebins, item.emin, item.emax,
                                                  particle,
                                                  item.printCoordinates,
                                                  item.printBins, toRound)                

            # PSF
            for i, item in enumerate(obj.penred_settings.talliesPSF):
                tallyName = f"{obj.name}_{i}_{item.name}"
                tallies.createTallyPSF(f, tallyName, outputPrefix,
                                       det, item.emin, item.emax,
                                       item.gamma, item.electron, item.positron)

            # Angular Detector
            for i, item in enumerate(obj.penred_settings.talliesAngularDetector):
                tallyName = f"{obj.name}_{i}_{item.name}"
                tallies.createTallyAngularDet(f, tallyName, outputPrefix, det,
                                              item.emin, item.emax, item.ebins,
                                              item.theta1, item.theta2,
                                              item.phi1, item.phi2,
                                              item.logScale)

def createVR(obj, simBodyName, f):

    if hasattr(obj, "penred_settings"):
        f.write(f"\n## Variance Reduction for body {simBodyName} (Blender Name: {obj.name}) ##\n\n")
        # Bremsstrahlung splitting
        if obj.penred_settings.enableBremssSplitting:
            bremssName = f"{obj.name}_bremss_splitting"
            vr.createBremssSplitting(f, bremssName,
                                     obj.penred_settings.bremssSplitting, (simBodyName,))

        # X-Ray splitting
        if obj.penred_settings.enableXRaySplitting:
            xrayName = f"{obj.name}_xray_splitting"
            vr.createXRaySplitting(f, xrayName,
                                   ((simBodyName, obj.penred_settings.xraySplitting),))
            
        # Interaction forcing
        for i, item in enumerate(obj.penred_settings.interactionForcing):
            IFName = f"{obj.name}_{i}_{item.name}"

            interaction = 0
            if item.particleType == "electron":
                interaction = int(item.electronInteraction)
            elif item.particleType == "gamma":
                interaction = int(item.gammaInteraction)
            elif item.particleType == "positron":
                interaction = int(item.positronInteraction)
                
            vr.createInteractionForcing(f, IFName, item.particleType, interaction,
                                        item.factorType, item.factor,
                                        item.weightWindow, (simBodyName,))
