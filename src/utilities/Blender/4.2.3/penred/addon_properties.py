
import bpy
from . import tracks

materialRows = 20
materialColumns = 20
maxMaterial = materialRows*materialColumns

def update_source_direction(self, context):
    # Force the 3D viewport to redraw
    for area in context.screen.areas:
        if area.type == 'VIEW_3D':
            area.tag_redraw()

# Tallies properties
#####################

## Object dependent tallies

# Cylindrical dose distrib
class tallyCylDoseDistrib(bpy.types.PropertyGroup):
    name : bpy.props.StringProperty(name = "Tally Name", default = "Cyl-Dose")
    show : bpy.props.BoolProperty(name = "Show Tally Properties", default = True)
    nr : bpy.props.IntProperty(name = "Radial Bins", min = 1, default = 1,
                               description="Number of radial bins")
    nz : bpy.props.IntProperty(name = "Z Bins", min = 1, default = 1,
                               description="Number of Z bins")
    nPhi : bpy.props.IntProperty(name = "Angular Bins", min = 1, default = 1,
                                 description="Number of angular bins")

    spatialBBFit : bpy.props.BoolProperty(
        name = "Spatial Bounding Box Fit",
        description = "When enabled, the spatial tally mesh will be sized to fit inside the object bounding box. If disabled, the mesh size is enlarged to include the entire bounding box",
        default = False)

# Impact detector
class tallyImpactDetector(bpy.types.PropertyGroup):
    name : bpy.props.StringProperty(name = "Tally Name", default = "Impact-Detector")
    show : bpy.props.BoolProperty(name = "Show Tally Properties", default = True)
    emin : bpy.props.FloatProperty(
        name = "Minimum Energy",
        default = 1.0e3,
        min = 0.0,
        description="Minimum energy, in eV, to be detected",
        update=lambda self, context: setattr(self, "eminEdit", False))
    eminEdit : bpy.props.BoolProperty(name = "Minimum Energy Edit", default = False)
    
    emax : bpy.props.FloatProperty(
        name = "Maximum Energy",
        default = 1.0e5,
        min = 0.0,
        description="Maximum energy, in eV, to be detected",
        update=lambda self, context: setattr(self, "emaxEdit", False))
    emaxEdit : bpy.props.BoolProperty(name = "Maximum Energy Edit", default = False)
    
    ebins : bpy.props.IntProperty(name = "Energy Bins", min = 1, default = 200,
                                   description="Number of energy bins")
    showEBox : bpy.props.BoolProperty(name = "Show Energy Properties", default = True)
    fluence : bpy.props.BoolProperty(name = "Fluence", default = False,
                                   description="Enable/disable fluence detection")
    fluenceLogScale : bpy.props.BoolProperty(name = "Fluence Log Scale", default = False,
                                             description="Enable/disable fluence logarithmic scale")
    spectrum : bpy.props.BoolProperty(name = "Spectrum", default = True,
                                      description="Enable/disable spectrum detection")
    spectrumLogScale : bpy.props.BoolProperty(name = "Spectrum Log Scale", default = False,
                                              description="Enable/disable spectrum logarithmic scale")
    edep : bpy.props.BoolProperty(name = "Energy Deposition", default = False,
                                  description="Enable/disable energy deposition detection")
    edepLogScale : bpy.props.BoolProperty(name = "Energy Deposition Log Scale", default = False,
                                          description="Enable/disable energy deposition logarithmic scale")
    age : bpy.props.BoolProperty(name = "Age", default = False,
                                 description="Enable/disable age detection")
    ageMin : bpy.props.FloatProperty(name = "Minimum Age", default = 0.0, min = 0.0,
                                     description="Minimum age, in seconds, to be detected")
    ageMax : bpy.props.FloatProperty(name = "Maximum Age", default = 30.0, min = 0.0,
                                     description="Maximum age, in seconds, to be detected")
    ageBins : bpy.props.IntProperty(name = "Age Bins", min = 1, default = 200,
                                    description="Number of age bins")
    ageLogScale : bpy.props.BoolProperty(name = "Age Log Scale", default = False,
                                         description="Enable/disable age logarithmic scale")
    showAgeBox : bpy.props.BoolProperty(name = "Show Age Properties", default = True)

# Spatial dose distribution
class tallySpatialDoseDistrib(bpy.types.PropertyGroup):
    name : bpy.props.StringProperty(name = "Tally Name", default = "Spatial-Dose")
    show : bpy.props.BoolProperty(name = "Show Tally Properties", default = True)
    nx : bpy.props.IntProperty(name = "X Bins", min = 1, default = 1,
                               description="Number of X axis bins")
    ny : bpy.props.IntProperty(name = "Y Bins", min = 1, default = 1,
                               description="Number of Y axis bins")
    nz : bpy.props.IntProperty(name = "Z Bins", min = 1, default = 1,
                               description="Number of Z axis bins")

# Spatial dose distribution
class tallySphericalDoseDistrib(bpy.types.PropertyGroup):
    name : bpy.props.StringProperty(name = "Tally Name", default = "Spherical-Dose")
    show : bpy.props.BoolProperty(name = "Show Tally Properties", default = True)
    nr : bpy.props.IntProperty(name = "Radial Bins", min = 1, default = 1,
                               description="Number of radial bins")
    ntheta : bpy.props.IntProperty(name = "Polar Bins", min = 1, default = 1,
                                   description="Number of polar bins")
    nphi : bpy.props.IntProperty(name = "Azimuth Bins", min = 1, default = 1,
                                 description="Number of azimuth bins")

    spatialBBFit : bpy.props.BoolProperty(
        name = "Spatial Bounding Box Fit",
        description = "When enabled, the spatial tally mesh will be sized to fit inside the object bounding box. If disabled, the mesh size is enlarged to include the entire bounding box",
        default = False)

# Phase Space File
class tallyPSF(bpy.types.PropertyGroup):
    name : bpy.props.StringProperty(name = "Tally Name", default = "PSF")
    show : bpy.props.BoolProperty(name = "Show Tally Properties", default = True)
    emin : bpy.props.FloatProperty(
        name = "Minimum Energy",
        default = 1.0e3,
        min = 0.0,
        description="Minimum energy, in eV, to be recorded",
        update=lambda self, context: setattr(self, "eminEdit", False))
    eminEdit : bpy.props.BoolProperty(name = "Minimum Energy Edit", default = False)
    
    emax : bpy.props.FloatProperty(
        name = "Maximum Energy",
        default = 1.0e6,
        min = 0.0,
        description="Maximum energy, in eV, to be recorded",
        update=lambda self, context: setattr(self, "emaxEdit", False))
    emaxEdit : bpy.props.BoolProperty(name = "Maximum Energy Edit", default = False)
    
    gamma : bpy.props.BoolProperty(name = "Record Gammas", default = True,
                                   description="Enable/disable gamma recording")
    electron : bpy.props.BoolProperty(name = "Record Electrons", default = True,
                                      description="Enable/disable electron recording")
    positron : bpy.props.BoolProperty(name = "Record Positrons", default = True,
                                      description="Enable/disable positron recording")

# Kerma track length
class tallyKerma(bpy.types.PropertyGroup):
    name : bpy.props.StringProperty(name = "Tally Name", default = "Kerma")
    show : bpy.props.BoolProperty(name = "Show Tally Properties", default = True)
    emin : bpy.props.FloatProperty(
        name = "Minimum Energy",
        default = 1.0e3,
        min = 0.0,
        description="Gamma minimum energy, in eV, to be recorded",
        update=lambda self, context: setattr(self, "eminEdit", False))
    eminEdit : bpy.props.BoolProperty(name = "Minimum Energy Edit", default = False)
    
    emax : bpy.props.FloatProperty(
        name = "Maximum Energy",
        default = 1.0e6,
        min = 0.0,
        description="Gamma maximum energy, in eV, to be recorded",
        update=lambda self, context: setattr(self, "emaxEdit", False))
    emaxEdit : bpy.props.BoolProperty(name = "Maximum Energy Edit", default = False)

    nx : bpy.props.IntProperty(name = "X Bins", min = 1, default = 1,
                               description="Number of X axis bins")
    ny : bpy.props.IntProperty(name = "Y Bins", min = 1, default = 1,
                               description="Number of Y axis bins")
    nz : bpy.props.IntProperty(name = "Z Bins", min = 1, default = 1,
                               description="Number of Z axis bins")
    
    nr : bpy.props.IntProperty(name = "Radial Bins", min = 1, default = 1,
                               description="Number of radial bins")
    nphi : bpy.props.IntProperty(name = "Azimuth Bins", min = 1, default = 1,
                               description="Number of azimuth bins")
    ntheta : bpy.props.IntProperty(name = "Polar Bins", min = 1, default = 1,
                               description="Number of polar bins")
    
    meshType : bpy.props.EnumProperty(
        name = "Mesh type",
        description = "Choose the mesh type to score the kerma",
        items = [
            ("MESH_CART" , "Cartesian", "Cartesian mesh"),
            ("MESH_CYL", "Cylindrical", "Cylindrical mesh"),
            ("MESH_SPH", "Spherical", "Spherical mesh"),
        ],
        default = "MESH_CART"
    )

    dataPath : bpy.props.StringProperty(name = "Data Prefix", default = "",
                                        description = "Prefix of the mu-en data files."
                                        " If empty, the data files will be created in the simulation folder")

    spatialBBFit : bpy.props.BoolProperty(
        name = "Spatial Bounding Box Fit",
        description = "When enabled, the spatial tally mesh will be sized to fit inside the object bounding box. If disabled, the mesh size is enlarged to include the entire bounding box",
        default = False)    

# Spatial distribution
class tallySpatialDistrib(bpy.types.PropertyGroup):
    name : bpy.props.StringProperty(name = "Tally Name", default = "Spatial-Distrib")
    show : bpy.props.BoolProperty(name = "Show Tally Properties", default = True)

    emin : bpy.props.FloatProperty(
        name = "Minimum Energy",
        default = 1.0e3,
        min = 0.0,
        description="Minimum energy, in eV, to be detected",
        update=lambda self, context: setattr(self, "eminEdit", False))
    eminEdit : bpy.props.BoolProperty(name = "Minimum Energy Edit", default = False)
    
    emax : bpy.props.FloatProperty(
        name = "Maximum Energy",
        default = 1.0e6,
        min = 0.0,
        description="Maximum energy, in eV, to be detected",
        update=lambda self, context: setattr(self, "emaxEdit", False))
    emaxEdit : bpy.props.BoolProperty(name = "Maximum Energy Edit", default = False)

    ebins : bpy.props.IntProperty(name = "Energy Bins", min = 1, default = 1,
                               description="Number of energy bins")

    
    nx : bpy.props.IntProperty(name = "X Bins", min = 1, default = 1,
                               description="Number of X axis bins")
    ny : bpy.props.IntProperty(name = "Y Bins", min = 1, default = 1,
                               description="Number of Y axis bins")
    nz : bpy.props.IntProperty(name = "Z Bins", min = 1, default = 1,
                               description="Number of Z axis bins")

    printCoordinates : bpy.props.BoolProperty(name = "Print Coordinates", default = True,
                                              description="Enable/disable coordinates print in results")

    printBins : bpy.props.BoolProperty(name = "Print Bins", default = True,
                                       description="Enable/disable bin numbers print in results")
    
    particleType : bpy.props.EnumProperty(
        name = "Detected Particle",
        description = "Choose the particle to be detected",
        items = [
            ("PART_GAMMA" , "Gamma", "Gamma"),
            ("PART_ELECTRON", "Electron", "Electron"),
            ("PART_POSITRON", "Positron", "Positron"),
        ],
        default = "PART_GAMMA"
    )

# Spatial distribution
class tallyCT(bpy.types.PropertyGroup):
    name : bpy.props.StringProperty(name = "Tally Name", default = "CT Sinogram")
    show : bpy.props.BoolProperty(name = "Show Tally Properties", default = True)

    emin : bpy.props.FloatProperty(
        name = "Minimum Energy",
        default = 1.0e3,
        min = 0.0,
        description="Minimum energy, in eV, to be detected",
        update=lambda self, context: setattr(self, "eminEdit", False))
    eminEdit : bpy.props.BoolProperty(name = "Minimum Energy Edit", default = False)
    
    emax : bpy.props.FloatProperty(
        name = "Maximum Energy",
        default = 1.0e6,
        min = 0.0,
        description="Maximum energy, in eV, to be detected",
        update=lambda self, context: setattr(self, "emaxEdit", False))
    emaxEdit : bpy.props.BoolProperty(name = "Maximum Energy Edit", default = False)

    nPixels : bpy.props.IntProperty(name = "Pixels", min = 1, default = 10,
                                    description="Number of detector pixels")

    pixelDepth : bpy.props.FloatProperty(name = "Pixel Depth",
                                         default = 1.0,
                                         min = 0.001,
                                         description = "Pixel depth in cm")

    aperture : bpy.props.FloatProperty(name = "Aperture",
                                       default = 10.0,
                                       min = 0.01,
                                       description = "Angle subtended by the detector with respect to the isocenter, in degrees")

    scatter : bpy.props.BoolProperty(name = "Scatter",
                                     default = True,
                                     description = "Enables/disables collecting scattered particles")

    particleType : bpy.props.EnumProperty(
        name = "Particle Type",
        description = "Choose the particle to detect",
        items = [
            ("gamma" , "Gamma", "Gamma"),
            ("electron", "Electron", "Electron"),
            ("positron", "Positron", "Positron"),
        ],
        default = "gamma"
    )
    
# Angular Detector
class tallyAngularDetector(bpy.types.PropertyGroup):
    name : bpy.props.StringProperty(name = "Tally Name", default = "Angular-Detector")
    show : bpy.props.BoolProperty(name = "Show Tally Properties", default = True)

    emin : bpy.props.FloatProperty(
        name = "Minimum Energy",
        default = 1.0e3,
        min = 0.0,
        description="Minimum energy, in eV, to be detected",
        update=lambda self, context: setattr(self, "eminEdit", False))
    eminEdit : bpy.props.BoolProperty(name = "Minimum Energy Edit", default = False)
    
    emax : bpy.props.FloatProperty(
        name = "Maximum Energy",
        default = 1.0e6,
        min = 0.0,
        description="Maximum energy, in eV, to be detected",
        update=lambda self, context: setattr(self, "emaxEdit", False))
    emaxEdit : bpy.props.BoolProperty(name = "Maximum Energy Edit", default = False)
    
    ebins : bpy.props.IntProperty(name = "Energy bins", min = 1, default = 1,
                                 description="Number of energy bins")
    
    theta1 : bpy.props.FloatProperty(name = "Minimum Polar Angle", default = 0.0, min = 0.0, max=180.0,
                                   description="Minimum polar angle, in degrees, to be detected")

    theta2 : bpy.props.FloatProperty(name = "Maximum Polar Angle", default = 180.0, min = 0.0, max = 180.0,
                                   description="Maximum polar angle, in degrees, to be detected")

    phi1 : bpy.props.FloatProperty(name = "Minimum Azimuthal Angle", default = 0.0, min = 0.0, max=360.0,
                                   description="Minimum polar angle, in degrees, to be detected")

    phi2 : bpy.props.FloatProperty(name = "Maximum Azimuthal Angle", default = 360.0, min = 0.0, max = 360.0,
                                   description="Maximum polar angle, in degrees, to be detected")
    
    logScale : bpy.props.BoolProperty(name = "Log Scale", default = False,
                                      description="Enable/disable energy logarithmic scale")

## World dependent tallies

# Emergint particle distribution
class tallyEmergingParticle(bpy.types.PropertyGroup):
    name : bpy.props.StringProperty(name = "Tally Name", default = "Emerging-Part")
    show : bpy.props.BoolProperty(name = "Show Tally Properties", default = True)

    emin : bpy.props.FloatProperty(
        name = "Minimum Energy",
        default = 1.0e3,
        min = 0.0,
        description="Minimum energy, in eV, to be scored",
        update=lambda self, context: setattr(self, "eminEdit", False))
    eminEdit : bpy.props.BoolProperty(name = "Minimum Energy Edit", default = False)
    
    emax : bpy.props.FloatProperty(
        name = "Maximum Energy",
        default = 1.0e6,
        min = 0.0,
        description="Maximum energy, in eV, to be detected",
        update=lambda self, context: setattr(self, "emaxEdit", False))
    emaxEdit : bpy.props.BoolProperty(name = "Maximum Energy Edit", default = False)
    
    ebins : bpy.props.IntProperty(name = "Energy bins", min = 1, default = 1,
                                 description="Number of energy bins")    

    nTheta : bpy.props.IntProperty(name = "Polar Bins", min = 1, default = 1,
                                 description="Number of polar bins")
    
    nPhi : bpy.props.IntProperty(name = "Azimuthal Bins", min = 1, default = 1,
                                 description="Number of azimuthal bins")

# Tracks
class tallyTracks(bpy.types.PropertyGroup):
    name : bpy.props.StringProperty(name = "Tally Name", default = "Tracks")
    show : bpy.props.BoolProperty(name = "Show Tally Properties", default = False)
    enable : bpy.props.BoolProperty(name = "Enable Track", default = False)

    nHists : bpy.props.IntProperty(name = "Histories", min = 1, default = 10,
                                   description="Number of histories to be tracked")
    
    
talliesPropsClasses = (
    tallyCylDoseDistrib,
    tallyImpactDetector,
    tallySpatialDoseDistrib,
    tallySphericalDoseDistrib,
    tallyPSF,
    tallyKerma,
    tallySpatialDistrib,
    tallyCT,
    tallyAngularDetector,
    tallyEmergingParticle,
    tallyTracks,
)

# Source properties
####################
class sourceProperties(bpy.types.PropertyGroup):
    
    # Eanble/disable using the object as a particle source
    enabled : bpy.props.BoolProperty(name = "Source",
                                     description = "Flag this object as particle source",
                                     default = False)

    # CT source
    ctEnable : bpy.props.BoolProperty(
        name = "CT Like Source",
        description = "Enables CT like source oriented on the Z axis. This source type moves arround the object, according to the specified parameters, drawing a circle on the XY plane.",
        default = False)
    
    ctSecondaries : bpy.props.IntProperty(
        name = "CT Generated Secondaries",
        min = 1,
        description = "On CT like sources, if no PSF is used, this value specify the sampled particles per history",
        default = 10)

    ctRad : bpy.props.FloatProperty(
        name = "CT Radius",
        default = 10.0,
        min = 0.0,
        description= "CT like source radius, in cm"
    )

    ctPhiInterval : bpy.props.FloatVectorProperty(
        name="CT Angular Interval",
        description="2D Vector with the initial and final azimuthal angle for the CT rotation, in degrees. A 0º angle coincides with the X axis.",
        size=2,
        default=(0.0, 360),
        update=lambda self, context: self.__setitem__(
            "ctPhiInterval",
            (self.ctPhiInterval[0], max(self.ctPhiInterval[0],
                                        min(self.ctPhiInterval[1],
                                            self.ctPhiInterval[0] + 360.0)))
        )
    )
    
    ctNSteps : bpy.props.IntProperty(
        name = "CT Projections",
        default = 10,
        min = 1,
        max = 7200,
        description= "CT number of projections"
    )
    
    ctTStart : bpy.props.FloatProperty(
        name = "CT Start Time",
        default = 0.0,
        min = 0.0,
        description= "Starting time for the first CT like source projection, in seconds"
    )

    ctDT : bpy.props.FloatProperty(
        name = "CT Time Between Projections",
        default = 10.0,
        min = 0.0,
        description= "Time interval between projections of the CT like source, in seconds"
    )
    
    # Generic source parameters #
    particleType : bpy.props.EnumProperty(
        name = "Source Type",
        description = "Choose the source particle type",
        items = [
            ("PART_GAMMA" , "Gamma", "Gamma"),
            ("PART_ELECTRON", "Electron", "Electron"),
            ("PART_POSITRON", "Positron", "Positron"),
            ("PART_PSF", "PSF", "Phase Space File"),
        ],
        default = "PART_GAMMA"
    )
    sourcePSF : bpy.props.StringProperty(name = "Source Phase Space File",
                                         description = "Path to the source Phase Space File",
                                         default = "data.psf")
    psfMaxE : bpy.props.FloatProperty(
        name = "PSF Maximum Energy",
        description = "Maximum energy, in eV, of particles recorded in the PSF",
        default = 1.0e6,
        min=50.0,
        max=1.0e9,
        update=lambda self, context: setattr(self, "psfMaxEEdit", False))
    psfMaxEEdit : bpy.props.BoolProperty(name = "Energy Edit", default = False)

    split : bpy.props.IntProperty(name = "PSF Split",
                                  min = 1,
                                  description = "Splitting factor applied to each PSF particle",
                                  default = 10)

    psfWindow : bpy.props.FloatVectorProperty(
        name="PSF Weight Window",
        description="PSF weight window. Russian Roulette is applied to particles below the minimum weight. Splitting is applied to particles over the maximum weight",
        size=2,
        default=(0.1, 0.5),
        update=lambda self, context: self.__setitem__(
            "psfWindow",
            (self.psfWindow[0], max(self.psfWindow[0], self.psfWindow[1]))
        )
    )


    nHists : bpy.props.FloatProperty(
        name = "Histories",
        description = "Number of histories to simulate",
        min = 1.0,
        max = 1.0e19,
        default = 1.0e6,
        update=lambda self, context: setattr(self, "nHistsEdit", False))
    nHistsEdit : bpy.props.BoolProperty(name = "History Number Edit", default = False)

    # Spatial parameters
    spatialType : bpy.props.EnumProperty(
        name = "Spatial Type",
        description = "Choose the source energy type",
        items = [
            ("SPATIAL_POINT" , "Point", "Point source"),
            ("SPATIAL_BOX", "Box", "Sampled inside bounding box"),
            ("SPATIAL_CYL", "Cylinder", "Sampled inside a cylindrical shape oriented in Z axis"),
        ],
        default = "SPATIAL_BOX"
    )
    spatialBBFit : bpy.props.BoolProperty(
        name = "Spatial Bounding Box Fit",
        description = "When enabled, the spatial source will be sized to fit inside the object bounding box. If disabled, the source size is enlarged to include the entire bounding box",
        default = False)
    
    # Energy parameters
    energyType : bpy.props.EnumProperty(
        name = "Energy Type",
        description = "Choose the source energy type",
        items = [
            ("ENERGY_MONO" , "Monoenergetic", "Monoenergetic"),
            ("ENERGY_SFILE", "Spectrum", "Spectrum file"),
        ],
        default = "ENERGY_MONO"
    )
    energy : bpy.props.FloatProperty(
        name = "Source Energy",
        description = "Source energy in eV",
        default = 1.0e3,
        min=50.0,
        max=1.0e9,
        update=lambda self, context: setattr(self, "energyEdit", False))
    energyEdit : bpy.props.BoolProperty(name = "Energy Edit", default = False)
    
    spcFile : bpy.props.StringProperty(name = "Source Spectrum File",
                                       description = "Source energy spectrum filename",
                                       default = "spectrum.spc")

    # Direction parameters
    aperture: bpy.props.FloatProperty(name = "Source Aperture",
                                      description = "Source aperture in degrees",
                                      default = 0.0,
                                      min=0.0,
                                      max=180.0)
    direction: bpy.props.FloatVectorProperty(name = "Direction",
                                             description = "Source direction vector",
                                             size = 3,
                                             default = (0.0,0.0,1.0),
                                             update=update_source_direction)

    # Time parameters
    timeRecord : bpy.props.BoolProperty(name = "Source Time Record",
                                        description = "Toggle record time during simulation",
                                        default = False)
    timeType : bpy.props.EnumProperty(
        name = "Time Type",
        description = "Choose the source time initialization type",
        items = [
            ("TIME_NOINIT" , "No initialization", "No initialization"),
            ("TIME_DECAY", "Decay", "Exponential decay"),
        ],
        default = "TIME_NOINIT"
    )
    decayHalf: bpy.props.FloatProperty(name = "Source Half Life",
                                       description = "Decay half life, in seconds",
                                       default = 1.0,
                                       min=0.0)

    timeWindow : bpy.props.FloatVectorProperty(
        name="Time Sampling Window",
        description="Particle time will be sampled inside this time window, specified in seconds",
        size=2,
        default=(0.0, 30.0),
        update=lambda self, context: self.__setitem__(
            "timeWindow",
            (self.timeWindow[0], max(self.timeWindow[0], self.timeWindow[1]))
        )
    )

    enableSourceMat : bpy.props.BoolProperty(
        name = "Enable Source Material",
        description = "Enable/disable restricted spatial sampling to a specified material",
        default = False)
    sourceMat : bpy.props.IntProperty(
        name = "Source Material",
        min = 1,
        description = "Discards particles sampled outside the specified material",
        default = 1)

# VR properties
####################    
class interactionForcingProperties(bpy.types.PropertyGroup):
    
    name : bpy.props.StringProperty(name = "Interaction Forcing Name", default = "IF")
    show : bpy.props.BoolProperty(name = "Show Interaction Forcing Properties",
                                  default = False)
    particleType : bpy.props.EnumProperty(
        name = "Particle Type",
        description = "Choose the particle to apply interaction forcing",
        items = [
            ("gamma" , "Gamma", "Gamma"),
            ("electron", "Electron", "Electron"),
            ("positron", "Positron", "Positron"),
        ],
        default = "electron"
    )

    electronInteraction : bpy.props.EnumProperty(
        name = "Electron Interaction",
        description = "Choose the electron interaction to force",
        items = [
            ("1" , "Inelastic", "Inelastic"),
            ("2", "Bremsstrahlung", "Bremsstrahlung"),
            ("3", "Inner Shell", "Inner Shell"),
        ],
        default = "2"
    )
    
    gammaInteraction : bpy.props.EnumProperty(
        name = "Gamma Interactions",
        description = "Choose the gamma interaction to force",
        items = [
            ("1", "Compton", "Compton"),
            ("2", "Photoelectric", "Photoelectric"),
            ("3", "Pair Production", "Pair Production"),
        ],
        default = "2"
    )

    positronInteraction : bpy.props.EnumProperty(
        name = "Positron Interaction",
        description = "Choose the positron interaction to force",
        items = [
            ("1" , "Inelastic", "Inelastic"),
            ("2", "Bremsstrahlung", "Bremsstrahlung"),
            ("3", "Inner Shell", "Inner Shell"),
        ],
        default = "2"
    )

    factorType : bpy.props.EnumProperty(
        name = "Factor Type",
        description = "Choose how to specify the amplification factor",
        items = [
            ("MULTIPLY" , "Multiply", "The interaction probability is multiplied by the factor"),
            ("AVERAGE" , "Average", "The factoris interpreted as the average number of interactions, of a particle with the maximum available energy, until rest (electrons and positrons) or across a mean free path (gammas)."),
        ],
        default = "MULTIPLY"
    )

    factor : bpy.props.IntProperty(
        name = "Factor",
        description="Defines the amplification factor for interaction forcing.",
        default=10,
        min = 1,
    )

    weightWindow : bpy.props.FloatVectorProperty(
        name="Weight Window",
        description="This interaction forcing will be applied only to those particles with a weight value inside the specified window.",
        size=2,
        default=(0.2, 2.0),
        update=lambda self, context: self.__setitem__(
            "weightWindow",
            (self.weightWindow[0], max(self.weightWindow[0], self.weightWindow[1]))
        )
    )
    
VRPropsClasses = (
    interactionForcingProperties,
)
    
# Object Properties group
###########################
class objectProperties(bpy.types.PropertyGroup):

    ## Generic object properties ##
    isMaterialObject : bpy.props.BoolProperty(
        name = "Flags if the object is a material object",
        description = "Sets this object as a detector",
        default = True)
    
    material : bpy.props.IntProperty(name = "Material Index",
                                     min = 0,
                                     max = maxMaterial,
                                     description = "Material index assigned to the object",
                                     default = 1)

    isDetector : bpy.props.BoolProperty(name = "Detector Toggle",
                                        description = "Sets this object as a detector",
                                        default = False)    
    detector : bpy.props.IntProperty(name = "Detector Index",
                                     min = 1,
                                     description = "Detector index assigned to the object",
                                     default = 1)

    dsmaxEnabled: bpy.props.BoolProperty(
        name = "Maximum Class II Distance Toggle",
        description=
        "Enable/disable limiting the distance between hard interactions in Class II scheme (dsmax)",
        default = False)
    dsmax : bpy.props.FloatProperty(
        name = "Maximum Class II Distance",
        default = 1.0,
        min = 0.0,
        description=
        "Limiting distance, in cm, for Class II scheme",
        update=lambda self, context: setattr(self, "dsmaxEdit", False)
    )
    dsmaxEdit: bpy.props.BoolProperty(
        name = "Maximum Class II Distance Edit",
        description=
        "Limiting distance, in cm, for Class II scheme",
        default = False)

    showInteractionForcing : bpy.props.BoolProperty(name = "Show Interaction Forcing",
                                                    default = False)
    interactionForcing: bpy.props.CollectionProperty(type=interactionForcingProperties)

    enableBremssSplitting : bpy.props.BoolProperty(
        name = "Bremsstrahlung Splitting Toggle",
        description = "Enables/disables splitting of generated Bremsstrahlung photons",
        default = False
    )
    bremssSplitting : bpy.props.IntProperty(
        name = "Bremsstrahlung Splitting",
        min = 2,
        description = "Times a Bremsstrahlung photon is splitted",
        default = 2
    )
    
    enableXRaySplitting : bpy.props.BoolProperty(
        name = "X-Ray Splitting Toggle",
        description = "Enables/disables splitting of generated x-ray",
        default = False
    )
    xraySplitting : bpy.props.IntProperty(
        name = "X-Ray Splitting",
        min = 2,
        description = "Times a x-ray is splitted",
        default = 2
    )
    
    ## Quadric object properties ##
    quadricType : bpy.props.StringProperty(name = "Quadric Type", default = "unknown")

    module : bpy.props.BoolProperty(name = "Module", default = False)
    r1 : bpy.props.FloatProperty(name = "r1", default = 0.0,
                                 min = 0.0,
                                 update=lambda self, context: (bpy.ops.remesh.remesh_operator(), None)[1])
    r2 : bpy.props.FloatProperty(name = "r2", default = 0.0,
                                 min = 0.0,
                                 update=lambda self, context: (bpy.ops.remesh.remesh_operator(), None)[1])

    topSize : bpy.props.FloatVectorProperty(
        name = "Top Size",
        description = "Top size",
        size = 2,
        default = (1.0,1.0),
        min = 0.0,
        update=lambda self, context: (bpy.ops.remesh.remesh_operator(), None)[1])
    botSize : bpy.props.FloatVectorProperty(
        name = "Bot Size",
        size = 2,
        default = (2.0, 2.0),
        min = 0.0,
        update=lambda self, context: (bpy.ops.remesh.remesh_operator(), None)[1])

    ## Source parameters ##
    source: bpy.props.PointerProperty(type=sourceProperties)
    
    ## Tally parameters ##
    showTalliesCylDose : bpy.props.BoolProperty(name = "Show Cylindrical Dose Tallies", default = True)
    talliesCylDose: bpy.props.CollectionProperty(type=tallyCylDoseDistrib)
    
    showTalliesImpactDet : bpy.props.BoolProperty(name = "Show Impact Detector Tallies", default = True)
    talliesImpactDetector : bpy.props.CollectionProperty(type=tallyImpactDetector)
    
    showTalliesSpatialDose : bpy.props.BoolProperty(name = "Show Spatial Dose Tallies", default = True)
    talliesSpatialDoseDistrib : bpy.props.CollectionProperty(type=tallySpatialDoseDistrib)
    
    showTalliesSphericalDose : bpy.props.BoolProperty(name = "Show Spherical Dose Tallies", default = True)
    talliesSphericalDoseDistrib : bpy.props.CollectionProperty(type=tallySphericalDoseDistrib)
    
    showTalliesPSF : bpy.props.BoolProperty(name = "Show PSF Tallies", default = True)
    talliesPSF : bpy.props.CollectionProperty(type=tallyPSF)
    
    showTalliesKerma : bpy.props.BoolProperty(name = "Show Kerma Tallies", default = True)
    talliesKerma : bpy.props.CollectionProperty(type=tallyKerma)
    
    showTalliesSpatialDistrib : bpy.props.BoolProperty(name="Show Spatial Distrib Tallies", default = True)
    talliesSpatialDistrib : bpy.props.CollectionProperty(type=tallySpatialDistrib)
    
    showTalliesAngularDet : bpy.props.BoolProperty(name = "Show Angular Detector Tallies", default = True)
    talliesAngularDetector : bpy.props.CollectionProperty(type=tallyAngularDetector)

    showTalliesCT : bpy.props.BoolProperty(name = "Show CT Tallies", default = True)
    talliesCT : bpy.props.CollectionProperty(type=tallyCT)
    
# Material properties groups
#############################
class elementProperties(bpy.types.PropertyGroup):
    z : bpy.props.IntProperty(name = "Z", min = 1, max = 99, default = 1)
    wFraction : bpy.props.FloatProperty(name = "Weight Factor", min = 0.0, default = 1.0)

class materialProperties(bpy.types.PropertyGroup):
    name : bpy.props.StringProperty(name = "Material Name", default = "material")
    show : bpy.props.BoolProperty(name = "Show Material Properties", default = True)
    index : bpy.props.IntProperty(name = "Index", min = 1, default = 1)
    composition : bpy.props.CollectionProperty(type=elementProperties)
    density : bpy.props.FloatProperty(
        name = "Density",
        default = 1.0,
        min = 1.0e-10,
        description="Material density in g/cm^3"
    )
    
    # Gamma
    gammaCutoffType : bpy.props.EnumProperty(
        name = "Gamma Cutoff Type",
        description = "Choose the cutoff type for gammas",
        items = [
            ("CUTOFF_EABS" , "Energy", "Energy absorption"),
            ("CUTOFF_RANGE", "Range", "Minimum range"),
        ],
        default = "CUTOFF_EABS"
    )
    gammaEabs : bpy.props.FloatProperty(
        name = "Gamma Energy Absorption",
        default = 1.0e3,
        min = 50.0,
        description="Absorption energy, in eV, for gammas",
        update=lambda self, context: setattr(self, "gammaEabsEdit", False)
    )
    gammaEabsEdit : bpy.props.BoolProperty(
        name = "Gamma Absorption Energy Edit",
        description="Absorption energy, in eV, for gammas",
        default = False)
    gammaRange : bpy.props.FloatProperty(
        name = "Gamma Minimum Range",
        default = 1.0,
        min = 0.0,
        description="Minimum range in terms of mean free path, in cm, for gammas"
    )

    # Electrons
    electronCutoffType : bpy.props.EnumProperty(
        name = "Electron Cutoff Type",
        description = "Choose the cutoff type for electrons",
        items = [
            ("CUTOFF_EABS" , "Energy", "Energy absorption"),
            ("CUTOFF_RANGE", "Range", "Minimum range"),
        ],
        default = "CUTOFF_EABS"
    )    
    electronEabs : bpy.props.FloatProperty(
        name = "Electron Energy Absorption",
        default = 1.0e3,
        min = 50.0,
        description="Absorption energy, in eV, for electrons",
        update=lambda self, context: setattr(self, "electronEabsEdit", False)
    )
    electronEabsEdit : bpy.props.BoolProperty(
        name = "Electron Absorption Energy Edit",
        description="Absorption energy, in eV, for electrons",
        default = False)
    electronRange : bpy.props.FloatProperty(
        name = "Electron Minimum Range",
        default = 1.0,
        min = 0.0,
        description=
        "Minimum range in terms of Continuous Slowing Down Approximation (CSDA), in cm, for electrons"
    )

    # Positrons
    positronCutoffType : bpy.props.EnumProperty(
        name = "Positron Cutoff Type",
        description = "Choose the cutoff type for positrons",
        items = [
            ("CUTOFF_EABS" , "Energy", "Energy absorption"),
            ("CUTOFF_RANGE", "Range", "Minimum range"),
        ],
        default = "CUTOFF_EABS"
    )
    positronEabs : bpy.props.FloatProperty(
        name = "Positron Energy Absorption",
        default = 1.0e3,
        min = 50.0,
        description="Absorption energy, in eV, for positrons",
        update=lambda self, context: setattr(self, "positronEabsEdit", False)
    )
    positronEabsEdit : bpy.props.BoolProperty(
        name = "Positron Absorption Energy Edit",
        description="Absorption energy, in eV, for positrons",        
        default = False)
    positronRange : bpy.props.FloatProperty(
        name = "Positron Minimum Range",
        default = 1.0,
        min = 0.0,
        description=
        "Minimum range in terms of Continuous Slowing Down Approximation (CSDA), in cm, for positrons"
    )

    # Advanced parameters
    enableAdvanced : bpy.props.BoolProperty(
        name = "Advanced Parameters Toggle",
        description = "Enable/disable material advanced parameters",
        default = False)

    showAdvanced : bpy.props.BoolProperty(
        name = "Show Advanced Properties",
        default = False)
    
    C1 : bpy.props.FloatProperty(
        name = "Average Angular Deflection",
        default = 0.05,
        min = 0.0,
        max = 0.2,
        description=
        "Average angular deflection produced by multiple elastic scattering along a path length equal to the mean free path between consecutive hard elastic events."
    )
        
    C2 : bpy.props.FloatProperty(
        name = "Maximum Average Fractional Energy Loss",
        default = 0.05,
        min = 0.0,
        max = 0.2,
        description=
        "Maximum average fractional energy loss between consecutive hard elastic events."
    )

    WCC : bpy.props.FloatProperty(
        name = "Inelastic Cutoff Energy Loss",
        default = 1.0e3,
        min = 0.0,
        max = 1.0e9,
        description=
        "Cutoff energy loss, in eV, for hard inelastic collisions.",
        update=lambda self, context: setattr(self, "WCCEdit", False)
    )
    WCCEdit: bpy.props.BoolProperty(
        name = "Inelastic Cutoff Energy Loss Edit",
        description=
        "Cutoff energy loss, in eV, for hard inelastic collisions.",        
        default = False)
        
    WCR : bpy.props.FloatProperty(
        name = "Bremsstrahlung Cutoff Energy Loss",
        default = 1.0e3,
        min = 0.0,
        max = 1.0e9,
        description=
        "Cutoff energy loss, in eV, for hard Bremsstrahlung collisions.",
        update=lambda self, context: setattr(self, "WCREdit", False)
    )
    WCREdit: bpy.props.BoolProperty(
        name = "Bremsstrahlung Cutoff Energy Loss Edit",
        description=
        "Cutoff energy loss, in eV, for hard Bremsstrahlung collisions",
        default = False)


# Simulation properties groups
#############################
class simulationProperties(bpy.types.PropertyGroup):
    
    enableDumps : bpy.props.BoolProperty(
        name = "Enable Dumps",
        description = "Enable/Disable write partial results dumps",
        default = False)
    dumpInterval : bpy.props.FloatProperty(
        name = "Dump Interval",
        default = 3600,
        min = 0.0,
        description=
        "Interval, in seconds, between results dump."
    )
    dumpWriteFile : bpy.props.StringProperty(
        name = "Write Dump Files",
        description = "Filename to write dumps. For each thread, a 'thN' prefix will be added, where N is the thread number",
        default = "dump.dat")

    readDumps : bpy.props.BoolProperty(
        name = "Read Dumps",
        description = "Enable/Disable resuming simulation from dump files",
        default = False)
    dumpReadFile : bpy.props.StringProperty(
        name = "Read Dump Files",
        description = "Filename to write dumps. For each thread, a 'thN' prefix will be added, where N is the thread number",
        default = "dump.dat")

    finalDump : bpy.props.BoolProperty(
        name = "Final Dumps",
        description = "Enable/Disable creating dumps at the simulation end",
        default = False)

    asciiResults : bpy.props.BoolProperty(
        name = "ASCII Results",
        description = "Enable/Disable creating results in ASCII format. If disabled, results will be writen to 'results.dump' in dump format",
        default = True)

    partialResults : bpy.props.BoolProperty(
        name = "ASCII Partial Results",
        description = "Enable/Disable creating partial results in ASCII format",
        default = False)
    
    
    outputPrefix : bpy.props.StringProperty(name = "Output Prefix",
                                         description = "Prefix for results files",
                                         default = "")
    

    # Threads
    threadSelType : bpy.props.EnumProperty(
        name = "Thread Selection Type",
        description = "Choose how the number of threads is selected",
        items = [
            ("AUTO" , "Auto", "All available threads"),
            ("MANUAL", "Manual", "Manual selected"),
        ],
        default = "AUTO"
    )
    nThreads : bpy.props.IntProperty(
        name = "Threads Number",
        min = 1,
        max = 1000,
        description = "Number of threads to be used during simulation",
        default = 4)
    seedPair : bpy.props.IntProperty(name = "Seed Pair",
                                     min = 0,
                                     max = 1000,
                                     description = "Initial seed pair to initialize random number generator",
                                    default = 0)

    limitSimTime : bpy.props.BoolProperty(
        name = "Toggle Maximum Simulation Time",
        description = "Enable/Disable the simulation time limit",
        default = False)
    maxSimTime : bpy.props.FloatProperty(
        name = "Maximum Simulation Time",
        default = 3600,
        min = 10.0,
        description=
        "Maximum simulation time in seconds"
    )

# Track properties group
#############################

def updateTrackERange(self, context):

    if self.trackERangeEdit:
        self.trackERangeEdit = False

    tm = tracks.getTrackManager()
    
    if self.trackAutoRange:
        # Get new ranges automaticatly
        if len(self.trackFiles) > 0:
            minE = min([tf.eRange[0] for tf in self.trackFiles])
            maxE = max([tf.eRange[1] for tf in self.trackFiles])
            
            # Update ranges if needed
            toUpdate = tm.setERange((minE, maxE))
            if toUpdate:
                tm.updateTracks(self.trackFiles)
                
    else:
        #Update ranges according to the property value
        toUpdate = tm.setERange(self.trackERange)
        if toUpdate:
            tm.updateTracks(self.trackFiles)

    # Update property values if needed
    if (self.trackERange[0] != tm.energyRange[0] or
        self.trackERange[1] != tm.energyRange[1]):
        self.trackERange = tm.energyRange

def updateTrackTRange(self, context):

    if self.trackTRangeEdit:
        self.trackTRangeEdit = False

    tm = tracks.getTrackManager()        
    
    if self.trackAutoRange:
        # Get new ranges automaticatly
        if len(self.trackFiles) > 0:
            minT = min([tf.tRange[0] for tf in self.trackFiles])
            maxT = max([tf.tRange[1] for tf in self.trackFiles])
            
            # Update ranges if needed
            toUpdate = tm.setTRange((minT, maxT))
            if toUpdate:
                tm.updateTracks(self.trackFiles)
                
    else:
        
        #Update ranges according to the property value
        toUpdate = tm.setTRange(self.trackTRange)
        if toUpdate:
            tm.updateTracks(self.trackFiles)

    # Update property values if needed
    if (self.trackTRange[0] != tm.timeRange[0] or
        self.trackTRange[1] != tm.timeRange[1]):
        self.trackTRange = tm.timeRange

def updateTrackRanges(self, context):
    updateTrackERange(self, context)
    updateTrackTRange(self, context)

def updateTracks(self, context):
    
    tm = tracks.getTrackManager()

    # Check if the logscale has been changed
    toUpdate = tm.setELogScale(self.trackERangeLog) 
    # Check if the out of range has been changed
    toUpdate = tm.setOutOfRange(self.trackDrawOutOfRange) or toUpdate

    if toUpdate:
        tm.updateTracks(self.trackFiles)

# Frame change handler
@bpy.app.handlers.persistent
def updateTracksFrame(scene):
    if hasattr(scene,"penred_settings"):
        updateTracks(scene.penred_settings, None)
        updateTrackRanges(scene.penred_settings, None)
    
class trackPointProperties(bpy.types.PropertyGroup):
    
    position : bpy.props.FloatVectorProperty(
        name = "Position",
        description = "Point position, in cm",
        size = 3,
        default = (0.0, 0.0, 0.0)
    )

    energy : bpy.props.FloatProperty(
        name = "Energy",
        default = 1.0,
        min = 0.0,
        description="Track energy at this point, in keV"
    )

    time : bpy.props.FloatProperty(
        name = "Time",
        default = 1.0,
        min = 0.0,
        description="Track time at this point, in seconds"
    )
    
    weight : bpy.props.FloatProperty(
        name = "Weight",
        default = 1.0,
        min = 0.0,
        description="Track weight at this point"
    )
    
class trackProperties(bpy.types.PropertyGroup):

    points: bpy.props.CollectionProperty(type=trackPointProperties)

    ILB : bpy.props.IntVectorProperty(
        name = "ILB",
        description = "ILB values",
        size = 5,
        default = (0, 0, 0, 0, 0)
    )

class trackFileProperties(bpy.types.PropertyGroup):

    tracks : bpy.props.CollectionProperty(type=trackProperties)
    filename : bpy.props.StringProperty(name = "Tracks File",
                                        description = "Tracks Filename",
                                        default = "")
    enabled : bpy.props.BoolProperty(
        name = "Enable Drawing",
        description = "Enable/Disable tracks drawing",
        default = False)

    eRange : bpy.props.FloatVectorProperty(
        name = "Tracks Energy Range",
        size = 2,
        default = (0.0, 1.0e6),
        description = "Energy range for this track file, in eV"
    )
    
    tRange : bpy.props.FloatVectorProperty(
        name = "Tracks Time Range",
        size = 2,
        default = (0.0, 30.0)
    )
    
tracksClasses = (
    trackPointProperties,
    trackProperties,
    trackFileProperties,
)

# World properties group
#############################
class worldProperties(bpy.types.PropertyGroup):

    ## Materials parameters ##
    materials : bpy.props.CollectionProperty(type=materialProperties)

    ## Simulation ##
    simulation : bpy.props.PointerProperty(type=simulationProperties)
    
    ## Tallies ##
    showTalliesEmergingPart : bpy.props.BoolProperty(name="Show Emerging Particles Tallies", default = True)
    talliesEmergingParticle : bpy.props.CollectionProperty(type=tallyEmergingParticle)
    tracksTally : bpy.props.PointerProperty(type=tallyTracks)

# Scene properties group
#############################
class penredSceneProperties(bpy.types.PropertyGroup):

    simulationState : bpy.props.EnumProperty(
        name = "Simulation Status",
        description = "Saves the current simulation status",
        items = [
            ("NONE" , "None", "No simulation running or on configuration"),
            ("EXPORTING" , "Exporting", "The simulation is exporting files"),
            ("EXPORTED" , "Exported", "The simulation files have been exported"),
            ("RUNNING", "Running", "Simulation is running"),
            ("CANCELLED", "Running", "Simulation cancelled"),
        ],
        default = "NONE"
    )

    simulationConfigPath : bpy.props.StringProperty(
        name = "Simulation Path",
        description = "Path to the simulation folder",
        default = "")

    # Tracks
    trackFiles : bpy.props.CollectionProperty(type=trackFileProperties)
    trackDrawOutOfRange : bpy.props.BoolProperty(
        name="Draw Out Of Range",
        description="Enable/disable drawing out of range track points",
        default = True,
        update = updateTracks
    )
    trackColorType : bpy.props.EnumProperty(
        name = "Track Color Type Selection",
        description = "Choose how color is assigned for tracks",
        items = [
            ("ENERGY" , "Energy", "By energy"),
            ("TIME", "Time", "By time"),
        ],
        default = "ENERGY"
    )
    
    trackAutoRange : bpy.props.BoolProperty(
        name = "Track Automatic Range",
        default = True,
        description = "Enable/disable automatic range selection for tracks",
        update = updateTrackRanges
    )
    
    trackERange : bpy.props.FloatVectorProperty(
        name = "Track Energy Range",
        size = 2,
        default = (0.0, 1.0e6),
        description = "Energy range for track drawing, in eV",
        update = updateTrackERange
    )
    trackERangeEdit : bpy.props.BoolProperty(name = "Track Energy Range Edit",
                                             default = False)
    trackERangeLog : bpy.props.BoolProperty(
        name = "Track Energy Range Logscale",
        default = False,
        description = "Enable/disable logscale for track energy range",
        update = updateTracks
    )

    trackTRange : bpy.props.FloatVectorProperty(
        name = "Track Time Range",
        size = 2,
        default = (0.0, 30.0),
        description = "Time range for track drawing, in seconds",
        update=updateTrackTRange
    )
    trackTRangeEdit : bpy.props.BoolProperty(name = "Track Time Range Edit",
                                             default = False)
def register():

    # Register source properties
    bpy.utils.register_class(sourceProperties)
    
    # Register tally properties
    for cls in talliesPropsClasses:
        bpy.utils.register_class(cls)

    # Register VR properties
    for cls in VRPropsClasses:
        bpy.utils.register_class(cls)

    #Register tracks property groups
    for cls in tracksClasses:
        bpy.utils.register_class(cls)        
        
    # Register globa properties
    bpy.utils.register_class(objectProperties)

    # Register element properties
    bpy.utils.register_class(elementProperties)

    # Register material properties
    bpy.utils.register_class(materialProperties)

    # Register simulation properties
    bpy.utils.register_class(simulationProperties)    
    
    # Register world properties
    bpy.utils.register_class(worldProperties)

    # Register scene properties
    bpy.utils.register_class(penredSceneProperties)

    # Add properties to objects
    bpy.types.Object.penred_settings = bpy.props.PointerProperty(type=objectProperties)

    # Add properties to world
    bpy.types.World.penred_settings = bpy.props.PointerProperty(type=worldProperties)

    # Add properties to scene
    bpy.types.Scene.penred_settings = bpy.props.PointerProperty(type=penredSceneProperties)

    # Register frame update for tracks
    bpy.app.handlers.frame_change_pre.append(updateTracksFrame)
        
def unregister():

    # Unregister source properties
    bpy.utils.unregister_class(sourceProperties)
    
    # Unregister tally properties
    for cls in talliesPropsClasses:
        bpy.utils.unregister_class(cls)

    # Unregister VR properties
    for cls in VRPropsClasses:
        bpy.utils.unregister_class(cls)        

    #Unregister tracks property groups
    for cls in tracksClasses:
        bpy.utils.unregister_class(cls)
        
    # Unregister element properties
    bpy.utils.unregister_class(elementProperties)

    # Unregister material properties
    bpy.utils.unregister_class(materialProperties)
    
    # Unregister simulation properties
    bpy.utils.unregister_class(simulationProperties)    
    
    del bpy.types.Object.penred_settings
    bpy.utils.unregister_class(objectProperties)    

    # Delete world penRed settings
    del bpy.types.World.penred_settings

    # Unregister world properties
    bpy.utils.unregister_class(worldProperties)

    # Delete scene penRed settings
    del bpy.types.Scene.penred_settings

    # Unregister frame update for tracks
    bpy.app.handlers.frame_change_pre.remove(updateTracksFrame)    
    
    # Unregister scene properties
    bpy.utils.unregister_class(penredSceneProperties)    
