
import bpy
from mathutils import Vector
from mathutils import Color
from . import addon_properties, utils
from math import pi, cos, sin

# Export menu function
def menu_func_penred_export(self, context):
    self.layout.operator("export_penred.data", text="PenRed")

# Object transforms functions
def menu_func_quadric_transform(self, context):
    self.layout.operator("transform.quadric_transform_putdown", text="Put Down")
    self.layout.operator("transform.quadric_transform_puttop", text="Put Top")
    self.layout.operator("transform.quadric_transform_putinsidetop", text="Put Inside Top")
    self.layout.operator("transform.quadric_transform_putinsidedown", text="Put Inside Down")
    self.layout.operator("transform.quadric_transform_putinsidecentered", text="Put Inside Centered")
    
# View Submenu
#####################
class OBJECT_MT_view_quadric_submenu(bpy.types.Menu):
    bl_idname = "OBJECT_MT_view_quadric_submenu"
    bl_label = "Quadric view"

    def draw(self, context):
        layout = self.layout

        layout.operator("view.quadric_material",text="Material")
        layout.operator("view.quadric_module",text="Module")

def menu_quadric_view_func(self, context):
    self.layout.menu("OBJECT_MT_view_quadric_submenu", text="Quadric")

# Add Quadrics Submenu class
##############################
class OBJECT_MT_quadric_submenu(bpy.types.Menu):
    bl_idname = "OBJECT_MT_quadric_submenu"
    bl_label = "Quadric"
    
    def draw(self, context):
        layout = self.layout
        
        layout.operator("mesh.add_cube_quadric",text="Cube",icon='CUBE')
        layout.operator("mesh.add_sphere_quadric",text="Sphere",icon='SPHERE')
        layout.operator("mesh.add_cone_quadric",text="Cone",icon='CONE')
        layout.operator("mesh.add_cylinder_quadric",text="Cylinder",icon='MESH_CYLINDER')
        layout.operator("mesh.add_plane_quadric",text="Plane",icon='MESH_PLANE')
        layout.operator("mesh.add_semisphere_quadric",text="Semi-Sphere",icon='SPHERE')

def menu_quadric_add_func(self, context):
    self.layout.menu("OBJECT_MT_quadric_submenu", text="Quadric", icon='CUBE')

# Generic body properties menu
##########################
class PenredBodyPropertiesPanel(bpy.types.Panel):
    bl_label = "Body properties"
    bl_idname = "OBJECT_PT_body"
    bl_space_type = 'PROPERTIES'
    bl_region_type = 'WINDOW'
    bl_context = "object"

    def draw(self, context):
        layout = self.layout
        obj = context.object

        if obj and obj.penred_settings and obj.type == "MESH" and obj.penred_settings.isMaterialObject:

            # Material
            row = layout.row()
            row.prop(obj.penred_settings, "material", text="Material Index")

            # Detector
            box = layout.box()
            box.label(text="Detector")
            row = box.row()
            row.prop(obj.penred_settings, "isDetector", text="Enabled")

            if obj.penred_settings.isDetector:
                row = box.row()
                row.prop(obj.penred_settings, "detector", text="Index")

            box = layout.box()
            box.label(text="Advanced Parameters")
            row = box.row()
            row.prop(obj.penred_settings, "dsmaxEnabled", text="Class II Limit Distance")
            if obj.penred_settings.dsmaxEnabled:
                row = box.row()
                if obj.penred_settings.dsmaxEdit:
                    row.prop(obj.penred_settings, "dsmax", text="")
                else:
                    row.label(text=f"DSMax :  {obj.penred_settings.dsmax:.3e} cm")
                    row.prop(obj.penred_settings, "dsmaxEdit",
                             text="", icon="GREASEPENCIL",
                             toggle=True)

            # X-Ray splitting
            ############################
            vrbox = box.box()
            vrbox.label(text="X-Ray Splitting")

            row = vrbox.row()
            row.prop(obj.penred_settings, "enableXRaySplitting", text="Enable")

            row = vrbox.row()
            row.enabled = obj.penred_settings.enableXRaySplitting
            row.prop(obj.penred_settings, "xraySplitting", text="Factor")

            # Bremss splitting
            ############################
            vrbox = box.box()
            vrbox.label(text="Bremsstrahlung Splitting")

            row = vrbox.row()
            row.prop(obj.penred_settings, "enableBremssSplitting", text="Enable")

            row = vrbox.row()
            row.enabled = obj.penred_settings.enableBremssSplitting
            row.prop(obj.penred_settings, "bremssSplitting", text="Factor")            
                    
            # Interaction Forcing
            ############################
            vrbox = box.box()
            row = vrbox.row()
            row.prop(obj.penred_settings, "showInteractionForcing",
                     text="Interaction Forcing", emboss=True,
                     icon="TRIA_DOWN" if obj.penred_settings.showInteractionForcing else "TRIA_RIGHT")
            
            row.operator("vr_if.add_item", text="", icon="ADD")
                    
            if obj.penred_settings.showInteractionForcing:
                for i, item in enumerate(obj.penred_settings.interactionForcing):

                    box = vrbox.box()
                    row = box.row()
                    row.prop(item, "show", text=item.name, emboss=True,
                             icon="TRIA_DOWN" if item.show else "TRIA_RIGHT")
                    row.operator("vr_if.remove_item", text="", icon="TRASH").index = i
                    if not item.show:
                        continue

                    # Name
                    row = box.row()
                    row.prop(item, "name", text="Name")

                    # Particle type
                    row = box.row()
                    row.prop(item, "particleType", text="Particle")

                    # Interaction
                    if item.particleType == "gamma":
                        row = box.row()
                        row.prop(item, "gammaInteraction", text="Interaction")
                    elif item.particleType == "electron":
                        row = box.row()
                        row.prop(item, "electronInteraction", text="Interaction")
                    elif item.particleType == "positron":
                        row = box.row()
                        row.prop(item, "positronInteraction", text="Interaction")

                    # Factor type
                    row = box.row()
                    row.prop(item, "factorType", text="Factor Type")
                    
                    # Factor
                    row = box.row()
                    row.prop(item, "factor", text="Factor")

                    # Factor
                    row = box.row()
                    row.prop(item, "weightWindow", text="Weight Window")
            

            # Check if it is a quadric object
            if obj.penred_settings.quadricType != "unknown":

                box = layout.box()
                box.label(text="Quadric Properties")
                row = box.row()
                row.label(text="Type: " + obj.penred_settings.quadricType)
        
                row = box.row()
                row.prop(obj.penred_settings, "module")
        
                if obj.penred_settings.quadricType == "CONE":
                    row = box.row()
                    row.prop(obj.penred_settings, "r2", text="Top Radius")
                    row = box.row()
                    row.prop(obj.penred_settings, "r1", text="Bottom Radius")

                # Add internal planes
                pbox = box.box()
                pbox.label(text="Cutting Planes")
                for mod in obj.modifiers:
                    if mod.type == "BOOLEAN":
                        bolObj = mod.object
                        if bolObj and hasattr(bolObj, "penred_settings"):
                            if bolObj.penred_settings.quadricType == "CUT_PLANE":
                                row = pbox.row()
                                row.operator("cut.select_operator", text=bolObj.name).cutName = bolObj.name
                                row.operator("quadric.remove_cut_plane",
                                             text="", icon="TRASH").modName = mod.name
                row = pbox.row()
                row.operator("mesh.add_cut_plane_quadric",
                             text="Add Cutting Plane",
                             icon="ADD")

# Source properties menu
##########################
class PenredSourcePropertiesPanel(bpy.types.Panel):
    bl_label = "Source properties"
    bl_idname = "OBJECT_PT_source"
    bl_space_type = 'PROPERTIES'
    bl_region_type = 'WINDOW'
    bl_context = "object"

    def draw(self, context):
        layout = self.layout
        obj = context.object

        if obj and obj.penred_settings:
            
            source = obj.penred_settings.source
            senabled = source.enabled

            row = layout.row()
            row.prop(source, "enabled")
            
            row = layout.row()
            row.enabled = senabled
            if source.nHistsEdit:
                row.prop(source, "nHists", text="")
            else:
                row.label(text=f"Histories :  {source.nHists:.5e}")
                row.prop(source, "nHistsEdit", text="", icon="GREASEPENCIL", toggle=True)

            # CT-like source
            if obj.type == "EMPTY":
                ctbox = layout.box()
                ctbox.enabled = senabled
                ctbox.label(text="CT-Like")
                row = ctbox.row()
                row.prop(source, "ctEnable", text="Enable")
                if source.ctEnable:
                    if source.particleType != "PART_PSF":
                        row = ctbox.row()
                        row.prop(source, "ctSecondaries", text="Secondaries to Sample")

                    row = ctbox.row()
                    row.label(text="Radius")
                    row.prop(source, "ctRad", text="")
                    row.label(text="cm")

                    row = ctbox.row()
                    row.label(text="Rotation Interval")
                    row = ctbox.row()
                    row.prop(source, "ctPhiInterval", text="")
                    row.label(text="Deg")

                    row = ctbox.row()
                    row.label(text="Projection Step")                
                    row.prop(source, "ctPhiStep", text="")
                    row.label(text="Deg")

                    row = ctbox.row()
                    row.label(text="Start Time")
                    row.prop(source, "ctTStart", text="")
                    row.label(text="s")

                    row = ctbox.row()
                    row.label(text="Projection Time")
                    row.prop(source, "ctDT", text="")
                    row.label(text="s")
                
            row = layout.row()
            row.enabled = senabled
            row.prop(source, "particleType", text="Type")

            if source.particleType == "PART_PSF":
                row = layout.row()
                row.enabled = senabled
                row.prop(source, "sourcePSF", text="PSF")

                row = layout.row()
                if source.psfMaxEEdit:
                    row.prop(source, "psfMaxE", text="")
                else:
                    row.label(text=f"Maximum Energy:  {source.psfMaxE:.3e} eV")
                    row.prop(source, "psfMaxEEdit", text="", icon="GREASEPENCIL",
                             toggle=True)

                box = layout.box()
                box.label(text="Splitting and Russian Roulette")
                row = box.row()
                row.label(text="Factor")
                row.prop(source, "split", text="")

                row = box.row()
                row.label(text="Weight Window")
                row.prop(source, "psfWindow", text="")

            else:
                # Spatial
                spatialBox = layout.box()
                spatialBox.label(text="Spatial")
                spatialBox.active = senabled

                row = spatialBox.row()
                row.prop(source, "spatialType", text="Type")
                if source.spatialType == "SPATIAL_CYL":
                    row = spatialBox.row()
                    row.prop(source, "spatialBBFit", text="Inside Bounding Box")
                
                # Energy
                energyBox = layout.box()
                energyBox.label(text="Energy")
                energyBox.active = senabled

                row = energyBox.row()
                row.prop(source, "energyType", text="Type")
        
                row = energyBox.row()
                if source.energyType == "ENERGY_MONO":
                    if source.energyEdit:
                        row.prop(source, "energy", text="")
                    else:
                        row.label(text=f"Energy :  {source.energy:.5e} eV")
                        row.prop(source, "energyEdit", text="", icon="GREASEPENCIL", toggle=True)
                elif source.energyType == "ENERGY_SFILE":
                    row.prop(source, "spcFile", text="File")

                # Direction
                dirBox = layout.box()
                dirBox.label(text="Direction")
                dirBox.active = senabled

                row = dirBox.row()
                row.prop(source, "direction", text="Direction")
        
                row = dirBox.row()
                row.label(text="Aperture")
                row.prop(source, "aperture", text="")
                row.label(text="Deg")

                
            # Time
            timeBox = layout.box()
            timeBox.label(text="Time")
            timeBox.active = senabled

            row = timeBox.row()
            row.prop(source, "timeRecord", text="Record time")

            if source.particleType != "PART_PSF" and not source.ctEnable:
        
                row = timeBox.row()
                row.prop(source, "timeType", text="Time initialization")
        
                if source.timeType == "TIME_DECAY":
                    row = timeBox.row()
                    row.prop(source, "decayHalf", text="Half life")
                    row = timeBox.row()
                    row.label(text="Time Window")
                    row.prop(source, "timeWindow", text="")
                    row.label(text="s")

# Tally properties menu
##########################

# Object talles
class PenredTallyPropertiesPanel(bpy.types.Panel):
    bl_label = "Tally properties"
    bl_idname = "OBJECT_PT_tallies"
    bl_space_type = 'PROPERTIES'
    bl_region_type = 'WINDOW'
    bl_context = "object"

    def draw(self, context):
        layout = self.layout
        obj = context.object

        if obj and obj.penred_settings:

            # Cylindrical dose
            ############################
            boxTallies = layout.box()
            row = boxTallies.row()
            row.prop(obj.penred_settings, "showTalliesCylDose",
                     text="Cylindrical Dose Tallies", emboss=True,
                     icon="TRIA_DOWN" if obj.penred_settings.showTalliesCylDose else "TRIA_RIGHT")
            
            row.operator("tallies_cyldose.add_item", text="", icon="ADD")

            if obj.penred_settings.showTalliesCylDose:
                for i, item in enumerate(obj.penred_settings.talliesCylDose):

                    box = boxTallies.box()
                    row = box.row()
                    row.prop(item, "show", text=item.name, emboss=True,
                             icon="TRIA_DOWN" if item.show else "TRIA_RIGHT")
                    row.operator("tallies_cyldose.remove_item", text="", icon="TRASH").index = i
                    if not item.show:
                        continue

                    # Name
                    row = box.row()
                    row.prop(item, "name", text="Name")

                    box.prop(item, "nr", text="Radial Bins")
                    box.prop(item, "nz", text="Z Bins")
                    box.prop(item, "nPhi", text="Angular Bins")

            # Spatial Dose Distribution
            ############################
            boxTallies = layout.box()
            row = boxTallies.row()
            row.prop(obj.penred_settings, "showTalliesSpatialDose",
                     text="Spatial Dose Tallies", emboss=True,
                     icon="TRIA_DOWN" if obj.penred_settings.showTalliesSpatialDose else "TRIA_RIGHT")
            
            row.operator("tallies_spatialdose.add_item", text="", icon="ADD")
            
            if obj.penred_settings.showTalliesSpatialDose:
                for i, item in enumerate(obj.penred_settings.talliesSpatialDoseDistrib):

                    box = boxTallies.box()
                    row = box.row()
                    row.prop(item, "show", text=item.name, emboss=True,
                             icon="TRIA_DOWN" if item.show else "TRIA_RIGHT")
                    row.operator("tallies_spatialdose.remove_item", text="", icon="TRASH").index = i
                    if not item.show:
                        continue

                    # Name
                    row = box.row()
                    row.prop(item, "name", text="Name")

                    # X
                    row = box.row()
                    row.prop(item, "nx", text="X Bins")

                    # Y
                    row = box.row()
                    row.prop(item, "ny", text="Y Bins")

                    # Z
                    row = box.row()
                    row.prop(item, "nz", text="Z Bins")
                            
            # Spherical dose distribution
            ############################
            boxTallies = layout.box()
            row = boxTallies.row()
            row.prop(obj.penred_settings, "showTalliesSphericalDose",
                     text="Spherical Dose Tallies", emboss=True,
                     icon="TRIA_DOWN" if obj.penred_settings.showTalliesSphericalDose else "TRIA_RIGHT")
            
            row.operator("tallies_sphdose.add_item", text="", icon="ADD")

            if obj.penred_settings.showTalliesSphericalDose:
                for i, item in enumerate(obj.penred_settings.talliesSphericalDoseDistrib):

                    box = boxTallies.box()
                    row = box.row()
                    row.prop(item, "show", text=item.name, emboss=True,
                             icon="TRIA_DOWN" if item.show else "TRIA_RIGHT")
                    row.operator("tallies_sphdose.remove_item", text="", icon="TRASH").index = i
                    if not item.show:
                        continue

                    # Name
                    row = box.row()
                    row.prop(item, "name", text="Name")

                    # R bins
                    row = box.row()
                    row.prop(item, "nr", text="Radial Bins")

                    # Polar
                    row = box.row()
                    row.prop(item, "ntheta", text="Polar Bins")

                    # Azimuth
                    row = box.row()
                    row.prop(item, "nphi", text="Azimuth Bins")

            # Kerma
            ############################
            boxTallies = layout.box()
            row = boxTallies.row()
            row.prop(obj.penred_settings, "showTalliesKerma",
                     text="Kerma Tallies", emboss=True,
                     icon="TRIA_DOWN" if obj.penred_settings.showTalliesKerma else "TRIA_RIGHT")
            
            row.operator("tallies_kerma.add_item", text="", icon="ADD")

            if obj.penred_settings.showTalliesKerma:            
                for i, item in enumerate(obj.penred_settings.talliesKerma):

                    box = boxTallies.box()
                    row = box.row()
                    row.prop(item, "show", text=item.name, emboss=True,
                             icon="TRIA_DOWN" if item.show else "TRIA_RIGHT")
                    row.operator("tallies_kerma.remove_item", text="", icon="TRASH").index = i
                    if not item.show:
                        continue

                    # Name
                    row = box.row()
                    row.prop(item, "name", text="Name")

                    # Energy
                    ebox = box.box()
                    ebox.label(text="Energy")

                    row = ebox.row()
                    if item.eminEdit:
                        row.prop(item, "emin", text="")
                    else:
                        row.label(text=f"Min :  {item.emin:.3e} eV")
                        row.prop(item, "eminEdit", text="", icon="GREASEPENCIL", toggle=True)

                    row = ebox.row()
                    if item.emaxEdit:
                        row.prop(item, "emax", text="")
                    else:
                        row.label(text=f"Max:  {item.emax:.3e} eV")
                        row.prop(item, "emaxEdit", text="", icon="GREASEPENCIL", toggle=True)

                    # Spatial
                    sbox = box.box()
                    sbox.label(text="Spatial")
                    row = sbox.row()
                    row.prop(item, "meshType", text="Mesh")

                    if item.meshType == "MESH_CART":
                        row = sbox.row()
                        row.prop(item, "nx", text="X Bins")
                        row = sbox.row()
                        row.prop(item, "ny", text="Y Bins")
                        row = sbox.row()
                        row.prop(item, "nz", text="Z Bins")
                    elif item.meshType == "MESH_CYL":
                        row = sbox.row()
                        row.prop(item, "nr", text="Radial Bins")
                        row = sbox.row()
                        row.prop(item, "nphi", text="Angular Bins")
                        row = sbox.row()
                        row.prop(item, "nz", text="Z Bins")
                    elif item.meshType == "MESH_SPH":
                        row = sbox.row()
                        row.prop(item, "nr", text="Radial Bins")
                        row = sbox.row()
                        row.prop(item, "ntheta", text="Polar Bins")
                        row = sbox.row()
                        row.prop(item, "nphi", text="Azimuth Bins")


                    row = box.row()
                    row.prop(item, "dataPath", text="Data Path")                
                    
            ## Detector based tallies ##

            if not obj.penred_settings.isDetector:
                return

            # Impact detector
            ############################
            boxTallies = layout.box()
            row = boxTallies.row()
            row.prop(obj.penred_settings, "showTalliesImpactDet",
                     text="Impact Detector Tallies", emboss=True,
                     icon="TRIA_DOWN" if obj.penred_settings.showTalliesImpactDet else "TRIA_RIGHT")
            
            row.operator("tallies_impactdetector.add_item", text="", icon="ADD")

            if obj.penred_settings.showTalliesImpactDet:            
                for i, item in enumerate(obj.penred_settings.talliesImpactDetector):

                    box = boxTallies.box()
                    row = box.row()
                    row.prop(item, "show", text=item.name, emboss=True,
                             icon="TRIA_DOWN" if item.show else "TRIA_RIGHT")
                    row.operator("tallies_impactdetector.remove_item", text="", icon="TRASH").index = i
                    if not item.show:
                        continue

                    # Name
                    row = box.row()
                    row.prop(item, "name", text="Name")

                    if item.fluence or item.spectrum or item.edep:
                        eEnabled = True
                    else:
                        eEnabled = False

                    # Fluence
                    fluenceBox = box.box()
                    fluenceBox.label(text="Fluence")
                    row = fluenceBox.row()
                    row.prop(item, "fluence", text="Enabled")
                    row = fluenceBox.row()
                    row.enabled = item.fluence
                    row.prop(item, "fluenceLogScale", text="Log Scale")

                    # Spectrum
                    spectrumBox = box.box()
                    spectrumBox.label(text="Spectrum")
                    row = spectrumBox.row()
                    row.prop(item, "spectrum", text="Enabled")
                    row = spectrumBox.row()
                    row.enabled = item.spectrum
                    row.prop(item, "spectrumLogScale", text="Log Scale")

                    # Edep
                    edepBox = box.box()
                    edepBox.label(text="Energy Deposition")
                    row = edepBox.row()
                    row.prop(item, "edep", text="Enabled")
                    row = edepBox.row()
                    row.enabled = item.edep
                    row.prop(item, "edepLogScale", text="Log Scale")

                    # Age
                    ageBox = box.box()
                    ageBox.label(text="Age")
                    row = ageBox.row()
                    row.prop(item, "age", text="Enabled")
                    row = ageBox.row()
                    row.enabled = item.age
                    row.prop(item, "ageLogScale", text="Log Scale")                

                    # Energy settings
                    ebox = box.box()
                    ebox.prop(item, "showEBox", text="Energy",emboss=True,
                              icon="TRIA_DOWN" if item.showEBox else "TRIA_RIGHT")

                    if item.showEBox:
                        row = ebox.row()
                        row.enabled = eEnabled
                        if item.eminEdit:
                            row.prop(item, "emin", text="")
                        else:
                            row.label(text=f"Min :  {item.emin:.3e} eV")
                            row.prop(item, "eminEdit", text="", icon="GREASEPENCIL", toggle=True)

                        row = ebox.row()
                        row.enabled = eEnabled
                        if item.emaxEdit:
                            row.prop(item, "emax", text="")
                        else:
                            row.label(text=f"Max:  {item.emax:.3e} eV")
                            row.prop(item, "emaxEdit", text="", icon="GREASEPENCIL", toggle=True)

                        row = ebox.row()
                        row.enabled = eEnabled
                        row.prop(item, "ebins", text="Bins")

                    # Age settings
                    agebox = box.box()
                    agebox.prop(item, "showAgeBox", text="Age",emboss=True,
                              icon="TRIA_DOWN" if item.showAgeBox else "TRIA_RIGHT")

                    if item.showAgeBox:
                        row = agebox.row()
                        row.enabled = item.age
                        row.label(text="Min")
                        row.prop(item, "ageMin", text="")
                        row.label(text="s")

                        row = agebox.row()
                        row.enabled = item.age
                        row.label(text="Max")
                        row.prop(item, "ageMax", text="")
                        row.label(text="s")

                        row = agebox.row()
                        row.enabled = item.age
                        row.prop(item, "ageBins", text="Bins")
                
            # Spatial distribution tally
            ############################
            boxTallies = layout.box()
            row = boxTallies.row()
            row.prop(obj.penred_settings, "showTalliesSpatialDistrib",
                     text="Spatial Distribution Tallies", emboss=True,
                     icon="TRIA_DOWN" if obj.penred_settings.showTalliesSpatialDistrib else "TRIA_RIGHT")
            
            row.operator("tallies_spatialdistrib.add_item", text="", icon="ADD")

            if obj.penred_settings.showTalliesSpatialDistrib:            
                for i, item in enumerate(obj.penred_settings.talliesSpatialDistrib):

                    box = boxTallies.box()
                    row = box.row()
                    row.prop(item, "show", text=item.name, emboss=True,
                             icon="TRIA_DOWN" if item.show else "TRIA_RIGHT")
                    row.operator("tallies_spatialdistrib.remove_item", text="", icon="TRASH").index = i
                    if not item.show:
                        continue

                    # Name
                    row = box.row()
                    row.prop(item, "name", text="Name")

                    # Energy box
                    ebox = box.box()
                    ebox.label(text="Energy")

                    row = ebox.row()
                    if item.eminEdit:
                        row.prop(item, "emin", text="")
                    else:
                        row.label(text=f"Min :  {item.emin:.3e} eV")
                        row.prop(item, "eminEdit", text="", icon="GREASEPENCIL", toggle=True)

                    row = ebox.row()
                    if item.emaxEdit:
                        row.prop(item, "emax", text="")
                    else:
                        row.label(text=f"Max:  {item.emax:.3e} eV")
                        row.prop(item, "emaxEdit", text="", icon="GREASEPENCIL", toggle=True)

                    row = ebox.row()
                    row.prop(item, "ebins", text="Bins")

                    # Spatial box
                    sbox = box.box()
                    sbox.label(text="Spatial")
                    row = sbox.row()
                    row.prop(item, "nx", text="X Bins")
                    row = sbox.row()
                    row.prop(item, "ny", text="Y Bins")
                    row = sbox.row()
                    row.prop(item, "nz", text="Z Bins")

                    row = box.row()
                    row.prop(item, "particleType", text="Particle")

                    row = box.row()
                    row.prop(item, "printCoordinates", text="Coordinates")

                    row = box.row()
                    row.prop(item, "printBins", text="Print Bins")

                    row = box.row()

            # PSF
            ############################
            boxTallies = layout.box()
            row = boxTallies.row()
            row.prop(obj.penred_settings, "showTalliesPSF",
                     text="PSF Tallies", emboss=True,
                     icon="TRIA_DOWN" if obj.penred_settings.showTalliesPSF else "TRIA_RIGHT")
            
            row.operator("tallies_psf.add_item", text="", icon="ADD")

            if obj.penred_settings.showTalliesPSF:            
                for i, item in enumerate(obj.penred_settings.talliesPSF):

                    box = boxTallies.box()
                    row = box.row()
                    row.prop(item, "show", text=item.name, emboss=True,
                             icon="TRIA_DOWN" if item.show else "TRIA_RIGHT")
                    row.operator("tallies_psf.remove_item", text="", icon="TRASH").index = i
                    if not item.show:
                        continue

                    # Name
                    row = box.row()
                    row.prop(item, "name", text="Name")

                    # Minimum and maximum energy
                    row = box.row()
                    if item.eminEdit:
                        row.prop(item, "emin", text="")
                    else:
                        row.label(text=f"Min :  {item.emin:.3e} eV")
                        row.prop(item, "eminEdit", text="", icon="GREASEPENCIL", toggle=True)

                    row = box.row()
                    if item.emaxEdit:
                        row.prop(item, "emax", text="")
                    else:
                        row.label(text=f"Max:  {item.emax:.3e} eV")
                        row.prop(item, "emaxEdit", text="", icon="GREASEPENCIL", toggle=True)

                    # Particles
                    row = box.row()
                    row.prop(item, "gamma", text="Gammas")
                    row = box.row()
                    row.prop(item, "electron", text="Electrons")
                    row = box.row()
                    row.prop(item, "positron", text="Positrons")
            
            # Angular Detector
            ############################
            boxTallies = layout.box()
            row = boxTallies.row()
            row.prop(obj.penred_settings, "showTalliesAngularDet",
                     text="Angular Detector Tallies", emboss=True,
                     icon="TRIA_DOWN" if obj.penred_settings.showTalliesAngularDet else "TRIA_RIGHT")
            
            row.operator("tallies_angdet.add_item", text="", icon="ADD")

            if obj.penred_settings.showTalliesAngularDet:            
                for i, item in enumerate(obj.penred_settings.talliesAngularDetector):

                    box = boxTallies.box()
                    row = box.row()
                    row.prop(item, "show", text=item.name, emboss=True,
                             icon="TRIA_DOWN" if item.show else "TRIA_RIGHT")
                    row.operator("tallies_angdet.remove_item", text="", icon="TRASH").index = i
                    if not item.show:
                        continue

                    # Name
                    row = box.row()
                    row.prop(item, "name", text="Name")

                    # Energy box
                    ebox = box.box()
                    ebox.label(text="Energy")

                    row = ebox.row()
                    if item.eminEdit:
                        row.prop(item, "emin", text="")
                    else:
                        row.label(text=f"Min :  {item.emin:.3e} eV")
                        row.prop(item, "eminEdit", text="", icon="GREASEPENCIL", toggle=True)

                    row = ebox.row()
                    if item.emaxEdit:
                        row.prop(item, "emax", text="")
                    else:
                        row.label(text=f"Max:  {item.emax:.3e} eV")
                        row.prop(item, "emaxEdit", text="", icon="GREASEPENCIL", toggle=True)

                    row = ebox.row()
                    row.prop(item, "ebins", text="Bins")


                    # Polar angle box
                    pbox = box.box()
                    pbox.label(text="Polar Angle")

                    row = pbox.row()
                    row.label(text="Min")
                    row.prop(item, "theta1", text="")
                    row.label(text="Deg")

                    row = pbox.row()
                    row.label(text="Max")
                    row.prop(item, "theta2", text="")
                    row.label(text="Deg")

                    # Azimuthal angle box
                    abox = box.box()
                    abox.label(text="Azimuthal Angle")

                    row = abox.row()
                    row.label(text="Min")
                    row.prop(item, "phi1", text="")
                    row.label(text="Deg")

                    row = abox.row()
                    row.label(text="Max")
                    row.prop(item, "phi2", text="")
                    row.label(text="Deg")

                    row = box.row()
                    row.prop(item, "logScale", text="Logarithmic scale")                

## World simulation properties
class PenredWorldSimulationPanel(bpy.types.Panel):
    bl_label = "Simulation properties"
    bl_idname = "WORLD_PT_simulation"
    bl_space_type = 'PROPERTIES'
    bl_region_type = 'WINDOW'
    bl_context = "world"

    def draw(self, context):
        layout = self.layout
        world = context.scene.world

        if world and hasattr(world, "penred_settings"):
            simulation = world.penred_settings.simulation

            # Dumps
            box = layout.box()
            box.label(text="Dumps")

            # Write
            row = box.row()
            row.prop(simulation, "enableDumps", text="Write")

            row = box.row()
            row.enabled = simulation.enableDumps
            row.label(text="Interval")
            row.prop(simulation, "dumpInterval", text="")
            row.label(text="s")

            row = box.row()
            row.enabled = simulation.enableDumps
            row.label(text="File")
            row.prop(simulation, "dumpWriteFile", text="")
            
            # Read
            row = box.row()
            row.prop(simulation, "readDumps", text="Read")
            
            row = box.row()
            row.enabled = simulation.readDumps
            row.label(text="File")
            row.prop(simulation, "dumpReadFile", text="")

            # Results
            box = layout.box()
            box.label(text="Results")
            
            # Output prefix
            row = box.row()
            row.prop(simulation, "outputPrefix", text="Prefix")
            
            # Final dump
            row = box.row()
            row.prop(simulation, "finalDump", text="Final Dump")

            # ASCII results
            row = box.row()
            row.prop(simulation, "asciiResults", text="ASCII Results")

            # ASCII partial results
            row = box.row()
            row.enabled = simulation.asciiResults
            row.prop(simulation, "asciiResults", text="ASCII Partial Results")

            # Threads
            box = layout.box()
            box.label(text="Threads")

            row = box.row()
            row.prop(simulation, "threadSelType", text="Selection")

            if simulation.threadSelType == "MANUAL":
                row = box.row()
                row.prop(simulation, "nThreads", text = "Number")

            row = box.row()
            row.prop(simulation, "seedPair", text="Seed Pair")

            # Sim Time
            box = layout.box()
            box.label(text="Time Limit")

            row = box.row()
            row.prop(simulation, "limitSimTime", text="Enabled")

            row = box.row()
            row.enabled = simulation.limitSimTime
            row.label(text="Limit")
            row.prop(simulation, "maxSimTime", text="")            
            row.label(text="s")
            
## World tallies
class PenredWorldTalliesPanel(bpy.types.Panel):
    bl_label = "Tallies properties"
    bl_idname = "WORLD_PT_tallies"
    bl_space_type = 'PROPERTIES'
    bl_region_type = 'WINDOW'
    bl_context = "world"

    def draw(self, context):
        layout = self.layout
        world = context.scene.world

        if world and hasattr(world, "penred_settings"):

            # Emerging particle tallies
            ############################
            boxTallies = layout.box()
            row = boxTallies.row()
            row.prop(world.penred_settings, "showTalliesEmergingPart",
                     text="Emerging Particles Tallies", emboss=True,
                     icon="TRIA_DOWN" if world.penred_settings.showTalliesEmergingPart else "TRIA_RIGHT")
            
            row.operator("tallies_emergingparticle.add_item", text="", icon="ADD")

            if world.penred_settings.showTalliesEmergingPart:
                for i, item in enumerate(world.penred_settings.talliesEmergingParticle):

                    box = boxTallies.box()
                    row = box.row()
                    row.prop(item, "show", text=item.name, emboss=True,
                             icon="TRIA_DOWN" if item.show else "TRIA_RIGHT")
                    row.operator("tallies_emergingparticle.remove_item", text="", icon="TRASH").index = i
                    if not item.show:
                        continue

                    # Name
                    row = box.row()
                    row.prop(item, "name", text="Name")

                    # Energy box
                    ebox = box.box()
                    ebox.label(text="Energy")

                    row = ebox.row()
                    if item.eminEdit:
                        row.prop(item, "emin", text="")
                    else:
                        row.label(text=f"Min :  {item.emin:.3e} eV")
                        row.prop(item, "eminEdit", text="", icon="GREASEPENCIL", toggle=True)

                    row = ebox.row()
                    if item.emaxEdit:
                        row.prop(item, "emax", text="")
                    else:
                        row.label(text=f"Max:  {item.emax:.3e} eV")
                        row.prop(item, "emaxEdit", text="", icon="GREASEPENCIL", toggle=True)

                    row = ebox.row()
                    row.prop(item, "ebins", text="Bins")

                    # Angular bins
                    sbox = box.box()
                    sbox.label(text="Spatial")
                    row = sbox.row()
                    row.prop(item, "nTheta", text="Polar bins")
                    row = sbox.row()
                    row.prop(item, "nPhi", text="Azimuthal bins")
                

        # Tracks tally
        tracksTally = world.penred_settings.tracksTally
        box = layout.box()
        box.prop(tracksTally, "show", text=tracksTally.name, emboss=True,
                 icon="TRIA_DOWN" if tracksTally.show else "TRIA_RIGHT")
        if tracksTally.show:
            row = box.row()
            row.prop(tracksTally, "enable", text="Enabled")
        
            row = box.row()
            row.enabled = tracksTally.enable
            row.prop(tracksTally, "nHists", text="Histories")
        

# Material properties menu
##########################
class PenredMaterialPropertiesPanel(bpy.types.Panel):
    bl_label = "Material properties"
    bl_idname = "WORLD_PT_materials"
    bl_space_type = 'PROPERTIES'
    bl_region_type = 'WINDOW'
    bl_context = "world"

    def draw(self, context):
        layout = self.layout
        world = context.scene.world

        if world and hasattr(world, "penred_settings"):

            #Materials
            for i, item in enumerate(world.penred_settings.materials):

                box = layout.box()
                box.prop(item, "show", text=item.name, emboss=True,
                         icon="TRIA_DOWN" if item.show else "TRIA_RIGHT")
                if not item.show:
                    continue

                # Name
                row = box.row()
                row.prop(item, "name", text="Name")

                # Index
                row = box.row()
                row.label(text=f"Index: {i+1}")
                
                # Cutoffs box
                cbox = box.box()
                cbox.label(text="Cutoffs")

                # Gamma
                gammaBox = cbox.box()
                gammaBox.label(text="Gammas")
                row = gammaBox.row()
                row.prop(item, "gammaCutoffType", text="Type")
                if item.gammaCutoffType == "CUTOFF_EABS":

                    row = gammaBox.row()
                    if item.gammaEabsEdit:
                        row.prop(item, "gammaEabs", text="")
                    else:
                        row.label(text=f"Energy :  {item.gammaEabs:.3e} eV")
                        row.prop(item, "gammaEabsEdit", text="", icon="GREASEPENCIL", toggle=True)
                elif item.gammaCutoffType == "CUTOFF_RANGE":

                    row = gammaBox.row()
                    row.label(text="Range: ")
                    row.prop(item, "gammaRange", text="")
                    row.label(text="cm")

                # Electron
                electronBox = cbox.box()
                electronBox.label(text="Electron")
                row = electronBox.row()
                row.prop(item, "electronCutoffType", text="Type")
                if item.electronCutoffType == "CUTOFF_EABS":

                    row = electronBox.row()
                    if item.electronEabsEdit:
                        row.prop(item, "electronEabs", text="")
                    else:
                        row.label(text=f"Energy :  {item.electronEabs:.3e} eV")
                        row.prop(item, "electronEabsEdit", text="", icon="GREASEPENCIL", toggle=True)
                elif item.electronCutoffType == "CUTOFF_RANGE":

                    row = electronBox.row()
                    row.label(text="Range: ")
                    row.prop(item, "electronRange", text="")
                    row.label(text="cm")

                # Positron
                positronBox = cbox.box()
                positronBox.label(text="Positrons")
                row = positronBox.row()
                row.prop(item, "positronCutoffType", text="Type")
                if item.positronCutoffType == "CUTOFF_EABS":

                    row = positronBox.row()
                    if item.positronEabsEdit:
                        row.prop(item, "positronEabs", text="")
                    else:
                        row.label(text=f"Energy :  {item.positronEabs:.3e} eV")
                        row.prop(item, "positronEabsEdit", text="", icon="GREASEPENCIL", toggle=True)
                elif item.positronCutoffType == "CUTOFF_RANGE":

                    row = positronBox.row()
                    row.label(text="Range: ")
                    row.prop(item, "positronRange", text="")
                    row.label(text="cm")

                # Advanced parameters                
                advancedBox = box.box()
                advancedBox.prop(item, "showAdvanced",
                                 text="Advanced Parameters",
                                 emboss=True,
                                 icon="TRIA_DOWN" if item.showAdvanced else "TRIA_RIGHT")

                if item.showAdvanced:
                    row = advancedBox.row()
                    row.label(text="Class II Parameters")
                    row = advancedBox.row()
                    row.prop(item, "enableAdvanced", text="Enable")

                    row = advancedBox.row()
                    row.enabled = item.enableAdvanced
                    row.prop(item, "C1", text="C1")
                
                    row = advancedBox.row()
                    row.enabled = item.enableAdvanced
                    row.prop(item, "C2", text="C2")

                    row = advancedBox.row()
                    row.enabled = item.enableAdvanced
                    if item.WCCEdit:
                        row.prop(item, "WCC", text="")
                    else:
                        row.label(text=f"WCC :  {item.WCC:.3e} eV")
                        row.prop(item, "WCCEdit",
                                 text="", icon="GREASEPENCIL",
                                 toggle=True)

                    row = advancedBox.row()
                    row.enabled = item.enableAdvanced
                    if item.WCREdit:
                        row.prop(item, "WCR", text="")
                    else:
                        row.label(text=f"WCR :  {item.WCR:.3e} eV")
                        row.prop(item, "WCREdit",
                                 text="", icon="GREASEPENCIL",
                                 toggle=True)
                
                # Material composition
                elementsBox = box.box()
                elementsBox.label(text="Composition")
                row = elementsBox.row()
                row.label(text="Density")
                row.prop(item, "density", text="")
                row.label(text="g/cm^3")
                for ie, element in enumerate(item.composition):
                    row = elementsBox.row()
                    split = row.split(factor=0.35)
                    col = split.column()
                    col.prop(element, "z", text="Z")
                    col = split.column()
                    col.prop(element, "wFraction", text="Weight fraction")

                # Operators to add and remove composition elements
                row = elementsBox.row()
                row.operator("materials_composition.add_item", text="Add Element").imat = i
                row = elementsBox.row()
                row.operator("materials_composition.remove_item", text="Remove Element").imat = i

            # Operators to add and remove composition elements
            row = layout.row()
            row.operator("materials_material.add_item", text="Add Material")
            row = layout.row()
            row.operator("materials_material.remove_item", text="Remove Material")
            
# Manual
##########
def add_object_manual_map():
    # TODO
    url_manual_prefix = "https://docs.blender.org/manual/en/latest/"
    url_manual_mapping = (
        ("bpy.ops.mesh.add_object", "scene_layout/object/types.html"),
    )
    return url_manual_prefix, url_manual_mapping

# Hints draw handler
def drawHintsHandler():

    # Define colors
    sourceColor = (1.0, 0.2, 0.2, 0.2)
    tallyColor = (0.2, 0.2, 1.0, 0.2)
    
    context = bpy.context
    for obj in context.scene.objects:
        if obj and obj.penred_settings:
            penSett = obj.penred_settings
            if penSett.quadricType == "CUT_PLANE":
                if obj.parent:
                    # Hide this objects if are not active
                    if obj == bpy.context.view_layer.objects.active:
                        if obj.hide_viewport:
                            obj.hide_viewport = False
                    elif not obj.hide_viewport:
                        obj.hide_viewport = True
                        
                else:
                    bpy.data.objects.remove(obj)
                
            if not obj.hide_get() and not obj.hide_viewport:

                if penSett.source and penSett.source.enabled:
                    source = penSett.source

                    # Print source direction
                    if source.particleType == "PART_PSF":
                        rotationMatrix = obj.matrix_world.to_3x3()
                        
                        localZ = Vector((0,0,1))
                        worldZ = rotationMatrix @ localZ

                        if source.ctEnable:
                            startPhi = source.ctPhiInterval[0]*pi/180.0
                            endPhi   = source.ctPhiInterval[1]*pi/180.0
                            dangle = source.ctPhiStep*pi/180.0
                            nProj = int((endPhi-startPhi)/dangle)+1

                            # Draw CT circle
                            utils.draw_zcircle(obj, source.ctRad,
                                               startPhi, endPhi,
                                               sourceColor)

                            #Draw projection arrows
                            iangle = 0
                            while iangle < nProj:
                                angle = startPhi + iangle*dangle
                                cAngle = cos(angle)
                                sAngle = sin(angle)
                                aux = (
                                    cAngle*worldZ[0] - sAngle*worldZ[1],
                                    sAngle*worldZ[0] - cAngle*worldZ[1],
                                    worldZ[2]
                                )
                                utils.draw_arrow(obj, aux, sourceColor,
                                                 (source.ctRad*cAngle,
                                                  source.ctRad*sAngle,
                                                  0.0))
                                iangle = iangle + 1
                                
                        else:                        
                            utils.draw_arrow(obj, worldZ, sourceColor)
                    else:
                        # Get the source direction
                        direction = source.direction

                        if source.ctEnable:
                            startPhi = source.ctPhiInterval[0]*pi/180.0
                            endPhi   = source.ctPhiInterval[1]*pi/180.0
                            dangle = source.ctPhiStep*pi/180.0
                            nProj = int((endPhi-startPhi)/dangle)+1

                            # Draw CT circle
                            utils.draw_zcircle(obj, source.ctRad,
                                               startPhi, endPhi,
                                               sourceColor)

                            #Draw projection arrows
                            iangle = 0
                            while iangle < nProj:
                                angle = startPhi + iangle*dangle
                                cAngle = cos(angle)
                                sAngle = sin(angle)
                                shift = (source.ctRad*cAngle, source.ctRad*sAngle, 0.0)
                                aux = (
                                    cAngle*direction[0] - sAngle*direction[1],
                                    sAngle*direction[0] - cAngle*direction[1],
                                    direction[2]
                                )
                                utils.draw_arrow(obj, aux, sourceColor, shift)
                                iangle = iangle + 1

                        else:
                            utils.draw_arrow(obj, direction, sourceColor)

                        if source.spatialType == "SPATIAL_BOX":
                            utils.draw_bbox(obj, sourceColor)
                        elif source.spatialType == "SPATIAL_CYL":
                            utils.draw_zcyl(obj, source.spatialBBFit, sourceColor)

                ## Tallies ##

                # Check kerma tallies
                kermaCart = False
                kermaCyl = False
                kermaSph = False
                for t in penSett.talliesKerma:
                    if t.meshType == "MESH_CART":
                        kermaCart = True
                    if t.meshType == "MESH_CYL":
                        kermaCyl = True
                    if t.meshType == "MESH_SPH":
                        kermaSph = True
                
                if kermaCyl or len(penSett.talliesCylDose) > 0:
                    utils.draw_zcyl(obj, False, tallyColor)
                if kermaCart or len(penSett.talliesSpatialDoseDistrib) + len(penSett.talliesSpatialDistrib) > 0:
                    utils.draw_bbox(obj, tallyColor)
                if kermaSph or len(penSett.talliesSphericalDoseDistrib) > 0:
                    utils.draw_sphere(obj, tallyColor)

                # Planes
                if penSett.quadricType == "PLANE" or penSett.quadricType == "CUT_PLANE":
                    # Get the object's local rotation matrix
                    rotationMatrix = obj.matrix_world.to_3x3()

                    localZ = Vector((0,0,1))
                    worldZ = rotationMatrix @ localZ

                    color = (0.0, 0.0, 0.0, 1.0)
                    utils.draw_arrow(obj, worldZ, color)

# Register functions
def register():    
    bpy.types.TOPBAR_MT_file_export.append(menu_func_penred_export)
    bpy.types.VIEW3D_MT_transform_object.append(menu_func_quadric_transform)

    #Register view submenu and the corresponding add function
    bpy.utils.register_class(OBJECT_MT_view_quadric_submenu)
    bpy.types.VIEW3D_MT_view.append(menu_quadric_view_func)        
        
    #Register quadric add menu and the corresponding add function
    bpy.utils.register_class(OBJECT_MT_quadric_submenu)
    bpy.types.VIEW3D_MT_add.append(menu_quadric_add_func)
    
    #Register object properties Panel
    bpy.utils.register_class(PenredBodyPropertiesPanel)

    #Register source properties Panel
    bpy.utils.register_class(PenredSourcePropertiesPanel)

    #Register objects based tallies Panel
    bpy.utils.register_class(PenredTallyPropertiesPanel)

    #Register materials properties Panel
    bpy.utils.register_class(PenredMaterialPropertiesPanel)

    #Register simulation parameters Panel
    bpy.utils.register_class(PenredWorldSimulationPanel)

    #Register world based tallies Panel
    bpy.utils.register_class(PenredWorldTalliesPanel)    

    #Register manual (To be done)
    bpy.utils.register_manual_map(add_object_manual_map)

    #Register hints draw handler
    bpy.types.SpaceView3D.draw_handler_add(drawHintsHandler, (), 'WINDOW', 'POST_VIEW')
    
def unregister():
    bpy.types.TOPBAR_MT_file_export.remove(menu_func_penred_export)
    bpy.types.VIEW3D_MT_transform_object.remove(menu_func_quadric_transform)
        
    #Unregister view submenu and the corresponding add function
    bpy.utils.unregister_class(OBJECT_MT_view_quadric_submenu)
    bpy.types.VIEW3D_MT_view.remove(menu_quadric_view_func)
    
    #Unregister quadric add menu and the corresponding add function
    bpy.utils.unregister_class(OBJECT_MT_quadric_submenu)
    bpy.types.VIEW3D_MT_add.remove(menu_quadric_add_func)
    
    #Unregister object properties Panel
    bpy.utils.unregister_class(PenredBodyPropertiesPanel)

    #Unregister source properties Panel
    bpy.utils.unregister_class(PenredSourcePropertiesPanel)    

    #Unregister tally properties Panel
    bpy.utils.unregister_class(PenredTallyPropertiesPanel)

    #Unregister materials properties Panel
    bpy.utils.unregister_class(PenredMaterialPropertiesPanel)

    #Unregister simulation parameters Panel
    bpy.utils.unregister_class(PenredWorldSimulationPanel)    

    #Unregister world based tallies Panel
    bpy.utils.unregister_class(PenredWorldTalliesPanel)    
    
    #Unregister manual (To be done)
    bpy.utils.unregister_manual_map(add_object_manual_map)    

    #Unregister hints draw handler
    bpy.types.SpaceView3D.draw_handler_remove(drawHintsHandler, 'WINDOW')
