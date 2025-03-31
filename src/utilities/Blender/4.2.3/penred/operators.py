
import bpy
import bmesh
import gpu
from gpu_extras.batch import batch_for_shader
import blf
from bpy_extras.io_utils import ExportHelper
from bpy.types import Operator
from bpy.props import FloatVectorProperty
from bpy_extras.object_utils import AddObjectHelper, object_data_add
from mathutils import Vector
from mathutils import Color
from math import cos, acos, sin, asin, tan, atan2, sqrt, pi
import os
import time
from . import surfaces, utils, addon_properties, conf

### Material view
class QUADRIC_OT_view_material(Operator, AddObjectHelper):
    """Set Quadric material view"""
    bl_idname = "view.quadric_material"
    bl_label = "Set Quadric Material View"

    def execute(self, context):
        for obj in context.scene.objects:
            if obj.penred_settings.quadricType != "unknown":
                #Is a quad type, assign a color according to material
                material = obj.penred_settings.material
                #Get alpha channel
                alpha = obj.color[3]
                
                H = float(material % addon_properties.materialColumns)/float(addon_properties.materialColumns)
                V = 1.0-float(material // addon_properties.materialRows)/float(addon_properties.materialRows)
                if V < 0.1:
                    V = 0.1
                
                c = Color()
                c.hsv = H, 1.0, V
                obj.color = c.r, c.g, c.b, alpha   
                
        return {'FINISHED'}    

### Module view
class QUADRIC_OT_view_module(Operator, AddObjectHelper):
    """Set Quadric module view"""
    bl_idname = "view.quadric_module"
    bl_label = "Set Quadric Inner View"

    def execute(self, context):
        for obj in context.scene.objects:
            if obj.penred_settings.quadricType != "unknown":
                #Is a quad type, assign a color according to module
                module = obj.penred_settings.module
                #Get alpha channel
                alpha = obj.color[3]
                
                if module:                
                    obj.color = 1.0, 1.0, 1.0, alpha   
                else:
                    obj.color = 0.0, 0.0, 0.0, alpha   
                
        return {'FINISHED'}

viewObjectClasses = (
    QUADRIC_OT_view_material,
    QUADRIC_OT_view_module,
)

### Tallies
#############

## Object based tallies ##

# Add Cylindrical dose
class TALLY_OT_addCylDoseTally(bpy.types.Operator):
    bl_idname = "tallies_cyldose.add_item"
    bl_label = "Add Item"
    bl_description = "Add a cylindrical dose tally"

    def execute(self, context):
        obj = context.object
        if obj and obj.penred_settings:
            obj.penred_settings.talliesCylDose.add()
        return {"FINISHED"}

# Remove Cylindrical dose
class TALLY_OT_removeCylDoseTally(bpy.types.Operator):
    bl_idname = "tallies_cyldose.remove_item"
    bl_label = "Remove Item"
    bl_description = "Remove the cylindrical dose tally"

    index: bpy.props.IntProperty()  # Index of the item to remove

    def execute(self, context):
        obj = context.object
        if obj and obj.penred_settings:
            obj.penred_settings.talliesCylDose.remove(self.index)
        return {"FINISHED"}

# Add Impact Detector
class TALLY_OT_addImpactDetectorTally(bpy.types.Operator):
    bl_idname = "tallies_impactdetector.add_item"
    bl_label = "Add Item"
    bl_description = "Add a impact detector tally"

    def execute(self, context):
        obj = context.object
        if obj and obj.penred_settings:
            obj.penred_settings.talliesImpactDetector.add()
        return {"FINISHED"}

# Remove Impact Detector
class TALLY_OT_removeImpactDetectorTally(bpy.types.Operator):
    bl_idname = "tallies_impactdetector.remove_item"
    bl_label = "Remove Item"
    bl_description = "Remove the impact detector tally"

    index: bpy.props.IntProperty()  # Index of the item to remove

    def execute(self, context):
        obj = context.object
        if obj and obj.penred_settings:
            obj.penred_settings.talliesImpactDetector.remove(self.index)
        return {"FINISHED"}

# Add Spatial dose
class TALLY_OT_addSpatialDoseTally(bpy.types.Operator):
    bl_idname = "tallies_spatialdose.add_item"
    bl_label = "Add Item"
    bl_description = "Add a spatial dose tally"

    def execute(self, context):
        obj = context.object
        if obj and obj.penred_settings:
            obj.penred_settings.talliesSpatialDoseDistrib.add()
        return {"FINISHED"}

# Remove Spatial dose
class TALLY_OT_removeSpatialDoseTally(bpy.types.Operator):
    bl_idname = "tallies_spatialdose.remove_item"
    bl_label = "Remove Item"
    bl_description = "Remove the spatial dose tally"

    index: bpy.props.IntProperty()  # Index of the item to remove

    def execute(self, context):
        obj = context.object
        if obj and obj.penred_settings:
            obj.penred_settings.talliesSpatialDoseDistrib.remove(self.index)
        return {"FINISHED"}

# Add Spherical dose
class TALLY_OT_addSphDoseTally(bpy.types.Operator):
    bl_idname = "tallies_sphdose.add_item"
    bl_label = "Add Item"
    bl_description = "Add a spherical dose tally"

    def execute(self, context):
        obj = context.object
        if obj and obj.penred_settings:
            obj.penred_settings.talliesSphericalDoseDistrib.add()
        return {"FINISHED"}

# Remove Spherical dose
class TALLY_OT_removeSphDoseTally(bpy.types.Operator):
    bl_idname = "tallies_sphdose.remove_item"
    bl_label = "Remove Item"
    bl_description = "Remove the spherical dose tally"

    index: bpy.props.IntProperty()  # Index of the item to remove

    def execute(self, context):
        obj = context.object
        if obj and obj.penred_settings:
            obj.penred_settings.talliesSphericalDoseDistrib.remove(self.index)
        return {"FINISHED"}    

# Add PSF
class TALLY_OT_addPSFTally(bpy.types.Operator):
    bl_idname = "tallies_psf.add_item"
    bl_label = "Add Item"
    bl_description = "Add a phase space file tally"

    def execute(self, context):
        obj = context.object
        if obj and obj.penred_settings:
            obj.penred_settings.talliesPSF.add()
        return {"FINISHED"}

# Remove PSF
class TALLY_OT_removePSFTally(bpy.types.Operator):
    bl_idname = "tallies_psf.remove_item"
    bl_label = "Remove Item"
    bl_description = "Remove the phase space file tally"

    index: bpy.props.IntProperty()  # Index of the item to remove

    def execute(self, context):
        obj = context.object
        if obj and obj.penred_settings:
            obj.penred_settings.talliesPSF.remove(self.index)
        return {"FINISHED"}

# Add kerma tally
class TALLY_OT_addKermaTally(bpy.types.Operator):
    bl_idname = "tallies_kerma.add_item"
    bl_label = "Add Item"
    bl_description = "Add a kerma track length tally"

    def execute(self, context):
        obj = context.object
        if obj and obj.penred_settings:
            obj.penred_settings.talliesKerma.add()
        return {"FINISHED"}

# Remove Kerma tally
class TALLY_OT_removeKermaTally(bpy.types.Operator):
    bl_idname = "tallies_kerma.remove_item"
    bl_label = "Remove Item"
    bl_description = "Remove the kerma track length tally"

    index: bpy.props.IntProperty()  # Index of the item to remove

    def execute(self, context):
        obj = context.object
        if obj and obj.penred_settings:
            obj.penred_settings.talliesKerma.remove(self.index)
        return {"FINISHED"}

# Add Spatial Distribution
class TALLY_OT_addSpatialDistribTally(bpy.types.Operator):
    bl_idname = "tallies_spatialdistrib.add_item"
    bl_label = "Add Item"
    bl_description = "Add a spatial distribution tally"

    def execute(self, context):
        obj = context.object
        if obj and obj.penred_settings:
            obj.penred_settings.talliesSpatialDistrib.add()
        return {"FINISHED"}

# Remove Spatial Distribution
class TALLY_OT_removeSpatialDistribTally(bpy.types.Operator):
    bl_idname = "tallies_spatialdistrib.remove_item"
    bl_label = "Remove Item"
    bl_description = "Remove the spatial distribution tally"

    index: bpy.props.IntProperty()  # Index of the item to remove

    def execute(self, context):
        obj = context.object
        if obj and obj.penred_settings:
            obj.penred_settings.talliesSpatialDistrib.remove(self.index)
        return {"FINISHED"}    

# Add Angular detector
class TALLY_OT_addAngDetTally(bpy.types.Operator):
    bl_idname = "tallies_angdet.add_item"
    bl_label = "Add Item"
    bl_description = "Add an angular detector tally"

    def execute(self, context):
        obj = context.object
        if obj and obj.penred_settings:
            obj.penred_settings.talliesAngularDetector.add()
        return {"FINISHED"}

# Remove Angular detector
class TALLY_OT_removeAngDetTally(bpy.types.Operator):
    bl_idname = "tallies_angdet.remove_item"
    bl_label = "Remove Item"
    bl_description = "Remove the angular detector tally"

    index: bpy.props.IntProperty()  # Index of the item to remove

    def execute(self, context):
        obj = context.object
        if obj and obj.penred_settings:
            obj.penred_settings.talliesAngularDetector.remove(self.index)
        return {"FINISHED"}

## World based tallies ##

# Add Cylindrical dose
class TALLY_OT_addEmergingParticleTally(bpy.types.Operator):
    bl_idname = "tallies_emergingparticle.add_item"
    bl_label = "Add Item"
    bl_description = "Add a emerging particle distribution tally"

    def execute(self, context):
        world = context.scene.world
        if world and world.penred_settings:
            world.penred_settings.talliesEmergingParticle.add()
        return {"FINISHED"}

# Remove Cylindrical dose
class TALLY_OT_EmergingParticleTally(bpy.types.Operator):
    bl_idname = "tallies_emergingparticle.remove_item"
    bl_label = "Remove Item"
    bl_description = "Remove the emerging particle distribution tally"

    index: bpy.props.IntProperty()  # Index of the item to remove

    def execute(self, context):
        world = context.scene.world
        if world and world.penred_settings:
            world.penred_settings.talliesEmergingParticle.remove(self.index)
        return {"FINISHED"}
    
talliesOperatorClasses = (
    TALLY_OT_addCylDoseTally,
    TALLY_OT_removeCylDoseTally,
    TALLY_OT_addImpactDetectorTally,
    TALLY_OT_removeImpactDetectorTally,
    TALLY_OT_addSpatialDoseTally,
    TALLY_OT_removeSpatialDoseTally,
    TALLY_OT_addSphDoseTally,
    TALLY_OT_removeSphDoseTally,
    TALLY_OT_addPSFTally,
    TALLY_OT_removePSFTally,
    TALLY_OT_addKermaTally,
    TALLY_OT_removeKermaTally,
    TALLY_OT_addSpatialDistribTally,
    TALLY_OT_removeSpatialDistribTally,
    TALLY_OT_addAngDetTally,
    TALLY_OT_removeAngDetTally,
    TALLY_OT_addEmergingParticleTally,
    TALLY_OT_EmergingParticleTally,
)

### Materials
###############

# Add Material
class MAT_OT_addMaterial(bpy.types.Operator):
    bl_idname = "materials_material.add_item"
    bl_label = "Add Material"
    bl_description = "Add a material"

    def execute(self, context):
        world = context.scene.world
        if world and world.penred_settings:
            world.penred_settings.materials.add()
            nMats = len(world.penred_settings.materials)
            world.penred_settings.materials[-1].name=f"material{nMats:03d}"
            world.penred_settings.materials[-1].composition.add()
        return {"FINISHED"}

# Remove Composition Element
class MAT_OT_removeMaterial(bpy.types.Operator):
    bl_idname = "materials_material.remove_item"
    bl_label = "Remove Material"
    bl_description = "Remove last material"

    def execute(self, context):
        world = context.scene.world
        if world and world.penred_settings:
            nMats = len(world.penred_settings.materials)
            if nMats > 1:
                world.penred_settings.materials.remove(nMats-1)
        return {"FINISHED"}

# Add Composition Element
class MAT_OT_addCompositionElement(bpy.types.Operator):
    bl_idname = "materials_composition.add_item"
    bl_label = "Add Element"
    bl_description = "Add a composition element for this material"

    imat: bpy.props.IntProperty()  # Material index

    def execute(self, context):
        world = context.scene.world
        if world and world.penred_settings:
            world.penred_settings.materials[self.imat].composition.add()
        return {"FINISHED"}

# Remove Composition Element
class MAT_OT_removeCompositionElement(bpy.types.Operator):
    bl_idname = "materials_composition.remove_item"
    bl_label = "Remove Element"
    bl_description = "Remove a composition element for this material"

    imat: bpy.props.IntProperty()  # Material index

    def execute(self, context):
        world = context.scene.world
        if world and world.penred_settings:
            nElements = len(world.penred_settings.materials[self.imat].composition)
            if nElements > 1:
                world.penred_settings.materials[self.imat].composition.remove(nElements-1)
        return {"FINISHED"}
    
materialsOperatorClasses = (
    MAT_OT_addMaterial,
    MAT_OT_removeMaterial,
    MAT_OT_addCompositionElement,
    MAT_OT_removeCompositionElement,
)    
    
### Put down tool
class QUADRIC_OT_transform_putdown(Operator, AddObjectHelper):
    bl_idname = "transform.quadric_transform_putdown"
    bl_label = "Tool to clip objects down of the active one"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        selected = context.selected_objects
        active = context.view_layer.objects.active
        
        if len(selected) == 0:
            return {'FINISHED'}
        
        #Get active object position
        transAct = active.matrix_world.to_translation()

        #Get Z size
        dza = active.dimensions.z
        
        #Get local Z direction
        dirZ = utils.getLocalZdir(active)

        for obj in selected:
            if obj != active:
                
                #Get Z size
                dz = obj.dimensions.z
                
                #Get position
                obj.matrix_world.translation = transAct - dirZ*(dza/2.0 + dz/2.0)
                
        return {'FINISHED'}


### Put top tool
class QUADRIC_OT_transform_puttop(Operator, AddObjectHelper):
    bl_idname = "transform.quadric_transform_puttop"
    bl_label = "Tool to clip objects top of the active one"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        selected = context.selected_objects
        active = context.view_layer.objects.active
        
        if len(selected) == 0:
            return {'FINISHED'}
        
        #Get active object position
        transAct = active.matrix_world.to_translation()

        #Get Z size
        dza = active.dimensions.z
        
        #Get local Z direction
        dirZ = utils.getLocalZdir(active)

        for obj in selected:
            if obj != active:
                
                #Get Z size
                dz = obj.dimensions.z
                
                #Get position
                obj.matrix_world.translation = transAct + dirZ*(dza/2.0 + dz/2.0)
                
        return {'FINISHED'}
    
### Put inside top tool
class QUADRIC_OT_transform_putinsidetop(Operator, AddObjectHelper):
    bl_idname = "transform.quadric_transform_putinsidetop"
    bl_label = "Tool to clip objects inside top of the active one"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        selected = context.selected_objects
        active = context.view_layer.objects.active
        
        if len(selected) == 0:
            return {'FINISHED'}
        
        #Get active object position
        transAct = active.matrix_world.to_translation()

        #Get Z size
        dza = active.dimensions.z
        
        #Get local Z direction
        dirZ = utils.getLocalZdir(active)

        for obj in selected:
            if obj != active:
                
                #Get Z size
                dz = obj.dimensions.z
                
                #Get position
                obj.matrix_world.translation = transAct + dirZ*(dza/2.0 - dz/2.0)
                
        return {'FINISHED'}

### Put inside down tool
class QUADRIC_OT_transform_putinsidedown(Operator, AddObjectHelper):
    bl_idname = "transform.quadric_transform_putinsidedown"
    bl_label = "Tool to clip objects inside down of the active one"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        selected = context.selected_objects
        active = context.view_layer.objects.active
        
        if len(selected) == 0:
            return {'FINISHED'}
        
        #Get active object position
        transAct = active.matrix_world.to_translation()

        #Get Z size
        dza = active.dimensions.z
        
        #Get local Z direction
        dirZ = utils.getLocalZdir(active)

        for obj in selected:
            if obj != active:
                
                #Get Z size
                dz = obj.dimensions.z
                
                #Get position
                obj.matrix_world.translation = transAct - dirZ*(dza/2.0 - dz/2.0)
                
        return {'FINISHED'}

### Put inside centered tool
class QUADRIC_OT_transform_putinsidecentered(Operator, AddObjectHelper):
    bl_idname = "transform.quadric_transform_putinsidecentered"
    bl_label = "Tool to set objects inside of the active one centered"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        selected = context.selected_objects
        active = context.view_layer.objects.active
        
        if len(selected) == 0:
            return {'FINISHED'}
        
        #Get active object position
        transAct = active.matrix_world.to_translation()

        #Get local Z direction
        dirZ = utils.getLocalZdir(active)

        for obj in selected:
            if obj != active:
                
                #Get position
                obj.matrix_world.translation = transAct
                
        return {'FINISHED'}
    
quadricTransformClasses = (
    QUADRIC_OT_transform_putdown,
    QUADRIC_OT_transform_puttop,
    QUADRIC_OT_transform_putinsidetop,
    QUADRIC_OT_transform_putinsidedown,
    QUADRIC_OT_transform_putinsidecentered
)

# Object operators
##########################
# Add Cylindrical dose
class VR_OT_addIF(bpy.types.Operator):
    bl_idname = "vr_if.add_item"
    bl_label = "Add Item"
    bl_description = "Add a interaction forcing"

    def execute(self, context):
        obj = context.object
        if obj and obj.penred_settings:
            obj.penred_settings.interactionForcing.add()
        return {"FINISHED"}

# Remove Cylindrical dose
class VR_OT_removeIF(bpy.types.Operator):
    bl_idname = "vr_if.remove_item"
    bl_label = "Remove Item"
    bl_description = "Remove the interaction forcing"

    index: bpy.props.IntProperty()  # Index of the item to remove

    def execute(self, context):
        obj = context.object
        if obj and obj.penred_settings:
            obj.penred_settings.interactionForcing.remove(self.index)
        return {"FINISHED"}

class OBJECT_OT_cut_select_operator(Operator):
    bl_idname = 'cut.select_operator'
    bl_label = 'Select Cut'
    bl_description = 'Select a cut plane'
    bl_options = {'REGISTER', 'UNDO'}

    cutName: bpy.props.StringProperty(name = "Cut Plane Name", default="")    
    
    def execute(self, context):
        obj = context.object
        if obj and obj.penred_settings and self.cutName:
            plane = None
            for child in obj.children:
                if child.name == self.cutName:
                    plane = child
                    break

            if plane and plane.penred_settings and plane.penred_settings.quadricType == "CUT_PLANE":
                # Resize cutting plane according to parents bounding box size
                x,y,z,dx,dy,dz,sx,sy,sz,omega,theta,phi,name,boxSize = utils.getObjInfo(obj)
                maxSize = Vector(boxSize).length
                plane.scale = Vector((1.01*maxSize,1.01*maxSize,1.0))

                # Set the plane as active object
                bpy.ops.object.select_all(action='DESELECT')
                bpy.context.view_layer.objects.active = plane
        return {'FINISHED'}

class REMESH_OT_remesh_operator(Operator):
    bl_idname = 'remesh.remesh_operator'
    bl_label = 'Remesh'
    bl_description = 'Remesh object'
    bl_options = {'REGISTER', 'UNDO'}
 
    def execute(self, context):
        obj = context.object
        if obj.penred_settings.quadricType == "CONE":
            
            #Get height
            h = obj.dimensions.z
            #Get radius
            r1 = obj.penred_settings.r1
            r2 = obj.penred_settings.r2
            
            r1 = max(1.0e-5,r1)
            r2 = max(1.0e-5,r2)
                        
            ### Create an empty mesh
            mesh = bpy.data.meshes.new(name=obj.data.name)
            
            #Create the new mesh
            bm = bmesh.new()
            #Despite the name, blender uses the diameter parameters as radius OLD
            if obj.penred_settings.quadricType == "CONE":
                bmesh.ops.create_cone(bm, cap_ends=True, segments=32, radius1=r1, radius2=r2, depth = h)
            else:
                bmesh.ops.create_cone(bm, cap_ends=False, segments=32, radius1=r1, radius2=r2, depth = h)
            
            bm.to_mesh(mesh)
            bm.free()
                        
            #Store old mesh
            oldMesh = obj.data
            
            #Assign new mesh
            obj.data = mesh
            
            #Remove old mesh
            if oldMesh.users == 0:
                oldMesh.user_clear()
                bpy.data.meshes.remove(oldMesh)
            
        return {'FINISHED'}

objectOperatorsClasses = (
    VR_OT_addIF,
    VR_OT_removeIF,
    OBJECT_OT_cut_select_operator,
    REMESH_OT_remesh_operator,
)

# Object add operators
##########################

# CUBE
class QUADRIC_OT_add_cube(Operator, AddObjectHelper):
    """Create a new Cube Quadric"""
    bl_idname = "mesh.add_cube_quadric"
    bl_label = "Add Cube Quadric"
    bl_options = {'REGISTER', 'UNDO'}

    scale: FloatVectorProperty(
        name="scale",
        default=(1.0, 1.0, 1.0),
        subtype='TRANSLATION',
        description="scaling",
    )

    def execute(self, context):
        utils.add_object(self, context, "Cube", "CUBE")
        return {'FINISHED'}

# SPHERE
class QUADRIC_OT_add_sphere(Operator, AddObjectHelper):
    """Create a new Spheric Quadric"""
    bl_idname = "mesh.add_sphere_quadric"
    bl_label = "Add Sphere Quadric"
    bl_options = {'REGISTER', 'UNDO'}

    scale: FloatVectorProperty(
        name="scale",
        default=(1.0, 1.0, 1.0),
        subtype='TRANSLATION',
        description="scaling",
    )

    def execute(self, context):
        utils.add_object(self, context, "Sphere", "SPHERE")
        return {'FINISHED'}

# CONE
class QUADRIC_OT_add_cone(Operator, AddObjectHelper):
    """Create a new Cone Quadric"""
    bl_idname = "mesh.add_cone_quadric"
    bl_label = "Add Cone Quadric"
    bl_options = {'REGISTER', 'UNDO'}

    scale: FloatVectorProperty(
        name="scale",
        default=(1.0, 1.0, 1.0),
        subtype='TRANSLATION',
        description="scaling",
    )

    def execute(self, context):
        utils.add_object(self, context, "Cone", "CONE")
        return {'FINISHED'}

# CYLINDER
class QUADRIC_OT_add_cylinder(Operator, AddObjectHelper):
    """Create a new Cylinder Quadric"""
    bl_idname = "mesh.add_cylinder_quadric"
    bl_label = "Add Cylinder Quadric"
    bl_options = {'REGISTER', 'UNDO'}

    scale: FloatVectorProperty(
        name="scale",
        default=(1.0, 1.0, 1.0),
        subtype='TRANSLATION',
        description="scaling",
    )

    def execute(self, context):
        utils.add_object(self, context, "Cylinder", "CYLINDER")
        return {'FINISHED'}

# PLANE
class QUADRIC_OT_add_plane(Operator, AddObjectHelper):
    """Create a new Plane Quadric"""
    bl_idname = "mesh.add_plane_quadric"
    bl_label = "Add Plane Quadric"
    bl_options = {'REGISTER', 'UNDO'}

    scale: FloatVectorProperty(
        name="scale",
        default=(1.0, 1.0, 1.0),
        subtype='TRANSLATION',
        description="scaling",
    )

    def execute(self, context):
        utils.add_object(self, context, "Plane", "PLANE")
        return {'FINISHED'}

# SEMI_SPHERE
class QUADRIC_OT_add_semiSphere(Operator, AddObjectHelper):
    """Create a new Tube Quadric"""
    bl_idname = "mesh.add_semisphere_quadric"
    bl_label = "Add SemiSphere Quadric"
    bl_options = {'REGISTER', 'UNDO'}

    scale: FloatVectorProperty(
        name="scale",
        default=(1.0, 1.0, 1.0),
        subtype='TRANSLATION',
        description="scaling",
    )

    def execute(self, context):
        utils.add_object(self, context, "SemiSphere", "SEMI_SPHERE")
        return {'FINISHED'}

# CUT PLANE
class QUADRIC_OT_add_cutplane(Operator, AddObjectHelper):
    """Create a new Cutting Plane Quadric"""
    bl_idname = "mesh.add_cut_plane_quadric"
    bl_label = "Add Cutting Plane Quadric"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):

        parent = context.object
        if parent and parent.penred_settings:
        
            # Create the cutting plane
            cutPlane = utils.add_object(self, context, "Plane", "CUT_PLANE")
            cutPlane.name = "Cut Plane"

            # Save the cut plane transform
            matrixWorld = cutPlane.matrix_world.copy()
            
            # Set parent
            cutPlane.parent = parent

            # Apply the parent's inverse transform
            cutPlane.matrix_parent_inverse = parent.matrix_world.inverted()
            
            # Restore the cut plane world transform
            cutPlane.matrix_world = matrixWorld

            # Add boolean modifier to the parent object
            mod = parent.modifiers.new(name="boolean", type="BOOLEAN")
            mod.operation = "DIFFERENCE"
            mod.object = cutPlane

            # Resize cutting plane according to parents bounding box size
            x,y,z,dx,dy,dz,sx,sy,sz,omega,theta,phi,name,boxSize = utils.getObjInfo(parent)
            maxSize = Vector(boxSize).length
            cutPlane.scale = Vector((1.01*maxSize,1.01*maxSize,1.0))
            cutPlane.location = Vector((x,y,z))

            # Set the plane as active object
            bpy.ops.object.select_all(action='DESELECT')
            cutPlane.select_set(True)
            bpy.context.view_layer.objects.active = cutPlane
                
        return {'FINISHED'}

class QUADRIC_OT_remove_cutplane(Operator):
    """Remove the Cutting Plane Quadric"""
    bl_idname = "quadric.remove_cut_plane"
    bl_label = "Remove Cutting Plane Quadric"
    bl_options = {'REGISTER', 'UNDO'}

    modName: bpy.props.StringProperty(name = "Modifier Name", default="")

    def execute(self, context):

        parent = context.object
        if parent and parent.penred_settings and self.modName and parent.modifiers:

            modifier = parent.modifiers.get(self.modName)

            if modifier and modifier.type == "BOOLEAN":

                # Get modifier object and remoth both, the modifier and the object
                bolObj = modifier.object

                # Remove modifier
                parent.modifiers.remove(modifier)

                # Remove object
                if bolObj:
                    bpy.data.objects.remove(bolObj)

        return {'FINISHED'}

quadricObjectAddOperators = (
    QUADRIC_OT_add_cube,
    QUADRIC_OT_add_sphere,
    QUADRIC_OT_add_cone,
    QUADRIC_OT_add_cylinder,
    QUADRIC_OT_add_plane,
    QUADRIC_OT_add_semiSphere,
    QUADRIC_OT_add_cutplane,
    QUADRIC_OT_remove_cutplane,
)

# Simulation operators
class SIMULATE_PENRED_OT_cancel(bpy.types.Operator):
    bl_idname = "world.cancel_penred_simulation"
    bl_label = "Cancel PenRed simulation"

    def execute(self, context):
        world = context.scene.world
        if world and world.penred_settings:
            if world.penred_settings.simulation.simulationState == "RUNNING":
                world.penred_settings.simulation.simulationState = "CANCELLED"
        return {'FINISHED'}

class SIMULATE_PENRED_OT_run(bpy.types.Operator):
    """Run a simulaiton using penred code"""
    bl_idname = "world.simulate_penred"
    bl_label = "Run PenRed simulation"
    bl_options = {'REGISTER', 'UNDO'}

    _timer = None
    _simu = None
    _ticks = 0
    _progress = 0
    _fade_alpha = 0.0
    _fade_in = True
    _draw_handler = None

    # Define the progress popup draw as a static method
    @staticmethod
    def draw_progress(self, context):
        region = context.region
        width, height = region.width, region.height

        # Get alpha
        alpha = getattr(self, "_fade_alpha", 1.0)

        # UI box size (relative to screen)
        bar_width = width // 3
        bar_height = max(30, bar_width // 10)
        border_size = max(1, bar_height // 10)
        corner_radius = 8

        x = (width - bar_width) // 2
        y = height // 5

        gpu.state.blend_set('ALPHA')

        # Colors (RGBA)
        bg_color = (0.0, 0.0, 0.0, 0.5 * alpha)
        fg_color = (0.2, 0.6, 1.0, 0.8 * alpha)
        border_color = (1, 1, 1, 0.2 * alpha)
        text_color = (1, 1, 1, alpha)
        
        # Draw background box
        shader = gpu.shader.from_builtin('UNIFORM_COLOR')


        def draw_rounded_rect(x, y, w, h, r, color):

            segments = 12  # More segments = smoother corners

            # Corner centers
            cx = x + r
            cy = y + r
            corners = [
                (cx, cy),                 # Bottom-left
                (x + w - r, cy),          # Bottom-right
                (x + w - r, y + h - r),   # Top-right
                (cx, y + h - r)           # Top-left
            ]

            # Angles for each corner (start angle, end angle)
            angles = [
                (pi, 1.5 * pi),
                (1.5 * pi, 2 * pi),
                (0, 0.5 * pi),
                (0.5 * pi, pi)
            ]

            verts = []
            indices = []

            # Center body quad (not covering corners)
            body_verts = [
                (x + r, y),
                (x + w - r, y),
                (x + w - r, y + h),
                (x + r, y + h),
            ]
            base_idx = len(verts)
            verts.extend(body_verts)
            indices.extend([
                (base_idx, base_idx + 1, base_idx + 2),
                (base_idx, base_idx + 2, base_idx + 3)
            ])

            # Side rects (excluding corners)
            side_verts = [
                (x, y + r), (x + r, y + r), (x + r, y + h - r), (x, y + h - r),  # Left
                (x + w - r, y + r), (x + w, y + r), (x + w, y + h - r), (x + w - r, y + h - r)  # Right
            ]
            side_indices = [
                (4, 5, 6), (4, 6, 7),  # Right side
                (0, 1, 2), (0, 2, 3)   # Left side
            ]
            for i, (v0, v1, v2) in enumerate(side_indices):
                base = len(verts)
                verts.extend(side_verts)
                indices.append((base + v0, base + v1, base + v2))

            # Top and bottom bars (excluding corners)
            top_verts = [
                (x + r, y + h - r), (x + w - r, y + h - r),
                (x + w - r, y + h), (x + r, y + h)
            ]
            bottom_verts = [
                (x + r, y), (x + w - r, y),
                (x + w - r, y + r), (x + r, y + r)
            ]
            for verts_set in [top_verts, bottom_verts]:
                base = len(verts)
                verts.extend(verts_set)
                indices.extend([
                    (base, base + 1, base + 2),
                    (base, base + 2, base + 3)
                ])

            # Draw each rounded corner
            for i, (cx, cy) in enumerate(corners):
                start_angle, end_angle = angles[i]
                arc = []
                for j in range(segments + 1):
                    t = start_angle + (end_angle - start_angle) * (j / segments)
                    arc.append((cx + r * cos(t), cy + r * sin(t)))

                base = len(verts)
                verts.append((cx, cy))  # Center of the fan
                verts.extend(arc)
                for j in range(len(arc)):
                    indices.append((base, base + j + 1, base + ((j + 1) % len(arc)) + 1))

            # Draw the batch
            batch = batch_for_shader(shader, 'TRIS', {"pos": verts}, indices=indices)
            shader.bind()
            shader.uniform_float("color", color)
            batch.draw(shader)


        def draw_rect(x, y, w, h, color):
            vertices = [
                (x,     y),
                (x + w, y),
                (x + w, y + h),
                (x,     y + h)
            ]
            batch = batch_for_shader(shader, 'TRI_FAN', {"pos": vertices})
            shader.bind()
            shader.uniform_float("color", color)
            batch.draw(shader)

        # Border
        draw_rounded_rect(x - border_size, y - border_size, bar_width + 2*border_size, bar_height + 2*border_size, corner_radius+border_size, border_color)

        # Background
        draw_rounded_rect(x, y, bar_width, bar_height, corner_radius, bg_color)
        
        # Foreground (progress fill)
        if self._progress > 2.5:
            progress_width = int(bar_width * (self._progress / 100.0))
            draw_rounded_rect(x, y, progress_width, bar_height, corner_radius, fg_color)

        # --- Text overlay ---
        percent_text = f"Simulation progress: {self._progress:.1f}%"
        font_id = 0
        blf.size(font_id, 16)
        text_width, text_height = blf.dimensions(font_id, percent_text)

        text_x = x + (bar_width - text_width) / 2
        text_y = y + (bar_height - text_height) / 2

        blf.position(font_id, text_x, text_y, 0)
        blf.color(font_id, 1.0, 1.0, 1.0, 1.0)
        blf.draw(font_id, percent_text)

        gpu.state.blend_set('NONE')

    def invoke(self, context, event):

        world = context.scene.world            
        if not (world and world.penred_settings):
            self.report({'ERROR'}, "Missing PenRed settings in World")
            return {'CANCELLED'}

        # Check if a previous simulation is on configuration or running:
        if world.penred_settings.simulation.simulationState != "NONE":
            self.report({'WARNING'}, "PenRed simulation already on configuration or running")
            return {'CANCELLED'}

        self.cleanSimu(context)
        
        # Set simulation to exporting state
        world.penred_settings.simulation.simulationState = "EXPORTING"

        # Invoke penred export UI
        bpy.ops.export_penred.data('INVOKE_DEFAULT',
                                   calledToSimulate=True)

        # Register modal function in window manager with no timer
        context.window_manager.modal_handler_add(self)

        # Trigger modal function
        return {'RUNNING_MODAL'}
    
    def modal(self, context, event):

        # Get simulation state
        world = context.scene.world            
        if not (world and world.penred_settings):
            self.report({'ERROR'}, "Missing PenRed settings in World")
            return self.cancel(context)

        # If no time is defined, the export is still running
        if not self._timer:
            
            # Check if export operator is still running
            if world.penred_settings.simulation.simulationState == "EXPORTING":
                # Still exporting
                return {'PASS_THROUGH'}  # Let export operator handle events

            elif world.penred_settings.simulation.simulationState == "CANCELLED":
                # Export cancelled
                return self.cancel(context)
                
            elif world.penred_settings.simulation.simulationState == "EXPORTED":

                # Files exported, set a timer
                self._timer = context.window_manager.event_timer_add(
                    0.1, # Ensure frequent calls to update the UI
                    window=context.window
                )
                
                # Start simulation
                self.setup_simulation(context)
                
                # Init progress
                self._progress = 0.0

                # Init fade
                self._fade_alpha = 0.0
                self._fade_in = True

                # Create the popup
                self._draw_handler = bpy.types.SpaceView3D.draw_handler_add(
                    self.__class__.draw_progress, (self, context), 'WINDOW', 'POST_PIXEL'
                )
                
                return {'RUNNING_MODAL'}
                        
        elif event.type == 'TIMER':

            if self._fade_in:
                self._fade_alpha = self._fade_alpha + 0.02
                if self._fade_alpha >= 1.0:
                    self._fade_in = False
            else:
                self._fade_alpha = self._fade_alpha - 0.02
                if self._fade_alpha <= 0.6:
                    self._fade_in = True

            # Check if the simulation should be cancelled
            if world.penred_settings.simulation.simulationState == "CANCELLED":
                # simulation cancelled
                return self.cancel(context)


            # Increase counter
            self._ticks = self._ticks + 1

            if self._ticks % 100 == 0: # Limit status requests rate
                self.update_progress(context)

                # Check if the simulation is still running
                if self._simu.isSimulating():
                    return {'RUNNING_MODAL'}
                else:
                    # Simulation finished
                    return self.finish(context)

            # Force redraw of windows in the viewport
            for window in context.window_manager.windows:
                screen = window.screen
                for area in screen.areas:
                    if area.type == 'VIEW_3D':
                        for region in area.regions:
                            if region.type == 'WINDOW':
                                region.tag_redraw()
            
        return {'PASS_THROUGH'}

    def setup_simulation(self, context):

        world = context.scene.world
        if not (world and world.penred_settings):
            self.report({'ERROR'}, "Missing PenRed settings in World")
            return self.cancel(context)

        # Try importing pyPenred and installing it if its not installed in
        # the blender environment
        try:
            import pyPenred
        except:
            try:
                self.report({'WARNING'}, f"pyPenred is not installed, trying to install it...")
                import sys
                import subprocess
                subprocess.call([sys.executable, "-m", "pip", "install", "pyPenred"])
                import pyPenred
            except Exception as e:            
                self.report({'ERROR'}, f"Setup failed: {str(e)}. Unable to install pyPenred. Please, install it manually in blender environment")
                return self.cancel(context)        
        try:
            
            paths = os.path.split(world.penred_settings.simulation.simulationConfigPath)

            # Verify export files
            if not os.path.exists(world.penred_settings.simulation.simulationConfigPath):
                self.report({'ERROR'}, f"Exported files ({world.penred_settings.simulation.simulationConfigPath}) not found!")
                return self.cancel(context)
                
            os.chdir(paths[0])
            
            pyPenred.simulation.setConfigurationLog("config.log")
            pyPenred.simulation.setSimulationLog("simulation.log")
            self._simu = pyPenred.simulation.create()
            
            if self._simu.configFromFile(paths[1]) != 0:
                self.report({'ERROR'}, "Invalid config file format. Please, report this error")
                return self.cancel(context)
                
            if self._simu.simulate(True) != 0:  # Async mode
                self.report({'ERROR'}, "Simulation failed to start. See config.log")
                return self.cancel(context)

            # Change the state to running
            world.penred_settings.simulation.simulationState = "RUNNING"
            
            return {'RUNNING_MODAL'}
                
        except Exception as e:

            self.report({'ERROR'}, f"Simulation setup failed: {str(e)}")
            return self.cancel(context)

    def update_progress(self, context):

        if not self._simu:
            return

        try:
            simulated = self._simu.simulated()
            self._progress = min((s[0]/s[1]*100.0 for s in simulated if s[1] > 0), default=0)
                
            self.report({'INFO'}, f"Simulation Progress: {self._progress:.2f}%")

        except Exception as e:
            self.report({'WARNING'}, f"Progress update failed: {str(e)}")
        
    def cleanSimu(self, context):
        if self._simu:
            if self._simu.isSimulating():
                self.report({'INFO'}, "Forcing simulation finish")
                self._simu.forceFinish()
                while self._simu.isSimulating(): # Wait until the simulation ends
                    time.sleep(5)
                    
            del self._simu
            _simu = None

        world = context.scene.world            
        if world and world.penred_settings:
            world.penred_settings.simulation.simulationState = "NONE"

        # Clean the popup
        if self._draw_handler:
            bpy.types.SpaceView3D.draw_handler_remove(self._draw_handler, 'WINDOW')
            self._draw_handler = None
            
    def finish(self, context):

        if self._timer:
            context.window_manager.event_timer_remove(self._timer)
            self._timer = None

        self.cleanSimu(context)
        
        self.report({'INFO'}, "Simulation finished")

        # Force redraw of windows in the viewport
        for window in context.window_manager.windows:
            screen = window.screen
            for area in screen.areas:
                if area.type == 'VIEW_3D':
                    for region in area.regions:
                        if region.type == 'WINDOW':
                            region.tag_redraw()        
        return {'FINISHED'}

    def cancel(self, context):

        if self._timer:
            context.window_manager.event_timer_remove(self._timer)
            self._timer = None
        
        self.cleanSimu(context)

        self.report({'WARNING'}, "Simulation cancelled")

        # Force redraw of windows in the viewport
        for window in context.window_manager.windows:
            screen = window.screen
            for area in screen.areas:
                if area.type == 'VIEW_3D':
                    for region in area.regions:
                        if region.type == 'WINDOW':
                            region.tag_redraw()        
        return {'CANCELLED'}
    
# Export operator
class export_penred(Operator, ExportHelper):
    """Export operator"""
    bl_idname = "export_penred.data"
    bl_label = "Export to PenRed"
    
    filename_ext = ".geo"

    # Add a flag to check if called to perform a simulation
    calledToSimulate: bpy.props.BoolProperty(
        name="Called To Simulate",
        description="Set to True if invoked programmatically to perform a simulation",
        default=False,
        options={'HIDDEN'}  # Hide from the UI
    )
    
    filter_glob: bpy.props.StringProperty(
        default = "*.quad;*.mesh;*.geo",
        options={'HIDDEN'},
        maxlen=255,
        )
        
    toRound: bpy.props.IntProperty(
        name="Rounding decimal",
        description="Decimal limit to round the export values (position, dimensions, etc)",
        default=5,
        min=3,
        max=10
        )

    
    exportType: bpy.props.EnumProperty(
        items=[("MESH","Mesh","Mesh based geometry",'',0), ("QUADRICS","Quadrics","Quadric based geometry",'',1)],
        name="Geometry type",
        description="PenRed geometry type to export",
        default=0,
        )

    onlyGeo: bpy.props.BoolProperty(
        name="Only geometry",
        description="Exports only the geometry, but not the configuration",
        default=False,
        )

    onlyActive: bpy.props.BoolProperty(
        name="Only active",
        description="Exports only the active object",
        default=False,
        )

    avoidHide: bpy.props.BoolProperty(
        name="Omit hide",
        description="Disable the export of hide objects and their children",
        default=False,
        )
    
    def getChildrens(self,context,parent):
        children = []
        for obj in context.scene.objects:
            if obj.parent == parent:
                children.append(obj)
        return children
    
    def getMeshParent(self,context,obj):
        parent = obj.parent
        while True:
            if parent:
                if parent.type == 'MESH':
                    return parent
                else:
                    parent = parent.parent
            else:
                return None
    
    def createTriangleMesh(self,f,context,obj,toRound,forceWorld,avoidHide,fconf):

        if obj.type != 'MESH':
            return
        
        #Check if hide objects must be avoided
        if avoidHide:
            if obj.hide_get():
                #Skip it
                return

        # Set object specific parameters
        if fconf:
            if obj.penred_settings.isDetector:
                fconf.write(f"geometry/kdet/{obj.name} {obj.penred_settings.detector}\n")
            if obj.penred_settings.dsmaxEnabled:
                fconf.write(f"geometry/dsmax/{obj.name} {obj.penred_settings.dsmax:.5e}\n")

            # Create variance reduction
            conf.createVR(obj, obj.name, fconf)
                
        #Get name
        name = obj.name
        name.replace(" ","_")
        if len(name) > 100:
            name = name[:100]
        
        #Get parent name
        if forceWorld:
            parentName = "void"
        else:
            parent = self.getMeshParent(context,obj)
            if parent:
                parentName = parent.name
                parentName.replace(" ","_")
                if len(parentName) > 100:
                    parentName = parentName[:100]
            else:
                parentName = "void"
        
        #Get number of vertex groups
        vgNames = obj.vertex_groups.keys()
        nVG = len(vgNames)

        #Get mesh data
        mesh = obj.data
        #Create triangles
        mesh.calc_loop_triangles()

        #Create lists with vertex belonging to each vertex group
        vgIndexLists = []
        for i in range(nVG):
            vgIndexLists.append([])

        for v in mesh.vertices:
            for g in v.groups:
                vgIndexLists[g.group].append(v.index)

        #Get the number of effective vertex groups (with at least 1 vertex)
        nVGEff = 0
        for i in range(nVG):
            if len(vgIndexLists[i]) > 0:
                nVGEff = nVGEff + 1
        
        f.write("# Object: %s\n" % name.replace(" ", "_"))
        f.write("#MAT      #NFACES     #NVERTEX     #NAME        #PARENT NAME    #N VERTEX GROUPS\n")
        f.write(" %03d      %07d     %08d     %s        %s   %04d\n" % (obj.penred_settings.material, len(mesh.loop_triangles),len(mesh.vertices), name.replace(" ", "_"), parentName.replace(" ", "_"), nVGEff))

        #Print vertex groups
        f.write("# VERTEX GROUPS\n")
        for i in range(nVG):

            #Skip empty groups
            if len(vgIndexLists[i]) == 0:
                continue

            f.write("#NAME  #NVERTEX\n")
            f.write(" %s   %04d\n" % (vgNames[i].replace(" ", "_"), len(vgIndexLists[i])))
            for index in vgIndexLists[i]:
                f.write(" %04d\n" % (index))

        #Print vertex
        f.write("# VERTEX LIST\n")
        f.write("# Index  (X Y Z)\n")
        
        for vertex in mesh.vertices:
            vertexWorld = obj.matrix_world @ vertex.co
            f.write("%04d %+.*E %+.*E %+.*E\n" % (vertex.index, toRound+3, round(vertexWorld[0],toRound), toRound+3, round(vertexWorld[1],toRound), toRound+3, round(vertexWorld[2],toRound)))
        
        f.write("# FACES(triangles)\n")
        for tri in mesh.loop_triangles:
            f.write(" %03d %03d %03d\n" % (tri.vertices[0], tri.vertices[1], tri.vertices[2]))
        f.write("#\n#\n")
        
    
    def createObject(self,f,context,obj,nSurf,nObj,toRound,createChilds,avoidHide,fconf):
        
        #Check the quadric type
        if obj.penred_settings.quadricType == "unknown" and obj.type != "EMPTY": 
            return [], nSurf, nObj

        # Avoid cut planes
        if obj.penred_settings.quadricType == "CUT_PLANE":
            return [], nSurf, nObj

        #Check if hide objects must be avoided
        if avoidHide:
            if obj.hide_get():
                #Skip it
                return [], nSurf, nObj

        #Get childrens
        childrens = []
        if createChilds:
            childrens = self.getChildrens(context,obj)
        
        #First, construct children objects
        tree = [] # Children tree information
        if len(childrens) > 0:
            for child in childrens:
                childTree, nSurf, nObj = self.createObject(f,context,child,nSurf,nObj,toRound,True,avoidHide,fconf)
                if len(childTree) > 0:
                    tree.extend(childTree)
                
        #If the object is an empty, no body must be created here
        if obj.type == "EMPTY":
            return tree,nSurf,nObj
        
        #Get object type
        quadType = obj.penred_settings.quadricType
                    
        #Get object properties
        x,y,z,dx,dy,dz,sx,sy,sz,omega,theta,phi,name,boxSize = utils.getObjInfo(obj)
                         
        #Save init surface
        initSurf = nSurf
        
        ### Create object surfaces
        if quadType == "CUBE":
            nSurf = surfaces.createCubeSurfaces(f,x,y,z,dx,dy,dz,omega,theta,phi,nSurf,name,toRound)
        elif quadType == "SPHERE":
            nSurf = surfaces.createSphereSurfaces(f,x,y,z,dx,dy,dz,nSurf,name,toRound)
        elif quadType == "CYLINDER":
            nSurf = surfaces.createCylinderSurfaces(f,x,y,z, dx,dy,dz,omega,theta,phi,nSurf,name,toRound)
        elif quadType == "CONE":
            r1 = obj.penred_settings.r1
            r2 = obj.penred_settings.r2
            if r1 == r2:
                nSurf = surfaces.createCylinderSurfaces(f,x,y,z, dx,dy,dz,omega,theta,phi,nSurf,name,toRound)
            else:
                nSurf = surfaces.createConeSurfaces(f,x,y,z,r1,r2,dz,sx,sy,omega,theta,phi,nSurf,name,toRound)
        elif quadType == "PLANE":
            nSurf = surfaces.createPlaneSurfaces(f,x,y,z,omega,theta,phi,nSurf,name,toRound)
        elif quadType == "SEMI_SPHERE":
            #Construct sphere. Take into account that the Z dimension
            #must not be multiplied by two to comepnsate the cut, because
            #we added an artificial vertex
            nSurf = surfaces.createSphereSurfaces(f,x,y,z,dx,dy,dz,nSurf,name,toRound)
            #Create limiting plane
            nSurf = surfaces.createPlaneSurfaces(f,x,y,z,omega,theta,phi,nSurf,name,toRound)

        # Create cutting plane surfaces, if defined
        nCuttingPlanes = 0
        for mod in obj.modifiers:
            if mod.type == "BOOLEAN":
                bolObj = mod.object
                if bolObj and hasattr(bolObj, "penred_settings"):
                    if bolObj.penred_settings.quadricType == "CUT_PLANE":
                        # Is a cutting plane, get its properties
                        xCut,yCut,zCut,_,_,_,_,_,_,omegaCut,thetaCut,phiCut,nameCut,_ = utils.getObjInfo(bolObj)
                        
                        nSurf = surfaces.createPlaneSurfaces(f,xCut,yCut,zCut,omegaCut,thetaCut,phiCut,nSurf,nameCut,toRound)
                        nCuttingPlanes = nCuttingPlanes+1
            
        ### Init body
        surfaces.initBody(f,nObj,name,obj.penred_settings.material,obj.penred_settings.module)

        # Set object specific parameters
        if fconf:
            if obj.penred_settings.isDetector:
                fconf.write(f"geometry/kdet/{nObj} {obj.penred_settings.detector}\n")
            if obj.penred_settings.dsmaxEnabled:
                fconf.write(f"geometry/dsmax/{nObj} {obj.penred_settings.dsmax:.5e}\n")

            # Create variance reduction
            conf.createVR(obj, nObj, fconf)
                
                
        ### Set object surfaces
        if quadType == "CUBE":
            initSurf = surfaces.setCubeSurfaces(f,initSurf, 1)
        elif quadType == "SPHERE":
            initSurf = surfaces.setSphereSurfaces(f,initSurf, 1)
        elif quadType == "CYLINDER" or quadType == "CONE":
            initSurf = surfaces.setCylinderConeSurfaces(f,initSurf, 1)
        elif quadType == "PLANE":
            initSurf = surfaces.setPlaneSurfaces(f,initSurf, 1)
        elif quadType == "SEMI_SPHERE":
            initSurf = surfaces.setSphereSurfaces(f,initSurf, 1)
            initSurf = surfaces.setPlaneSurfaces(f,initSurf, 1)

        # Set cutting planes surfaces
        for i in range(nCuttingPlanes):
            initSurf = surfaces.setPlaneSurfaces(f,initSurf, 1)            
            
        #Set childrens
        if len(tree) > 0:
            surfaces.addChilds(f,tree)
        
        #Add current body or module to tree
        if obj.penred_settings.module:
            #Module childs can't be repeated by other bodies/modules
            tree = [(nObj,True,name)] 
        else:
            tree.append((nObj,False,name))
        nObj += 1
            
        return tree,nSurf,nObj
                 
    def countNonHideMeshes(self,context,obj):

        #Number of non hide meshes
        nMeshes = 0

        #Check this mesh status
        if not obj.hide_get():
            #Count this object if it is a mesh
            if obj.type == 'MESH':
                nMeshes = 1

            #Check childrens
            childrens = self.getChildrens(context,obj)
            for child in childrens:
                #Add contribution of all children
                nMeshes += self.countNonHideMeshes(context,child)

        #Return the number of non hide meshes
        return nMeshes

    def invoke(self, context, event):
        return super().invoke(context, event)

    def cancel(self, context):
        if self.calledToSimulate:
            world = context.scene.world
            if world and world.penred_settings:
                world.penred_settings.simulation.simulationState = "CANCELLED"
        return {'CANCELLED'}    
    
    def execute(self, context):
        
        #Open output file
        f = open(self.filepath,'w',encoding='utf-8')
        f.write("# Geometry file created with PenRed blender plugin v.2.0\n")

        #If required, open output configuration file
        fconf = None
        if not self.onlyGeo:
            filenameConf = os.path.splitext(self.filepath)[0] + ".in"
            fconf = open(filenameConf,'w',encoding='utf-8')
            fconf.write("# Configuration file created with PenRed blender plugin v.2.0\n\n")

            # Create simulation configuration
            conf.createSim(context, fconf)

            # Create materials configuration
            conf.createMaterials(context, fconf)

            # Create source configuration
            conf.createSources(context, fconf, self.toRound)

            # Create tallies configuration
            conf.createTallies(context, fconf, self.toRound)

            # Create geometry configuration
            fconf.write("\n## Geometry configuration ##\n")
            if self.exportType == 'QUADRICS':
                fconf.write( "geometry/type \"PEN_QUADRIC\"\n")
                fconf.write(f"geometry/input-file \"{self.filepath}\"\n")
                fconf.write("geometry/processed-geo-file \"report.geo\"\n\n")
            else:
                fconf.write( "geometry/type \"MESH_BODY\"\n")
                fconf.write(f"geometry/input-file \"{self.filepath}\"\n")                

        if self.exportType == 'QUADRICS':
            #Create an array for object names
            objNames = []
            
            nSurf = 1 #Number of surface to be created
            nObj = 1 #Number of object to be created

            if self.onlyActive:
                self.createObject(f,context,bpy.context.active_object,nSurf,nObj,self.toRound,False,False,fconf)
            else:
                #Find objects with no parents
                for obj in context.scene.objects:
                    if not obj.parent:
                        #This object has no parent, create it
                        nSurf,nObj = self.createObject(f,context,obj,nSurf,nObj,self.toRound,True,self.avoidHide,fconf)[1:]

            surfaces.endFile(f)
        elif self.exportType == 'MESH':
            
            if self.onlyActive:
                #Print number of objects
                f.write("# Number of objects:\n 1\n")

                self.createTriangleMesh(f,context,bpy.context.active_object,self.toRound,True,False,fconf)
            else:

                #Count number of meshes
                nMeshes = 0

                if self.avoidHide:
                    for obj in context.scene.objects:
                        if not obj.parent: #Begin from objects with no parent
                            nMeshes += self.countNonHideMeshes(context,obj)
                else:
                    #Count all mesh objects
                    for obj in context.scene.objects:
                        if obj.type == 'MESH':
                            nMeshes = nMeshes + 1

                #Print number of objects
                f.write("# Number of objects:\n %d\n" % (nMeshes))

                for obj in context.scene.objects:
                    self.createTriangleMesh(f,context,obj,self.toRound,False,self.avoidHide,fconf)
        else:
            f.write("# Unknown export format. Please, report this issue\n")
                
        f.close()

        if self.calledToSimulate:
            world = context.scene.world
            if world and world.penred_settings:
                # Save simulation path
                world.penred_settings.simulation.simulationConfigPath = os.path.splitext(self.filepath)[0] + ".in"
                # Flag the simulation state to exported
                world.penred_settings.simulation.simulationState = "EXPORTED"
                
        return {'FINISHED'}

def register():

    #Register View classes
    for cls in viewObjectClasses:
        bpy.utils.register_class(cls)
    
    #Register transformation classes
    for cls in quadricTransformClasses:
        bpy.utils.register_class(cls)

    #Register object operators
    for cls in objectOperatorsClasses:
        bpy.utils.register_class(cls)
        
    #Register object add operators
    for cls in quadricObjectAddOperators:
        bpy.utils.register_class(cls)

    #Register tally operators
    for cls in talliesOperatorClasses:
        bpy.utils.register_class(cls)

    #Register material operators
    for cls in materialsOperatorClasses:
        bpy.utils.register_class(cls)

    bpy.utils.register_class(export_penred)
    
    #Register simulation operators
    bpy.utils.register_class(SIMULATE_PENRED_OT_run)
    bpy.utils.register_class(SIMULATE_PENRED_OT_cancel)
def unregister():

    #Unregister view classes
    for cls in viewObjectClasses:
        bpy.utils.unregister_class(cls)
    
    #Unregister transformation classes
    for cls in quadricTransformClasses:
        bpy.utils.unregister_class(cls)

    #Unregister object operators
    for cls in objectOperatorsClasses:
        bpy.utils.unregister_class(cls)
    
    #Unregister object add operators
    for cls in quadricObjectAddOperators:
        bpy.utils.unregister_class(cls)

    #Unregister tally operators
    for cls in talliesOperatorClasses:
        bpy.utils.unregister_class(cls)
        
    #Unregister material operators
    for cls in materialsOperatorClasses:
        bpy.utils.unregister_class(cls)

    bpy.utils.unregister_class(export_penred)

    #Unregister simulation operator
    bpy.utils.unregister_class(SIMULATE_PENRED_OT_run)
    bpy.utils.unregister_class(SIMULATE_PENRED_OT_cancel)
