
import bpy
import bmesh
from bpy_extras.io_utils import ExportHelper
from bpy.types import Operator
from bpy.props import FloatVectorProperty
from bpy_extras.object_utils import AddObjectHelper, object_data_add
from mathutils import Vector
from mathutils import Color
from math import cos, acos, sin, asin, tan, atan2, sqrt, pi
import os

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

# Export operator
class export_penred(Operator, ExportHelper):
    """Export operator"""
    bl_idname = "export_penred.data"
    bl_label = "Export to PenRed"
    
    filename_ext = ".geo"
    
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
            r1 = obj.penred_settings.r3
            r2 = obj.penred_settings.r4
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
