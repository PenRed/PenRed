
import bpy
import gpu
from gpu_extras.batch import batch_for_shader
import bmesh
from bpy_extras.io_utils import ExportHelper
from bpy.types import Operator
from bpy.props import FloatVectorProperty
from bpy_extras.object_utils import AddObjectHelper, object_data_add
from mathutils import Vector
from mathutils import Color
from math import cos, acos, sin, asin, tan, atan2, sqrt, pi

def getObjInfo(context,obj):
    #Set object origin to origin geometry
    context.view_layer.objects.active = obj
    #bpy.ops.object.origin_set( type = 'ORIGIN_GEOMETRY' )

    #Get name
    name = obj.name
        
    #Get scale
    sx,sy,sz = obj.matrix_world.to_scale()
        
    if obj.type == 'MESH':

        # Get the bounding box corners in world coordinates
        bounding_box = [obj.matrix_world @ Vector(corner) for corner in obj.bound_box]

        # Calculate the min and max coordinates of the bounding box
        min_coords = [min([v[i] for v in bounding_box]) for i in range(3)]
        max_coords = [max([v[i] for v in bounding_box]) for i in range(3)]

        # Calculate the size of the bounding box
        size = [max_coords[i] - min_coords[i] for i in range(3)]

        # Calculate the center position of the bounding box
        center = [(max_coords[i] + min_coords[i]) / 2.0 for i in range(3)]
        
        #Get sizes
        dx = obj.dimensions.x
        dy = obj.dimensions.y
        dz = obj.dimensions.z
            
    else:
        #Set bounding box size
        size = (sx, sy, sz)
            
        #Get position
        center = obj.matrix_world.to_translation()

        #Set scale and dimensions
        dx = 1.0
        dy = 1.0
        dz = 1.0
                            
    #Get rotation in "ZYZ" euler angles
    rotMatrix = obj.matrix_world.to_quaternion().to_matrix()
                    
    cosTheta = rotMatrix[2][2]
    sinTheta = sqrt(1.0-cosTheta*cosTheta)
    if sinTheta > 1.0e-15:
        theta = acos(cosTheta)
        phi = atan2(rotMatrix[1][2]/sinTheta,rotMatrix[0][2]/sinTheta)
        omega = atan2(rotMatrix[2][1]/sinTheta,-rotMatrix[2][0]/sinTheta)
    else:
        omega = 0.0
        if cosTheta > 0.0:
            theta = 0.0
            phi = atan2(rotMatrix[1][0],rotMatrix[0][0])
        else:
            theta = pi
            phi = atan2(-rotMatrix[1][0],-rotMatrix[0][0])
    #Return values
    return center[0],center[1],center[2],dx,dy,dz,sx,sy,sz,omega,theta,phi,name,size

# Check bounding boxes overlap
def boundingBoxOverlap(obj1, obj2):
    # Get the world-space bounding box corners for both objects
    bbox1 = [obj1.matrix_world @ Vector(corner) for corner in obj1.bound_box]
    bbox2 = [obj2.matrix_world @ Vector(corner) for corner in obj2.bound_box]

    # Check for overlap in all axes
    for axis in range(3):  # 0 = x, 1 = y, 2 = z
        min1 = min(corner[axis] for corner in bbox1)
        max1 = max(corner[axis] for corner in bbox1)
        min2 = min(corner[axis] for corner in bbox2)
        max2 = max(corner[axis] for corner in bbox2)

        if max1 < min2 or min1 > max2:
            return False  # No overlap on this axis
    return True  # Overlap on all axes

# Get overlapping materials
def overlappingMats(sourceObj):
    mats = set()
    for obj in bpy.context.scene.objects:        
        if hasattr(obj, "penred_settings") and obj.type == 'MESH' and boundingBoxOverlap(sourceObj, obj):
            mats.add(obj.penred_settings.material)
    return mats

# Draw an arrow
def draw_arrow(context, obj):
    if not obj.penred_settings.source or not obj.penred_settings.source.enabled:
        return  # Skip drawing if the source flag is disabled

    source = obj.penred_settings.source
    
    # Get the object's location and direction
    location = obj.location
    direction = Vector(source.direction).normalized()

    # Define the arrow geometry
    arrow_length = 1.0 + obj.dimensions.length
    arrow_tip_length = arrow_length/4.0
    arrow_tip_width = arrow_length/10.0

    # Arrow line
    start = location
    end = location + direction * arrow_length

    # Arrow tip (a small triangle)
    tip_base = end - direction * arrow_tip_length
    perpendicular = direction.orthogonal().normalized() * arrow_tip_width
    tip_left = tip_base + perpendicular
    tip_right = tip_base - perpendicular

    # Combine vertices and indices for the arrow
    vertices = [start, end, tip_left, tip_right]
    indices = [(0, 1), (1, 2), (1, 3), (2, 3)]

    # Create a batch for drawing
    shader = gpu.shader.from_builtin('UNIFORM_COLOR')
    batch = batch_for_shader(shader, 'LINES', {"pos": vertices}, indices=indices)

    # Draw the arrow
    shader.bind()
    shader.uniform_float("color", (1.0, 0.0, 0.0, 1.0))  # Red color
    batch.draw(shader)

def getLocalZdir(obj):
    
        nz = Vector((0.0,0.0,1.0))
        inv = obj.matrix_world.copy()
        inv.invert()
        
        dirZ = nz @ inv
        dirZ.normalize()
        return dirZ

def add_object(self, context, meshType, quadType):

    coneDefaultR1 = 1.0
    coneDefaultR2 = 1.0e-5

    coneShellDefaultR1 = 1.0
    coneShellDefaultR2 = 0.8

    ### Create an empty mesh
    mesh = bpy.data.meshes.new(name=meshType)

    ### Construct the primary bmesh and assign it to the blender mesh.
    bm = bmesh.new()
    # Notice that, although blender parameters are named "diameter",
    # blender uses radius instead of diameter (bug?) OLD VERSIONS
    if(quadType == "CUBE"):
        bmesh.ops.create_cube(bm, size=2.0)
    elif(quadType == "SPHERE"):
        bmesh.ops.create_uvsphere(bm, u_segments=32, v_segments=16, radius=0.5)
    elif(quadType == "CONE"):
        bmesh.ops.create_cone(bm, cap_ends=True, segments=32, radius1=coneDefaultR1, radius2=coneDefaultR2, depth = 1.0)
    elif(quadType == "CYLINDER"):
        bmesh.ops.create_cone(bm, cap_ends=True, segments=32, radius1=1.0, radius2=1.0, depth = 2.0)
    elif(quadType == "PLANE"):
        bmesh.ops.create_circle(bm, cap_ends=True, segments=32, radius=1.0)
    elif(quadType == "TUBE"):
        #Create outer cylinder
        bmesh.ops.create_cone(bm, cap_ends=True, segments=32, radius1=1.0, radius2=1.0, depth = 2.0)
    elif(quadType == "SEMI_SPHERE" or quadType == "SEMI_SPHERE_SHELL"):
        bmesh.ops.create_uvsphere(bm, u_segments=32, v_segments=16, radius=0.5)

        elementsBot = [e for e in bm.edges if e.verts[1].co[2] < -1.0e-3]
        bmesh.ops.delete(bm,geom=elementsBot,context="EDGES_FACES")
        
        if quadType == "SEMI_SPHERE":
            verts = [v for v in bm.verts if v.co[2] < 1.0e-3 and v.co[2] > -1.0e-3]
            bm.faces.new(verts)
    elif(quadType == "CONE_SHELL"):
        bmesh.ops.create_cone(bm, cap_ends=False, segments=32, radius1=coneShellDefaultR1, radius2=coneShellDefaultR2, depth = 1.0)        
            
        
    bm.to_mesh(mesh)
    bm.free()

    ### Create the base object
    object_data_add(context, mesh, operator=self)
    context.object.penred_settings.quadricType=quadType

    ### Create secondary objects, if required
    if(quadType == "TUBE"):
        #Create inner cylinder to create the hole
        innerMesh = bpy.data.meshes.new(name="Inner")
        bmInner = bmesh.new()
        bmesh.ops.create_cone(bmInner, cap_ends=True, segments=32, radius1=0.7, radius2=0.7, depth = 2.1)
        bmInner.to_mesh(innerMesh)
        bmInner.free()
        
        innerObj=bpy.data.objects.new(name="Inner", object_data=innerMesh)
        #Set the object in the scene
        bpy.context.collection.objects.link(innerObj)        
        #Set the hole as child of tube object
        innerObj.parent = context.object
        #Make hole invisible
        innerObj.hide_set(1)
        innerObj.show_instancer_for_render = 0
        innerObj.show_instancer_for_viewport = 0
        
        #As active object is the being created (New "TUBE" object), add a boolean modifier
        bpy.ops.object.modifier_add(type="BOOLEAN")
        #Set the operation to "Difference"
        bpy.context.object.modifiers["Boolean"].operation = "DIFFERENCE"
        #Set the "inner" object to do the operation 
        bpy.context.object.modifiers["Boolean"].object = innerObj
        
    elif(quadType == "PLANE"):
        #Add normal line
        normalMesh = bpy.data.meshes.new(name="Normal")
        bmNormal = bmesh.new()
        
        v0 = bmNormal.verts.new((0,0,0))
        v1 = bmNormal.verts.new((0,0,1))
        
        bmNormal.edges.new((v0,v1))
        bmNormal.to_mesh(normalMesh)
        bmNormal.free()
        
        normalObj=bpy.data.objects.new(name="Normal", object_data=normalMesh)
        #Set the normal in the scene
        bpy.context.collection.objects.link(normalObj)        
        #Set the normal object as child of plane object
        normalObj.parent = context.object
    elif(quadType == "CONE"):
        #Add defualt radius
        context.object.penred_settings.r1 = coneDefaultR1
        context.object.penred_settings.r2 = coneDefaultR2
        context.object.penred_settings.r3 = coneDefaultR1
        context.object.penred_settings.r4 = coneDefaultR2
    elif quadType == "SEMI_SPHERE_SHELL" or quadType == "CONE_SHELL":
        #Add solidify modifier
        bpy.ops.object.modifier_add(type="SOLIDIFY")
        #Set the operation to "Difference"
        bpy.context.object.modifiers['Solidify'].solidify_mode = "NON_MANIFOLD"
        bpy.context.object.modifiers['Solidify'].nonmanifold_boundary_mode = "FLAT"
        bpy.context.object.modifiers['Solidify'].thickness = 0.3
        
        if quadType == "CONE_SHELL":
            #Add defualt radius
            context.object.penred_settings.r1 = coneShellDefaultR1
            context.object.penred_settings.r2 = coneShellDefaultR2
            context.object.penred_settings.r3 = coneShellDefaultR1
            context.object.penred_settings.r4 = coneShellDefaultR2
