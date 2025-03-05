
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

def getObjPosSize(obj):

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

    return center,dx,dy,dz,sx,sy,sz,size
    
def getObjInfo(obj):
    #Set object origin to origin geometry
    #context.view_layer.objects.active = obj
    #bpy.ops.object.origin_set( type = 'ORIGIN_GEOMETRY' )

    #Get name
    name = obj.name
        
    center,dx,dy,dz,sx,sy,sz,size = getObjPosSize(obj)
                            
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

# Draw a object bounding box
def draw_bbox(obj, color):

    # Get object spatial properties
    center,dx,dy,dz,sx,sy,sz,bsize = getObjPosSize(obj)

    width05 = bsize[0]/2.0
    height05 = bsize[1]/2.0
    depth05 = bsize[2]/2.0

    xmin = center[0] - width05
    xmax = center[0] + width05

    ymin = center[1] - height05
    ymax = center[1] + height05

    zmin = center[2] - depth05
    zmax = center[2] + depth05
    
    # Define the 8 vertices of the box
    vertices = [
        Vector(( xmin, ymin, zmin)),
        Vector(( xmax, ymin, zmin)),
        Vector(( xmax, ymax, zmin)),
        Vector(( xmin, ymax, zmin)),
        Vector(( xmin, ymin, zmax)),
        Vector(( xmax, ymin, zmax)),
        Vector(( xmax, ymax, zmax)),
        Vector(( xmin, ymax, zmax)),
    ]

    # Define the 12 edges (lines) of the box
    indices = [
        # Bottom face
        (0, 1), (1, 2), (2, 3), (3, 0),
        # Top face
        (4, 5), (5, 6), (6, 7), (7, 4),
        # Vertical edges
        (0, 4), (1, 5), (2, 6), (3, 7),
    ]

    # Create a shader for 3D drawing
    shader = gpu.shader.from_builtin('UNIFORM_COLOR')
    batch = batch_for_shader(shader, 'LINES', {"pos": vertices}, indices=indices)

    # Set the color
    shader.bind()
    shader.uniform_float("color", color)

    # Draw the box as lines
    batch.draw(shader)    

# Draw a sphere
def draw_sphere(obj, color, segments=16, rings=16):

    # Get object spatial properties
    center,dx,dy,dz,sx,sy,sz,bsize = getObjPosSize(obj)
    
    vertices = []
    indices = []

    radius = sqrt(bsize[0]*bsize[0] + bsize[1]*bsize[1] + bsize[2]*bsize[2])/2.0

    # Generate vertices for the sphere
    for i in range(rings + 1):
        theta = i * pi / rings  # Polar angle (0 to π)
        for j in range(segments):
            phi = j * 2 * pi / segments  # Azimuthal angle (0 to 2π)
            x = center[0] + radius * sin(theta) * cos(phi)
            y = center[1] + radius * sin(theta) * sin(phi)
            z = center[2] + radius * cos(theta)
            vertices.append(Vector((x, y, z)))

    # Generate indices for the sphere's wireframe
    for i in range(rings):
        for j in range(segments):
            # Horizontal lines (rings)
            indices.append((i * segments + j, i * segments + (j + 1) % segments))
            # Vertical lines (segments)
            indices.append((i * segments + j, (i + 1) * segments + j))


    # Create a shader for 3D drawing
    shader = gpu.shader.from_builtin('UNIFORM_COLOR')
    batch = batch_for_shader(shader, 'LINES', {"pos": vertices}, indices=indices)

    # Set the color
    shader.bind()
    shader.uniform_float("color", color)

    # Draw the sphere as lines
    batch.draw(shader)            
    
# Draw a cylinder along Z axis
def draw_zcyl(obj, inbox, color, segments = 32):

    # Get object spatial properties
    center,dx,dy,dz,sx,sy,sz,bsize = getObjPosSize(obj)

    if inbox:
        radius = min((bsize[0], bsize[1]))/2.0
    else:
        radius = sqrt(bsize[0]*bsize[0] + bsize[1]*bsize[1])/2.0

    height05 = bsize[2]/2.0

    vertices = []
    indices = []

    # Create the top and bottom circles
    angleStep = (2 * pi ) / segments
    for i in range(segments):
        angle = angleStep * i
        x = center[0] + radius * cos(angle)
        y = center[1] + radius * sin(angle)

        # Bottom circle
        vertices.append((x, y, center[2]-height05))
        # Top circle
        vertices.append((x, y, center[2]+height05))

    # Create indices to connect the circles
    for i in range(segments):
        next_i = (i + 1) % segments
        # Tris index
        # Bottom face (TRIS)
        #indices.append((2 * i, 2 * next_i, 2 * i + 1))
        #indices.append((2 * next_i, 2 * next_i + 1, 2 * i + 1))
        # Side faces (TRIS)
        #indices.append((2 * i, 2 * next_i, 2 * i + 1))
        #indices.append((2 * next_i, 2 * next_i + 1, 2 * i + 1))

        # Lines index
        # Caps (LINES)
        indices.append((2 * i, 2 * next_i))
        indices.append((2 * next_i, 2 * next_i + 1))
        # Side faces (LINES)
        indices.append((2 * i, 2 * i + 1))
        
    # Add vertices for the center of the top and bottom caps
    bottom_center = len(vertices)
    vertices.append((center[0], center[1], center[2]-height05))  # Bottom center
    top_center = len(vertices)
    vertices.append((center[0], center[1], center[2]+height05))  # Top center

    # Create indices for the bottom cap
    for i in range(segments):
        next_i = (i + 1) % segments
        # TRIS
        #indices.append((2 * i, 2 * next_i, bottom_center))
        # LINES
        indices.append((2 * i, bottom_center))

    # Create indices for the top cap
    for i in range(segments):
        next_i = (i + 1) % segments
        # TRIS
        #indices.append((2 * i + 1, 2 * next_i + 1, top_center))
        # LINES
        indices.append((2 * i + 1, top_center))
        
    # Create a shader for 3D drawing
    shader = gpu.shader.from_builtin('UNIFORM_COLOR')
    batch = batch_for_shader(shader, 'LINES', {"pos": vertices}, indices=indices)

    # Set the color
    shader.bind()
    shader.uniform_float("color", color)

    # Draw the cylinder
    batch.draw(shader)

# Draw an arrow
def draw_arrow(obj, direction, color):

    # Get object spatial properties
    center,dx,dy,dz,sx,sy,sz,bsize = getObjPosSize(obj)

    # Normalize the direction
    direction = Vector(direction)
    direction = direction.normalized()

    location = Vector(center)

    # Define the arrow geometry
    arrow_length = 1.0 + max(bsize)*1.5
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
    shader.uniform_float("color", color)  # Red color
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

    ### Create an empty mesh
    mesh = bpy.data.meshes.new(name=meshType)

    ### Construct the primary bmesh and assign it to the blender mesh.
    bm = bmesh.new()
    # Notice that, although blender parameters are named "diameter",
    # blender uses radius instead of diameter (bug?) OLD VERSIONS
    if quadType == "CUBE":
        bmesh.ops.create_cube(bm, size=2.0)
    elif quadType == "SPHERE":
        bmesh.ops.create_uvsphere(bm, u_segments=32, v_segments=16, radius=0.5)
    elif quadType == "CONE":
        bmesh.ops.create_cone(bm, cap_ends=True, segments=32, radius1=coneDefaultR1, radius2=coneDefaultR2, depth = 1.0)
    elif quadType == "CYLINDER":
        bmesh.ops.create_cone(bm, cap_ends=True, segments=32, radius1=1.0, radius2=1.0, depth = 2.0)
    elif quadType == "PLANE":
        bmesh.ops.create_circle(bm, cap_ends=True, segments=32, radius=1.0)
    elif quadType == "CUT_PLANE":
        bmesh.ops.create_circle(bm, cap_ends=True, segments=32, radius=1.0)
    elif quadType == "SEMI_SPHERE":
        bmesh.ops.create_uvsphere(bm, u_segments=32, v_segments=16, radius=0.5)

        elementsBot = [e for e in bm.edges if e.verts[1].co[2] < -1.0e-3]
        bmesh.ops.delete(bm,geom=elementsBot,context="EDGES_FACES")
        
        verts = [v for v in bm.verts if v.co[2] < 1.0e-3 and v.co[2] > -1.0e-3]
        bm.faces.new(verts)
        
    bm.to_mesh(mesh)
    bm.free()

    ### Create the base object
    obj = object_data_add(context, mesh, operator=self)
    obj.penred_settings.quadricType=quadType

    ### Set other parameters, if needed        
    if quadType == "CONE":
        #Add defualt radius
        obj.penred_settings.r1 = coneDefaultR1
        obj.penred_settings.r2 = coneDefaultR2
    elif quadType == "CUT_PLANE":
        obj.penred_settings.isMaterialObject = False

    return obj
