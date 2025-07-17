
import bpy
from . import penelopeMatDB
from . import ICRP_AF_DB
from . import ICRP_AM_DB

class PenRedMaterial(bpy.types.PropertyGroup):
    name: bpy.props.StringProperty(name="Material Name")

class PenRedMaterialCategory(bpy.types.PropertyGroup):
    name: bpy.props.StringProperty(name="Category Name")
    materials: bpy.props.CollectionProperty(type=PenRedMaterial)

class PenRedMaterialDB(bpy.types.PropertyGroup):
    categories: bpy.props.CollectionProperty(type=PenRedMaterialCategory)

categories = [
    ("PENELOPE" , "Penelope", "Predefined PENELOPE materials list"),
    ("ICRP_AF", "ICRP AF", "ICRP materials for adult female phantom"),
    ("ICRP_AM", "ICRP AM", "ICRP materials for adult male phantom"),
]
def getCategoryIndex(name : str):
    i = 0
    for c in categories:
        if name == c[0]:
            return i
        i = i + 1
    return -1 # not found

def getCategory(name : str):
    iDB = getCategoryIndex(name)
    if iDB < len(categories):
        return bpy.context.scene.penred_material_db.categories[iDB]
    return None

def getMatName(categoryName : str, index : int):
    category = getCategory(categoryName)
    if category:
        return category.materials[index].name
    return ""
    
def initDB(scene):
    db = scene.penred_material_db

    db.categories.clear()

    # Add PENELOPE category
    penelope = db.categories.add()
    penelope.name = "PENELOPE"

    # Add materials to the penelope category
    for mat in penelopeMatDB.materials:
        penelope.materials.add().name = mat

    # Add ICRP categories

    # Adult Female
    icrp_af = db.categories.add()
    icrp_af.name = "ICRP AF"

    for mat in ICRP_AF_DB.materials:
        icrp_af.materials.add().name = mat

    # Adult Male
    icrp_am = db.categories.add()
    icrp_am.name = "ICRP AM"

    for mat in ICRP_AM_DB.materials:
        icrp_am.materials.add().name = mat

# --- Handler Registration ---
@bpy.app.handlers.persistent  # Survives file loads
def on_blender_loaded(dummy):
    if hasattr(bpy.types.Scene, "penred_material_db"):
        print("Loading materials DB")
        initDB(bpy.context.scene)
    else:
        print("Loading materials DB failed. No scene provided")
        
def register():
    bpy.utils.register_class(PenRedMaterial)
    bpy.utils.register_class(PenRedMaterialCategory)
    bpy.utils.register_class(PenRedMaterialDB)

    bpy.types.Scene.penred_material_db = bpy.props.PointerProperty(type=PenRedMaterialDB)

    # Add handler (runs after file load/startup)
    bpy.app.handlers.load_post.append(on_blender_loaded)
        
def unregister():

    bpy.app.handlers.load_post.remove(on_blender_loaded)
    
    del bpy.types.Scene.penred_material_db
    
    bpy.utils.unregister_class(PenRedMaterialCategory)
    bpy.utils.unregister_class(PenRedMaterialDB)
    bpy.utils.unregister_class(PenRedMaterial)
