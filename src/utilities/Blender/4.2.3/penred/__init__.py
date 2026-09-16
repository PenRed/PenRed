#
#
#    Copyright (C) 2022-2025 Universitat de València - UV
#    Copyright (C) 2022-2025 Universitat Politècnica de València - UPV
#    Copyright (C) 2024-2025 Vicent Giménez Alventosa
#
#    This file is part of PenRed: Parallel Engine for Radiation Energy Deposition.
#
#    PenRed is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Affero General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    PenRed is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Affero General Public License for more details.
#
#    You should have received a copy of the GNU Affero General Public License
#    along with PenRed.  If not, see <https://www.gnu.org/licenses/>. 
#
#    contact emails:
#
#        vicent.gimenez.alventosa@gmail.com (Vicent Giménez Alventosa)
#        sanolgi@upvnet.upv.es              (Sandra Oliver Gil)
#        vicente.gimenez@uv.es              (Vicent Giménez Gómez)
#

import bpy
import importlib
from . import addon_properties, operators, ui, penred_import

importlib.reload(penred_import)
importlib.reload(addon_properties)
importlib.reload(operators)
importlib.reload(ui)

# --- Add-on Preferences UI ---

class PYPENRED_Preferences(bpy.types.AddonPreferences):
    bl_idname = __package__

    def draw(self, context):
        layout = self.layout
        layout.label(text="PenRed for Blender", icon="PREFERENCES")
        layout.label(text="pyPenred and dependencies are bundled with this extension.")

addonClasses = (PYPENRED_Preferences,)

# Register all modules
def register():

    for c in addonClasses:
        bpy.utils.register_class(c)
    
    addon_properties.register()
    operators.register()
    ui.register()
    
# Unregister all modules
def unregister():
    
    ui.unregister()
    operators.unregister()
    addon_properties.unregister()

    for c in addonClasses:
        bpy.utils.unregister_class(c)    

if __name__ == "__main__":
    register()
