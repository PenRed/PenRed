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

bl_info = {
    "name": "PenRed for blender",
    "author": "PenRed colaboration",
    "version": (2, 0),
    "blender": (4, 5, 0),
    "location": "File > Export",
    "description": "Adds geometry construction capabilities and simulation configuration for PenRed simulations",
    "warning": "",
    "doc_url": "https://penred.github.io/PenRed/Blender",
    "category": "Import-Export",
}

import bpy
import importlib
from . import addon_properties, operators, ui, dependency_manager

importlib.reload(addon_properties)
importlib.reload(dependency_manager)
importlib.reload(operators)
importlib.reload(ui)

# --- Operators for Install / Update / Uninstall pyPenred ---

class PYPENRED_OT_install_dependency(bpy.types.Operator):
    bl_idname = "pypenred.install_dependency"
    bl_label = "Install pyPenred"
    bl_description = "Downloads and installs the required pyPenred package"

    def execute(self, context):
        try:
            dependency_manager.install_package()
            self.report({"INFO"}, "pyPenred installed successfully!")
        except Exception as e:
            self.report({"ERROR"}, f"Installation failed: {e}")
            return {"CANCELLED"}
        return {"FINISHED"}


class PYPENRED_OT_update_dependency(bpy.types.Operator):
    bl_idname = "pypenred.update_dependency"
    bl_label = "Update pyPenred"
    bl_description = "Upgrades pyPenred to the latest version"

    def execute(self, context):
        try:
            dependency_manager.update_package()
            self.report({"INFO"}, "pyPenred updated successfully!")
        except Exception as e:
            self.report({"ERROR"}, f"Update failed: {e}")
            return {"CANCELLED"}
        return {"FINISHED"}


class PYPENRED_OT_uninstall_dependency(bpy.types.Operator):
    bl_idname = "pypenred.uninstall_dependency"
    bl_label = "Remove pyPenred"
    bl_description = "Uninstalls the pyPenred package"

    def execute(self, context):
        try:
            dependency_manager.uninstall_package()
            self.report({"INFO"}, "pyPenred removed.")
        except Exception as e:
            self.report({"ERROR"}, f"Removal failed: {e}")
            return {"CANCELLED"}
        return {"FINISHED"}

# --- Add-on Preferences UI ---

class PYPENRED_Preferences(bpy.types.AddonPreferences):
    bl_idname = __package__

    def draw(self, context):
        layout = self.layout
        installed = dependency_manager.is_installed()

        box = layout.box()
        box.label(text="Dependency Management", icon="PREFERENCES")

        status_row = box.row()
        if installed:
            status_row.label(text="Status: pyPenred is Installed", icon="CHECKMARK")
            row = box.row(align=True)
            row.operator(
                PYPENRED_OT_update_dependency.bl_idname, icon="FILE_REFRESH"
            )
            row.operator(
                PYPENRED_OT_uninstall_dependency.bl_idname, icon="TRASH"
            )
        else:
            status_row.label(
                text="Status: pyPenred is NOT installed", icon="ERROR"
            )
            box.operator(
                PYPENRED_OT_install_dependency.bl_idname, icon="IMPORT"
            )

class PYPENRED_OT_open_preferences(bpy.types.Operator):
    bl_idname = "pypenred.open_preferences"
    bl_label = "Open pyPenred Preferences"
    bl_description = "Opens preferences to manage pyPenred dependencies"

    def execute(self, context):
        bpy.ops.preferences.addon_show(module=__package__)
        return {'FINISHED'}

addonClasses = (PYPENRED_OT_install_dependency,
                PYPENRED_OT_update_dependency,
                PYPENRED_OT_uninstall_dependency,
                PYPENRED_Preferences,
                PYPENRED_OT_open_preferences)

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
