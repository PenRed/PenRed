#
#
#    Copyright (C) 2022-2023 Universitat de València - UV
#    Copyright (C) 2022-2023 Universitat Politècnica de València - UPV
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
    "blender": (4, 2, 3),
    "location": "File > Export",
    "description": "Adds geometry construction capabilities and simulation configuration for PenRed simulations",
    "warning": "",
    "doc_url": "https://penred.github.io/PenRed/",
    "category": "Import-Export",
}

import bpy
from . import addon_properties, operators, ui

# Register all modules
def register():
    addon_properties.register()
    operators.register()
    ui.register()
    
# Unregister all modules
def unregister():
    ui.unregister()
    operators.unregister()
    addon_properties.unregister()

if __name__ == "__main__":
    register()
