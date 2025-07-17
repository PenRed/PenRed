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

def createMaterial(f, name, index,
                   gammaCutoffType, gammaCutoffValue,
                   electronCutoffType, electronCutoffValue,
                   positronCutoffType, positronCutoffValue,
                   definitionType,
                   db, dbmat,
                   density, composition,
                   advanced,
                   C1, C2, WCC, WCR):
    f.write(f"# Material '{name}' configuration\n")    
    f.write(f"materials/{name}/number {index}\n")
    f.write(f"materials/{name}/filename \"{name}.mat\"\n")
    
    # Gamma
    if gammaCutoffType == "CUTOFF_EABS":
        f.write(f"materials/{name}/eabs/gamma {gammaCutoffValue:.3e}\n")
    else:
        f.write(f"materials/{name}/range/gamma {gammaCutoffValue:.3e}\n")

    # Electron
    if electronCutoffType == "CUTOFF_EABS":
        f.write(f"materials/{name}/eabs/electron {electronCutoffValue:.3e}\n")
    else:
        f.write(f"materials/{name}/range/electron {electronCutoffValue:.3e}\n")

    # Positron
    if positronCutoffType == "CUTOFF_EABS":
        f.write(f"materials/{name}/eabs/positron {positronCutoffValue:.3e}\n")
    else:
        f.write(f"materials/{name}/range/positron {positronCutoffValue:.3e}\n")

    if advanced:
        f.write(f"materials/{name}/C1 {C1:.2f}\n")
        f.write(f"materials/{name}/C2 {C2:.2f}\n")
        f.write(f"materials/{name}/WCC {WCC:.3e}\n")
        f.write(f"materials/{name}/WCR {WCR:.3e}\n")


    if definitionType == "COMPOSITION":
        # Define material by weight fraction
        f.write(f"materials/{name}/density {density:.4e}\n")
        for e in composition:
            f.write(f"materials/{name}/elements/{e[0]} {e[1]:.4e}\n")
    else:
        # Create a predefined DB material
        f.write(f"materials/{name}/DB \"{db}\"\n")
        f.write(f"materials/{name}/DB-material \"{dbmat}\"\n")        
        
    f.write("\n")
