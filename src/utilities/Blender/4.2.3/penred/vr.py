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

def createInteractionForcing(f, name, partType, interactionIndex, factorType, factor, weightWindow, bodies):
    f.write(f"# Interaction Forcing configuration for '{name}'\n")
    f.write(f"VR/IForcing/{name}/particle \"{partType}\"\n")
    f.write(f"VR/IForcing/{name}/interaction {interactionIndex}\n")

    if factorType == "AVERAGE":
        f.write(f"VR/IForcing/{name}/factor -{factor}\n")
    else:
        f.write(f"VR/IForcing/{name}/factor  {factor}\n")

    f.write(f"VR/IForcing/{name}/min-weight {weightWindow[0]:.4e}\n")
    f.write(f"VR/IForcing/{name}/max-weight {weightWindow[1]:.4e}\n")

    for body in bodies:
        f.write(f"VR/IForcing/{name}/bodies/{body} true\n")
    
    f.write("\n")

def createBremssSplitting(f, name, factor, bodies):
    f.write(f"# Bremsstrahlung splitting configuration for '{name}'\n")
    f.write(f"VR/bremss/{name}/splitting {factor}\n")
    for body in bodies:
        f.write(f"VR/bremss/{name}/bodies/{body} true\n")    
    f.write("\n")

def createXRaySplitting(f, name, bodies):
    f.write(f"# X-Ray splitting configuration for '{name}'\n")
    f.write(f"VR/photon/{name}/type \"XRAY_SPLITTING\"\n")
    for body in bodies:
        f.write(f"VR/photon/{name}/bodies/{body[0]}/splitting {body[1]}\n")    
    f.write("\n")
