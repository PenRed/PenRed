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

def createHeaderCT(f, name, phiIni, phiEnd, nProj,
                   center, rad, tmin, dt, nSecond, toRound):
    f.write(f"# CT Source configuration for '{name}'\n")

    prefix = f"sources/generic/{name}/specific"

    f.write(f"{prefix}/type \"CT\"\n")
    f.write(f"{prefix}/phi-ini {phiIni:.3f}\n")
    f.write(f"{prefix}/phi-end {phiEnd:.3f}\n")
    f.write(f"{prefix}/nProjections {nProj}\n\n")

    f.write(f"{prefix}/CTx0 {round(center[0],toRound)}\n")
    f.write(f"{prefix}/CTy0 {round(center[1],toRound)}\n")
    f.write(f"{prefix}/CTz0 {round(center[2],toRound)}\n")
    f.write(f"{prefix}/rad {round(rad,toRound)}\n\n")
    
    f.write(f"{prefix}/tmin {tmin:.4}\n")
    f.write(f"{prefix}/dt {dt:.4}\n\n")

    if nSecond > 0:
        f.write(f"{prefix}/nSecondaries {nSecond}\n\n")

def createPSF(f, name, nhist, isCT, psfFile, psfMaxE,
              split, window, shift, orientation, toRound):
    f.write(f"# PSF Source configuration for '{name}'\n")

    prefix = f"sources/generic/{name}"
    
    # Write number of histories
    f.write(f"{prefix}/nhist {nhist:.5e}\n")
    
    if not isCT:
        prefix = f"{prefix}/specific"
        f.write(f"{prefix}/type \"PSF\"\n")
    else:
        prefix = f"{prefix}/specific/psf"
        
    f.write(f"{prefix}/filename \"{psfFile}\"\n")
    f.write(f"{prefix}/Emax {psfMaxE:.3e}\n")

    f.write(f"{prefix}/nsplit {split}\n")
    f.write(f"{prefix}/wght-window [{window[0]:.4e}, {window[1]:.4e}]\n\n")

    if shift:
        f.write(f"{prefix}/translation/dx {round(shift[0], toRound)}\n")
        f.write(f"{prefix}/translation/dy {round(shift[1], toRound)}\n")
        f.write(f"{prefix}/translation/dz {round(shift[2], toRound)}\n")

    if orientation:
        f.write(f"{prefix}/rotation/omega {round(orientation[0], toRound)}\n")
        f.write(f"{prefix}/rotation/theta {round(orientation[1], toRound)}\n")
        f.write(f"{prefix}/rotation/phi {round(orientation[2], toRound)}\n")
    

def createGeneric(f, name, nhist, particleType,
                  energyType, energy, spcFile,
                  aperture, direction,
                  spatialType, position, size,
                  sourceMat, toRound):
    f.write(f"\n# Generic Source configuration for '{name}'\n")

    prefix = f"sources/generic/{name}"
    
    # Write number of histories
    f.write(f"{prefix}/nhist {nhist:.5e}\n")

    # Write particle type
    if particleType == "PART_GAMMA":
        f.write(f"{prefix}/kpar \"gamma\"\n")
    elif particleType == "PART_ELECTRON":
        f.write(f"{prefix}/kpar \"electron\"\n")
    elif particleType == "PART_POSITRON":
        f.write(f"{prefix}/kpar \"positron\"\n")
    else:
        f.write(f"{prefix}/kpar \"unknown\"\n")

    # Write source material
    if sourceMat > 0:
        f.write(f"{prefix}/source-material {sourceMat}\n")        

    # Energy
    f.write(f"\n# - Energy parameters\n")
    if energyType == "ENERGY_MONO":
        f.write(f"{prefix}/energy/type \"MONOENERGETIC\"\n")
        f.write(f"{prefix}/energy/energy {energy:.3e}\n")
    else:
        f.write(f"{prefix}/energy/type \"FILE_SPECTRUM\"\n")
        f.write(f"{prefix}/energy/filename \"{spcFile}\"\n")

    # Direction
    f.write(f"\n# - Direction parameters\n")
    
    f.write(f"{prefix}/direction/type \"SOLID_ANGLE\"\n\n")
    f.write(f"{prefix}/direction/u {direction[0]:.3f}\n")
    f.write(f"{prefix}/direction/v {direction[1]:.3f}\n")
    f.write(f"{prefix}/direction/w {direction[2]:.3f}\n\n")

    f.write(f"{prefix}/direction/theta0 0.0\n")
    f.write(f"{prefix}/direction/theta1 { aperture/2.0:.3f}\n\n")
    
    f.write(f"{prefix}/direction/phi0 0.0\n")
    f.write(f"{prefix}/direction/dphi 360.0\n\n")

    # Spatial
    f.write(f"\n# - Spatial parameters\n")

    if spatialType == "SPATIAL_POINT":
        f.write(f"{prefix}/spatial/type \"POINT\"\n\n")
        
    elif spatialType == "SPATIAL_BOX":
        f.write(f"{prefix}/spatial/type \"BOX\"\n\n")
        
        f.write(f"{prefix}/spatial/size/dx {round(size[0],toRound)}\n")
        f.write(f"{prefix}/spatial/size/dy {round(size[1],toRound)}\n")
        f.write(f"{prefix}/spatial/size/dz {round(size[2],toRound)}\n\n")
        
    elif spatialType == "SPATIAL_CYL":
        f.write(f"{prefix}/spatial/type \"CYLINDER\"\n\n")

        f.write(f"{prefix}/spatial/size/rmin {round(size[0],toRound)}\n")
        f.write(f"{prefix}/spatial/size/rmax {round(size[1],toRound)}\n")
        f.write(f"{prefix}/spatial/size/dzmin {round(size[2],toRound)}\n\n")        
        f.write(f"{prefix}/spatial/size/dzmax {round(size[3],toRound)}\n\n")        
        
    f.write(f"{prefix}/spatial/position/x {round(position[0],toRound)}\n")
    f.write(f"{prefix}/spatial/position/y {round(position[1],toRound)}\n")
    f.write(f"{prefix}/spatial/position/z {round(position[2],toRound)}\n\n")

def createTime(f, name, timeType, halfLife, timeWindow):

    f.write(f"\n# - Time parameters\n")

    prefix = f"sources/generic/{name}"
    f.write(f"{prefix}/record-time true\n\n")

    if timeType == "TIME_DECAY":
        f.write(f"{prefix}/time/type \"DECAY\"\n\n")
        f.write(f"{prefix}/time/halfLife {halfLife:.5f}\n")
        f.write(f"{prefix}/time/time/time0 {timeWindow[0]:.3f}\n")
        f.write(f"{prefix}/time/time/time1 {timeWindow[1]:.3f}\n\n")
    
    
