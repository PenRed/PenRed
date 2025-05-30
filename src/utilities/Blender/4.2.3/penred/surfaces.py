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

# Construct quadric file functions
####################################
from math import cos
from math import sin

def initBody(f, nObj, name, material, module):
    ### Create the object
    f.write("0000000000000000000000000000000000000000000000000000000000000000\n")
    if module:
        f.write("MODULE  (%4i) %s\n" % (nObj,name))
    else:
        f.write("BODY    (%4i) %s\n" % (nObj,name))
    f.write("MATERIAL(%4i)\n" % (material))    

def addChilds(f, childs):
    for nObj,module,name in childs:
        if module:
            f.write("MODULE  (%4i) %s\n" % (nObj,name))
        else:
            f.write("BODY    (%4i) %s\n" % (nObj,name))
        
def endFile(f):
    f.write("0000000000000000000000000000000000000000000000000000000000000000\n")
    f.write("END      0000000000000000000000000000000000000000000000000000000\n")

### CUBE
#############

def createCubeSurfaces(f, x, y, z, dx, dy, dz, omega, theta, phi, nSurf, name, toRound):
    
    ### Plane -Z
    
    if omega == 0.0 and theta == 0.0 and phi == 0.0:
        xRot = x
        yRot = y
        zRot = z - dz/2.0
    else:
        xRot=-dz/2.0*cos(phi)*sin(theta)+x
        yRot=-dz/2.0*sin(theta)*sin(phi)+y
        zRot=-dz/2.0*cos(theta)+z
    
    f.write("0000000000000000000000000000000000000000000000000000000000000000\n")
    f.write("SURFACE (%4i) %s -Z plane\n" % (nSurf,name))
    f.write("INDICES=( 0, 0, 0, 1, 0)\n")
    f.write("X-SHIFT=(%+.15E,   0)\n" % (round(xRot,toRound)))
    f.write("Y-SHIFT=(%+.15E,   0)\n" % (round(yRot,toRound)))
    f.write("Z-SHIFT=(%+.15E,   0)\n" % (round(zRot,toRound)))
    f.write("  OMEGA=(%+.15E,   0) RAD\n" % (round(omega,toRound)))
    f.write("  THETA=(%+.15E,   0) RAD\n" % (round(theta,toRound)))
    f.write("    PHI=(%+.15E,   0) RAD\n" % (round(phi,toRound)))

    #Increment surface number
    nSurf += 1
    
    ### Plane +Z
    
    if omega == 0.0 and theta == 0.0 and phi == 0.0:
        xRot = x
        yRot = y
        zRot = z + dz/2.0
    else:
        xRot=+dz/2.0*cos(phi)*sin(theta)+x
        yRot=+dz/2.0*sin(theta)*sin(phi)+y
        zRot=+dz/2.0*cos(theta)+z
    
    f.write("0000000000000000000000000000000000000000000000000000000000000000\n")
    f.write("SURFACE (%4i) %s +Z plane\n" % (nSurf,name))
    f.write("INDICES=( 0, 0, 0, 1, 0)\n")
    f.write("X-SHIFT=(%+.15E,   0)\n" % (round(xRot,toRound)))
    f.write("Y-SHIFT=(%+.15E,   0)\n" % (round(yRot,toRound)))
    f.write("Z-SHIFT=(%+.15E,   0)\n" % (round(zRot,toRound)))
    f.write("  OMEGA=(%+.15E,   0) RAD\n" % (round(omega,toRound)))
    f.write("  THETA=(%+.15E,   0) RAD\n" % (round(theta,toRound)))
    f.write("    PHI=(%+.15E,   0) RAD\n" % (round(phi,toRound)))

    #Increment surface number
    nSurf += 1
    
    ### Plane -X
    
    if omega == 0.0 and theta == 0.0 and phi == 0.0:
        xRot = x - dx/2.0
        yRot = y
        zRot = z
    else:
        xRot=-dx/2.0*(cos(theta)*cos(phi)*cos(omega)-sin(phi)*sin(omega))+x
        yRot=-dx/2.0*(cos(theta)*cos(omega)*sin(phi)+cos(phi)*sin(omega))+y
        zRot=+dx/2.0*cos(omega)*sin(theta)+z
    
    f.write("0000000000000000000000000000000000000000000000000000000000000000\n")
    f.write("SURFACE (%4i) %s -X plane\n" % (nSurf,name))
    f.write("INDICES=( 0, 0, 0, 0, 0)\n")
    f.write("     AX=(+1.000000000000000E+00,   0)\n")
    f.write("1111111111111111111111111111111111111111111111111111111111111111\n")
    f.write("X-SHIFT=(%+.15E,   0)\n" % (round(xRot,toRound)))
    f.write("Y-SHIFT=(%+.15E,   0)\n" % (round(yRot,toRound)))
    f.write("Z-SHIFT=(%+.15E,   0)\n" % (round(zRot,toRound)))
    f.write("  OMEGA=(%+.15E,   0) RAD\n" % (round(omega,toRound)))
    f.write("  THETA=(%+.15E,   0) RAD\n" % (round(theta,toRound)))
    f.write("    PHI=(%+.15E,   0) RAD\n" % (round(phi,toRound)))
    
    #Increment surface number
    nSurf += 1
    
    ### Plane +X
    
    if omega == 0.0 and theta == 0.0 and phi == 0.0:
        xRot = x + dx/2.0
        yRot = y
        zRot = z
    else:
	    xRot=dx/2.0*(cos(theta)*cos(phi)*cos(omega)-sin(phi)*sin(omega))+x
	    yRot=dx/2.0*(cos(theta)*cos(omega)*sin(phi)+cos(phi)*sin(omega))+y
	    zRot=-dx/2.0*cos(omega)*sin(theta)+z
    
    f.write("0000000000000000000000000000000000000000000000000000000000000000\n")
    f.write("SURFACE (%4i) %s +X plane\n" % (nSurf,name))
    f.write("INDICES=( 0, 0, 0, 0, 0)\n")
    f.write("     AX=(+1.000000000000000E+00,   0)\n")
    f.write("1111111111111111111111111111111111111111111111111111111111111111\n")
    f.write("X-SHIFT=(%+.15E,   0)\n" % (round(xRot,toRound)))
    f.write("Y-SHIFT=(%+.15E,   0)\n" % (round(yRot,toRound)))
    f.write("Z-SHIFT=(%+.15E,   0)\n" % (round(zRot,toRound)))
    f.write("  OMEGA=(%+.15E,   0) RAD\n" % (round(omega,toRound)))
    f.write("  THETA=(%+.15E,   0) RAD\n" % (round(theta,toRound)))
    f.write("    PHI=(%+.15E,   0) RAD\n" % (round(phi,toRound)))
    
    #Increment surface number
    nSurf += 1

    ### Plane -Y
    
    if omega == 0.0 and theta == 0.0 and phi == 0.0:
        xRot = x
        yRot = y - dy/2.0
        zRot = z
    else:
	    xRot=-dy/2.0*(-cos(omega)*sin(phi)-cos(theta)*cos(phi)*sin(omega))+x
	    yRot=-dy/2.0*(cos(phi)*cos(omega)-cos(theta)*sin(phi)*sin(omega))+y
	    zRot=-dy/2.0*sin(theta)*sin(omega)+z
    
    f.write("0000000000000000000000000000000000000000000000000000000000000000\n")
    f.write("SURFACE (%4i) %s -Y plane\n" % (nSurf,name))
    f.write("INDICES=( 0, 0, 0, 0, 0)\n")
    f.write("     AY=(+1.000000000000000E+00,   0)\n")
    f.write("1111111111111111111111111111111111111111111111111111111111111111\n")
    f.write("X-SHIFT=(%+.15E,   0)\n" % (round(xRot,toRound)))
    f.write("Y-SHIFT=(%+.15E,   0)\n" % (round(yRot,toRound)))
    f.write("Z-SHIFT=(%+.15E,   0)\n" % (round(zRot,toRound)))
    f.write("  OMEGA=(%+.15E,   0) RAD\n" % (round(omega,toRound)))
    f.write("  THETA=(%+.15E,   0) RAD\n" % (round(theta,toRound)))
    f.write("    PHI=(%+.15E,   0) RAD\n" % (round(phi,toRound)))
    
    #Increment surface number
    nSurf += 1

    ### Plane +Y
    
    if omega == 0.0 and theta == 0.0 and phi == 0.0:
        xRot = x
        yRot = y + dy/2.0
        zRot = z
    else:
	    xRot=dy/2.0*(-cos(omega)*sin(phi)-cos(theta)*cos(phi)*sin(omega))+x
	    yRot=dy/2.0*(cos(phi)*cos(omega)-cos(theta)*sin(phi)*sin(omega))+y
	    zRot=dy/2.0*sin(theta)*sin(omega)+z
    
    f.write("0000000000000000000000000000000000000000000000000000000000000000\n")
    f.write("SURFACE (%4i) %s +Y plane\n" % (nSurf,name))
    f.write("INDICES=( 0, 0, 0, 0, 0)\n")
    f.write("     AY=(+1.000000000000000E+00,   0)\n")
    f.write("1111111111111111111111111111111111111111111111111111111111111111\n")
    f.write("X-SHIFT=(%+.15E,   0)\n" % (round(xRot,toRound)))
    f.write("Y-SHIFT=(%+.15E,   0)\n" % (round(yRot,toRound)))
    f.write("Z-SHIFT=(%+.15E,   0)\n" % (round(zRot,toRound)))
    f.write("  OMEGA=(%+.15E,   0) RAD\n" % (round(omega,toRound)))
    f.write("  THETA=(%+.15E,   0) RAD\n" % (round(theta,toRound)))
    f.write("    PHI=(%+.15E,   0) RAD\n" % (round(phi,toRound)))
    
    #Increment surface number
    nSurf += 1
    
    return nSurf

def setCubeSurfaces(f, initSurf, sign): 
    
    ###  Set the sufraces
    # -Z
    f.write("SURFACE (%4i), SIDE POINTER=(%+2i)\n" % (initSurf,sign))
    initSurf += 1
    # +Z
    f.write("SURFACE (%4i), SIDE POINTER=(%+2i)\n" % (initSurf,-sign))
    initSurf += 1

    # -X
    f.write("SURFACE (%4i), SIDE POINTER=(%+2i)\n" % (initSurf,sign))
    initSurf += 1
    # +X
    f.write("SURFACE (%4i), SIDE POINTER=(%+2i)\n" % (initSurf,-sign))
    initSurf += 1
    
    # -Y
    f.write("SURFACE (%4i), SIDE POINTER=(%+2i)\n" % (initSurf,sign))
    initSurf += 1
    # +Y
    f.write("SURFACE (%4i), SIDE POINTER=(%+2i)\n" % (initSurf,-sign))
    initSurf += 1 

    return initSurf
    
### TRAPEZOID
##############

def getLineEq(x1,y1,x2,y2):
    
    m = (y1-y2)/(x1-x2)
    b = y1 - m*x1

    return m, b

def getSlope(x1,y1,x2,y2):
    
    m = (y1-y2)/(x1-x2)
    return m

def createTrapezoidSurfaces(f, x, y, z,
                            dxBot, dyBot, dxTop, dyTop, dz,
                            omega, theta, phi, nSurf, name, toRound):

    # Calculate the dx and dy as the mid point between top and bot x
    # and y positions to define X and Y planes cutting the origin (0,0)
    dx = (dxBot + dxTop)/2.0
    dy = (dyBot + dyTop)/2.0
    
    ### Plane -Z
    
    if omega == 0.0 and theta == 0.0 and phi == 0.0:
        xRot = x
        yRot = y
        zRot = z - dz/2.0
    else:
        xRot=-dz/2.0*cos(phi)*sin(theta)+x
        yRot=-dz/2.0*sin(theta)*sin(phi)+y
        zRot=-dz/2.0*cos(theta)+z
    
    f.write("0000000000000000000000000000000000000000000000000000000000000000\n")
    f.write("SURFACE (%4i) %s -Z plane\n" % (nSurf,name))
    f.write("INDICES=( 0, 0, 0, 1, 0)\n")
    f.write("X-SHIFT=(%+.15E,   0)\n" % (round(xRot,toRound)))
    f.write("Y-SHIFT=(%+.15E,   0)\n" % (round(yRot,toRound)))
    f.write("Z-SHIFT=(%+.15E,   0)\n" % (round(zRot,toRound)))
    f.write("  OMEGA=(%+.15E,   0) RAD\n" % (round(omega,toRound)))
    f.write("  THETA=(%+.15E,   0) RAD\n" % (round(theta,toRound)))
    f.write("    PHI=(%+.15E,   0) RAD\n" % (round(phi,toRound)))

    #Increment surface number
    nSurf += 1
    
    ### Plane +Z
    
    if omega == 0.0 and theta == 0.0 and phi == 0.0:
        xRot = x
        yRot = y
        zRot = z + dz/2.0
    else:
        xRot=+dz/2.0*cos(phi)*sin(theta)+x
        yRot=+dz/2.0*sin(theta)*sin(phi)+y
        zRot=+dz/2.0*cos(theta)+z
    
    f.write("0000000000000000000000000000000000000000000000000000000000000000\n")
    f.write("SURFACE (%4i) %s +Z plane\n" % (nSurf,name))
    f.write("INDICES=( 0, 0, 0, 1, 0)\n")
    f.write("X-SHIFT=(%+.15E,   0)\n" % (round(xRot,toRound)))
    f.write("Y-SHIFT=(%+.15E,   0)\n" % (round(yRot,toRound)))
    f.write("Z-SHIFT=(%+.15E,   0)\n" % (round(zRot,toRound)))
    f.write("  OMEGA=(%+.15E,   0) RAD\n" % (round(omega,toRound)))
    f.write("  THETA=(%+.15E,   0) RAD\n" % (round(theta,toRound)))
    f.write("    PHI=(%+.15E,   0) RAD\n" % (round(phi,toRound)))

    #Increment surface number
    nSurf += 1
    
    ### Plane -X

    # Calculate plane parameters
    if dxBot != dxTop:
        z1 =  dz/2.0
        z2 = -dz/2.0
        # Defined the plane to cut the origin (0,0)
        x1 = (dx-dxTop)/2.0
        x2 = (dx-dxBot)/2.0
        m = getSlope(x1,z1,x2,z2)
        Ax = m
        Az = -1
        A0 = 0.0
    else:
        Ax = 1.0
        Az = 0.0
        A0 = 0.0
    
    if omega == 0.0 and theta == 0.0 and phi == 0.0:
        xRot = x - dx/2.0
        yRot = y
        zRot = z
    else:
        xRot=-dx/2.0*(cos(theta)*cos(phi)*cos(omega)-sin(phi)*sin(omega))+x
        yRot=-dx/2.0*(cos(theta)*cos(omega)*sin(phi)+cos(phi)*sin(omega))+y
        zRot=+dx/2.0*cos(omega)*sin(theta)+z
    
    f.write("0000000000000000000000000000000000000000000000000000000000000000\n")
    f.write("SURFACE (%4i) %s -X plane\n" % (nSurf,name))
    f.write("INDICES=( 0, 0, 0, 0, 0)\n")
    f.write("     AX=(%+.15E,   0)\n" % (round(Ax,toRound)))
    f.write("     AZ=(%+.15E,   0)\n" % (round(Az,toRound)))
    f.write("     A0=(%+.15E,   0)\n" % (round(A0,toRound)))
    f.write("1111111111111111111111111111111111111111111111111111111111111111\n")
    f.write("X-SHIFT=(%+.15E,   0)\n" % (round(xRot,toRound)))
    f.write("Y-SHIFT=(%+.15E,   0)\n" % (round(yRot,toRound)))
    f.write("Z-SHIFT=(%+.15E,   0)\n" % (round(zRot,toRound)))
    f.write("  OMEGA=(%+.15E,   0) RAD\n" % (round(omega,toRound)))
    f.write("  THETA=(%+.15E,   0) RAD\n" % (round(theta,toRound)))
    f.write("    PHI=(%+.15E,   0) RAD\n" % (round(phi,toRound)))
    
    #Increment surface number
    nSurf += 1
    
    ### Plane +X
    
    # Calculate plane parameters
    if dxBot != dxTop:
        z1 =  dz/2.0
        z2 = -dz/2.0
        # Defined the plane to cut the origin (0,0)
        x1 = (dxTop-dx)/2.0
        x2 = (dxBot-dx)/2.0
        m = getSlope(x1,z1,x2,z2)
        Ax = m
        Az = -1
        A0 = 0.0
    else:
        Ax = 1.0
        Az = 0.0
        A0 = 0.0
    
    if omega == 0.0 and theta == 0.0 and phi == 0.0:
        xRot = x + dx/2.0
        yRot = y
        zRot = z
    else:
	    xRot=dx/2.0*(cos(theta)*cos(phi)*cos(omega)-sin(phi)*sin(omega))+x
	    yRot=dx/2.0*(cos(theta)*cos(omega)*sin(phi)+cos(phi)*sin(omega))+y
	    zRot=-dx/2.0*cos(omega)*sin(theta)+z
    
    f.write("0000000000000000000000000000000000000000000000000000000000000000\n")
    f.write("SURFACE (%4i) %s +X plane\n" % (nSurf,name))
    f.write("INDICES=( 0, 0, 0, 0, 0)\n")
    f.write("     AX=(%+.15E,   0)\n" % (round(Ax,toRound)))
    f.write("     AZ=(%+.15E,   0)\n" % (round(Az,toRound)))
    f.write("     A0=(%+.15E,   0)\n" % (round(A0,toRound)))
    f.write("1111111111111111111111111111111111111111111111111111111111111111\n")
    f.write("X-SHIFT=(%+.15E,   0)\n" % (round(xRot,toRound)))
    f.write("Y-SHIFT=(%+.15E,   0)\n" % (round(yRot,toRound)))
    f.write("Z-SHIFT=(%+.15E,   0)\n" % (round(zRot,toRound)))
    f.write("  OMEGA=(%+.15E,   0) RAD\n" % (round(omega,toRound)))
    f.write("  THETA=(%+.15E,   0) RAD\n" % (round(theta,toRound)))
    f.write("    PHI=(%+.15E,   0) RAD\n" % (round(phi,toRound)))
    
    #Increment surface number
    nSurf += 1

    ### Plane -Y

    # Calculate plane parameters
    if dyBot != dyTop:
        z1 =  dz/2.0
        z2 = -dz/2.0
        # Defined the plane to cut the origin (0,0)
        y1 = (dy-dyTop)/2.0
        y2 = (dy-dyBot)/2.0
        m = getSlope(y1,z1,y2,z2)
        Ay = m
        Az = -1
        A0 = 0.0
    else:
        Ay = 1.0
        Az = 0.0
        A0 = 0.0    
    
    if omega == 0.0 and theta == 0.0 and phi == 0.0:
        xRot = x
        yRot = y - dy/2.0
        zRot = z
    else:
	    xRot=-dy/2.0*(-cos(omega)*sin(phi)-cos(theta)*cos(phi)*sin(omega))+x
	    yRot=-dy/2.0*(cos(phi)*cos(omega)-cos(theta)*sin(phi)*sin(omega))+y
	    zRot=-dy/2.0*sin(theta)*sin(omega)+z
    
    f.write("0000000000000000000000000000000000000000000000000000000000000000\n")
    f.write("SURFACE (%4i) %s -Y plane\n" % (nSurf,name))
    f.write("INDICES=( 0, 0, 0, 0, 0)\n")
    f.write("     AY=(%+.15E,   0)\n" % (round(Ay,toRound)))
    f.write("     AZ=(%+.15E,   0)\n" % (round(Az,toRound)))
    f.write("     A0=(%+.15E,   0)\n" % (round(A0,toRound)))
    f.write("1111111111111111111111111111111111111111111111111111111111111111\n")
    f.write("X-SHIFT=(%+.15E,   0)\n" % (round(xRot,toRound)))
    f.write("Y-SHIFT=(%+.15E,   0)\n" % (round(yRot,toRound)))
    f.write("Z-SHIFT=(%+.15E,   0)\n" % (round(zRot,toRound)))
    f.write("  OMEGA=(%+.15E,   0) RAD\n" % (round(omega,toRound)))
    f.write("  THETA=(%+.15E,   0) RAD\n" % (round(theta,toRound)))
    f.write("    PHI=(%+.15E,   0) RAD\n" % (round(phi,toRound)))
    
    #Increment surface number
    nSurf += 1

    ### Plane +Y

    # Calculate plane parameters
    if dyBot != dyTop:
        z1 =  dz/2.0
        z2 = -dz/2.0
        # Defined the plane to cut the origin (0,0)
        y1 = (dyTop-dy)/2.0
        y2 = (dyBot-dy)/2.0
        m = getSlope(y1,z1,y2,z2)
        Ay = m
        Az = -1
        A0 = 0.0
    else:
        Ay = 1.0
        Az = 0.0
        A0 = 0.0    
    
    if omega == 0.0 and theta == 0.0 and phi == 0.0:
        xRot = x
        yRot = y + dy/2.0
        zRot = z
    else:
	    xRot=dy/2.0*(-cos(omega)*sin(phi)-cos(theta)*cos(phi)*sin(omega))+x
	    yRot=dy/2.0*(cos(phi)*cos(omega)-cos(theta)*sin(phi)*sin(omega))+y
	    zRot=dy/2.0*sin(theta)*sin(omega)+z
    
    f.write("0000000000000000000000000000000000000000000000000000000000000000\n")
    f.write("SURFACE (%4i) %s +Y plane\n" % (nSurf,name))
    f.write("INDICES=( 0, 0, 0, 0, 0)\n")
    f.write("     AY=(%+.15E,   0)\n" % (round(Ay,toRound)))
    f.write("     AZ=(%+.15E,   0)\n" % (round(Az,toRound)))
    f.write("     A0=(%+.15E,   0)\n" % (round(A0,toRound)))
    f.write("1111111111111111111111111111111111111111111111111111111111111111\n")
    f.write("X-SHIFT=(%+.15E,   0)\n" % (round(xRot,toRound)))
    f.write("Y-SHIFT=(%+.15E,   0)\n" % (round(yRot,toRound)))
    f.write("Z-SHIFT=(%+.15E,   0)\n" % (round(zRot,toRound)))
    f.write("  OMEGA=(%+.15E,   0) RAD\n" % (round(omega,toRound)))
    f.write("  THETA=(%+.15E,   0) RAD\n" % (round(theta,toRound)))
    f.write("    PHI=(%+.15E,   0) RAD\n" % (round(phi,toRound)))
    
    #Increment surface number
    nSurf += 1
    
    return nSurf

def setTrapezoidSurfaces(f, dxBot, dyBot, dxTop, dyTop, initSurf, sign): 
    
    ###  Set the sufraces
    # -Z
    f.write("SURFACE (%4i), SIDE POINTER=(%+2i)\n" % (initSurf,sign))
    initSurf += 1
    # +Z
    f.write("SURFACE (%4i), SIDE POINTER=(%+2i)\n" % (initSurf,-sign))
    initSurf += 1

    if dxTop == dxBot:
        # -X
        f.write("SURFACE (%4i), SIDE POINTER=(%+2i)\n" % (initSurf,sign))
        initSurf += 1
        # +X
        f.write("SURFACE (%4i), SIDE POINTER=(%+2i)\n" % (initSurf,-sign))
        initSurf += 1
    elif dxTop > dxBot:
        # -X
        f.write("SURFACE (%4i), SIDE POINTER=(%+2i)\n" % (initSurf,-sign))
        initSurf += 1
        # +X
        f.write("SURFACE (%4i), SIDE POINTER=(%+2i)\n" % (initSurf,-sign))
        initSurf += 1
    else:
        # -X
        f.write("SURFACE (%4i), SIDE POINTER=(%+2i)\n" % (initSurf,sign))
        initSurf += 1
        # +X
        f.write("SURFACE (%4i), SIDE POINTER=(%+2i)\n" % (initSurf,sign))
        initSurf += 1

    if dyTop == dyBot:
        # -Y
        f.write("SURFACE (%4i), SIDE POINTER=(%+2i)\n" % (initSurf,sign))
        initSurf += 1
        # +Y
        f.write("SURFACE (%4i), SIDE POINTER=(%+2i)\n" % (initSurf,-sign))
        initSurf += 1
    elif dyTop > dyBot:
        # -Y
        f.write("SURFACE (%4i), SIDE POINTER=(%+2i)\n" % (initSurf,-sign))
        initSurf += 1
        # +Y
        f.write("SURFACE (%4i), SIDE POINTER=(%+2i)\n" % (initSurf,-sign))
        initSurf += 1
    else:
        # -Y
        f.write("SURFACE (%4i), SIDE POINTER=(%+2i)\n" % (initSurf,sign))
        initSurf += 1
        # +Y
        f.write("SURFACE (%4i), SIDE POINTER=(%+2i)\n" % (initSurf,sign))
        initSurf += 1
        
    return initSurf

### SPHERE
#############

def createSphereSurfaces(f, x, y, z, dx, dy, dz, nSurf, name, toRound): 

    f.write("0000000000000000000000000000000000000000000000000000000000000000\n")
    f.write("SURFACE (%4i) %s sphere surface\n" % (nSurf,name))
    f.write("INDICES=( 1, 1, 1, 0,-1)\n")
    f.write("X-SCALE=(%+.15E,   0)\n" % (round(dx/2.0,toRound)))
    f.write("Y-SCALE=(%+.15E,   0)\n" % (round(dy/2.0,toRound)))
    f.write("Z-SCALE=(%+.15E,   0)\n" % (round(dz/2.0,toRound)))	
    f.write("X-SHIFT=(%+.15E,   0)\n" % (round(x,toRound)))
    f.write("Y-SHIFT=(%+.15E,   0)\n" % (round(y,toRound)))
    f.write("Z-SHIFT=(%+.15E,   0)\n" % (round(z,toRound)))
    
    nSurf += 1
    
    return nSurf
    
def setSphereSurfaces(f, initSurf, sign): 
    
    ###  Set sphere surface
    f.write("SURFACE (%4i), SIDE POINTER=(%+2i)\n" % (initSurf,-sign))
    initSurf += 1

    return initSurf

### CYLINDER
#############

def printCylinderQuad(f,x,y,z,sx,sy,omega,theta,phi,nSurf,name,toRound):
    f.write("0000000000000000000000000000000000000000000000000000000000000000\n")
    f.write("SURFACE (%4i) %s cylinder surface\n" % (nSurf,name))
    f.write("INDICES=( 1, 1, 0, 0,-1)\n")
    f.write("X-SCALE=(%+.15E,   0)\n" % (round(sx,toRound)))
    f.write("Y-SCALE=(%+.15E,   0)\n" % (round(sy,toRound)))
    f.write("X-SHIFT=(%+.15E,   0)\n" % (round(x,toRound)))
    f.write("Y-SHIFT=(%+.15E,   0)\n" % (round(y,toRound)))
    f.write("Z-SHIFT=(%+.15E,   0)\n" % (round(z,toRound)))
    f.write("  OMEGA=(%+.15E,   0) RAD\n" % (round(omega,toRound)))
    f.write("  THETA=(%+.15E,   0) RAD\n" % (round(theta,toRound)))
    f.write("    PHI=(%+.15E,   0) RAD\n" % (round(phi,toRound)))    
    
    nSurf += 1
    return nSurf

def createCylinderSurfaces(f, x, y, z, dx, dy, dz, omega, theta, phi, nSurf, name, toRound):

    ### Cylinder surface
    
    nSurf = printCylinderQuad(f,x,y,z,dx/2.0,dy/2.0,omega,theta,phi,nSurf,name,toRound)
    
    ### -Z plane surface
    
    if omega == 0.0 and theta == 0.0 and phi == 0.0:
        xRot = x
        yRot = y
        zRot = z - dz/2.0
    else:
        xRot=-dz/2.0*cos(phi)*sin(theta)+x
        yRot=-dz/2.0*sin(theta)*sin(phi)+y
        zRot=-dz/2.0*cos(theta)+z    
    
    f.write("0000000000000000000000000000000000000000000000000000000000000000\n")
    f.write("SURFACE (%4i) %s -Z plane\n" % (nSurf,name))
    f.write("INDICES=( 0, 0, 0, 1, 0)\n")
    f.write("X-SHIFT=(%+.15E,   0)\n" % (round(xRot,toRound)))
    f.write("Y-SHIFT=(%+.15E,   0)\n" % (round(yRot,toRound)))
    f.write("Z-SHIFT=(%+.15E,   0)\n" % (round(zRot,toRound)))
    f.write("  OMEGA=(%+.15E,   0) RAD\n" % (round(omega,toRound)))
    f.write("  THETA=(%+.15E,   0) RAD\n" % (round(theta,toRound)))
    f.write("    PHI=(%+.15E,   0) RAD\n" % (round(phi,toRound)))
    
    nSurf += 1
    
    ### +Z plane surface    
    
    if omega == 0.0 and theta == 0.0 and phi == 0.0:
        xRot = x
        yRot = y
        zRot = z + dz/2.0
    else:
        xRot=+dz/2.0*cos(phi)*sin(theta)+x
        yRot=+dz/2.0*sin(theta)*sin(phi)+y
        zRot=+dz/2.0*cos(theta)+z
    
    f.write("0000000000000000000000000000000000000000000000000000000000000000\n")
    f.write("SURFACE (%4i) %s +Z plane\n" % (nSurf,name))
    f.write("INDICES=( 0, 0, 0, 1, 0)\n")
    f.write("X-SHIFT=(%+.15E,   0)\n" % (round(xRot,toRound)))
    f.write("Y-SHIFT=(%+.15E,   0)\n" % (round(yRot,toRound)))
    f.write("Z-SHIFT=(%+.15E,   0)\n" % (round(zRot,toRound)))
    f.write("  OMEGA=(%+.15E,   0) RAD\n" % (round(omega,toRound)))
    f.write("  THETA=(%+.15E,   0) RAD\n" % (round(theta,toRound)))
    f.write("    PHI=(%+.15E,   0) RAD\n" % (round(phi,toRound)))

    #Increment surface number
    nSurf += 1
    
    return nSurf


### CONE
#############

def printConeQuad(f,x,y,z,dz,r,h,s1,s2,omega,theta,phi,nSurf,name,toRound):
    f.write("0000000000000000000000000000000000000000000000000000000000000000\n")
    f.write("SURFACE (%4i) %s cone surface\n" % (nSurf,name))
    f.write("INDICES=( 1, 1,-1, 0, 0)\n")
    f.write("X-SCALE=(%+.15E,   0)\n" % (round(r*s1,toRound)))
    f.write("Y-SCALE=(%+.15E,   0)\n" % (round(r*s2,toRound)))
    f.write("Z-SCALE=(%+.15E,   0)\n" % (round(h+dz,toRound)))
    f.write("X-SHIFT=(%+.15E,   0)\n" % (round(x,toRound)))
    f.write("Y-SHIFT=(%+.15E,   0)\n" % (round(y,toRound)))
    f.write("Z-SHIFT=(%+.15E,   0)\n" % (round(z,toRound)))
    f.write("  OMEGA=(%+.15E,   0) RAD\n" % (round(omega,toRound)))
    f.write("  THETA=(%+.15E,   0) RAD\n" % (round(theta,toRound)))
    f.write("    PHI=(%+.15E,   0) RAD\n" % (round(phi,toRound)))
    
    nSurf += 1
    return nSurf

def createConeSurfaces(f, x, y, z, r1, r2, h, s1, s2, omega, theta, phi, nSurf, name, toRound):

    if r1 > r2: #Top radius larger than bottom
        # Calculate the distance between the cone
        # quadric center and the closest plane
        dz=r2*h/(r1-r2)

        ### Cone surface

        if omega == 0.0 and theta == 0.0 and phi == 0.0:
            xRot = x
            yRot = y
            zRot = z + h/2.0+dz
        else:
            xRot=(h/2.0+dz)*cos(phi)*sin(theta)+x
            yRot=(h/2.0+dz)*sin(theta)*sin(phi)+y
            zRot=(h/2.0+dz)*cos(theta)+z
        
        nSurf = printConeQuad(f,xRot,yRot,zRot,dz,r1,h,s1,s2,omega,theta,phi,nSurf,name,toRound)        
        
    else: #Bottom radius larger than top
        # Calculate the distance between the cone
        # quadric center and the closest plane        
        dz=r1*h/(r2-r1)
        
        ### Cone surface

        if omega == 0.0 and theta == 0.0 and phi == 0.0:
            xRot = x
            yRot = y
            zRot = z - h/2.0-dz
        else:
            
            xRot=(-h/2.0-dz)*cos(phi)*sin(theta)+x;
            yRot=(-h/2.0-dz)*sin(theta)*sin(phi)+y;
            zRot=(-h/2.0-dz)*cos(theta)+z;

        nSurf = printConeQuad(f,xRot,yRot,zRot,dz,r2,h,s1,s2,omega,theta,phi,nSurf,name,toRound)
            
    ### bottom Z plane

    if omega == 0.0 and theta == 0.0 and phi == 0.0:
        xRot = x
        yRot = y
        zRot = z - h/2.0
    else:
        xRot=-h/2.0*cos(phi)*sin(theta)+x;
        yRot=-h/2.0*sin(theta)*sin(phi)+y;
        zRot=-h/2.0*cos(theta)+z;
    
    f.write("0000000000000000000000000000000000000000000000000000000000000000\n")
    f.write("SURFACE (%4i) %s cone bottom Z plane\n" % (nSurf,name))
    f.write("INDICES=( 0, 0, 0, 1, 0)\n")
    f.write("X-SHIFT=(%+.15E,   0)\n" % (round(xRot,toRound)))
    f.write("Y-SHIFT=(%+.15E,   0)\n" % (round(yRot,toRound)))
    f.write("Z-SHIFT=(%+.15E,   0)\n" % (round(zRot,toRound)))
    f.write("  OMEGA=(%+.15E,   0) RAD\n" % (round(omega,toRound)))
    f.write("  THETA=(%+.15E,   0) RAD\n" % (round(theta,toRound)))
    f.write("    PHI=(%+.15E,   0) RAD\n" % (round(phi,toRound)))        

    #Increment surface number
    nSurf += 1

    ### Top Z plane

    if omega == 0.0 and theta == 0.0 and phi == 0.0:
        xRot = x
        yRot = y
        zRot = z + h/2.0
    else:
        xRot=h/2.0*cos(phi)*sin(theta)+x;
        yRot=h/2.0*sin(theta)*sin(phi)+y;
        zRot=h/2.0*cos(theta)+z;
    
    f.write("0000000000000000000000000000000000000000000000000000000000000000\n")
    f.write("SURFACE (%4i) %s cone top Z plane\n" % (nSurf,name))
    f.write("INDICES=( 0, 0, 0, 1, 0)\n")
    f.write("X-SHIFT=(%+.15E,   0)\n" % (round(xRot,toRound)))
    f.write("Y-SHIFT=(%+.15E,   0)\n" % (round(yRot,toRound)))
    f.write("Z-SHIFT=(%+.15E,   0)\n" % (round(zRot,toRound)))
    f.write("  OMEGA=(%+.15E,   0) RAD\n" % (round(omega,toRound)))
    f.write("  THETA=(%+.15E,   0) RAD\n" % (round(theta,toRound)))
    f.write("    PHI=(%+.15E,   0) RAD\n" % (round(phi,toRound)))        

    #Increment surface number
    nSurf += 1
        
    return nSurf

def setCylinderConeSurfaces(f, initSurf, sign): 
    
    ###  Set cylinder/cone surface
    f.write("SURFACE (%4i), SIDE POINTER=(%+2i)\n" % (initSurf,-sign))
    initSurf += 1

    ###  Set -Z plane surface
    f.write("SURFACE (%4i), SIDE POINTER=(%+2i)\n" % (initSurf,sign))
    initSurf += 1
    
    ###  Set +Z plane surface
    f.write("SURFACE (%4i), SIDE POINTER=(%+2i)\n" % (initSurf,-sign))
    initSurf += 1
    
    return initSurf

### TUBE
#############

def createTubeSurfaces(f, x, y, z, rOut, rInner, sx, sy, dz, omega, theta, phi, nSurf, name, toRound):

    ### Outer Cylinder surface
    nSurf = printCylinderQuad(f,x,y,z,sx*rOut,sy*rOut,omega,theta,phi,nSurf,name,toRound)

    ### Inner Cylinder surface
    nSurf = printCylinderQuad(f,x,y,z,sx*rInner,sy*rInner,omega,theta,phi,nSurf,name,toRound)
    
    ### -Z plane surface
    
    if omega == 0.0 and theta == 0.0 and phi == 0.0:
        xRot = x
        yRot = y
        zRot = z - dz/2.0
    else:
        xRot=-dz/2.0*cos(phi)*sin(theta)+x
        yRot=-dz/2.0*sin(theta)*sin(phi)+y
        zRot=-dz/2.0*cos(theta)+z    
    
    f.write("0000000000000000000000000000000000000000000000000000000000000000\n")
    f.write("SURFACE (%4i) %s -Z plane\n" % (nSurf,name))
    f.write("INDICES=( 0, 0, 0, 1, 0)\n")
    f.write("X-SHIFT=(%+.15E,   0)\n" % (round(xRot,toRound)))
    f.write("Y-SHIFT=(%+.15E,   0)\n" % (round(yRot,toRound)))
    f.write("Z-SHIFT=(%+.15E,   0)\n" % (round(zRot,toRound)))
    f.write("  OMEGA=(%+.15E,   0) RAD\n" % (round(omega,toRound)))
    f.write("  THETA=(%+.15E,   0) RAD\n" % (round(theta,toRound)))
    f.write("    PHI=(%+.15E,   0) RAD\n" % (round(phi,toRound)))
    
    nSurf += 1
    
    ### +Z plane surface    
    
    if omega == 0.0 and theta == 0.0 and phi == 0.0:
        xRot = x
        yRot = y
        zRot = z + dz/2.0
    else:
        xRot=+dz/2.0*cos(phi)*sin(theta)+x
        yRot=+dz/2.0*sin(theta)*sin(phi)+y
        zRot=+dz/2.0*cos(theta)+z
    
    f.write("0000000000000000000000000000000000000000000000000000000000000000\n")
    f.write("SURFACE (%4i) %s +Z plane\n" % (nSurf,name))
    f.write("INDICES=( 0, 0, 0, 1, 0)\n")
    f.write("X-SHIFT=(%+.15E,   0)\n" % (round(xRot,toRound)))
    f.write("Y-SHIFT=(%+.15E,   0)\n" % (round(yRot,toRound)))
    f.write("Z-SHIFT=(%+.15E,   0)\n" % (round(zRot,toRound)))
    f.write("  OMEGA=(%+.15E,   0) RAD\n" % (round(omega,toRound)))
    f.write("  THETA=(%+.15E,   0) RAD\n" % (round(theta,toRound)))
    f.write("    PHI=(%+.15E,   0) RAD\n" % (round(phi,toRound)))

    #Increment surface number
    nSurf += 1
    
    return nSurf

def setTubeSurfaces(f, initSurf, sign): 
    
    ###  Set outer cylinder surface
    f.write("SURFACE (%4i), SIDE POINTER=(%+2i)\n" % (initSurf,-sign))
    initSurf += 1

    ###  Set inner cylinder surface
    f.write("SURFACE (%4i), SIDE POINTER=(%+2i)\n" % (initSurf, sign))
    initSurf += 1
    
    ###  Set -Z plane surface
    f.write("SURFACE (%4i), SIDE POINTER=(%+2i)\n" % (initSurf,sign))
    initSurf += 1
    
    ###  Set +Z plane surface
    f.write("SURFACE (%4i), SIDE POINTER=(%+2i)\n" % (initSurf,-sign))
    initSurf += 1
    
    return initSurf

### PLANE
#############

def createPlaneSurfaces(f, x, y, z, omega, theta, phi, nSurf, name, toRound):

    f.write("0000000000000000000000000000000000000000000000000000000000000000\n")
    f.write("SURFACE (%4i) %s infinite plane\n" % (nSurf,name))
    f.write("INDICES=( 0, 0, 0, 1, 0)\n")
    f.write("X-SHIFT=(%+.15E,   0)\n" % (round(x,toRound)))
    f.write("Y-SHIFT=(%+.15E,   0)\n" % (round(y,toRound)))
    f.write("Z-SHIFT=(%+.15E,   0)\n" % (round(z,toRound)))
    f.write("  OMEGA=(%+.15E,   0) RAD\n" % (round(omega,toRound)))
    f.write("  THETA=(%+.15E,   0) RAD\n" % (round(theta,toRound)))
    f.write("    PHI=(%+.15E,   0) RAD\n" % (round(phi,toRound)))

    #Increment surface number
    nSurf += 1
    
    return nSurf

def setPlaneSurfaces(f, initSurf, sign): 
    
    ###  Set sphere surface
    f.write("SURFACE (%4i), SIDE POINTER=(%+2i)\n" % (initSurf,sign))
    initSurf += 1

    return initSurf

def setSurface(f, initSurf, sign): 
    
    ###  Set sphere surface
    f.write("SURFACE (%4i), SIDE POINTER=(%+2i)\n" % (initSurf,sign))
    initSurf += 1

    return initSurf
