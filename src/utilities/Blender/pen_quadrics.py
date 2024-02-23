#
#
#    Copyright (C) 2022-2023 Universitat de València - UV
#    Copyright (C) 2022-2023 Universitat Politècnica de València - UPV
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
    "name": "PenRed Geometry builder",
    "author": "Vicent Giménez Alventosa",
    "version": (1, 0),
    "blender": (2, 93, 0),
    "location": "View3D > Add > Quadric",
    "description": "Adds quadric objects to construct PenRed geometries and export meshes to PenRed format",
    "warning": "",
    "doc_url": "",
    "category": "Add Mesh",
}

##############################################################################
##############################################################################
####       EXPORT PRINT FUNCTIONS
##############################################################################
##############################################################################

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

def createConeSurfaces(f, x, y, z, r1, r2, h, s1, s2, omega, theta, phi, nSurf, name, skinSize, toRound):

    shell = False
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
        
        if skinSize > 0.0 and skinSize < r2:
            #Print the inner cone to create the shell 
            nSurf = printConeQuad(f,xRot,yRot,zRot,dz,r1-skinSize,h,s1,s2,omega,theta,phi,nSurf,name,toRound)
            shell = True
        
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

        if skinSize > 0.0 and skinSize < r1:
            #Print the inner cone to create the shell 
            nSurf = printConeQuad(f,xRot,yRot,zRot,dz,r2-skinSize,h,s1,s2,omega,theta,phi,nSurf,name,toRound)
            shell = True    
            
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
        
    return nSurf, shell

def setCylinderConeSurfaces(f, initSurf, sign, shell = False): 
    
    ###  Set cylinder/cone surface
    f.write("SURFACE (%4i), SIDE POINTER=(%+2i)\n" % (initSurf,-sign))
    initSurf += 1

    if shell:
        f.write("SURFACE (%4i), SIDE POINTER=(%+2i)\n" % (initSurf,sign))
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

##############################################################################
##############################################################################
##############################################################################

# Internal Blender menus, objects, etc
########################################

import bpy
import bmesh
from bpy_extras.io_utils import ExportHelper
from bpy.types import Operator
from bpy.props import FloatVectorProperty
from bpy_extras.object_utils import AddObjectHelper, object_data_add
from mathutils import Vector
from mathutils import Color
from math import cos, acos, sin, asin, tan, atan2, sqrt, pi

##############################################################################
##############################################################################
####       EXTRA PROPERTIES
##############################################################################
##############################################################################

materialRows = 20
materialColumns = 20
maxMaterial = materialRows*materialColumns

# Properties group
###################
class PenredProperties(bpy.types.PropertyGroup):    
    quadricType : bpy.props.StringProperty(name = "Quadric Type", default = "unknown")
    material : bpy.props.IntProperty(name = "Material", min = 0, max = maxMaterial, default = 1)
    module : bpy.props.BoolProperty(name = "Module", default = False)
    r1 : bpy.props.FloatProperty(name = "r1", default = 0.0)
    r2 : bpy.props.FloatProperty(name = "r2", default = 0.0)
    r3 : bpy.props.FloatProperty(name = "r3", default = 0.0)
    r4 : bpy.props.FloatProperty(name = "r4", default = 0.0)


##############################################################################
##############################################################################
####       EXPORT BLENDER MENUS AND CLASSES
##############################################################################
##############################################################################

# Export class
################
class export_penred(Operator, ExportHelper):
    """Export operator"""
    bl_idname = "export_penred.data"
    bl_label = "Export to PenRed"
    
    filename_ext = ".geo"
    
    filter_glob: bpy.props.StringProperty(
        default = "*.quad;*.mesh;*.geo",
        options={'HIDDEN'},
        maxlen=255,
        )
        
    toRound: bpy.props.IntProperty(
        name="Rounding decimal",
        description="Decimal limit to round the export values (position, dimensions, etc)",
        default=6,
        )

    
    exportType: bpy.props.EnumProperty(
        items=[("QUADRICS","Quadrics","Quadric based geometry",'',0), ("MESH","Mesh","Mesh based geometry",'',1)],
        name="Geometry type",
        description="PenRed geometry type to export",
        default=0,
        )

    onlyActive: bpy.props.BoolProperty(
        name="Only active",
        description="Exports only the active object",
        default=False,
        )

    avoidHide: bpy.props.BoolProperty(
        name="Omit hide",
        description="Disable the export of hide objects and their children",
        default=False,
        )
    
    def getObjInfo(self,context,obj):
        #Set object origin to origin geometry
        context.view_layer.objects.active = obj
        bpy.ops.object.origin_set( type = 'ORIGIN_GEOMETRY' )

        #Get name
        name = obj.name

        #Get position
        trans = obj.matrix_world.to_translation()
        x = trans[0]
        y = trans[1]
        z = trans[2]
        
        #Get scale
        sx,sy,sz = obj.matrix_world.to_scale()
                    
        #Get sizes
        dx = obj.dimensions.x
        dy = obj.dimensions.y
        dz = obj.dimensions.z
                    
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
        return x,y,z,dx,dy,dz,sx,sy,sz,omega,theta,phi,name
    
    def getChildrens(self,context,parent):
        children = []
        for obj in context.scene.objects:
            if obj.parent == parent:
                children.append(obj)
        return children
    
    def getMeshParent(self,context,obj):
        parent = obj.parent
        while True:
            if parent:
                if parent.type == 'MESH':
                    return parent
                else:
                    parent = parent.parent
            else:
                return None
    
    def createTriangleMesh(self,f,context,obj,toRound,forceWorld,avoidHide):

        if obj.type != 'MESH':
            return

        #Check if hide objects must be avoided
        if avoidHide:
            if obj.hide_get():
                #Skip it
                return

        #Get name
        name = obj.name
        name.replace(" ","_")
        if len(name) > 100:
            name = name[:100]
        
        #Get parent name
        if forceWorld:
            parentName = "void"
        else:
            parent = self.getMeshParent(context,obj)
            if parent:
                parentName = parent.name
                parentName.replace(" ","_")
                if len(parentName) > 100:
                    parentName = parentName[:100]
            else:
                parentName = "void"
        
        #Get number of vertex groups
        vgNames = obj.vertex_groups.keys()
        nVG = len(vgNames)

        #Get mesh data
        mesh = obj.data
        #Create triangles
        mesh.calc_loop_triangles()
        
        f.write("# Object: %s\n" % name)
        f.write("#MAT      #NFACES     #NVERTEX     #NAME        #PARENT NAME    #N VERTEX GROUPS\n")
        f.write(" %03d      %07d     %08d     %s        %s   %04d\n" % (obj.penred_settings.material, len(mesh.loop_triangles),len(mesh.vertices), name, parentName, nVG))

        #Create lists with vertex belonging to each vertex group
        vgIndexLists = []
        for i in range(nVG):
            vgIndexLists.append([])

        for v in mesh.vertices:
            for g in v.groups:
                vgIndexLists[g.group].append(v.index)

        #Print vertex groups
        f.write("# VERTEX GROUPS\n")
        for i in range(nVG):
            f.write("#NAME  #NVERTEX\n")
            f.write(" %s   %04d\n" % (vgNames[i], len(vgIndexLists[i])))
            for index in vgIndexLists[i]:
                f.write(" %04d\n" % (index))

        #Print vertex
        f.write("# VERTEX LIST\n")
        f.write("# Index  (X Y Z)\n")
        
        for vertex in mesh.vertices:
            vertexWorld = obj.matrix_world @ vertex.co
            f.write("%04d %+.*E %+.*E %+.*E\n" % (vertex.index, toRound+3, round(vertexWorld[0],toRound), toRound+3, round(vertexWorld[1],toRound), toRound+3, round(vertexWorld[2],toRound)))
        
        f.write("# FACES(triangles)\n")
        for tri in mesh.loop_triangles:
            f.write(" %03d %03d %03d\n" % (tri.vertices[0], tri.vertices[1], tri.vertices[2]))
        f.write("#\n#\n")
        
    
    def createObject(self,f,context,obj,nSurf,nObj,toRound,createChilds,avoidHide):
        
        #Check the quadric type
        if obj.penred_settings.quadricType == "unknown" and obj.type != "EMPTY": 
            return [], nSurf, nObj

        #Check if hide objects must be avoided
        if avoidHide:
            if obj.hide_get():
                #Skip it
                return [], nSurf, nObj

        #Get childrens
        childrens = []
        if createChilds:
            childrens = self.getChildrens(context,obj)
        
        #First, construct children objects
        tree = [] # Children tree information
        if len(childrens) > 0:
            for child in childrens:
                childTree, nSurf, nObj = self.createObject(f,context,child,nSurf,nObj,toRound,True,avoidHide)
                if len(childTree) > 0:
                    tree.extend(childTree)
                
        #If the object is an empty, no body must be created here
        if obj.type == "EMPTY":
            return tree,nSurf,nObj
        
        #Get object type
        quadType = obj.penred_settings.quadricType
                    
        #Get object properties
        x,y,z,dx,dy,dz,sx,sy,sz,omega,theta,phi,name = self.getObjInfo(context,obj)
                         
        #Save init surface
        initSurf = nSurf
        
        ### Create object surfaces
        if quadType == "CUBE":
            nSurf = createCubeSurfaces(f,x,y,z,dx,dy,dz,omega,theta,phi,nSurf,name,toRound)
        elif quadType == "SPHERE":
            nSurf = createSphereSurfaces(f,x,y,z,dx,dy,dz,nSurf,name,toRound)
        elif quadType == "CYLINDER":
            nSurf = createCylinderSurfaces(f,x,y,z, dx,dy,dz,omega,theta,phi,nSurf,name,toRound)
        elif quadType == "CONE":
            r1 = obj.penred_settings.r3
            r2 = obj.penred_settings.r4
            if r1 == r2:
                nSurf = createCylinderSurfaces(f,x,y,z, dx,dy,dz,omega,theta,phi,nSurf,name,toRound)
            else:
                nSurf, shell = createConeSurfaces(f,x,y,z,r1,r2,dz,sx,sy,omega,theta,phi,nSurf,name,-1.0,toRound)
        elif quadType == "PLANE":
            nSurf = createPlaneSurfaces(f,x,y,z,omega,theta,phi,nSurf,name,toRound)
        elif quadType == "TUBE":
            innerCreated = False
            #Construct outer cylinder
            nSurf = createCylinderSurfaces(f,x,y,z, dx,dy,dz,omega,theta,phi,nSurf,name,toRound)            
            #Check if the object has the boolean modifier
            if "Boolean" in obj.modifiers:
                #Get limiting cylinder
                innerCyl = obj.modifiers['Boolean'].object
                x,y,z,dx,dy,dz,sx,sy,sz,omega,theta,phi,name2 = self.getObjInfo(context,innerCyl)
                #Construct inner cylinder
                nSurf = printCylinderQuad(f,x,y,z,dx/2.0,dy/2.0,omega,theta,phi,nSurf,name2,toRound)
                #Flag inner quad
                innerCreated = True
        elif quadType == "SEMI_SPHERE":
            #Construct sphere. Take into account that the Z dimension
            #must not be multiplied by two to comepnsate the cut, because
            #we added an artificial vertex
            nSurf = createSphereSurfaces(f,x,y,z,dx,dy,dz,nSurf,name,toRound)
            #Create limiting plane
            nSurf = createPlaneSurfaces(f,x,y,z,omega,theta,phi,nSurf,name,toRound)            
        elif quadType == "SEMI_SPHERE_SHELL":
            #Construct sphere. Take into account that the Z dimension
            #must be multiplied by two to comepnsate the cut, because
            #the modifier remove the effect of artificial vertex
            innerCreated = False
            nSurf = createSphereSurfaces(f,x,y,z,dx,dy,2.0*dz,nSurf,name,toRound)
            #Get solidify modifier information
            if "Solidify" in obj.modifiers:
                #thickness
                thickness = obj.modifiers['Solidify'].thickness
                t2 = 2.0*thickness
                if t2 < dx and t2 < dy and t2 < dz:
                    nameInt = "Internal " + name
                    nSurf = createSphereSurfaces(f,x,y,z,dx-t2,dy-t2,2.0*dz-t2,nSurf,nameInt,toRound)
                    innerCreated = True
            #Create limiting planen
            nSurf = createPlaneSurfaces(f,x,y,z,omega,theta,phi,nSurf,name,toRound)            
        elif quadType == "CONE_SHELL":
            innerCreated = False
            r1 = obj.penred_settings.r3
            r2 = obj.penred_settings.r4
            #Get solidify modifier information
            if "Solidify" in obj.modifiers:
                #thickness
                thickness = obj.modifiers['Solidify'].thickness
            else:
                thickness = -1.0
            
            if r1 == r2:
                nSurf = createCylinderSurfaces(f,x,y,z, dx,dy,dz,omega,theta,phi,nSurf,name,toRound)               
                shell = False
            else:
                nSurf, shell = createConeSurfaces(f,x,y,z,r1,r2,dz,sx,sy,omega,theta,phi,nSurf,name,thickness,toRound)
            
        ### Init body
        initBody(f,nObj,name,obj.penred_settings.material,obj.penred_settings.module)
        
        ### Set object surfaces
        if quadType == "CUBE":
            initSurf = setCubeSurfaces(f,initSurf, 1)
        elif quadType == "SPHERE":
            initSurf = setSphereSurfaces(f,initSurf, 1)
        elif quadType == "CYLINDER" or quadType == "CONE":
            initSurf = setCylinderConeSurfaces(f,initSurf, 1)
        elif quadType == "PLANE":
            initSurf = setPlaneSurfaces(f,initSurf, 1)
        elif quadType == "TUBE":
            #Set cylinder surface
            initSurf = setCylinderConeSurfaces(f,initSurf, 1)
            #Set inner surface
            if innerCreated:
                initSurf = setSurface(f,initSurf, 1)            
        elif quadType == "SEMI_SPHERE":
            initSurf = setSphereSurfaces(f,initSurf, 1)
            initSurf = setPlaneSurfaces(f,initSurf, 1)
        elif quadType == "SEMI_SPHERE_SHELL":
            #External sphere
            initSurf = setSphereSurfaces(f,initSurf, 1)
            if innerCreated:
                #Internal sphere
                initSurf = setSurface(f,initSurf, 1)            
            #Plane
            initSurf = setPlaneSurfaces(f,initSurf, 1)
        elif quadType == "CONE_SHELL":
            initSurf = setCylinderConeSurfaces(f,initSurf, 1, shell)
            
        #Set childrens
        if len(tree) > 0:
            addChilds(f,tree)
        
        #Add current body or module to tree
        if obj.penred_settings.module:
            #Module childs can't be repeated by other bodies/modules
            tree = [(nObj,True,name)] 
        else:
            tree.append((nObj,False,name))
        nObj += 1
            
        return tree,nSurf,nObj
                 
    def countNonHideMeshes(self,context,obj):

        #Number of non hide meshes
        nMeshes = 0

        #Check this mesh status
        if not obj.hide_get():
            #Count this object if it is a mesh
            if obj.type == 'MESH':
                nMeshes = 1

            #Check childrens
            childrens = self.getChildrens(context,obj)
            for child in childrens:
                #Add contribution of all children
                nMeshes += self.countNonHideMeshes(context,child)

        #Return the number of non hide meshes
        return nMeshes

    def execute(self, context):
        
        #Open output file
        f = open(self.filepath,'w',encoding='utf-8')
        f.write("# Geometry file created with PenRed blender plugin v.0.2\n")
        
        if self.exportType == 'QUADRICS':
            #Create an array for object names
            objNames = []
            
            nSurf = 1 #Number of surface to be created
            nObj = 1 #Number of object to be created

            if self.onlyActive:
                self.createObject(f,context,bpy.context.active_object,nSurf,nObj,self.toRound,False,False)
            else:
                #Find objects with no parents
                for obj in context.scene.objects:
                    if not obj.parent:
                        #This object has no parent, create it
                        nSurf,nObj = self.createObject(f,context,obj,nSurf,nObj,self.toRound,True,self.avoidHide)[1:]

            endFile(f)
        elif self.exportType == 'MESH':
            
            if self.onlyActive:
                #Print number of objects
                f.write("# Number of objects:\n 1\n")

                self.createTriangleMesh(f,context,bpy.context.active_object,self.toRound,True,False)
            else:

                #Count number of meshes
                nMeshes = 0

                if self.avoidHide:
                    for obj in context.scene.objects:
                        if not obj.parent: #Begin from objects with no parent
                            nMeshes += self.countNonHideMeshes(context,obj)
                else:
                    #Count all mesh objects
                    for obj in context.scene.objects:
                        if obj.type == 'MESH':
                            nMeshes = nMeshes + 1

                #Print number of objects
                f.write("# Number of objects:\n %d\n" % (nMeshes))

                for obj in context.scene.objects:
                    self.createTriangleMesh(f,context,obj,self.toRound,False,self.avoidHide)
        else:
            f.write("# Unknown export format. Please, report this issue\n")
                
        f.close()
        
        return {'FINISHED'}
    
def menu_func_penred_export(self, context):
    self.layout.operator(export_penred.bl_idname, text="PenRed")

##############################################################################
##############################################################################
####       TRANSFORM TOOLS
##############################################################################
##############################################################################

###################################
### Quadric Transform options  ####
###################################            
    
def getLocalZdir(obj):
    
        nz = Vector((0.0,0.0,1.0))
        inv = obj.matrix_world.copy()
        inv.invert()
        
        dirZ = nz @ inv
        dirZ.normalize()
        return dirZ
    
    
### Put down tool
class QUADRIC_OT_transform_putdown(Operator, AddObjectHelper):
    """Set Quadric material view"""
    bl_idname = "transform.quadric_transform_putdown"
    bl_label = "Tool to clip objects down of the active one"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        selected = context.selected_objects
        active = context.view_layer.objects.active
        
        if len(selected) == 0:
            return {'FINISHED'}
        
        #Get active object position
        transAct = active.matrix_world.to_translation()

        #Get Z size
        dza = active.dimensions.z
        
        #Get local Z direction
        dirZ = getLocalZdir(active)

        for obj in selected:
            if obj != active:
                
                #Get Z size
                dz = obj.dimensions.z
                
                #Get position
                obj.matrix_world.translation = transAct - dirZ*(dza/2.0 + dz/2.0)
                
        return {'FINISHED'}

### Put top tool
class QUADRIC_OT_transform_puttop(Operator, AddObjectHelper):
    """Set Quadric material view"""
    bl_idname = "transform.quadric_transform_puttop"
    bl_label = "Tool to clip objects top of the active one"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        selected = context.selected_objects
        active = context.view_layer.objects.active
        
        if len(selected) == 0:
            return {'FINISHED'}
        
        #Get active object position
        transAct = active.matrix_world.to_translation()

        #Get Z size
        dza = active.dimensions.z
        
        #Get local Z direction
        dirZ = getLocalZdir(active)

        for obj in selected:
            if obj != active:
                
                #Get Z size
                dz = obj.dimensions.z
                
                #Get position
                obj.matrix_world.translation = transAct + dirZ*(dza/2.0 + dz/2.0)
                
        return {'FINISHED'}

### Put inside top tool
class QUADRIC_OT_transform_putinsidetop(Operator, AddObjectHelper):
    """Set Quadric material view"""
    bl_idname = "transform.quadric_transform_putinsidetop"
    bl_label = "Tool to clip objects inside top of the active one"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        selected = context.selected_objects
        active = context.view_layer.objects.active
        
        if len(selected) == 0:
            return {'FINISHED'}
        
        #Get active object position
        transAct = active.matrix_world.to_translation()

        #Get Z size
        dza = active.dimensions.z
        
        #Get local Z direction
        dirZ = getLocalZdir(active)

        for obj in selected:
            if obj != active:
                
                #Get Z size
                dz = obj.dimensions.z
                
                #Get position
                obj.matrix_world.translation = transAct + dirZ*(dza/2.0 - dz/2.0)
                
        return {'FINISHED'}

### Put inside down tool
class QUADRIC_OT_transform_putinsidedown(Operator, AddObjectHelper):
    """Set Quadric material view"""
    bl_idname = "transform.quadric_transform_putinsidedown"
    bl_label = "Tool to clip objects inside down of the active one"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        selected = context.selected_objects
        active = context.view_layer.objects.active
        
        if len(selected) == 0:
            return {'FINISHED'}
        
        #Get active object position
        transAct = active.matrix_world.to_translation()

        #Get Z size
        dza = active.dimensions.z
        
        #Get local Z direction
        dirZ = getLocalZdir(active)

        for obj in selected:
            if obj != active:
                
                #Get Z size
                dz = obj.dimensions.z
                
                #Get position
                obj.matrix_world.translation = transAct - dirZ*(dza/2.0 - dz/2.0)
                
        return {'FINISHED'}

### Put inside centered tool
class QUADRIC_OT_transform_putinsidecentered(Operator, AddObjectHelper):
    """Set Quadric material view"""
    bl_idname = "transform.quadric_transform_putinsidecentered"
    bl_label = "Tool to set objects inside of the active one centered"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        selected = context.selected_objects
        active = context.view_layer.objects.active
        
        if len(selected) == 0:
            return {'FINISHED'}
        
        #Get active object position
        transAct = active.matrix_world.to_translation()

        #Get local Z direction
        dirZ = getLocalZdir(active)

        for obj in selected:
            if obj != active:
                
                #Get position
                obj.matrix_world.translation = transAct
                
        return {'FINISHED'}


quadricTransformClasses = (
    QUADRIC_OT_transform_putdown,
    QUADRIC_OT_transform_puttop,
    QUADRIC_OT_transform_putinsidetop,
    QUADRIC_OT_transform_putinsidedown,
    QUADRIC_OT_transform_putinsidecentered
)

def menu_func_quadric_transform(self, context):
    self.layout.operator(QUADRIC_OT_transform_putdown.bl_idname, text="Put Down")
    self.layout.operator(QUADRIC_OT_transform_puttop.bl_idname, text="Put Top")
    self.layout.operator(QUADRIC_OT_transform_putinsidetop.bl_idname, text="Put Inside Top")
    self.layout.operator(QUADRIC_OT_transform_putinsidedown.bl_idname, text="Put Inside Down")
    self.layout.operator(QUADRIC_OT_transform_putinsidecentered.bl_idname, text="Put Inside Centered")

##############################################################################
##############################################################################
####       VIEW OPTIONS
##############################################################################
##############################################################################

# View Submenu class
#####################
class OBJECT_MT_view_quadric_submenu(bpy.types.Menu):
    bl_idname = "OBJECT_MT_view_quadric_submenu"
    bl_label = "Quadric view"
    bl_options = {'REGISTER', 'UNDO'}
    
    def draw(self, context):
        layout = self.layout

        layout.operator("view.quadric_material",text="Material")
        layout.operator("view.quadric_module",text="Module")

def menu_quadric_view_func(self, context):
    self.layout.menu("OBJECT_MT_view_quadric_submenu", text="Quadric")        

### Material view
class QUADRIC_OT_view_material(Operator, AddObjectHelper):
    """Set Quadric material view"""
    bl_idname = "view.quadric_material"
    bl_label = "Set Quadric Material View"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        for obj in context.scene.objects:
            if obj.penred_settings.quadricType != "unknown":
                #Is a quad type, assign a colour according to material
                material = obj.penred_settings.material
                #Get alpha channel
                alpha = obj.color[3]
                
                H = float(material % materialColumns)/float(materialColumns)
                V = 1.0-float(material // materialRows)/float(materialRows)
                if V < 0.1:
                    V = 0.1
                
                c = Color()
                c.hsv = H, 1.0, V
                obj.color = c.r, c.g, c.b, alpha   
                
        return {'FINISHED'}

### Module view
class QUADRIC_OT_view_module(Operator, AddObjectHelper):
    """Set Quadric module view"""
    bl_idname = "view.quadric_module"
    bl_label = "Set Quadric Inner View"
    bl_options = {'REGISTER', 'UNDO'}

    def execute(self, context):
        for obj in context.scene.objects:
            if obj.penred_settings.quadricType != "unknown":
                #Is a quad type, assign a colour according to module
                module = obj.penred_settings.module
                #Get alpha channel
                alpha = obj.color[3]
                
                if module:                
                    obj.color = 1.0, 1.0, 1.0, alpha   
                else:
                    obj.color = 0.0, 0.0, 0.0, alpha   
                
        return {'FINISHED'}

viewObjectClasses = (
    OBJECT_MT_view_quadric_submenu,
    QUADRIC_OT_view_material,
    QUADRIC_OT_view_module,
)

##############################################################################
##############################################################################
####       ADD OBJECTS CLASSES
##############################################################################
##############################################################################

# Add Submenu class
####################
class OBJECT_MT_quadric_submenu(bpy.types.Menu):
    bl_idname = "OBJECT_MT_quadric_submenu"
    bl_label = "Quadric"
    bl_options = {'REGISTER', 'UNDO'}
    
    def draw(self, context):
        layout = self.layout
        
        layout.operator("mesh.add_cube_quadric",text="Cube",icon='CUBE')
        layout.operator("mesh.add_sphere_quadric",text="Sphere",icon='SPHERE')
        layout.operator("mesh.add_cone_quadric",text="Cone",icon='CONE')
        layout.operator("mesh.add_cylinder_quadric",text="Cylinder",icon='MESH_CYLINDER')
        layout.operator("mesh.add_plane_quadric",text="Plane",icon='MESH_PLANE')
        layout.operator("mesh.add_tube_quadric",text="Tube",icon='MESH_TORUS')
        layout.operator("mesh.add_semisphere_quadric",text="Semi-Sphere",icon='SPHERE')
        layout.operator("mesh.add_semisphereshell_quadric",text="Semi-Sphere Shell",icon='SPHERE')
        layout.operator("mesh.add_coneshell_quadric",text="Cone Shell",icon='CONE')

def menu_quadric_add_func(self, context):
    self.layout.menu("OBJECT_MT_quadric_submenu", text="Quadric", icon='CUBE')        

# Object remesh operator
##########################
class REMESH_OT_remesh_operator(Operator):
    bl_idname = 'remesh.remesh_operator'
    bl_label = 'Remesh'
    bl_description = 'Remesh object'
    bl_options = {'REGISTER', 'UNDO'}
 
    def execute(self, context):
        obj = context.object
        if obj.penred_settings.quadricType == "CONE" or obj.penred_settings.quadricType == "CONE_SHELL":
            
            #Get height
            h = obj.dimensions.z
            #Get radius
            r1 = obj.penred_settings.r1
            r2 = obj.penred_settings.r2
            
            r1 = max(1.0e-5,r1)
            r2 = max(1.0e-5,r2)
            
            #Set scale to identity
            obj.scale[0] = 1.0
            obj.scale[1] = 1.0
            obj.scale[2] = 1.0
            
            ### Create an empty mesh
            mesh = bpy.data.meshes.new(name=obj.data.name)
            
            #Create the new mesh
            bm = bmesh.new()
            #Despite the name, blender uses the diameter parameters as radius OLD
            if obj.penred_settings.quadricType == "CONE":
                bmesh.ops.create_cone(bm, cap_ends=True, segments=32, radius1=r1, radius2=r2, depth = h)
            else:
                bmesh.ops.create_cone(bm, cap_ends=False, segments=32, radius1=r1, radius2=r2, depth = h)
            
            bm.to_mesh(mesh)
            bm.free()
            
            obj.penred_settings.r3 = r1
            obj.penred_settings.r4 = r2
            
            #Store old mesh
            oldMesh = obj.data
            
            #Assign new mesh
            obj.data = mesh
            
            #Remove old mesh
            if oldMesh.users == 0:
                oldMesh.user_clear()
                bpy.data.meshes.remove(oldMesh)
            
        return {'FINISHED'}

# Object read operator
##########################
class READ_OT_read_operator(Operator):
    bl_idname = 'read.read_operator'
    bl_label = 'Read'
    bl_description = 'Read object specific information'
    bl_options = {'REGISTER', 'UNDO'}
 
    def execute(self, context):
        obj = context.object
        if obj.penred_settings.quadricType == "CONE" or obj.penred_settings.quadricType == "CONE_SHELL":
            obj.penred_settings.r1 = obj.penred_settings.r3
            obj.penred_settings.r2 = obj.penred_settings.r4            
            
        return {'FINISHED'}

# Object properties menu
##########################
class PenredPropertiesPanel(bpy.types.Panel):
    """Creates a Panel in the Object properties window"""
    bl_label = "Quadric properties"
    bl_idname = "OBJECT_PT_quad"
    bl_space_type = 'PROPERTIES'
    bl_region_type = 'WINDOW'
    bl_context = "object"

    def draw(self, context):
        layout = self.layout

        obj = context.object

        row = layout.row()
        row.label(text="Quadric Type: " + obj.penred_settings.quadricType)

        row = layout.row()
        row.prop(obj.penred_settings, "material")
        
        row = layout.row()
        row.prop(obj.penred_settings, "module")
        
        if obj.penred_settings.quadricType == "CONE" or obj.penred_settings.quadricType == "CONE_SHELL":
            row = layout.row()
            row.prop(obj.penred_settings, "r2")
            row = layout.row()
            row.prop(obj.penred_settings, "r1")
            row = layout.row()
            row.operator("remesh.remesh_operator",text="Remesh")
            row.operator("read.read_operator",text="Reset")
            
        

# Generic add object function
###############################

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

#####################
#  Objects classes  # 
#####################

# CUBE
class QUADRIC_OT_add_cube(Operator, AddObjectHelper):
    """Create a new Cube Quadric"""
    bl_idname = "mesh.add_cube_quadric"
    bl_label = "Add Cube Quadric"
    bl_options = {'REGISTER', 'UNDO'}

    scale: FloatVectorProperty(
        name="scale",
        default=(1.0, 1.0, 1.0),
        subtype='TRANSLATION',
        description="scaling",
    )

    def execute(self, context):
        add_object(self, context, "Cube", "CUBE")
        return {'FINISHED'}

# SPHERE
class QUADRIC_OT_add_sphere(Operator, AddObjectHelper):
    """Create a new Spheric Quadric"""
    bl_idname = "mesh.add_sphere_quadric"
    bl_label = "Add Sphere Quadric"
    bl_options = {'REGISTER', 'UNDO'}

    scale: FloatVectorProperty(
        name="scale",
        default=(1.0, 1.0, 1.0),
        subtype='TRANSLATION',
        description="scaling",
    )

    def execute(self, context):
        add_object(self, context, "Sphere", "SPHERE")
        return {'FINISHED'}

# CONE
class QUADRIC_OT_add_cone(Operator, AddObjectHelper):
    """Create a new Cone Quadric"""
    bl_idname = "mesh.add_cone_quadric"
    bl_label = "Add Cone Quadric"
    bl_options = {'REGISTER', 'UNDO'}

    scale: FloatVectorProperty(
        name="scale",
        default=(1.0, 1.0, 1.0),
        subtype='TRANSLATION',
        description="scaling",
    )

    def execute(self, context):
        add_object(self, context, "Cone", "CONE")
        return {'FINISHED'}

# CYLINDER
class QUADRIC_OT_add_cylinder(Operator, AddObjectHelper):
    """Create a new Cylinder Quadric"""
    bl_idname = "mesh.add_cylinder_quadric"
    bl_label = "Add Cylinder Quadric"
    bl_options = {'REGISTER', 'UNDO'}

    scale: FloatVectorProperty(
        name="scale",
        default=(1.0, 1.0, 1.0),
        subtype='TRANSLATION',
        description="scaling",
    )

    def execute(self, context):
        add_object(self, context, "Cylinder", "CYLINDER")
        return {'FINISHED'}

# PLANE
class QUADRIC_OT_add_plane(Operator, AddObjectHelper):
    """Create a new Plane Quadric"""
    bl_idname = "mesh.add_plane_quadric"
    bl_label = "Add Plane Quadric"
    bl_options = {'REGISTER', 'UNDO'}

    scale: FloatVectorProperty(
        name="scale",
        default=(1.0, 1.0, 1.0),
        subtype='TRANSLATION',
        description="scaling",
    )

    def execute(self, context):
        add_object(self, context, "Plane", "PLANE")
        return {'FINISHED'}

# TUBE
class QUADRIC_OT_add_tube(Operator, AddObjectHelper):
    """Create a new Tube Quadric"""
    bl_idname = "mesh.add_tube_quadric"
    bl_label = "Add Tube Quadric"
    bl_options = {'REGISTER', 'UNDO'}

    scale: FloatVectorProperty(
        name="scale",
        default=(1.0, 1.0, 1.0),
        subtype='TRANSLATION',
        description="scaling",
    )

    def execute(self, context):
        add_object(self, context, "Tube", "TUBE")
        return {'FINISHED'}

# SEMI_SPHERE
class QUADRIC_OT_add_semiSphere(Operator, AddObjectHelper):
    """Create a new Tube Quadric"""
    bl_idname = "mesh.add_semisphere_quadric"
    bl_label = "Add SemiSphere Quadric"
    bl_options = {'REGISTER', 'UNDO'}

    scale: FloatVectorProperty(
        name="scale",
        default=(1.0, 1.0, 1.0),
        subtype='TRANSLATION',
        description="scaling",
    )

    def execute(self, context):
        add_object(self, context, "SemiSphere", "SEMI_SPHERE")
        return {'FINISHED'}

# SEMI_SPHERE_SHELL
class QUADRIC_OT_add_semiSphereShell(Operator, AddObjectHelper):
    """Create a new Tube Quadric"""
    bl_idname = "mesh.add_semisphereshell_quadric"
    bl_label = "Add SemiSphere Shell Quadric"
    bl_options = {'REGISTER', 'UNDO'}

    scale: FloatVectorProperty(
        name="scale",
        default=(1.0, 1.0, 1.0),
        subtype='TRANSLATION',
        description="scaling",
    )

    def execute(self, context):
        add_object(self, context, "SemiSphere Shell", "SEMI_SPHERE_SHELL")
        return {'FINISHED'}

# CONE_SHELL
class QUADRIC_OT_add_coneShell(Operator, AddObjectHelper):
    """Create a new Cone Quadric"""
    bl_idname = "mesh.add_coneshell_quadric"
    bl_label = "Add Cone Shell Quadric"
    bl_options = {'REGISTER', 'UNDO'}

    scale: FloatVectorProperty(
        name="scale",
        default=(1.0, 1.0, 1.0),
        subtype='TRANSLATION',
        description="scaling",
    )

    def execute(self, context):
        add_object(self, context, "Cone Shell", "CONE_SHELL")
        return {'FINISHED'}

quadricObjectClasses = (
    OBJECT_MT_quadric_submenu,
    QUADRIC_OT_add_cube,
    QUADRIC_OT_add_sphere,
    QUADRIC_OT_add_cone,
    QUADRIC_OT_add_cylinder,
    QUADRIC_OT_add_plane,
    QUADRIC_OT_add_tube,
    QUADRIC_OT_add_semiSphere,
    QUADRIC_OT_add_semiSphereShell,
    QUADRIC_OT_add_coneShell,
)

# Manual
##########
def add_object_manual_map():
    url_manual_prefix = "https://docs.blender.org/manual/en/latest/"
    url_manual_mapping = (
        ("bpy.ops.mesh.add_object", "scene_layout/object/types.html"),
    )
    return url_manual_prefix, url_manual_mapping


##############################################################################
##############################################################################
####       REGISTER/UNREGISTER FUNCTIONS
##############################################################################
##############################################################################

# Register function
#####################
def register():
    #Register properties class
    bpy.utils.register_class(PenredProperties)
    bpy.types.Object.penred_settings = bpy.props.PointerProperty(type=PenredProperties)
    
    #Register object properties UI
    bpy.utils.register_class(PenredPropertiesPanel)

    #Register View classes
    for cls in viewObjectClasses:
        bpy.utils.register_class(cls)    
    
    #Append view menu
    bpy.types.VIEW3D_MT_view.append(menu_quadric_view_func)
    
    #Register object classes
    for cls in quadricObjectClasses:
        bpy.utils.register_class(cls)    

    #Append "add" menu
    bpy.types.VIEW3D_MT_add.append(menu_quadric_add_func)
    
    #Register export class
    bpy.utils.register_class(export_penred)
    
    #Append export menu
    bpy.types.TOPBAR_MT_file_export.append(menu_func_penred_export)
    
    #Register remesh class
    bpy.utils.register_class(REMESH_OT_remesh_operator)
    
    #Register read object information class
    bpy.utils.register_class(READ_OT_read_operator)
    
    #Register transformation classes
    for cls in quadricTransformClasses:
        bpy.utils.register_class(cls)
            
    #Append object transformation button menu
    bpy.types.VIEW3D_MT_transform_object.append(menu_func_quadric_transform)
    
    bpy.utils.register_manual_map(add_object_manual_map)


# Unregister function
######################
def unregister():
    #Unregister properties class
    bpy.utils.unregister_class(PenredProperties)

    #Unregister object properties UI
    bpy.utils.unregister_class(PenredPropertiesPanel)

    #Unregister view classes
    for cls in viewObjectClasses:
        bpy.utils.unregister_class(cls)    

    #Remove view menu
    bpy.types.VIEW3D_MT_view.remove(menu_quadric_view_func)

    #Unregister object classes
    for cls in quadricObjectClasses:
        bpy.utils.unregister_class(cls)    

    #Remove "add" menu
    bpy.types.VIEW3D_MT_add.remove(menu_quadric_add_func)

    #Unregister export class
    bpy.utils.unregister_class(export_penred)

    #Remove export menu
    bpy.types.TOPBAR_MT_file_export.remove(menu_func_penred_export)

    #Unregister remesh class
    bpy.utils.unregister_class(REMESH_OT_remesh_operator)

    #Unregister read object information class
    bpy.utils.unregister_class(READ_OT_read_operator)

    #Unregister transformation classes
    for cls in quadricTransformClasses:
        bpy.utils.unregister_class(cls)
        
    #Remove object transformation button menu
    bpy.types.VIEW3D_MT_transform_object.remove(menu_func_quadric_transform)

    bpy.utils.unregister_manual_map(add_object_manual_map)


if __name__ == "__main__":
    register()
