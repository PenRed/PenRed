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
