                                                                       
A right prism (do NOT use double planes!)                              
                                                                       
0000000000000000000000000000000000000000000000000000000000000000
SURFACE (   1)   Plane Z=-2.0                                       
INDICES=( 0, 0, 0, 1, 0)
Z-SHIFT=(-2.000000000000000E+00,   0)       
0000000000000000000000000000000000000000000000000000000000000000
SURFACE (   2)   Plane Z=+2.0                                       
INDICES=( 0, 0, 0, 1, 0)
Z-SHIFT=( 2.000000000000000E+00,   0)       
0000000000000000000000000000000000000000000000000000000000000000
SURFACE (   3)   Plane X=-2.0                                       
INDICES=( 0, 0, 0, 0, 0)
     AX=( 1.000000000000000E+00,   0)       
     A0=( 2.000000000000000E+00,   0)       
0000000000000000000000000000000000000000000000000000000000000000
SURFACE (   4)   Plane X=+2.0                                       
INDICES=( 0, 0, 0, 0, 0)
     AX=( 1.000000000000000E+00,   0)       
     A0=(-2.000000000000000E+00,   0)       
0000000000000000000000000000000000000000000000000000000000000000
SURFACE (   5)   Plane Y=-2.0                                       
INDICES=( 0, 0, 0, 0, 0)
     AY=( 1.000000000000000E+00,   0)       
     A0=( 2.000000000000000E+00,   0)       
0000000000000000000000000000000000000000000000000000000000000000
SURFACE (   6)   Plane Y=+2.0                                       
INDICES=( 0, 0, 0, 0, 0)
     AY=( 1.000000000000000E+00,   0)       
     A0=(-2.000000000000000E+00,   0)       
0000000000000000000000000000000000000000000000000000000000000000
MODULE  (   1)   Solid box                                          
MATERIAL(   1)
SURFACE (   1), SIDE POINTER=( 1)
SURFACE (   2), SIDE POINTER=(-1)
SURFACE (   3), SIDE POINTER=( 1)
SURFACE (   4), SIDE POINTER=(-1)
SURFACE (   5), SIDE POINTER=( 1)
SURFACE (   6), SIDE POINTER=(-1)
0000000000000000000000000000000000000000000000000000000000000000
END      0000000000000000000000000000000000000000000000000000000



*****************************************
****     PENGEOM (version 2014)      ****
****  Constructive Quadric Geometry  ****
*****************************************



************  Genealogical tree. 


*** MODULE =    1,  KMOTH =    0,  MAT =  1
KDGHT =    1
KSURF =    1    2    3    4    5    6
KFLAG =    2    1    2    1    2    1

The module    1 is the enclosure.


************  Adequacy of the geometry definition.

The largest number of bodies in a module or
     bodies limiting a single body is ............    1

The largest number of limiting surfaces for
     a single body or module is ..................    6


The simulation of this geometry will be relatively fast,
     no further optimization seems to be required.

************  The end.
