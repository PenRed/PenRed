>>> GEOMETRY FILE FOR THE PENELOPE SYSTEM >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

Manufacturer   : AAPM (virtual source)
Isotope        : Ir-192
Source         : MBDCA Source model
Consensus data : 

Author         : Javier Vijande (javier.vijande@uv.es)

Material:

1 - air
2 - Ir: Density 22.42 g/cm^3
3 - capsule: Stainless Steel 316L density 8.02 g/cm^3
4 - cable: Stainless Steel 316L density 5.00 g/cm^3

Disclaimer and Copyright

This geometry is released “as is”. No warranties, express or implied, that this geometry is free of error, or is consistent with any particular standard of merchantability, 
or that it will meet your requirements for any particular application, is made. No responsibility for any mathematical or technical limitations of the procedures and functions is accepted. 
This geometry should not be relied on for solving a problem whose incorrect solution could result in injury to a person or loss of property. The authors of this geometry shall not in any event 
be liable for any damages, whether direct or indirect, special or general, consequential or incidental, arising from its use. Your use of this software is entirely at your own risk. 
Permission to use this geometry for any purpose is hereby granted without fee. This geometry, or any part of it, cannot be sold, modified or re-distributed unless a written consent from the author is obtained.

Anyone using this seed for research purposes is kindly requested to cite the web site of the Medical Physic Group (www.uv.es/radiofisica).

>>> END OF HEADER >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

0000000000000000000000000000000000000000000000000000000000000000
C
C Active part
C
SURFACE (   1) Plane limiting upper active part
INDICES=( 0, 0, 0, 1, 0)
Z-SHIFT=(+0.175000000000000E+00,   0)              (DEFAULT=0.0)
0000000000000000000000000000000000000000000000000000000000000000
SURFACE (   2) Plane limiting upper active part
INDICES=( 0, 0, 0, 1, 0)
Z-SHIFT=(-0.175000000000000E+00,   0)              (DEFAULT=0.0)
0000000000000000000000000000000000000000000000000000000000000000
SURFACE (   3)  Cylinder core
INDICES=( 1, 1, 0, 0,-1)
X-SCALE=(+0.030000000000000E+00,   0)              (DEFAULT=1.0)
Y-SCALE=(+0.030000000000000E+00,   0)              (DEFAULT=1.0)
0000000000000000000000000000000000000000000000000000000000000000
BODY    (   1)   Active part.
MATERIAL(   2)
SURFACE (   1), SIDE POINTER=(-1)
SURFACE (   2), SIDE POINTER=(+1)
SURFACE (   3), SIDE POINTER=(-1)
0000000000000000000000000000000000000000000000000000000000000000
C
C Capsule
C
SURFACE (   4)   Sphere core 
INDICES=( 1, 1, 1, 0,-1)
X-SCALE=(+0.050000000000000E+00,   0)              (DEFAULT=1.0)
Y-SCALE=(+0.050000000000000E+00,   0)              (DEFAULT=1.0)
Z-SCALE=(+0.050000000000000E+00,   0)              (DEFAULT=1.0)
Z-SHIFT=(+0.185000000000000E+00,   0)              (DEFAULT=0.0)
0000000000000000000000000000000000000000000000000000000000000000
SURFACE (   5) Sphere base
INDICES=( 0, 0, 0, 1, 0)
Z-SHIFT=(+0.185000000000000E+00,   0)              (DEFAULT=0.0)
0000000000000000000000000000000000000000000000000000000000000000
SURFACE (   6) Cable top
INDICES=( 0, 0, 0, 1, 0)
Z-SHIFT=(-0.265000000000000E+00,   0)              (DEFAULT=0.0)
0000000000000000000000000000000000000000000000000000000000000000
SURFACE (   7)  Cylinder capsule
INDICES=( 1, 1, 0, 0,-1)
X-SCALE=(+0.050000000000000E+00,   0)              (DEFAULT=1.0)
Y-SCALE=(+0.050000000000000E+00,   0)              (DEFAULT=1.0)
0000000000000000000000000000000000000000000000000000000000000000
MODULE  (   2)   Capsule
MATERIAL(   3)
SURFACE (   5), SIDE POINTER=(-1)
SURFACE (   6), SIDE POINTER=(+1)
SURFACE (   7), SIDE POINTER=(-1)
BODY    (   1)
0000000000000000000000000000000000000000000000000000000000000000
BODY    (   3)   Capsule
MATERIAL(   3)
SURFACE (   4), SIDE POINTER=(-1)
SURFACE (   5), SIDE POINTER=(+1)
0000000000000000000000000000000000000000000000000000000000000000
C
C Cable
C
SURFACE (   8) Cable base
INDICES=( 0, 0, 0, 1, 0)
Z-SHIFT=(-0.465000000000000E+00,   0)              (DEFAULT=0.0)
0000000000000000000000000000000000000000000000000000000000000000
BODY    (   4)   Capsule
MATERIAL(   4)
SURFACE (   6), SIDE POINTER=(-1)
SURFACE (   7), SIDE POINTER=(-1)
SURFACE (   8), SIDE POINTER=(+1)
0000000000000000000000000000000000000000000000000000000000000000
C
C  air enclosure to obtain PSF
SURFACE (  18)   sphere, R=0.75
INDICES=( 1, 1, 1, 0,-1)
X-SCALE=( 0.5000000000000000E+00,   0)
Y-SCALE=( 0.5000000000000000E+00,   0)
Z-SCALE=( 0.5000000000000000E+00,   0)
Z-SHIFT=(-0.2000000000000000E+00,   0)
0000000000000000000000000000000000000000000000000000000000000000
MODULE  (   5)   sphere
MATERIAL(   1)
SURFACE (  18), SIDE POINTER=(-1)
MODULE  (   2)
BODY    (   3)   
BODY    (   4)   
0000000000000000000000000000000000000000000000000000000000000000
END      0000000000000000000000000000000000000000000000000000000
