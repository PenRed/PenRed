0000000000000000000000000000000000000000000000000000000000000000
SURFACE (   1)   Plane Z=1.00
INDICES=( 0, 0, 0, 1, 0)
Z-SHIFT=(+1.0000000000000000000,   0)
0000000000000000000000000000000000000000000000000000000000000000
SURFACE (   2)   Plane Z=1.005
INDICES=( 0, 0, 0, 1, 0)
Z-SHIFT=(+1.1000000000000000000,   0)
0000000000000000000000000000000000000000000000000000000000000000
SURFACE (   3)   Cylinder R=0.1
INDICES=( 1, 1, 0, 0,-1)
X-SCALE=(+0.500000000000000E+00,   0)
Y-SCALE=(+0.800000000000000E+00,   0)
0000000000000000000000000000000000000000000000000000000000000000
MODULE  (   1)   Solid cylinder
MATERIAL(   1)
SURFACE (   1), SIDE POINTER=(+1)
SURFACE (   2), SIDE POINTER=(-1)
SURFACE (   3), SIDE POINTER=(-1)
0000000000000000000000000000000000000000000000000000000000000000
SURFACE (   4)   Plane Z=-10.0
INDICES=( 0, 0, 0, 1, 0)
Z-SHIFT=(-1.000000000000000E+01,   0)              (DEFAULT=0.0)
0000000000000000000000000000000000000000000000000000000000000000
MODULE  (   2)   detector
MATERIAL(   1)
SURFACE (   4), SIDE POINTER=(-1)
0000000000000000000000000000000000000000000000000000000000000000
END      0000000000000000000000000000000000000000000000000000000
