XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
   A simple model of electron accelerator head and a water
   phantom. Illustrates the use of the program 'penmain'.

   The phase-space file at the exit of the accelerator head
   is generated using a planar impact detector.

   Materials: 1.- Tungsten
              2.- Tungsten
              3.- Air
              4.- Water

C
C  ****  Accelerator head.
C
0000000000000000000000000000000000000000000000000000000000000000
SURFACE (   1)   Plane z=0
INDICES=( 0, 0, 0, 1, 0)
0000000000000000000000000000000000000000000000000000000000000000
SURFACE (   2)   Plane z=0.2
INDICES=( 0, 0, 0, 1, 0)
Z-SHIFT=( 2.500000000000000E-01,   0)
0000000000000000000000000000000000000000000000000000000000000000
SURFACE (   3)   Cylinder, r=0.75
INDICES=( 1, 1, 0, 0,-1)
X-SCALE=( 7.500000000000000E-01,   0)
Y-SCALE=( 7.500000000000000E-01,   0)
Z-SCALE=( 7.500000000000000E-01,   0)
0000000000000000000000000000000000000000000000000000000000000000
MODULE  (   1)   Tungsten slab (target)
MATERIAL(   1)
SURFACE (   1), SIDE POINTER=( 1)
SURFACE (   2), SIDE POINTER=(-1)
SURFACE (   3), SIDE POINTER=(-1)
0000000000000000000000000000000000000000000000000000000000000000
SURFACE (   4)   Plane z=0.5
INDICES=( 0, 0, 0, 1, 0)
Z-SHIFT=( 5.000000000000000E-01,   0)
0000000000000000000000000000000000000000000000000000000000000000
SURFACE (   5)   Plane z=5
INDICES=( 0, 0, 0, 1, 0)
Z-SHIFT=( 5.000000000000000E+00,   0)
0000000000000000000000000000000000000000000000000000000000000000
SURFACE (   6)   Cone
INDICES=( 1, 1,-1, 0, 0)
X-SCALE=( 2.000000000000000E-01,   0)
Y-SCALE=( 2.000000000000000E-01,   0)
Z-SHIFT=(-6.000000000000000E-01,   0)
0000000000000000000000000000000000000000000000000000000000000000
SURFACE (   7)   Cylinder, r=8.0
INDICES=( 1, 1, 0, 0,-1)
X-SCALE=( 8.000000000000000E+00,   0)
Y-SCALE=( 8.000000000000000E+00,   0)
0000000000000000000000000000000000000000000000000000000000000000
MODULE  (   2)   Tungsten collimator
MATERIAL(   2)
SURFACE (   4), SIDE POINTER=( 1)
SURFACE (   5), SIDE POINTER=(-1)
SURFACE (   6), SIDE POINTER=( 1)
SURFACE (   7), SIDE POINTER=(-1)
0000000000000000000000000000000000000000000000000000000000000000
C
C  ****  Detector.
C
SURFACE (   8)   Plane z=6
INDICES=( 0, 0, 0, 1, 0)
Z-SHIFT=( 6.000000000000000E+00,   0)
0000000000000000000000000000000000000000000000000000000000000000
SURFACE (   9)   Plane z=6.1
INDICES=( 0, 0, 0, 1, 0)
Z-SHIFT=( 6.100000000000000E+00,   0)
0000000000000000000000000000000000000000000000000000000000000000
SURFACE (  10)   Cylinder, r=15.0
INDICES=( 1, 1, 0, 0,-1)
X-SCALE=( 1.500000000000000E+01,   0)
Y-SCALE=( 1.500000000000000E+01,   0)
0000000000000000000000000000000000000000000000000000000000000000
MODULE  (   3)   Detector volume
MATERIAL(   3)
SURFACE (   8), SIDE POINTER=( 1)
SURFACE (   9), SIDE POINTER=(-1)
SURFACE (  10), SIDE POINTER=(-1)
0000000000000000000000000000000000000000000000000000000000000000
C
C  ****  Water phantom.
C
SURFACE (  11)   Plane Z=10
INDICES=( 0, 0, 0, 1, 0)
Z-SHIFT=( 1.000000000000000E+01,   0)
0000000000000000000000000000000000000000000000000000000000000000
SURFACE (  12)   Plane Z=12
INDICES=( 0, 0, 0, 1, 0)
Z-SHIFT=( 1.200000000000000E+01,   0)
0000000000000000000000000000000000000000000000000000000000000000
MODULE  (   4)   Water phantom slice
MATERIAL(   4)
SURFACE (  11), SIDE POINTER=( 1)
SURFACE (  12), SIDE POINTER=(-1)
SURFACE (  10), SIDE POINTER=(-1)
0000000000000000000000000000000000000000000000000000000000000000
SURFACE (  13)   Plane Z=14
INDICES=( 0, 0, 0, 1, 0)
Z-SHIFT=( 1.400000000000000E+01,   0)
0000000000000000000000000000000000000000000000000000000000000000
MODULE  (   5)   Water phantom slice
MATERIAL(   4)
SURFACE (  12), SIDE POINTER=( 1)
SURFACE (  13), SIDE POINTER=(-1)
SURFACE (  10), SIDE POINTER=(-1)
0000000000000000000000000000000000000000000000000000000000000000
SURFACE (  14)   Plane Z=16
INDICES=( 0, 0, 0, 1, 0)
Z-SHIFT=( 1.600000000000000E+01,   0)
0000000000000000000000000000000000000000000000000000000000000000
MODULE  (   6)   Water phantom slice
MATERIAL(   4)
SURFACE (  13), SIDE POINTER=( 1)
SURFACE (  14), SIDE POINTER=(-1)
SURFACE (  10), SIDE POINTER=(-1)
0000000000000000000000000000000000000000000000000000000000000000
SURFACE (  15)   Plane Z=18
INDICES=( 0, 0, 0, 1, 0)
Z-SHIFT=( 1.800000000000000E+01,   0)
0000000000000000000000000000000000000000000000000000000000000000
MODULE  (   7)   Water phantom slice
MATERIAL(   4)
SURFACE (  14), SIDE POINTER=( 1)
SURFACE (  15), SIDE POINTER=(-1)
SURFACE (  10), SIDE POINTER=(-1)
0000000000000000000000000000000000000000000000000000000000000000
SURFACE (  16)   Plane Z=20
INDICES=( 0, 0, 0, 1, 0)
Z-SHIFT=( 2.000000000000000E+01,   0)
0000000000000000000000000000000000000000000000000000000000000000
MODULE  (   8)   Water phantom slice
MATERIAL(   4)
SURFACE (  15), SIDE POINTER=( 1)
SURFACE (  16), SIDE POINTER=(-1)
SURFACE (  10), SIDE POINTER=(-1)
0000000000000000000000000000000000000000000000000000000000000000
SURFACE (  17)   Plane Z=22
INDICES=( 0, 0, 0, 1, 0)
Z-SHIFT=( 2.200000000000000E+01,   0)
0000000000000000000000000000000000000000000000000000000000000000
MODULE  (   9)   Water phantom slice
MATERIAL(   4)
SURFACE (  16), SIDE POINTER=( 1)
SURFACE (  17), SIDE POINTER=(-1)
SURFACE (  10), SIDE POINTER=(-1)
0000000000000000000000000000000000000000000000000000000000000000
SURFACE (  18)   Plane Z=24
INDICES=( 0, 0, 0, 1, 0)
Z-SHIFT=( 2.400000000000000E+01,   0)
0000000000000000000000000000000000000000000000000000000000000000
MODULE  (  10)   Water phantom slice
MATERIAL(   4)
SURFACE (  17), SIDE POINTER=( 1)
SURFACE (  18), SIDE POINTER=(-1)
SURFACE (  10), SIDE POINTER=(-1)
0000000000000000000000000000000000000000000000000000000000000000
SURFACE (  19)   Plane Z=26
INDICES=( 0, 0, 0, 1, 0)
Z-SHIFT=( 2.600000000000000E+01,   0)
0000000000000000000000000000000000000000000000000000000000000000
MODULE  (  11)   Water phantom slice
MATERIAL(   4)
SURFACE (  18), SIDE POINTER=( 1)
SURFACE (  19), SIDE POINTER=(-1)
SURFACE (  10), SIDE POINTER=(-1)
0000000000000000000000000000000000000000000000000000000000000000
SURFACE (  20)   Plane Z=28
INDICES=( 0, 0, 0, 1, 0)
Z-SHIFT=( 2.800000000000000E+01,   0)
0000000000000000000000000000000000000000000000000000000000000000
MODULE  (  12)   Water phantom slice
MATERIAL(   4)
SURFACE (  19), SIDE POINTER=( 1)
SURFACE (  20), SIDE POINTER=(-1)
SURFACE (  10), SIDE POINTER=(-1)
0000000000000000000000000000000000000000000000000000000000000000
SURFACE (  21)   Plane Z=30
INDICES=( 0, 0, 0, 1, 0)
Z-SHIFT=( 3.000000000000000E+01,   0)
0000000000000000000000000000000000000000000000000000000000000000
MODULE  (  13)   Water phantom slice
MATERIAL(   4)
SURFACE (  20), SIDE POINTER=( 1)
SURFACE (  21), SIDE POINTER=(-1)
SURFACE (  10), SIDE POINTER=(-1)
0000000000000000000000000000000000000000000000000000000000000000
C
C  The complete system is defined as a single module, which is
C  filled with air.
C
SURFACE (  22)   Plane Z=40
INDICES=( 0, 0, 0, 1, 0)
Z-SHIFT=( 4.000000000000000E+01,   0)
0000000000000000000000000000000000000000000000000000000000000000
MODULE  (  14)   Enclosure
MATERIAL(   3)
SURFACE (   1), SIDE POINTER=( 1)
SURFACE (  22), SIDE POINTER=(-1)
SURFACE (  10), SIDE POINTER=(-1)
MODULE  (   1)
MODULE  (   2)
MODULE  (   3)
MODULE  (   4)
MODULE  (   5)
MODULE  (   6)
MODULE  (   7)
MODULE  (   8)
MODULE  (   9)
MODULE  (  10)
MODULE  (  11)
MODULE  (  12)
MODULE  (  13)
1111111111111111111111111111111111111111111111111111111111111111
  OMEGA=( 2.000000000000000E+01,   0) DEG
  THETA=( 4.500000000000000E+01,   0) DEG
    PHI=( 2.000000000000000E+01,   0) DEG
X-SHIFT=( 0.000000000000000E+00,   0)
Y-SHIFT=( 1.000000000000000E+01,   0)
Z-SHIFT=( 0.000000000000000E+00,   0)
0000000000000000000000000000000000000000000000000000000000000000
END      0000000000000000000000000000000000000000000000000000000