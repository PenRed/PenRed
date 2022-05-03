import sys
import numpy as np

if len ( sys.argv ) != 4:
    print( "Usage: %s spatialdosefile headerdatafile outputfile" % ( sys.argv[0] ) )
    sys.exit ( 1 )

dosefilename=str(sys.argv[1])
headerdatafile=str(sys.argv[2])
outputfile=str(sys.argv[3])

# Read image parameters from headerdatafile
(Ox,Oy,Oz) = np.loadtxt(headerdatafile,skiprows=6,max_rows=1,usecols=(2,3,4),dtype='float')
(dx,dy,dz) = np.loadtxt(headerdatafile,skiprows=9,max_rows=1,usecols=(2,3,4),dtype='float')
(nx,ny,nz) = np.loadtxt(headerdatafile,skiprows=10,max_rows=1,usecols=(2,3,4),dtype='int')
nbins=nx*ny*nz

print("Processing ",dosefilename,"with image parameters from",headerdatafile)

with open(dosefilename,'r') as f:
  for line in f.readlines():
    li = line.lstrip()
    if li.startswith("#") or li.startswith("\n"):
       continue
    else:
       break
  f.close()
ncols=len(line.split())
print(ncols," columns found")
if ncols==2:
   readcol=0
elif ncols==5:
   readcol=3
elif ncols==11:
   readcol=9
print("Assuming the dose data is in column ",readcol+1)

factorz=dx*dy/100.
factorx=dy*dz/100.
factory=dx*dz/100.

print("Loading spatial dose ...")
print("   Processing ",dosefilename)
sum  = np.loadtxt(dosefilename,usecols=readcol)
sume = np.loadtxt(dosefilename,usecols=readcol+1)
sume = sume/2.0
sum  = sum.reshape(nz,ny,nx)
sume = sume.reshape(nz,ny,nx)
print(" ... done")

print("Processing depth dose on X axis ...")
with open(outputfile+'-x.dat', 'w') as f:
     f.write("Dose distribution information:\n")
     f.write("range x = {}\n".format(sum.shape[2]))
     f.write("\n")
     f.write("======================================\n")
     f.write("IndexX    Dose (eV/g/cm^2 per history)\n")
     f.write("======================================\n")
     sum_x  = np.sum(sum,axis=(0,1))*factorx
     sume_x = np.sqrt(np.sum(np.square(sume),axis=(0,1)))*factorx
     for xslice in range(nx):
         doseSum  = sum_x[xslice]
         edoseSum = sume_x[xslice]
         f.write("{}       {:.6e}  {:.6e}\n".format(xslice,doseSum,edoseSum))
     f.close()
print(" ... done")

print("Processing depth dose on Y axis ...")
with open(outputfile+'-y.dat', 'w') as f:
     f.write("Dose distribution information:\n")
     f.write("range y = {}\n".format(sum.shape[1]))
     f.write("\n")
     f.write("======================================\n")
     f.write("IndexY    Dose (eV/g/cm^2 per history)\n")
     f.write("======================================\n")
     sum_y  = np.sum(sum,axis=(0,2))*factory
     sume_y = np.sqrt(np.sum(np.square(sume),axis=(0,2)))*factory
     for yslice in range(ny):
         doseSum  = sum_y[yslice]
         edoseSum = sume_y[yslice]
         f.write("{}       {:.6e}  {:.6e}\n".format(yslice,doseSum,edoseSum))
     f.close()
print(" ... done")

print("Processing depth dose on Z axis ...")
with open(outputfile+'-z.dat', 'w') as f:   
     f.write("Dose distribution information:\n")
     f.write("range z = {}\n".format(sum.shape[0]))
     f.write("\n")
     f.write("======================================\n")
     f.write("IndexZ    Dose (eV/g/cm^2 per history)\n")
     f.write("======================================\n")
     sum_z  = np.sum(sum,axis=(1,2))*factorz
     sume_z = np.sqrt(np.sum(np.square(sume),axis=(1,2)))*factorz
     for zslice in range(nz):
         doseSum  = sum_z[zslice]
         edoseSum = sume_z[zslice]
         f.write("{}       {:.6e}  {:.6e}\n".format(zslice,doseSum,edoseSum))
     f.close()  
print(" ... done")
print("Goodbye")

