import SimpleITK as sitk
import sys
import numpy as np

#nz=209; ny=128; nx=128
#dx=dy=dz=4.41806
#Ox=-281.857
#Oy=-279.547
#Oz=-460.547

if len ( sys.argv ) != 4 and len ( sys.argv ) != 5:
    print( "Usage: %s sourceactivityfile headerdatafile outputfile [asciioutfile(0:1)]" % ( sys.argv[0] ) )
    sys.exit ( 1 )

filename=str(sys.argv[1])
headerdatafile=str(sys.argv[2])
outputfile=str(sys.argv[3])

# Read image parameters from headerdatafile
(Ox,Oy,Oz) = np.loadtxt(headerdatafile,skiprows=6,max_rows=1,usecols=(2,3,4),dtype='float')
(dx,dy,dz) = np.loadtxt(headerdatafile,skiprows=9,max_rows=1,usecols=(2,3,4),dtype='float')
(nx,ny,nz) = np.loadtxt(headerdatafile,skiprows=10,max_rows=1,usecols=(2,3,4),dtype='int')

isaveres=0
if len(sys.argv) == 5:
   isaveres=int(sys.argv[4])

dosearray2=np.zeros((nz,ny,nx))

print("Processing ",filename,"with image parameters from",headerdatafile)
dosearray  = np.loadtxt(filename,usecols=0)
voxz  = np.loadtxt(filename,usecols=1,dtype='float').astype(int)
voxy  = np.loadtxt(filename,usecols=2,dtype='float').astype(int)
voxx  = np.loadtxt(filename,usecols=3,dtype='float').astype(int)
print("Array voxz: ",voxz)
print("Dosearray ndim : ",dosearray.ndim)
print("Dosearray shape: ",dosearray.shape)
print("Reshaping ...")
for ivox in range(dosearray.shape[0]):
    dosearray2[voxz[ivox],voxy[ivox],voxx[ivox]]=dosearray[ivox]
print("Dosearray ndim : ",dosearray2.ndim)
print("Dosearray shape: ",dosearray2.shape)
if isaveres==1:
 print("Creating image ascii file ",outputfile+'.dat')
 with open(outputfile+'.dat', 'w') as f:
     nbin=0
     for iz in range(nz):
        for iy in range(ny):
           for ix in range(nx):
                f.write("   {}     {}     {}     {:12.6E} \n".format(ix,iy,iz,dosearray2[iz,iy,ix]))
                nbin=nbin+1
     f.close()
print("Creating image mhd/raw file ",outputfile+'.mhd')
print("Getting image from numpy array ...")
img = sitk.GetImageFromArray(dosearray2)
img.SetSpacing([dx,dy,dz])
img.SetOrigin([Ox,Oy,Oz])
print("Image dim      : ",img.GetDimension())
print("Image size     : ",img.GetSize())
print("Image origin   : ",img.GetOrigin())
print("Image spacing  : ",img.GetSpacing())
print("Image direction: ",img.GetDirection())
sitk.WriteImage(img,outputfile+".mhd")
