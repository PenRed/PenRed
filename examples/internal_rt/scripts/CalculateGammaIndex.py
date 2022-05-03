import SimpleITK as sitk
import sys
import numpy as np


if len ( sys.argv ) != 6:
    print( "Usage: %s GammaInputImage dosedistance dosedistanceunit distancetoagreement threshold" % ( sys.argv[0] ) )
    sys.exit ( 1 )

filename=str(sys.argv[1])
dd=str(sys.argv[2])
ddunit=str(sys.argv[3])
dta=str(sys.argv[4])
T=str(sys.argv[5])

# Reads the image using SimpleITK
itkimage = sitk.ReadImage(filename,sitk.sitkFloat32)

# Convert the image to a  numpy array 
ct_scans = sitk.GetArrayFromImage(itkimage)

f=sys.stdout
f.write("Image information:\n")
f.write("range= %d\n" % ct_scans.shape[2])
size = np.array(list(itkimage.GetSize()))
f.write("Size= {}\n".format(size))
# Read the origin of the ct_scan, will be used to convert the coordinates from world to voxel and vice versa.
origin = np.array(list(itkimage.GetOrigin()))
f.write("Origin= {}\n".format(origin))
# Read the spacing along each dimension
spacing = np.array(list(itkimage.GetSpacing()))
f.write("Spacing= {}\n".format(spacing))
f.write("\n")

f.write("======================================\n")
f.write("Gamma index parameters:\n")
f.write("(see [Daniel Low, 1998] for details\n")
f.write(" or/and read help gt_gamma_index -h)\n")
f.write("======================================\n")
f.write("Dose distance        : {}\n".format(dd))
f.write("Dose distance units  : {}\n".format(ddunit))
f.write("Distance to agreement: {}\n".format(dta))
f.write("Threshold            : {}\n".format(T))
f.write("\n")

f.write("======================================\n")
f.write("Gamma Index results\n")
f.write("======================================\n")
passed=0;notpassed=0;underthreshold=0
for z in range(ct_scans.shape[0]):#
    for y in range(ct_scans.shape[1]):
        for x in range(ct_scans.shape[2]):
            if ct_scans[z,y,x] >= 1.0:
                notpassed=notpassed+1
            elif ct_scans[z,y,x] >= 0 and ct_scans[z,y,x] < 1.0:
                passed=passed+1
            else:
                underthreshold=underthreshold+1

f.write("Passed         : {:8d} {:10.4f}%\n".format(passed,passed/(passed+notpassed)*100))
f.write("Not passed     : {:8d} {:10.4f}%\n".format(notpassed,notpassed/(passed+notpassed)*100))
f.write("Under threshold: {:8d} {:10.4f}%\n".format(underthreshold,underthreshold/(passed+notpassed+underthreshold)*100))
f.close()
