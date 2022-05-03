import SimpleITK as sitk
import sys
import numpy as np

if len ( sys.argv ) != 4:
    print( "Usage: %s inputImage1 inputImage2 outputText" % ( sys.argv[0] ) )
    sys.exit ( 1 )

filename1=str(sys.argv[1])
filename2=str(sys.argv[2])
outfile=str(sys.argv[3])

# Reads the image using SimpleITK
itkimage1 = sitk.ReadImage(filename1,sitk.sitkFloat32)
itkimage2 = sitk.ReadImage(filename2,sitk.sitkFloat32)

# Convert the image to a  numpy array 
ct_scans1 = sitk.GetArrayFromImage(itkimage1)
ct_scans2 = sitk.GetArrayFromImage(itkimage2)

with open(outfile, 'w') as f:
     f.write("# Images information:\n")
     f.write("# range= {:3d}  {:3d}\n".format(ct_scans1.shape[2],ct_scans2.shape[2]))
     size1 = np.array(list(itkimage1.GetSize()))
     size2 = np.array(list(itkimage2.GetSize()))
     f.write("# Size= {}  {}\n".format(size1,size2))
     # Read the origin of the ct_scan, will be used to convert the coordinates from world to voxel and vice versa.
     origin1 = np.array(list(itkimage1.GetOrigin()))
     origin2 = np.array(list(itkimage2.GetOrigin()))
     f.write("# Origin= {}  {}\n".format(origin1,origin2))
     # Read the spacing along each dimension
     spacing1 = np.array(list(itkimage1.GetSpacing()))
     spacing2 = np.array(list(itkimage2.GetSpacing()))
     f.write("# Spacing= {}  {}\n".format(spacing1,spacing2))
     f.write("# \n")
     f.write("# ============================================\n")
     f.write("# Ix  Iy  Iz    Pixel value 1   Pixel value 2\n")
     f.write("# ============================================\n")
     nbin=0
     for z in range(ct_scans1.shape[0]):
         for y in range(ct_scans1.shape[1]):
             for x in range(ct_scans1.shape[2]):
                 nbin=nbin+1
                 f.write(" {:3d} {:3d} {:3d}    {:10.6e}    {:10.6e}\n".format(x,y,z,ct_scans1[z,y,x],ct_scans2[z,y,x]))
     f.close()
