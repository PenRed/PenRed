
//
//
//    Copyright (C) 2019 Universitat de València - UV
//    Copyright (C) 2019 Universitat Politècnica de València - UPV
//
//    This file is part of PenRed: Parallel Engine for Radiation Energy Deposition.
//
//    PenRed is free software: you can redistribute it and/or modify
//    it under the terms of the GNU Affero General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    PenRed is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU Affero General Public License for more details.
//
//    You should have received a copy of the GNU Affero General Public License
//    along with PenRed.  If not, see <https://www.gnu.org/licenses/>. 
//
//    contact emails:
//
//        vicent.gimenez.alventosa@gmail.com
//        vicente.gimenez@uv.es
//    
//


#include <cstdio>
#include "pen_dicom.hh"


int main(int argc, char** argv){

  if(argc < 2){
    printf("usage: %s dicomFolder\n",argv[0]);
    return 1;
  }
  
  pen_dicom dicom;

  //Load dicom directory
  int err = dicom.loadDicom(argv[1],3);
  if(err != PEN_DICOM_SUCCESS){
    printf("Error loading dicoms.\n");
    printf("   Error code: %d\n",err);
    return -1;
  }

  //Assign contours without priorities
  err = dicom.assignContours();
  if(err != PEN_DICOM_SUCCESS){
    printf("Error assigning contours.\n");
    printf("   Error code: %d\n",err);
    return -2;
  }

  //Print contour information
  printf("Printing contours...\n");
  dicom.printContours("contour.dat");

  //Print seed information
  printf("Printing seeds...\n");
  dicom.printSeeds("seeds.dat");

  //Print image
  printf("Printing image...\n");
  dicom.printImage("image.dat");

  //Print contour assign
  printf("Printing contour assign...\n");
  dicom.printContourVox("imageContour.dat");

  printf("Done!\n");
  return 0;
}
