

//
//
//    Copyright (C) 2022 Universitat de València - UV
//    Copyright (C) 2022 Universitat Politècnica de València - UPV
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
//        vicent.gimenez.alventosa@gmail.com  (Vicent Giménez Alventosa)
//        sanolgi@upvnet.upv.es (Sandra Oliver Gil)
//        vicente.gimenez@uv.es (Vicente Giménez Gómez)
//    
//

#include "pen_image.hh"

// ** Constructors for images without asociated error

pen_imageExporter::pen_imageExporter(std::function<float(unsigned long long, size_t)> fin) :
  nDimensions(1),
  nElements{0},
  elementSizes{0},
  origin{0},
  fbytes(std::bind(pen_imageExporter::conv2byte<float>,
		   std::placeholders::_1,
		   std::placeholders::_2,
		   std::placeholders::_3,
		   std::placeholders::_4,
		   std::placeholders::_5,
		   std::placeholders::_6, fin)),
  fASCII(std::bind(pen_imageExporter::conv2ASCII<float>,
		   std::placeholders::_1,
		   std::placeholders::_2,
		   std::placeholders::_3,
		   std::placeholders::_4,
		   std::placeholders::_5,fin)),    
  typeSize(sizeof(float)),
  integer(false),
  withSign(true),
  baseName("none")    
{}
  
pen_imageExporter::pen_imageExporter(std::function<double(unsigned long long, size_t)> fin) :
  nDimensions(1),
  nElements{0},
  elementSizes{0},
  origin{0},
  fbytes(std::bind(pen_imageExporter::conv2byte<double>,
		   std::placeholders::_1,
		   std::placeholders::_2,
		   std::placeholders::_3,
		   std::placeholders::_4,
		   std::placeholders::_5,
		   std::placeholders::_6, fin)),
  fASCII(std::bind(pen_imageExporter::conv2ASCII<double>,
		   std::placeholders::_1,
		   std::placeholders::_2,
		   std::placeholders::_3,
		   std::placeholders::_4,
		   std::placeholders::_5, fin)),    
  typeSize(sizeof(double)),
  integer(false),
  withSign(true),
  withErrors(false),
  baseName("none")
{}
  
pen_imageExporter::pen_imageExporter(std::function<std::int8_t(unsigned long long, size_t)> fin) :
  nDimensions(1),
  nElements{0},
  elementSizes{0},
  origin{0},
  fbytes(std::bind(pen_imageExporter::conv2byte<std::int8_t>,
		   std::placeholders::_1,
		   std::placeholders::_2,
		   std::placeholders::_3,
		   std::placeholders::_4,
		   std::placeholders::_5,
		   std::placeholders::_6, fin)),
  fASCII(std::bind(pen_imageExporter::conv2ASCII<std::int8_t>,
		   std::placeholders::_1,
		   std::placeholders::_2,
		   std::placeholders::_3,
		   std::placeholders::_4,
		   std::placeholders::_5, fin)),    
  typeSize(sizeof(std::int8_t)),
  integer(true),
  withSign(true),
  withErrors(false),
  baseName("none")
{}

pen_imageExporter::pen_imageExporter(std::function<std::uint8_t(unsigned long long, size_t)> fin) :
  nDimensions(1),
  nElements{0},
  elementSizes{0},
  origin{0},
  fbytes(std::bind(pen_imageExporter::conv2byte<std::uint8_t>,
		   std::placeholders::_1,
		   std::placeholders::_2,
		   std::placeholders::_3,
		   std::placeholders::_4,
		   std::placeholders::_5,
		   std::placeholders::_6, fin)),
  fASCII(std::bind(pen_imageExporter::conv2ASCII<std::uint8_t>,
		   std::placeholders::_1,
		   std::placeholders::_2,
		   std::placeholders::_3,
		   std::placeholders::_4,
		   std::placeholders::_5, fin)),
  typeSize(sizeof(std::uint8_t)),
  integer(true),
  withSign(false),
  withErrors(false),
  baseName("none")
{}
  
pen_imageExporter::pen_imageExporter(std::function<std::int16_t(unsigned long long, size_t)> fin) :
  nDimensions(1),
  nElements{0},
  elementSizes{0},
  origin{0},
  fbytes(std::bind(pen_imageExporter::conv2byte<std::int16_t>,
		   std::placeholders::_1,
		   std::placeholders::_2,
		   std::placeholders::_3,
		   std::placeholders::_4,
		   std::placeholders::_5,
		   std::placeholders::_6, fin)),
  fASCII(std::bind(pen_imageExporter::conv2ASCII<std::int16_t>,
		   std::placeholders::_1,
		   std::placeholders::_2,
		   std::placeholders::_3,
		   std::placeholders::_4,
		   std::placeholders::_5, fin)),
  typeSize(sizeof(std::int16_t)),
  integer(true),
  withSign(true),
  withErrors(false),
  baseName("none")
{}

pen_imageExporter::pen_imageExporter(std::function<std::uint16_t(unsigned long long, size_t)> fin) :
  nDimensions(1),
  nElements{0},
  elementSizes{0},
  origin{0},
  fbytes(std::bind(pen_imageExporter::conv2byte<std::uint16_t>,
		   std::placeholders::_1,
		   std::placeholders::_2,
		   std::placeholders::_3,
		   std::placeholders::_4,
		   std::placeholders::_5,
		   std::placeholders::_6, fin)),
  fASCII(std::bind(pen_imageExporter::conv2ASCII<std::uint16_t>,
		   std::placeholders::_1,
		   std::placeholders::_2,
		   std::placeholders::_3,
		   std::placeholders::_4,
		   std::placeholders::_5, fin)),
  typeSize(sizeof(std::uint16_t)),
  integer(true),
  withSign(false),
  withErrors(false),
  baseName("none")
{}

pen_imageExporter::pen_imageExporter(std::function<std::int32_t(unsigned long long, size_t)> fin) :
  nDimensions(1),
  nElements{0},
  elementSizes{0},
  origin{0},
  fbytes(std::bind(pen_imageExporter::conv2byte<std::int32_t>,
		   std::placeholders::_1,
		   std::placeholders::_2,
		   std::placeholders::_3,
		   std::placeholders::_4,
		   std::placeholders::_5,
		   std::placeholders::_6, fin)),
  fASCII(std::bind(pen_imageExporter::conv2ASCII<std::int32_t>,
		   std::placeholders::_1,
		   std::placeholders::_2,
		   std::placeholders::_3,
		   std::placeholders::_4,
		   std::placeholders::_5, fin)),
  typeSize(sizeof(std::int32_t)),
  integer(true),
  withSign(true),
  withErrors(false),
  baseName("none")
{}

pen_imageExporter::pen_imageExporter(std::function<std::uint32_t(unsigned long long, size_t)> fin) :
  nDimensions(1),
  nElements{0},
  elementSizes{0},
  origin{0},
  fbytes(std::bind(pen_imageExporter::conv2byte<std::uint32_t>,
		   std::placeholders::_1,
		   std::placeholders::_2,
		   std::placeholders::_3,
		   std::placeholders::_4,
		   std::placeholders::_5,
		   std::placeholders::_6, fin)),
  fASCII(std::bind(pen_imageExporter::conv2ASCII<std::uint32_t>,
		   std::placeholders::_1,
		   std::placeholders::_2,
		   std::placeholders::_3,
		   std::placeholders::_4,
		   std::placeholders::_5, fin)),
  typeSize(sizeof(std::uint32_t)),
  integer(true),
  withSign(false),
  withErrors(false),
  baseName("none")
{}

pen_imageExporter::pen_imageExporter(std::function<std::int64_t(unsigned long long, size_t)> fin) :
  nDimensions(1),
  nElements{0},
  elementSizes{0},
  origin{0},
  fbytes(std::bind(pen_imageExporter::conv2byte<std::int64_t>,
		   std::placeholders::_1,
		   std::placeholders::_2,
		   std::placeholders::_3,
		   std::placeholders::_4,
		   std::placeholders::_5,
		   std::placeholders::_6, fin)),
  fASCII(std::bind(pen_imageExporter::conv2ASCII<std::int64_t>,
		   std::placeholders::_1,
		   std::placeholders::_2,
		   std::placeholders::_3,
		   std::placeholders::_4,
		   std::placeholders::_5, fin)),
  typeSize(sizeof(std::int64_t)),
  integer(true),
  withSign(true),
  withErrors(false),
  baseName("none")
{}

pen_imageExporter::pen_imageExporter(std::function<std::uint64_t(unsigned long long, size_t)> fin) :
  nDimensions(1),
  nElements{0},
  elementSizes{0},
  origin{0},
  fbytes(std::bind(pen_imageExporter::conv2byte<std::uint64_t>,
		   std::placeholders::_1,
		   std::placeholders::_2,
		   std::placeholders::_3,
		   std::placeholders::_4,
		   std::placeholders::_5,
		   std::placeholders::_6, fin)),
  fASCII(std::bind(pen_imageExporter::conv2ASCII<std::uint64_t>,
		   std::placeholders::_1,
		   std::placeholders::_2,
		   std::placeholders::_3,
		   std::placeholders::_4,
		   std::placeholders::_5, fin)),
  typeSize(sizeof(std::uint64_t)),
  integer(true),
  withSign(true),
  withErrors(false),
  baseName("none")
{}

//------------------------------------------------------

// ** Constructors for images with asociated error

pen_imageExporter::pen_imageExporter(std::function<float(unsigned long long, size_t, float&)> fin) :
  nDimensions(1),
  nElements{0},
  elementSizes{0},
  origin{0},
  fbytes(std::bind(pen_imageExporter::conv2byteErr<float>,
		   std::placeholders::_1,
		   std::placeholders::_2,
		   std::placeholders::_3,
		   std::placeholders::_4,
		   std::placeholders::_5,
		   std::placeholders::_6, fin)),
  fASCII(std::bind(pen_imageExporter::conv2ASCIIErr<float>,
		   std::placeholders::_1,
		   std::placeholders::_2,
		   std::placeholders::_3,
		   std::placeholders::_4,
		   std::placeholders::_5, fin)),    
  typeSize(sizeof(float)),
  integer(false),
  withSign(true),
  withErrors(true),
  baseName("none")    
{}
  
pen_imageExporter::pen_imageExporter(std::function<double(unsigned long long, size_t, double&)> fin) :
  nDimensions(1),
  nElements{0},
  elementSizes{0},
  origin{0},
  fbytes(std::bind(pen_imageExporter::conv2byteErr<double>,
		   std::placeholders::_1,
		   std::placeholders::_2,
		   std::placeholders::_3,
		   std::placeholders::_4,
		   std::placeholders::_5,
		   std::placeholders::_6, fin)),
  fASCII(std::bind(pen_imageExporter::conv2ASCIIErr<double>,
		   std::placeholders::_1,
		   std::placeholders::_2,
		   std::placeholders::_3,
		   std::placeholders::_4,
		   std::placeholders::_5, fin)),    
  typeSize(sizeof(double)),
  integer(false),
  withSign(true),
  withErrors(true),
  baseName("none")
{}
  
pen_imageExporter::pen_imageExporter(std::function<std::int8_t(unsigned long long, size_t, std::int8_t&)> fin) :
  nDimensions(1),
  nElements{0},
  elementSizes{0},
  origin{0},
  fbytes(std::bind(pen_imageExporter::conv2byteErr<std::int8_t>,
		   std::placeholders::_1,
		   std::placeholders::_2,
		   std::placeholders::_3,
		   std::placeholders::_4,
		   std::placeholders::_5,
		   std::placeholders::_6, fin)),
  fASCII(std::bind(pen_imageExporter::conv2ASCIIErr<std::int8_t>,
		   std::placeholders::_1,
		   std::placeholders::_2,
		   std::placeholders::_3,
		   std::placeholders::_4,
		   std::placeholders::_5, fin)),    
  typeSize(sizeof(std::int8_t)),
  integer(true),
  withSign(true),
  withErrors(true),
  baseName("none")
{}

pen_imageExporter::pen_imageExporter(std::function<std::uint8_t(unsigned long long, size_t, std::uint8_t&)> fin) :
  nDimensions(1),
  nElements{0},
  elementSizes{0},
  origin{0},
  fbytes(std::bind(pen_imageExporter::conv2byteErr<std::uint8_t>,
		   std::placeholders::_1,
		   std::placeholders::_2,
		   std::placeholders::_3,
		   std::placeholders::_4,
		   std::placeholders::_5,
		   std::placeholders::_6, fin)),
  fASCII(std::bind(pen_imageExporter::conv2ASCIIErr<std::uint8_t>,
		   std::placeholders::_1,
		   std::placeholders::_2,
		   std::placeholders::_3,
		   std::placeholders::_4,
		   std::placeholders::_5, fin)),
  typeSize(sizeof(std::uint8_t)),
  integer(true),
  withSign(false),
  withErrors(true),
  baseName("none")
{}
  
pen_imageExporter::pen_imageExporter(std::function<std::int16_t(unsigned long long, size_t, std::int16_t&)> fin) :
  nDimensions(1),
  nElements{0},
  elementSizes{0},
  origin{0},
  fbytes(std::bind(pen_imageExporter::conv2byteErr<std::int16_t>,
		   std::placeholders::_1,
		   std::placeholders::_2,
		   std::placeholders::_3,
		   std::placeholders::_4,
		   std::placeholders::_5,
		   std::placeholders::_6, fin)),
  fASCII(std::bind(pen_imageExporter::conv2ASCIIErr<std::int16_t>,
		   std::placeholders::_1,
		   std::placeholders::_2,
		   std::placeholders::_3,
		   std::placeholders::_4,
		   std::placeholders::_5, fin)),
  typeSize(sizeof(std::int16_t)),
  integer(true),
  withSign(true),
  withErrors(true),
  baseName("none")
{}

pen_imageExporter::pen_imageExporter(std::function<std::uint16_t(unsigned long long, size_t, std::uint16_t&)> fin) :
  nDimensions(1),
  nElements{0},
  elementSizes{0},
  origin{0},
  fbytes(std::bind(pen_imageExporter::conv2byteErr<std::uint16_t>,
		   std::placeholders::_1,
		   std::placeholders::_2,
		   std::placeholders::_3,
		   std::placeholders::_4,
		   std::placeholders::_5,
		   std::placeholders::_6, fin)),
  fASCII(std::bind(pen_imageExporter::conv2ASCIIErr<std::uint16_t>,
		   std::placeholders::_1,
		   std::placeholders::_2,
		   std::placeholders::_3,
		   std::placeholders::_4,
		   std::placeholders::_5, fin)),
  typeSize(sizeof(std::uint16_t)),
  integer(true),
  withSign(false),
  withErrors(true),
  baseName("none")
{}

pen_imageExporter::pen_imageExporter(std::function<std::int32_t(unsigned long long, size_t, std::int32_t&)> fin) :
  nDimensions(1),
  nElements{0},
  elementSizes{0},
  origin{0},
  fbytes(std::bind(pen_imageExporter::conv2byteErr<std::int32_t>,
		   std::placeholders::_1,
		   std::placeholders::_2,
		   std::placeholders::_3,
		   std::placeholders::_4,
		   std::placeholders::_5,
		   std::placeholders::_6, fin)),
  fASCII(std::bind(pen_imageExporter::conv2ASCIIErr<std::int32_t>,
		   std::placeholders::_1,
		   std::placeholders::_2,
		   std::placeholders::_3,
		   std::placeholders::_4,
		   std::placeholders::_5, fin)),
  typeSize(sizeof(std::int32_t)),
  integer(true),
  withSign(true),
  withErrors(true),
  baseName("none")
{}

pen_imageExporter::pen_imageExporter(std::function<std::uint32_t(unsigned long long, size_t, std::uint32_t&)> fin) :
  nDimensions(1),
  nElements{0},
  elementSizes{0},
  origin{0},
  fbytes(std::bind(pen_imageExporter::conv2byteErr<std::uint32_t>,
		   std::placeholders::_1,
		   std::placeholders::_2,
		   std::placeholders::_3,
		   std::placeholders::_4,
		   std::placeholders::_5,
		   std::placeholders::_6, fin)),
  fASCII(std::bind(pen_imageExporter::conv2ASCIIErr<std::uint32_t>,
		   std::placeholders::_1,
		   std::placeholders::_2,
		   std::placeholders::_3,
		   std::placeholders::_4,
		   std::placeholders::_5, fin)),
  typeSize(sizeof(std::uint32_t)),
  integer(true),
  withSign(false),
  withErrors(true),
  baseName("none")
{}

pen_imageExporter::pen_imageExporter(std::function<std::int64_t(unsigned long long, size_t, std::int64_t&)> fin) :
  nDimensions(1),
  nElements{0},
  elementSizes{0},
  origin{0},
  fbytes(std::bind(pen_imageExporter::conv2byteErr<std::int64_t>,
		   std::placeholders::_1,
		   std::placeholders::_2,
		   std::placeholders::_3,
		   std::placeholders::_4,
		   std::placeholders::_5,
		   std::placeholders::_6, fin)),
  fASCII(std::bind(pen_imageExporter::conv2ASCIIErr<std::int64_t>,
		   std::placeholders::_1,
		   std::placeholders::_2,
		   std::placeholders::_3,
		   std::placeholders::_4,
		   std::placeholders::_5, fin)),
  typeSize(sizeof(std::int64_t)),
  integer(true),
  withSign(true),
  withErrors(true),
  baseName("none")
{}

pen_imageExporter::pen_imageExporter(std::function<std::uint64_t(unsigned long long, size_t, std::uint64_t&)> fin) :
  nDimensions(1),
  nElements{0},
  elementSizes{0},
  origin{0},
  fbytes(std::bind(pen_imageExporter::conv2byteErr<std::uint64_t>,
		   std::placeholders::_1,
		   std::placeholders::_2,
		   std::placeholders::_3,
		   std::placeholders::_4,
		   std::placeholders::_5,
		   std::placeholders::_6, fin)),
  fASCII(std::bind(pen_imageExporter::conv2ASCIIErr<std::uint64_t>,
		   std::placeholders::_1,
		   std::placeholders::_2,
		   std::placeholders::_3,
		   std::placeholders::_4,
		   std::placeholders::_5, fin)),
  typeSize(sizeof(std::uint64_t)),
  integer(true),
  withSign(true),
  withErrors(true),
  baseName("none")
{}

void pen_imageExporter::writeMHDheader(FILE* fout,
				       const std::string& rawFilename) const{

  fprintf(fout,"ObjectType = Image\n");
  fprintf(fout,"NDims = %u\n",nDimensions);
  fprintf(fout,"BinaryData = True\n");
  fprintf(fout,"BinaryDataByteOrderMSB = False\n");
  fprintf(fout,"CompressedData = False\n");
  fprintf(fout,"TransformMatrix = 1 0 0 0 1 0 0 0 1\n");

  fprintf(fout,"Offset =");
  for(unsigned i = 0; i < nDimensions; ++i){
    fprintf(fout," %12.6f",(origin[i] + elementSizes[i]/2.0)*10.0);
  }
  fprintf(fout,"\n");
    
  fprintf(fout,"CenterOfRotation = 0 0 0\n");
  fprintf(fout,"AnatomicalOrientation = RAI\n");

  fprintf(fout,"ElementSpacing = ");
  for(unsigned i = 0; i < nDimensions; ++i){
    fprintf(fout," %10.6f",elementSizes[i]*10.0);
  }
  fprintf(fout,"\n");
      
  fprintf(fout,"DimSize = ");
  for(unsigned i = 0; i < nDimensions; ++i){
    fprintf(fout," %u",nElements[i]);
  }
  fprintf(fout,"\n");

  if(!integer){
    //Non integer types
    if(typeSize == sizeof(float))
      fprintf(fout,"ElementType = MET_FLOAT\n");
    else
      fprintf(fout,"ElementType = MET_DOUBLE\n");
  }else{
    //Integer types

    if(withSign){
      //Signed integers
      if(typeSize == sizeof(int8_t))
	fprintf(fout,"ElementType = MET_CHAR\n");
      else if(typeSize == sizeof(int16_t))
	fprintf(fout,"ElementType = MET_SHORT\n");
      else if(typeSize == sizeof(int32_t))
	fprintf(fout,"ElementType = MET_INT\n");
      else if(typeSize == sizeof(int64_t))
	fprintf(fout,"ElementType = MET_LONG\n");
    }else{
      //Unsigned integers
      if(typeSize == sizeof(uint8_t))
	fprintf(fout,"ElementType = MET_UCHAR\n");
      else if(typeSize == sizeof(uint16_t))
	fprintf(fout,"ElementType = MET_USHORT\n");
      else if(typeSize == sizeof(uint32_t))
	fprintf(fout,"ElementType = MET_UINT\n");
      else if(typeSize == sizeof(uint64_t))
	fprintf(fout,"ElementType = MET_ULONG\n");	  
    }
  }

  std::string filename;
  auto const pos = rawFilename.find_last_of('/');
  if(pos != std::string::npos){
    filename = rawFilename.substr(pos+1);
  }else{
    auto const pos2 = rawFilename.find_last_of('\\');      
    if(pos2 != std::string::npos){
      filename = rawFilename.substr(pos2+1);
    }else{
      filename = rawFilename;
    }
  }

  fprintf(fout,"ElementDataFile = %-s\n",filename.c_str());
  
}

void pen_imageExporter::writeGNUheader(FILE* fout) const{

  fprintf(fout,"# NDims = %u\n",nDimensions);
  fprintf(fout,"# Offset =");
  for(unsigned i = 0; i < nDimensions; ++i){
    fprintf(fout," %.5E",origin[i]);
  }
  fprintf(fout," cm\n");
    
  fprintf(fout,"# Element Spacing = ");
  for(unsigned i = 0; i < nDimensions; ++i){
    fprintf(fout," %10.6f",elementSizes[i]);
  }
  fprintf(fout," cm\n");
      
  fprintf(fout,"# Number of elements = ");
  for(unsigned i = 0; i < nDimensions; ++i){
    fprintf(fout," %u",nElements[i]);
  }
  fprintf(fout,"\n\n");  
}


void pen_imageExporter::exportImage(const unsigned long long nhists,
				    const formatTypes format) const{

  size_t totalElements = 1;
  for(unsigned i = 0; i < nDimensions; ++i){
    totalElements *= nElements[i];
  }
    
  switch(format){
  case formatTypes::MHD:{

    //Create file names
    std::string headerFilename(baseName);
    headerFilename.append(".mhd");

    std::string rawFilename(baseName);
    rawFilename.append(".raw");

    std::string headerErrFilename(baseName);
    headerErrFilename.append("-Uncertainty.mhd");

    std::string rawErrFilename(baseName);
    rawErrFilename.append("-Uncertainty.raw");

    //Create the mhd header file
    FILE* fout = nullptr;
    fout = fopen(headerFilename.c_str(), "w");
    if(fout ==  nullptr)
      return;

    //Write and close the header
    writeMHDheader(fout, rawFilename);
    fclose(fout);

    //Check if the values have an asociated error
    if(withErrors){
      //Do the same procedure for errors mhd header

      fout = nullptr;
      fout = fopen(headerErrFilename.c_str(), "w");
      if(fout ==  nullptr)
	return;

      //Write and close the header
      writeMHDheader(fout, rawErrFilename);
      fclose(fout);
      
    }
    
    FILE* fraw = nullptr;
    fraw = fopen(rawFilename.c_str(), "wb");
    if(fraw ==  nullptr)
      return;

    FILE* frawErr = nullptr;
    if(withErrors){
      frawErr = fopen(rawErrFilename.c_str(), "wb");
      if(frawErr ==  nullptr)
	return;
    }
    
    //Create a buffer to read a chunk of elements
    const size_t readBufferSize = sizeof(double)*100;
    unsigned char readBuffer [readBufferSize];
    unsigned char readBufferErr [readBufferSize];
    size_t elementsPerRead = readBufferSize/typeSize;
    size_t next2read = 0;
    while(next2read < totalElements){
      //Obtain the number of elements to read
      const size_t toRead = std::min(totalElements-next2read, elementsPerRead);

      //Read the elements from the specified function
      size_t offset = 0;
      fbytes(nhists, readBuffer, readBufferErr, offset, next2read, toRead);
      next2read += toRead;
	
      //Write elements to raw file
      fwrite(static_cast<void*>(readBuffer), typeSize, toRead, fraw);
      if(withErrors)
	fwrite(static_cast<void*>(readBufferErr), typeSize, toRead, frawErr);
	
    }

    fclose(fraw);
    if(withErrors)
      fclose(frawErr);
      
    break;
  }
  case formatTypes::GNU:{

    const size_t maxCharsPerValue = 22;
    const size_t maxValuesInBuffer = 1000;
    const size_t readBufferSize = maxCharsPerValue*maxValuesInBuffer;
    char readBuffer[readBufferSize];
    char readBufferErr[readBufferSize];

    FILE* fout = nullptr;
    std::string filename(baseName);
    filename.append(".gnu");
    fout = fopen(filename.c_str(), "w");
    if(fout == nullptr)
      return;

    writeGNUheader(fout);

    //Check if associated errors must be printed
    FILE* foutErr = nullptr;
    if(withErrors){
      filename.assign(baseName);
      filename.append("-Uncertainty.gnu");
      foutErr = fopen(filename.c_str(), "w");
      if(foutErr == nullptr)
	return;      

      writeGNUheader(foutErr);
    }
    
    size_t printed = 0;
    size_t dimCounters[maxDims] = {0};
    //Check if a row of the dimension 0 fits in the buffer

    if(nElements[0] <= maxValuesInBuffer){

      while(printed < totalElements){
	//Read a 0 dimension row
	fASCII(nhists, readBuffer, readBufferErr, printed, nElements[0]);
	printed += nElements[0];
	dimCounters[1] += 1;

	//Write data
	fprintf(fout,"%s\n",readBuffer);

	if(withErrors)
	  fprintf(foutErr,"%s\n",readBufferErr);
	  
	
	//Check if some dimension has been ended the row
	for(size_t idim = 1; idim < nDimensions-1; ++idim){
	  if(dimCounters[idim] >= nElements[idim]){
	    //This dimension ends its row, go to next
	    //dimension and print two blank lines
	    dimCounters[idim] = 0;
	    ++dimCounters[idim+1];
	    fprintf(fout,"\n\n");
	    if(withErrors)
	      fprintf(foutErr,"\n\n");
	  }else{
	    //This dimension row remains incomplete, stop
	    break;
	  }
	}
      }
    }
    else{

      while(printed < totalElements){
	//Read a 0 dimension row
	size_t rowPos = 0;
	while(rowPos < nElements[0]){

	  //Calculate elements to read
	  size_t toRead = std::min(maxValuesInBuffer,nElements[0]-rowPos);
	  //Read values
	  fASCII(nhists, readBuffer, readBufferErr, printed, toRead);
	  //Write data
	  fprintf(fout,"%s",readBuffer);
	  if(withErrors)
	    fprintf(foutErr,"%s",readBufferErr);
	  
	  //Update counters
	  printed += toRead;
	  rowPos += toRead;
	}
	
	//End of row, print a end of line and increase row counter
	fprintf(fout,"\n");
	if(withErrors)
	  fprintf(foutErr,"\n");
	dimCounters[1] += 1;


	//Check if some dimension has been ended the row
	for(size_t idim = 1; idim < nDimensions-1; ++idim){
	  if(dimCounters[idim] >= nElements[idim]){
	    //This dimension ends its row, go to next
	    //dimension and print a blank line
	    dimCounters[idim] = 0;
	    ++dimCounters[idim+1];
	    fprintf(fout,"\n\n");
	    if(withErrors)
	      fprintf(foutErr,"\n\n");
	  }else{
	    //This dimension row remains incomplete, stop
	    break;
	  }
	}
      }
      
    }
    
    fclose(fout);
    if(withErrors)
      fclose(foutErr);
    
    break;
  }
  case formatTypes::NONE:{
    break;
  }
  default:{
    break;
  }
  }
    
}
