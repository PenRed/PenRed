 

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


#ifndef __PEN_IMAGE_DUMP__
#define __PEN_IMAGE_DUMP__

#include <cstdlib>
#include <functional>
#include <cstdint>
#include <string>
#include <cstring>
#include <array>

struct pen_imageExporter{

  enum formatTypes{
    MHD,
    GNU,
    NONE
  };

  static const constexpr std::array<const char*, formatTypes::NONE> formatNames{
    "MHD",
    "GNU"
  };  
  
  static inline formatTypes toFormat(const char* name){
    if(strcmp(name,"MHD") == 0)
      return formatTypes::MHD;
    if(strcmp(name,"GNU") == 0)
      return formatTypes::GNU;
    
    return formatTypes::NONE;
  }

  static inline const char* toString(const formatTypes format){
    switch(format){
    case formatTypes::MHD:
      return "MHD";
    case formatTypes::GNU:
      return "GNU";
    default:
      return nullptr;
    };
    return nullptr;
  }
  
  static inline bool isFormat(const char* name){
    if(strcmp(name,"MHD") == 0)
      return true;
    if(strcmp(name,"GNU") == 0)
      return true;

    return false;    
  }

  static constexpr size_t nFormats(){return formatTypes::NONE;}
  
  static unsigned const maxDims = 10;
  
private:

  // With no error asociated
  template<class T>
  static void conv2byte(unsigned long long nhists,
			unsigned char* buffer, unsigned char* /*errBuffer*/,
			size_t& offset,
			size_t init, size_t n,
			std::function<T(unsigned long long, size_t)> f){

    size_t end = init + n;
    
    for(size_t i = init; i < end; ++i){
      T value = f(nhists,i);
      memcpy(buffer+offset, static_cast<const void*>(&value), sizeof(T));
      offset += sizeof(T);
    }
    
  }

  // With error asociated
  template<class T>
  static void conv2byteErr(unsigned long long nhists,
			   unsigned char* buffer, unsigned char* errBuffer,
			   size_t& offset,
			   size_t init, size_t n,
			   std::function<T(unsigned long long, size_t, T&)> f){

    size_t end = init + n;
    
    for(size_t i = init; i < end; ++i){
      T err;
      T value = f(nhists,i,err);
      memcpy(buffer+offset   , static_cast<const void*>(&value), sizeof(T));
      memcpy(errBuffer+offset, static_cast<const void*>(&err)  , sizeof(T));
      offset += sizeof(T);
    }
    
  }  

  // With no error asociated
  template<class T>
  static void conv2ASCII(unsigned long long nhists,
			 char* buffer, char* /*errBuffer*/,
			 size_t init, size_t n,
			 std::function<T(unsigned long long, size_t)> f){

    size_t end = init + n;
    int pos = 0;

    if(std::is_floating_point<T>::value){

      for(size_t i = init; i < end; ++i){
	T value = f(nhists,i);
	int printed;
	sprintf(buffer+pos," %15.5E%n",
		static_cast<double>(value),&printed);
	pos += printed;
      }
	
    }else{
      if(std::is_signed<T>::value){
	  
	for(size_t i = init; i < end; ++i){
	  T value = f(nhists,i);	
	  int printed;
	  sprintf(buffer+pos," %lld%n",
		  static_cast<long long int>(value),&printed);
	  pos += printed;
	}
      }
      else

	for(size_t i = init; i < end; ++i){
	  T value = f(nhists,i);	
	  int printed;
	  sprintf(buffer+pos," %llu%n",
		  static_cast<long long unsigned>(value),&printed);
	  pos += printed;
	}
    }
    
  }

  // With error asociated
  template<class T>
  static void conv2ASCIIErr(unsigned long long nhists,
			    char* buffer, char* errBuffer,
			    size_t init, size_t n,
			    std::function<T(unsigned long long, size_t, T&)> f){

    size_t end = init + n;
    int pos = 0;
    int errPos = 0;

    if(std::is_floating_point<T>::value){

      for(size_t i = init; i < end; ++i){
	T err;
	T value = f(nhists,i,err);
	int printed;
	sprintf(buffer+pos," %15.5E%n",
		static_cast<double>(value),&printed);
	pos += printed;

	sprintf(errBuffer+errPos," %15.5E%n",
		static_cast<double>(err),&printed);
	errPos += printed;
      }
	
    }else{
      if(std::is_signed<T>::value){
	  
	for(size_t i = init; i < end; ++i){
	  T err;
	  T value = f(nhists,i,err);	
	  int printed;
	  sprintf(buffer+pos," %lld%n",
		  static_cast<long long int>(value),&printed);
	  pos += printed;

	  sprintf(errBuffer+errPos," %lld%n",
		  static_cast<long long int>(err),&printed);
	  errPos += printed;
	  
	}
      }
      else

	for(size_t i = init; i < end; ++i){
	  T err;
	  T value = f(nhists,i,err);	
	  int printed;
	  sprintf(buffer+pos," %llu%n",
		  static_cast<long long unsigned>(value),&printed);
	  pos += printed;

	  sprintf(errBuffer+errPos," %llu%n",
		  static_cast<long long unsigned>(err),&printed);
	  errPos += printed;
	  
	}
    }
    
  }
  
  unsigned nDimensions;
  unsigned nElements[maxDims];
  float elementSizes[maxDims];
  float origin[maxDims];

  std::function<void(unsigned long long,
		     unsigned char* , unsigned char*,
		     size_t&, size_t, size_t)> fbytes;

  std::function<void(unsigned long long,
		     char* , char* , size_t, size_t)> fASCII;
  
  size_t typeSize;
  bool integer;
  bool withSign;
  bool withErrors;
public:

  std::string baseName;

  // ** Constructors for images without associated uncertainty
  pen_imageExporter(std::function<float(unsigned long long, size_t)> fin); 
  pen_imageExporter(std::function<double(unsigned long long, size_t)> fin);
  pen_imageExporter(std::function<std::int8_t(unsigned long long, size_t)> fin);
  pen_imageExporter(std::function<std::uint8_t(unsigned long long, size_t)> fin);
  pen_imageExporter(std::function<std::int16_t(unsigned long long, size_t)> fin);
  pen_imageExporter(std::function<std::uint16_t(unsigned long long, size_t)> fin);
  pen_imageExporter(std::function<std::int32_t(unsigned long long, size_t)> fin);
  pen_imageExporter(std::function<std::uint32_t(unsigned long long, size_t)> fin);
  pen_imageExporter(std::function<std::int64_t(unsigned long long, size_t)> fin);
  pen_imageExporter(std::function<std::uint64_t(unsigned long long, size_t)> fin);

  // ** Constructors for images associated uncertainty
  pen_imageExporter(std::function<float(unsigned long long, size_t, float&)> fin); 
  pen_imageExporter(std::function<double(unsigned long long, size_t, double&)> fin);
  pen_imageExporter(std::function<std::int8_t(unsigned long long, size_t, std::int8_t&)> fin);
  pen_imageExporter(std::function<std::uint8_t(unsigned long long, size_t, std::uint8_t&)> fin);
  pen_imageExporter(std::function<std::int16_t(unsigned long long, size_t, std::int16_t&)> fin);
  pen_imageExporter(std::function<std::uint16_t(unsigned long long, size_t, std::uint16_t&)> fin);
  pen_imageExporter(std::function<std::int32_t(unsigned long long, size_t, std::int32_t&)> fin);
  pen_imageExporter(std::function<std::uint32_t(unsigned long long, size_t, std::uint32_t&)> fin);
  pen_imageExporter(std::function<std::int64_t(unsigned long long, size_t, std::int64_t&)> fin);
  pen_imageExporter(std::function<std::uint64_t(unsigned long long, size_t, std::uint64_t&)> fin);
  
  void setDimensions(const unsigned nDim,
		     const unsigned* elements,
		     const float* delements);

  void setOrigin(const double* newOrigin);
  
  void exportImage(const unsigned long long nhists,
		   const formatTypes format) const;

  unsigned getDim() const{return nDimensions;}

  void writeMHDheader(FILE*, const std::string&) const;
  void writeGNUheader(FILE*) const;
};

inline void pen_imageExporter::setDimensions(const unsigned nDim,
					     const unsigned* elements,
					     const float* delements){

  if(nDim > 0 && nDim <= maxDims){
    nDimensions = nDim;
    for(unsigned i = 0; i < nDimensions; ++i){
      nElements[i] = elements[i];
      elementSizes[i] = delements[i];
    }
  }
}

inline void pen_imageExporter::setOrigin(const double* newOrigin){

  for(unsigned i = 0; i < nDimensions; ++i){
    origin[i] = newOrigin[i];
  }
}

#endif
