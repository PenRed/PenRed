
//
//
//    Copyright (C) 2019-2023 Universitat de València - UV
//    Copyright (C) 2019-2023 Universitat Politècnica de València - UPV
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


#ifndef __PEN_BINARY_DUMP__
#define __PEN_BINARY_DUMP__

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <vector>
#include <cmath>
#include <algorithm>
#include <climits>

template<class T>
inline size_t sizeofBits(){
  return sizeof(T)*CHAR_BIT;
}

typedef int16_t dump_Short;
typedef int16_t   dump_Int;
typedef int32_t  dump_Long;
typedef int64_t dump_Llong;

typedef uint16_t dump_uShort;
typedef uint16_t   dump_uInt;
typedef uint32_t  dump_uLong;
typedef uint64_t dump_uLlong;

//Constants with minimum size in bits of each integer type
enum class intTypes {SHORT, INT, LONG, LLONG, UNKNOWN};
enum class minIntSizes_b {SHORTb = 16,
			  INTb = 16,
			  LONGb = 32,
			  LLONGb = 64,
			  INVALIDb = 0};
enum class minIntSizes_B {SHORTB = sizeof(dump_Short),
			  INTB = sizeof(dump_Int),
			  LONGB = sizeof(dump_Long),
			  LLONGB = sizeof(dump_Llong),
			  INVALIDB = 0,};


enum dumpState{
  PEN_DUMP_SUCCESS = 0,
  PEN_DUMP_NULL_POINTER,
  PEN_DUMP_ZERO_ELEMENTS,
  PEN_DUMP_ELEMENT_NOT_FOUND,
  PEN_DUMP_BAD_ALLOCATION,
  PEN_DUMP_INSUFFICIENT_MEMORY,
  PEN_DUMP_INT_OUT_OF_RANGE,
  PEN_DUMP_UNSIGNED_OUT_OF_RANGE,
  PEN_DUMP_ERROR_DOUBLE_DUMP,  
  PEN_DUMP_ERROR_INT_DUMP,  
  PEN_DUMP_ERROR_UNSIGNED_DUMP,  
  PEN_DUMP_ERROR_CHAR_DUMP,
  PEN_DUMP_INCORRECT_DATA_SIZE,
  PEN_DUMP_NARRAY_NOT_MATCH,
  PEN_DUMP_ELEMENT_SIZE_NOT_MATCH,
  PEN_DUMP_ELEMENT_NUMBER_NOT_MATCH,
  PEN_DUMP_UNABLE_TO_READ_DOUBLE_ARRAYS,
  PEN_DUMP_UNABLE_TO_READ_INT_ARRAYS,
  PEN_DUMP_UNABLE_TO_READ_UNSIGNED_ARRAYS,
  PEN_DUMP_UNABLE_TO_READ_CHAR_ARRAYS,
  PEN_DUMP_INTEGER_LARGER_THAN_MAXIMUM,
  PEN_DUMP_INVALID_TYPE
};


template<class T>
inline minIntSizes_b minContainer(size_t& nbytes){
  
  const size_t nbits = sizeofBits<T>();

  if(nbits <= static_cast<size_t>(minIntSizes_b::SHORTb)){
    nbytes = static_cast<size_t>(minIntSizes_B::SHORTB);
    return minIntSizes_b::SHORTb;
  }
  if(nbits <= static_cast<size_t>(minIntSizes_b::INTb)){
    nbytes = static_cast<size_t>(minIntSizes_B::INTB);
    return minIntSizes_b::INTb;
  }  
  if(nbits <= static_cast<size_t>(minIntSizes_b::LONGb)){
    nbytes = static_cast<size_t>(minIntSizes_B::LONGB);
    return minIntSizes_b::LONGb;
  }
  if(nbits <= static_cast<size_t>(minIntSizes_b::LLONGb)){
    nbytes = static_cast<size_t>(minIntSizes_B::LLONGB);
    return minIntSizes_b::LLONGb;
  }

  nbytes = static_cast<size_t>(minIntSizes_B::INVALIDB);
  return minIntSizes_b::INVALIDb;
  
}

template<class T>
inline int dumpSInt(unsigned char* pout,
		    size_t& pos,
		    const T* p,
		    const size_t nElements,
		    const minIntSizes_b containerBits){
 
  if(containerBits == minIntSizes_b::LONGb){
    //This data requires a container of "LONGb" bits 
    const size_t byteSize = static_cast<size_t>(minIntSizes_B::LONGB);
    for(size_t i = 0; i < nElements; ++i){
      dump_Long value = dump_Long(p[i]);
      memcpy(&pout[pos],&value,byteSize);
      pos += byteSize;
    }
    return PEN_DUMP_SUCCESS;
  }
  
  if(containerBits == minIntSizes_b::LLONGb){
    //This data requires a container of "LLONGb" bits 
    const size_t byteSize = static_cast<size_t>(minIntSizes_B::LLONGB);
    for(size_t i = 0; i < nElements; ++i){
      dump_Llong value = dump_Llong(p[i]);
      memcpy(&pout[pos],&value,byteSize);
      pos += byteSize;
    }
    return PEN_DUMP_SUCCESS;
  }
  
  if(containerBits == minIntSizes_b::SHORTb){
    //This data requires a container of "SHORTb" bits 
    const size_t byteSize = static_cast<size_t>(minIntSizes_B::SHORTB);
    for(size_t i = 0; i < nElements; ++i){
      dump_Short value = dump_Short(p[i]);
      memcpy(&pout[pos],&value,byteSize);
      pos += byteSize;
    }
    return PEN_DUMP_SUCCESS;
  }  

  if(containerBits == minIntSizes_b::INTb){
    //This data requires a container of "INTb" bits
    const size_t byteSize = static_cast<size_t>(minIntSizes_B::INTB);
    for(size_t i = 0; i < nElements; ++i){
      dump_Int value = dump_Int(p[i]);
      memcpy(&pout[pos],&value,byteSize);
      pos += byteSize;
    }
    return PEN_DUMP_SUCCESS;
  }
  
  //Invalid container size
  return PEN_DUMP_INVALID_TYPE;
}

template<class T>
inline int dumpUInt(unsigned char* pout,
		    size_t& pos,
		    const T* p,
		    const size_t nElements,
		    const minIntSizes_b containerBits){

  if(containerBits == minIntSizes_b::LONGb){
    //This data requires a container of "LONGb" bits 
    const size_t byteSize = static_cast<size_t>(minIntSizes_B::LONGB);
    for(size_t i = 0; i < nElements; ++i){
      dump_uLong value = dump_uLong(p[i]);
      memcpy(&pout[pos],&value,byteSize);
      pos += byteSize;
    }
    return PEN_DUMP_SUCCESS;
  }
  
  if(containerBits == minIntSizes_b::LLONGb){
    //This data requires a container of "LLONGb" bits 
    const size_t byteSize = static_cast<size_t>(minIntSizes_B::LLONGB);
    for(size_t i = 0; i < nElements; ++i){
      dump_uLlong value = dump_uLlong(p[i]);
      memcpy(&pout[pos],&value,byteSize);
      pos += byteSize;
    }
    return PEN_DUMP_SUCCESS;
  }
  
  if(containerBits == minIntSizes_b::SHORTb){
    //This data requires a container of "SHORTb" bits 
    const size_t byteSize = static_cast<size_t>(minIntSizes_B::SHORTB);
    for(size_t i = 0; i < nElements; ++i){
      dump_uShort value = dump_uShort(p[i]);
      memcpy(&pout[pos],&value,byteSize);
      pos += byteSize;
    }
    return PEN_DUMP_SUCCESS;
  }  

  if(containerBits == minIntSizes_b::INTb){
    //This data requires a container of "INTb" bits
    const size_t byteSize = static_cast<size_t>(minIntSizes_B::INTB);
    for(size_t i = 0; i < nElements; ++i){
      dump_uInt value = dump_uInt(p[i]);
      memcpy(&pout[pos],&value,byteSize);
      pos += byteSize;
    }
    return PEN_DUMP_SUCCESS;
  }
  
  //Invalid container size
  return PEN_DUMP_INVALID_TYPE;
}

template<class T>
inline int readSInt(const unsigned char* pin,
		    size_t& pos,
		    T* p,
		    const size_t nElements,
		    const size_t containerBits){

  if(containerBits == static_cast<size_t>(minIntSizes_b::LONGb)){
    //This data requires a container of "LONGb" bits 
    const size_t byteSize = static_cast<size_t>(minIntSizes_B::LONGB);
    for(size_t i = 0; i < nElements; ++i){
      dump_Long value;
      memcpy(&value,&pin[pos],byteSize);
      pos += byteSize;
      //Covert value
      p[i] = static_cast<T>(value);
    }
    return PEN_DUMP_SUCCESS;
  }
  
  if(containerBits == static_cast<size_t>(minIntSizes_b::LLONGb)){
    //This data requires a container of "LLONGb" bits 
    const size_t byteSize = static_cast<size_t>(minIntSizes_B::LLONGB);
    for(size_t i = 0; i < nElements; ++i){
      dump_Llong value;
      memcpy(&value,&pin[pos],byteSize);
      pos += byteSize;
      //Covert value
      p[i] = static_cast<T>(value);
    }
    return PEN_DUMP_SUCCESS;
  }
  
  if(containerBits == static_cast<size_t>(minIntSizes_b::SHORTb)){
    //This data requires a container of "SHORTb" bits 
    const size_t byteSize = static_cast<size_t>(minIntSizes_B::SHORTB);
    for(size_t i = 0; i < nElements; ++i){
      dump_Short value;
      memcpy(&value,&pin[pos],byteSize);
      pos += byteSize;
      //Covert value
      p[i] = static_cast<T>(value);
    }
    return PEN_DUMP_SUCCESS;
  }  

  if(containerBits == static_cast<size_t>(minIntSizes_b::INTb)){
    //This data requires a container of "INTb" bits
    const size_t byteSize = static_cast<size_t>(minIntSizes_B::INTB);
    for(size_t i = 0; i < nElements; ++i){
      dump_Int value;
      memcpy(&value,&pin[pos],byteSize);
      pos += byteSize;
      //Covert value
      p[i] = static_cast<T>(value);
    }
    return PEN_DUMP_SUCCESS;
  }
  
  //Invalid container size
  return PEN_DUMP_INVALID_TYPE;
}

template<class T>
inline int readUInt(const unsigned char* pin,
		    size_t& pos,
		    T* p,
		    const size_t nElements,
		    const size_t containerBits){

  if(containerBits == static_cast<size_t>(minIntSizes_b::LONGb)){
    //This data requires a container of "LONGb" bits 
    const size_t byteSize = static_cast<size_t>(minIntSizes_B::LONGB);
    for(size_t i = 0; i < nElements; ++i){
      dump_uLong value;
      memcpy(&value,&pin[pos],byteSize);
      pos += byteSize;
      //Covert value
      p[i] = static_cast<T>(value);
    }
    return PEN_DUMP_SUCCESS;
  }
  
  if(containerBits == static_cast<size_t>(minIntSizes_b::LLONGb)){
    //This data requires a container of "LLONGb" bits 
    const size_t byteSize = static_cast<size_t>(minIntSizes_B::LLONGB);
    for(size_t i = 0; i < nElements; ++i){
      dump_uLlong value;
      memcpy(&value,&pin[pos],byteSize);
      pos += byteSize;
      //Covert value
      p[i] = static_cast<T>(value);
    }
    return PEN_DUMP_SUCCESS;
  }
  
  if(containerBits == static_cast<size_t>(minIntSizes_b::SHORTb)){
    //This data requires a container of "SHORTb" bits 
    const size_t byteSize = static_cast<size_t>(minIntSizes_B::SHORTB);
    for(size_t i = 0; i < nElements; ++i){
      dump_uShort value;
      memcpy(&value,&pin[pos],byteSize);
      pos += byteSize;
      //Covert value
      p[i] = static_cast<T>(value);
    }
    return PEN_DUMP_SUCCESS;
  }  

  if(containerBits == static_cast<size_t>(minIntSizes_b::INTb)){
    //This data requires a container of "INTb" bits
    const size_t byteSize = static_cast<size_t>(minIntSizes_B::INTB);
    for(size_t i = 0; i < nElements; ++i){
      dump_uInt value;
      memcpy(&value,&pin[pos],byteSize);
      pos += byteSize;
      //Covert value
      p[i] = static_cast<T>(value);
    }
    return PEN_DUMP_SUCCESS;
  }
  
  //Invalid container size
  return PEN_DUMP_INVALID_TYPE;
}


struct dArray{

  friend class pen_dump;

private:
  double* p;
  size_t n;

  size_t dumpBits;
public:

  static const size_t doubMem = sizeof(int64_t) + sizeof(int16_t);
  static const uint16_t doubleBits = 64 + 16;
  
  dArray(double* pin, const size_t nin) : p(pin), n(nin){
    dumpBits = doubleBits*n;
  }

  inline const double* getPointer(){return p;}
  
  inline void updateN(const size_t nin){
    n = nin;
    dumpBits = doubleBits*n;
  }

  inline bool isStored(const double* pin) const {
    return pin == p;
  }
  
};

struct iArray{

  friend class pen_dump;
  
private:
  intTypes type;
  union{
    short int*      ps;
    int*            pi;
    long int*       pl;
    long long int*  pll;
  };
  size_t n;
  minIntSizes_b containerBits;
  size_t containerByte;
  size_t dumpBits;
  size_t typeBits;
  
public:

  iArray(short int*      pin, const size_t nin) : type(intTypes::SHORT), ps(pin), n(nin){
    containerBits = minContainer<short int>(containerByte);
    dumpBits = static_cast<size_t>(containerBits)*n;
    typeBits = sizeofBits<short int>();
  }
  iArray(int*            pin, const size_t nin) : type(intTypes::INT)  , pi(pin), n(nin){
    containerBits = minContainer<int>(containerByte);
    dumpBits = static_cast<size_t>(containerBits)*n;
    typeBits = sizeofBits<int>();
  }
  iArray(long int*       pin, const size_t nin) : type(intTypes::LONG) , pl(pin), n(nin){
    containerBits = minContainer<long int>(containerByte);
    dumpBits = static_cast<size_t>(containerBits)*n;
    typeBits = sizeofBits<long int>();
  }
  iArray(long long int*  pin, const size_t nin) : type(intTypes::LLONG), pll(pin), n(nin){
    containerBits = minContainer<long long int>(containerByte);
    dumpBits = static_cast<size_t>(containerBits)*n;
    typeBits = sizeofBits<long long int>();
  }

  inline intTypes getType() const {return type;}
  inline const void* getPointer(){
    switch(type){
    case intTypes::SHORT: return (void*)ps;
    case intTypes::INT: return (void*)pi;
    case intTypes::LONG: return (void*)pl;
    case intTypes::LLONG: return (void*)pll;
    default: return nullptr;
    };
  }
  
  inline void updateN(const size_t nin){
    n = nin;
    dumpBits = static_cast<size_t>(containerBits)*n;
  }  

  inline bool isStored(const short int* pin) const {
    if(type != intTypes::SHORT)
      return false;
    return pin == ps;
  }
  inline bool isStored(const int* pin) const {
    if(type != intTypes::INT)
      return false;
    return pin == pi;
  }
  inline bool isStored(const long int* pin) const {
    if(type != intTypes::LONG)
      return false;
    return pin == pl;
  }
  inline bool isStored(const long long int* pin) const {
    if(type != intTypes::LLONG)
      return false;
    return pin == pll;
  }
  
};

struct uiArray{

  friend class pen_dump;

private:
  
  intTypes type;
  union{
    short unsigned int*      ps;
    unsigned int*            pi;
    long unsigned int*       pl;
    long long unsigned int*  pll;
  };
  
  size_t n;
  minIntSizes_b containerBits;
  size_t containerByte;
  size_t dumpBits;
  size_t typeBits;
public:
  
  uiArray(short unsigned int*      pin, const size_t nin) : type(intTypes::SHORT), ps(pin), n(nin){
    containerBits = minContainer<short unsigned int>(containerByte);
    dumpBits = static_cast<size_t>(containerBits)*n;
    typeBits = sizeofBits<short unsigned int>();
  }
  uiArray(unsigned int*            pin, const size_t nin) : type(intTypes::INT), pi(pin), n(nin){
    containerBits = minContainer<unsigned int>(containerByte);
    dumpBits = static_cast<size_t>(containerBits)*n;
    typeBits = sizeofBits<unsigned int>();
  }
  uiArray(long unsigned int*       pin, const size_t nin) : type(intTypes::LONG), pl(pin), n(nin){
    containerBits = minContainer<long unsigned int>(containerByte);
    dumpBits = static_cast<size_t>(containerBits)*n;
    typeBits = sizeofBits<long unsigned int>();
  }
  uiArray(long long unsigned int*  pin, const size_t nin) : type(intTypes::LLONG), pll(pin), n(nin){
    containerBits = minContainer<long long unsigned int>(containerByte);
    dumpBits = static_cast<size_t>(containerBits)*n;
    typeBits = sizeofBits<long long unsigned int>();
  }

  inline intTypes getType() const {return type;}
  inline const void* getPointer(){
    switch(type){
    case intTypes::SHORT: return (void*)ps;
    case intTypes::INT: return (void*)pi;
    case intTypes::LONG: return (void*)pl;
    case intTypes::LLONG: return (void*)pll;
    default: return nullptr;
    };
  }
  
  inline void updateN(const size_t nin){
    n = nin;
    dumpBits = static_cast<size_t>(containerBits)*n;
  }
  
  inline bool isStored(const short unsigned int* pin) const {
    if(type != intTypes::SHORT)
      return false;
    return pin == ps;
  }
  inline bool isStored(const unsigned int* pin) const {
    if(type != intTypes::INT)
      return false;
    return pin == pi;
  }
  inline bool isStored(const long unsigned int* pin) const {
    if(type != intTypes::LONG)
      return false;
    return pin == pl;
  }
  inline bool isStored(const long long unsigned int* pin) const {
    if(type != intTypes::LLONG)
      return false;
    return pin == pll;
  }  
  
};

struct ucArray{

  friend class pen_dump;
  
private:
  unsigned char* p;
  size_t n;

  size_t dumpBits;
  
public:

  static const size_t charMem = sizeof(unsigned char);  
  static const uint16_t charBits = CHAR_BIT;
  
  ucArray(unsigned char* pin, const size_t nin) : p(pin), n(nin){
    dumpBits = n*charBits;
  }

  inline const unsigned char* getPointer(){return p;}

  inline void updateN(const size_t nin){
    n = nin;
    dumpBits = n*charBits;
  }

  inline bool isStored(const unsigned char* pin) const {
    return pin == p;
  }
};

class pen_dump{

  typedef std::vector<dArray >::iterator iteratorD;
  typedef std::vector<iArray >::iterator iteratorI;
  typedef std::vector<uiArray>::iterator iteratorUI;
  typedef std::vector<ucArray>::iterator iteratorUC;
  
private:
  std::vector<dArray> pdouble;
  std::vector<iArray> pint;
  std::vector<uiArray> punsigned;
  std::vector<ucArray> puchar;

  std::vector<pen_dump*> subDumps;
  
  size_t dataBits; //Store dumped data size (bits)
  
  int dumpDouble(unsigned char* pout, size_t& pos) const;
  int dumpInt(unsigned char* pout, size_t& pos) const;
  int dumpUnsigned(unsigned char* pout, size_t& pos) const;
  int dumpChar(unsigned char* pout, size_t& pos) const;
  int dumpSubDumps(unsigned char* pout,
		   size_t& pos,
		   const size_t outputSize,
		   const unsigned verbose) const;

  int readDouble(const unsigned char* const pin, size_t& pos, const unsigned verbose = 0);
  int readInt(const unsigned char* const pin, size_t& pos, const unsigned verbose = 0);
  int readUnsigned(const unsigned char* const pin, size_t& pos, const unsigned verbose = 0);
  int readChar(const unsigned char* const pin, size_t& pos, const unsigned verbose = 0);
  int readSubDumps(const unsigned char* const pin, size_t& pos, const unsigned verbose);  

public:
  
  static const size_t doubMem = dArray::doubMem;
  static const size_t charMem = ucArray::charMem;

  static const uint16_t doubleBits = dArray::doubleBits;
  static const uint16_t charBits   = ucArray::charBits;

  //Storage size (in bits) for number of arrais and elements
  static const size_t metadataNelem = 32;
  //Storage size (in bits) for element size
  static const size_t metadataESize = 16;
  //Storage size (in bits) for elements used bits 
  static const size_t metadataEUsed = 16;
  
  static const uint64_t precision = 100000000000000000; //1e17  
  
  pen_dump();
  
  inline size_t nDoubles() const {return pdouble.size();}
  inline size_t nInts() const {return pint.size();}
  inline size_t nUnsigneds() const {return punsigned.size();}
  inline size_t nChars() const {return puchar.size();}
  inline size_t nSubDumps() const{return subDumps.size();}

  inline size_t nRegistered() const {
    return nDoubles()+nInts()+nUnsigneds()+nChars()+nSubDumps();
  }

  inline size_t memory() const {
    size_t mem = dataBits/charBits;
    for(const pen_dump* p : subDumps)
      mem += p->memory();
    return mem;
  }
  inline size_t bits() const {
    size_t b = dataBits;
    for(const pen_dump* p : subDumps)
      b += p->bits();
    return b;
  }

  inline int toDump(pen_dump& subDump){

    //Check if it is already stored
    auto it = std::find(subDumps.begin(), subDumps.end(), &subDump);
    if(it == subDumps.end()){
      subDumps.push_back(&subDump);
    }

    return PEN_DUMP_SUCCESS;    
  }
  
  int toDump(double*        p, const size_t n);

  // toDump for signed integrals:
  template <class signedT>

  // Enable only overloads with signed types using return type
  typename std::enable_if<std::is_signed<signedT>::value,int>::type //(int)  
  toDump(signedT*       p, const size_t n){

    if(p == nullptr)
      return PEN_DUMP_NULL_POINTER;

    if(n == 0)
      return PEN_DUMP_ZERO_ELEMENTS;

    //Check if this pointer has already been registered
    for(iteratorI i = pint.begin(); i != pint.end(); ++i){

      if(i->isStored(p)){
	
	dataBits -= i->dumpBits;
	i->updateN(n);
	dataBits += i->dumpBits;

	return PEN_DUMP_SUCCESS;
      }
    }

    //Add new element to signed integer array 
    pint.push_back(iArray(p,n));

    iArray* plast = &(pint.back());

    //Check if this integer type fits in 64 bits element
    if(plast->containerByte == static_cast<size_t>(minIntSizes_B::INVALIDB)){
      pint.pop_back(); //Remove element
      return PEN_DUMP_INTEGER_LARGER_THAN_MAXIMUM;
    }
    
    //Add required memory for all elements
    dataBits += plast->dumpBits;
    //Add required memory to store elements size
    dataBits += metadataESize;
    //Add used bits per element
    dataBits += metadataEUsed;
    //Add required memory to store the number of elements
    dataBits += metadataNelem;

    return PEN_DUMP_SUCCESS;    
    
  }

  // toDump for unsigned integrals:
  template <class unsignedT>
  // Enable only overloads with unsigned types using return type
  typename std::enable_if<std::is_unsigned<unsignedT>::value,int>::type //(int)
  toDump(unsignedT*     p, const size_t n){

    if(p == nullptr)
      return PEN_DUMP_NULL_POINTER;

    if(n == 0)
      return PEN_DUMP_ZERO_ELEMENTS;

    //Check if this pointer has already been registered
    for(iteratorUI i = punsigned.begin(); i != punsigned.end(); ++i){

      if(i->isStored(p)){
	
	dataBits -= i->dumpBits;
	i->updateN(n);
	dataBits += i->dumpBits;

	return PEN_DUMP_SUCCESS;
      }
    }

    //Add new element to unsigned integer array 
    punsigned.push_back(uiArray(p,n));

    uiArray* plast = &(punsigned.back());

    //Check if this integer type fits in 64 bits element
    if(plast->containerByte == static_cast<size_t>(minIntSizes_B::INVALIDB)){
      punsigned.pop_back(); //Remove element
      return PEN_DUMP_INTEGER_LARGER_THAN_MAXIMUM;
    }
    
    //Add required memory for all elements
    dataBits += plast->dumpBits;
    //Add required memory to store elements size
    dataBits += metadataESize;
    //Add used bits per element
    dataBits += metadataEUsed;
    //Add required memory to store the number of elements
    dataBits += metadataNelem;

    return PEN_DUMP_SUCCESS;    
    
  }
  
  int toDump(unsigned char* p, const size_t n);

  inline int remove(const pen_dump* subDump){

    if(subDump == nullptr)
      return PEN_DUMP_NULL_POINTER;

    auto it = std::find(subDumps.begin(), subDumps.end(), subDump);
    if(it != subDumps.end())
      subDumps.erase(it);
    
    return PEN_DUMP_SUCCESS;
  }
  
  int remove(const double*        p);

  template <class signedT>
  // Enable only overloads with signed types using return type
  typename std::enable_if<std::is_signed<signedT>::value,int>::type //(int)  
  remove(const signedT* p){

    if(p == nullptr)
      return PEN_DUMP_NULL_POINTER;
    
    //Check if this pointer has already been registered
    std::vector<iArray>::iterator it;
    for(it = pint.begin(); it != pint.end(); ++it){
      
      if(it->isStored(p)){
	dataBits -= it->dumpBits;
	dataBits -= metadataESize;
	dataBits -= metadataEUsed;
	dataBits -= metadataNelem;
	pint.erase(it);
	return PEN_DUMP_SUCCESS;  
      }
    }
    
    return PEN_DUMP_ELEMENT_NOT_FOUND;    
  }
  
  template <class unsignedT>
  // Enable only overloads with unsigned types using return type
  typename std::enable_if<!(std::is_signed<unsignedT>::value),int>::type //(int)
  remove(const unsignedT* p){

    if(p == nullptr)
      return PEN_DUMP_NULL_POINTER;
    
    //Check if this pointer has already been registered
    std::vector<uiArray>::iterator it;
    for(it = punsigned.begin(); it != punsigned.end(); ++it){
      
      if(it->isStored(p)){
	dataBits -= it->dumpBits;
	dataBits -= metadataESize;
	dataBits -= metadataEUsed;
	dataBits -= metadataNelem;
	punsigned.erase(it);
	return PEN_DUMP_SUCCESS;  
      }
    }
    
    return PEN_DUMP_ELEMENT_NOT_FOUND;    
  }
    
  int remove(const unsigned char* p);
  
  int dump(unsigned char*& pout,
	   size_t& written,
	   const size_t outputSize,
	   const unsigned verbose) const;
  
  int read(const unsigned char* const pin,
	   size_t& pos,
	   const unsigned verbose = 0);
  
  inline void clear(){
    pdouble.clear();
    pint.clear();
    punsigned.clear();
    puchar.clear();
    subDumps.clear();
    dataBits = 0;
    //Add global metadata for double arrays (num arrays and element bits)
    dataBits += metadataNelem + metadataESize; 
    //Add global metadata for integer arrays
    dataBits += metadataNelem;   //num arrays 
    //Add global metadata for unsigned arrays
    dataBits += metadataNelem;   //num arrays
    //Add global metadata for char arrays (num arrays and element bits)
    dataBits += metadataNelem + metadataESize;
    //Add global metadata for subdumps (num subdumps)
    dataBits += metadataNelem;   //num sub dumps
  }

  ~pen_dump();
};


#endif
