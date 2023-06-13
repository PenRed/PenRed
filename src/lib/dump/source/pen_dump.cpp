
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



#include "pen_dump.hh"

const uint16_t pen_dump::doubleBits;
const uint16_t pen_dump::charBits;

pen_dump::pen_dump()
{
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

int pen_dump::toDump(double*        p, const size_t n){

  if(p == nullptr)
    return PEN_DUMP_NULL_POINTER;

  if(n == 0)
    return PEN_DUMP_ZERO_ELEMENTS;

  //Check if this pointer has already been registered
  for(iteratorD i = pdouble.begin(); i != pdouble.end(); ++i){

      if(i->isStored(p)){
	
	dataBits -= i->dumpBits;
	i->updateN(n);
	dataBits += i->dumpBits;

	return PEN_DUMP_SUCCESS;
      }
  }

  //Add new element to double array 
  pdouble.push_back(dArray(p,n));

  dArray* plast = &(pdouble.back());
  
  //Add required memory for all elements
  dataBits += plast->dumpBits;
  //Add required memory to store the number of elements
  dataBits += metadataNelem;
  
  return PEN_DUMP_SUCCESS;
}
int pen_dump::toDump(unsigned char* p, const size_t n){

  if(p == nullptr)
    return PEN_DUMP_NULL_POINTER;

  if(n == 0)
    return PEN_DUMP_ZERO_ELEMENTS;

  //Check if this pointer has already been registered
  for(iteratorUC i = puchar.begin(); i != puchar.end(); ++i){

      if(i->isStored(p)){
	
	dataBits -= i->dumpBits;
	i->updateN(n);
	dataBits += i->dumpBits;

	return PEN_DUMP_SUCCESS;
      }
  }

  //Add new element to unsigned integer array 
  puchar.push_back(ucArray(p,n));

  ucArray* plast = &(puchar.back());
  
  //Add required memory for all elements
  dataBits += plast->dumpBits;
  //Add required memory to store the number of elements
  dataBits += metadataNelem;

  return PEN_DUMP_SUCCESS;
}

int pen_dump::remove(const double*        p){

  if(p == nullptr)
    return PEN_DUMP_NULL_POINTER;
  
  //Check if this pointer has already been registered
  std::vector<dArray>::iterator it;
  for(it = pdouble.begin(); it != pdouble.end(); ++it){
  
    if(it->isStored(p)){
      dataBits -= it->dumpBits;
      dataBits -= metadataNelem;
      pdouble.erase(it);
      return PEN_DUMP_SUCCESS;
    }
  }
  
  return PEN_DUMP_ELEMENT_NOT_FOUND;
}
int pen_dump::remove(const unsigned char* p){

  if(p == nullptr)
    return PEN_DUMP_NULL_POINTER;
  
  //Check if this pointer has already been registered
  std::vector<ucArray>::iterator it;
  for(it = puchar.begin(); it != puchar.end(); ++it){
  
    if(it->isStored(p)){
      dataBits -= it->dumpBits;
      dataBits -= metadataNelem;
      puchar.erase(it);
      return PEN_DUMP_SUCCESS;      
    }
  }
  return PEN_DUMP_ELEMENT_NOT_FOUND;
}

int pen_dump::dumpDouble(unsigned char* pout, size_t& pos) const{
  
  const uint32_t nArrays = pdouble.size();
  //Store number of arrays
  memcpy(&pout[pos],&nArrays,sizeof(uint32_t));
  pos += sizeof(uint32_t);
  //Store bits per element
  memcpy(&pout[pos],&doubleBits,sizeof(uint16_t));
  pos += sizeof(uint16_t);
    
  //Iterate over double arrays
  std::vector<dArray>::const_iterator it;
  for(it = pdouble.begin(); it != pdouble.end(); ++it){

    double* p = it->p;
    //Save number of elements in array
    uint32_t nElements = it->n;
    memcpy(&pout[pos],&nElements,sizeof(uint32_t));
    pos += sizeof(uint32_t);

    //Save elements
    for(size_t i = 0; i < nElements; i++){

      //Split number in significand and exponent part
      // num = signif * 2^exp
      int exponent;
      double significand;
      
      significand = frexp(p[i],&exponent);
      
      //Convert value to ints
      int64_t sig64 = static_cast<int64_t>(significand*precision);
      int16_t exp16 = static_cast<int16_t>(exponent);

      //Save values
      memcpy(&pout[pos],&sig64,sizeof(int64_t));
      pos += sizeof(int64_t);
      memcpy(&pout[pos],&exp16,sizeof(int16_t));
      pos += sizeof(int16_t);
    }
  }
  return PEN_DUMP_SUCCESS;
}

int pen_dump::dumpInt(unsigned char* pout, size_t& pos) const{
  
  uint32_t nArrays = pint.size();
  memcpy(&pout[pos],&nArrays,sizeof(uint32_t));
  pos += sizeof(uint32_t);
  
  //Iterate over integer arrays
  std::vector<iArray>::const_iterator it;
  for(it = pint.begin(); it != pint.end(); ++it){
    
    //Store element memory size
    uint16_t elementMem = static_cast<uint16_t>(it->containerBits);
    memcpy(&pout[pos],&elementMem,sizeof(uint16_t));
    pos += sizeof(uint16_t);

    //Store element used bits
    uint16_t usedBits = it->typeBits;
    memcpy(&pout[pos],&usedBits,sizeof(uint16_t));
    pos += sizeof(uint16_t);
    
    //Save number of elements
    uint32_t nElements = it->n;
    memcpy(&pout[pos],&nElements,sizeof(uint32_t));
    pos += sizeof(uint32_t);    
    
    //Check data type and save it
    if(it->getType() == intTypes::LONG){
      dumpSInt(pout,pos,it->pl,nElements,it->containerBits);
    }
    else if(it->getType() == intTypes::LLONG){
      dumpSInt(pout,pos,it->pll,nElements,it->containerBits);
    }
    else if(it->getType() == intTypes::SHORT){
      dumpSInt(pout,pos,it->ps,nElements,it->containerBits);      
    }
    else if(it->getType() == intTypes::INT){
      dumpSInt(pout,pos,it->pi,nElements,it->containerBits);
    }
    else{
      return PEN_DUMP_INVALID_TYPE;
    }
    
  }
    
  return PEN_DUMP_SUCCESS;
}

int pen_dump::dumpUnsigned(unsigned char* pout, size_t& pos) const{
  
  uint32_t nArrays = punsigned.size();
  memcpy(&pout[pos],&nArrays,sizeof(uint32_t));
  pos += sizeof(uint32_t);
  
  //Iterate over integer arrays
  std::vector<uiArray>::const_iterator it;
  for(it = punsigned.begin(); it != punsigned.end(); ++it){

    //Store element memory size
    uint16_t elementMem = static_cast<uint16_t>(it->containerBits);
    memcpy(&pout[pos],&elementMem,sizeof(uint16_t));
    pos += sizeof(uint16_t);

    //Store element used bits
    uint16_t usedBits = it->typeBits;
    memcpy(&pout[pos],&usedBits,sizeof(uint16_t));
    pos += sizeof(uint16_t);

    //Save number of elements
    uint32_t nElements = it->n;    
    memcpy(&pout[pos],&nElements,sizeof(uint32_t));
    pos += sizeof(uint32_t);

    //Check data type and save it
    if(it->getType() == intTypes::LONG){
      dumpUInt(pout,pos,it->pl,nElements,it->containerBits);
    }
    else if(it->getType() == intTypes::LLONG){
      dumpUInt(pout,pos,it->pll,nElements,it->containerBits);
    }
    else if(it->getType() == intTypes::SHORT){
      dumpUInt(pout,pos,it->ps,nElements,it->containerBits);      
    }
    else if(it->getType() == intTypes::INT){
      dumpUInt(pout,pos,it->pi,nElements,it->containerBits);
    }
    else{
      return PEN_DUMP_INVALID_TYPE;
    }
    
    
  }
    
  return PEN_DUMP_SUCCESS;
}

int pen_dump::dumpChar(unsigned char* pout, size_t& pos) const{
  
  uint32_t nArrays = puchar.size();
  memcpy(&pout[pos],&nArrays,sizeof(uint32_t));
  pos += sizeof(uint32_t);
  //Store bits per element
  memcpy(&pout[pos],&charBits,sizeof(uint16_t));
  pos += sizeof(uint16_t);
  
  //Iterate over integer arrays
  std::vector<ucArray>::const_iterator it;
  for(it = puchar.begin(); it != puchar.end(); ++it){
    unsigned char* p = it->p;

    //Save number of elements in array and elements size
    uint32_t nElements = it->n;    
    memcpy(&pout[pos],&nElements,sizeof(uint32_t));
    pos += sizeof(uint32_t);

    //Save data
    memcpy(&pout[pos],p,charMem*nElements);
    pos += charMem*nElements;
  }
    
  return PEN_DUMP_SUCCESS;
}

int pen_dump::dumpSubDumps(unsigned char* pout,
			   size_t& pos,
			   const size_t outputSize,
			   const unsigned verbose) const{
  
  const uint32_t nArrays = subDumps.size();
  //Write number of sub dumps to store 
  memcpy(&pout[pos],&nArrays,sizeof(uint32_t));
  pos += sizeof(uint32_t);

  //Write sub dumps
  for(const pen_dump* p : subDumps){
    size_t subWritten;
    unsigned char* nextp = pout+pos;
    int err = p->dump(nextp,subWritten,outputSize,verbose);
    pos += subWritten;
    if(err != PEN_DUMP_SUCCESS){
      if(verbose > 0){
	auto it = find(subDumps.begin(), subDumps.end(), p);
	printf("dumpSubDumps:Error: Error dumping sub dump %ld.\n", it-subDumps.begin());
	printf("             Error code: %d\n",err);
      }
      return err;
    }
  }

  return PEN_DUMP_SUCCESS;
}

int pen_dump::dump(unsigned char*& pout,
		   size_t& written,
		   const size_t outputSize,
		   const unsigned verbose) const{

  size_t finalOutSize = outputSize;
  if(outputSize == 0){
    //Allocate memmory for all elements
    finalOutSize = memory();
    pout = nullptr;
    pout = (unsigned char*) malloc(finalOutSize);
    if(pout == nullptr)
      return PEN_DUMP_BAD_ALLOCATION;
  }
  else{
    //Check if the provided array is not null
    if(pout == nullptr)
      return PEN_DUMP_NULL_POINTER;
    
    //Check if output array has the required dimension
    if(memory() > outputSize)
      return PEN_DUMP_INSUFFICIENT_MEMORY;
  }
  
  written = 0;
  int err = 0;
  //Dump double arrays
  err = dumpDouble(pout,written);
  if(err != PEN_DUMP_SUCCESS){
    if(verbose > 0){
      printf("dump:Error: Error dumping double arrays.\n");
      printf("            Error code: %d\n",err);
    }
    free(pout);
    pout = nullptr;
    return PEN_DUMP_ERROR_DOUBLE_DUMP;
  }

  //Dump integer arrays
  err = dumpInt(pout,written);
  if(err != PEN_DUMP_SUCCESS){
    if(verbose > 0){
      printf("dump:Error: Error dumping integer arrays.\n");
      printf("            Error code: %d\n",err);
    }
    free(pout);
    pout = nullptr;
    return PEN_DUMP_ERROR_INT_DUMP;
  }

  //Dump unsigned arrays
  err = dumpUnsigned(pout,written);
  if(err != PEN_DUMP_SUCCESS){
    if(verbose > 0){
      printf("dump:Error: Error dumping unsigned arrays.\n");
      printf("            Error code: %d\n",err);
    }
    free(pout);
    pout = nullptr;
    return PEN_DUMP_ERROR_UNSIGNED_DUMP;
  }

  //Dump char arrays
  err = dumpChar(pout,written);
  if(err != PEN_DUMP_SUCCESS){
    if(verbose > 0){
      printf("dump:Error: Error dumping char arrays.\n");
      printf("            Error code: %d\n",err);
    }
    free(pout);
    pout = nullptr;
    return PEN_DUMP_ERROR_CHAR_DUMP;
  }

  //Dump sub dumps
  err = dumpSubDumps(pout,written,finalOutSize,verbose);
  if(err != PEN_DUMP_SUCCESS){
    return err;
  }
  
  //Check written data
  if(written != memory()){
    if(verbose > 0){
      printf("dump:Error: written data bytes mismatch with expected.\n");
      printf("            written: %lu\n",written);
      printf("            expected: %lu\n",memory());
    }
    free(pout);
    pout = nullptr;
    return PEN_DUMP_INCORRECT_DATA_SIZE;
  }
  
  return PEN_DUMP_SUCCESS;
}

int pen_dump::readDouble(const unsigned char* const pin,
			 size_t& pos,
			 const unsigned verbose){

  uint32_t nArrays;
  uint16_t elementMem; //In bits
  memcpy(&nArrays,&pin[pos],sizeof(uint32_t));
  pos += sizeof(uint32_t);
  memcpy(&elementMem,&pin[pos],sizeof(uint16_t));
  pos += sizeof(uint16_t);

  //Check array number
  if(nArrays != pdouble.size()){
    if(verbose > 0){
      printf("pen_dump:readDouble: Error: Number of arrays mismatch.\n");
      printf("                   Read: %u\n",nArrays);
      printf("               Expected: %lu\n",pdouble.size());
    }
    return PEN_DUMP_NARRAY_NOT_MATCH;
  }

  //Check element memory size
  if(elementMem != doubleBits){
    if(verbose > 0){
      printf("pen_dump:readDouble: Error: Element memory size mismatch.\n");
      printf("                   Read: %u bits\n",elementMem);
      printf("               Expected: %u bits\n",doubleBits);
      printf("  Has this dumped file generated with an older version of this lib?\n");
    }    
    return PEN_DUMP_ELEMENT_SIZE_NOT_MATCH;
  }

  //Read all arrays
  for(size_t i = 0; i < nArrays; i++){
    
    double* p = pdouble[i].p;
    //Read number of elements in array
    uint32_t nElements;
    memcpy(&nElements,&pin[pos],sizeof(uint32_t));
    pos += sizeof(uint32_t);

    if(nElements != pdouble[i].n){
      if(verbose > 0){
	printf("pen_dump:readDouble: Error: Number of elements in array %lu mismatch.\n",i);
	printf("                   Read: %u\n",nElements);
	printf("               Expected: %lu\n",pdouble[i].n);
      }
      return PEN_DUMP_ELEMENT_NUMBER_NOT_MATCH;
    }

    double Dprecision = double(precision);

    //Read elements
    for(size_t j = 0; j < nElements; ++j){

      int64_t sig64;
      int16_t exp16;

      //Get number parts
      memcpy(&sig64,&pin[pos],sizeof(int64_t));
      pos += sizeof(int64_t);
      memcpy(&exp16,&pin[pos],sizeof(int16_t));
      pos += sizeof(int16_t);
      
      //Join number using significand and exponent part
      // num = signif * 2^exp
      p[j] = ldexp((double)sig64/Dprecision,int(exp16));
    }
  }
  return PEN_DUMP_SUCCESS;
}

int pen_dump::readInt(const unsigned char* const pin,
		      size_t& pos,
		      const unsigned verbose){

  uint32_t nArrays;
  memcpy(&nArrays,&pin[pos],sizeof(uint32_t));
  pos += sizeof(uint32_t);

  //Check array number
  if(nArrays != pint.size()){
    if(verbose > 0){
      printf("pen_dump:readInt: Error: Number of arrays mismatch.\n");
      printf("                   Read: %u\n",nArrays);
      printf("               Expected: %lu\n",pint.size());
    }    
    return PEN_DUMP_NARRAY_NOT_MATCH;
  }

  //Iterate over integer arrays
  std::vector<iArray>::const_iterator it;
  for(it = pint.begin(); it != pint.end(); ++it){

    //Read element size in bits
    uint16_t elementMem;
    memcpy(&elementMem,&pin[pos],sizeof(uint16_t));
    pos += sizeof(uint16_t);

    //Read used bits per element
    uint16_t usedBits;
    memcpy(&usedBits,&pin[pos],sizeof(uint16_t));
    pos += sizeof(uint16_t);

    //Check if elements to read fits in the specified type
    if(usedBits > it->typeBits){
      if(verbose > 0){
	printf("pen_dump:readInt: Error: Element memory size doesn't fit in specified type.\n");
	printf("                   Read: %u bits\n",usedBits);
	printf("               Expected: %lu bits\n",it->typeBits);
      }    
      return PEN_DUMP_ELEMENT_SIZE_NOT_MATCH;
    }
    
    //Read number of elements in array
    uint32_t nElements;
    memcpy(&nElements,&pin[pos],sizeof(uint32_t));
    pos += sizeof(uint32_t);

    if(nElements != it->n){
      if(verbose > 0){
	printf("pen_dump:readInt: Error: Number of elements in array %lu mismatch.\n",size_t(it - pint.begin()));
	printf("                   Read: %u\n",nElements);
	printf("               Expected: %lu\n",it->n);
      }      
      return PEN_DUMP_ELEMENT_NUMBER_NOT_MATCH;
    }

    //Check data type and read it
    if(it->getType() == intTypes::LONG){
      readSInt(pin,pos,it->pl,nElements,elementMem);
    }
    else if(it->getType() == intTypes::LLONG){
      readSInt(pin,pos,it->pll,nElements,elementMem);
    }
    else if(it->getType() == intTypes::SHORT){
      readSInt(pin,pos,it->ps,nElements,elementMem);      
    }
    else if(it->getType() == intTypes::INT){
      readSInt(pin,pos,it->pi,nElements,elementMem);
    }
    else{
      return PEN_DUMP_INVALID_TYPE;
    }
  }
  return PEN_DUMP_SUCCESS;
}

int pen_dump::readUnsigned(const unsigned char* const pin,
			   size_t& pos,
			   const unsigned verbose){

  uint32_t nArrays;
  memcpy(&nArrays,&pin[pos],sizeof(uint32_t));
  pos += sizeof(uint32_t);

  //Check array number
  if(nArrays != punsigned.size()){
    if(verbose > 0){
      printf("pen_dump:readUnsigned: Error: Number of arrays mismatch.\n");
      printf("                   Read: %u\n",nArrays);
      printf("               Expected: %lu\n",punsigned.size());
    }    
    return PEN_DUMP_NARRAY_NOT_MATCH;
  }

  //Iterate over unsigned integer arrays
  std::vector<uiArray>::const_iterator it;
  for(it = punsigned.begin(); it != punsigned.end(); ++it){

    //Read element size in bits
    uint16_t elementMem;
    memcpy(&elementMem,&pin[pos],sizeof(uint16_t));
    pos += sizeof(uint16_t);

    //Read used bits per element
    uint16_t usedBits;
    memcpy(&usedBits,&pin[pos],sizeof(uint16_t));
    pos += sizeof(uint16_t);    

    //Check if elements to read fits in the specified type
    if(usedBits > it->typeBits){
      if(verbose > 0){
	printf("pen_dump:readInt: Error: Element memory size doesn't fit in specified type.\n");
	printf("                   Read: %u bits\n",usedBits);
	printf("               Expected: %lu bits\n",it->typeBits);
      }    
      return PEN_DUMP_ELEMENT_SIZE_NOT_MATCH;
    }
    
    //Read number of elements in array
    uint32_t nElements;
    memcpy(&nElements,&pin[pos],sizeof(uint32_t));
    pos += sizeof(uint32_t);

    if(nElements != it->n){
      if(verbose > 0){
	printf("pen_dump:readInt: Error: Number of elements in array %lu mismatch.\n",size_t(it - punsigned.begin()));
	printf("                   Read: %u\n",nElements);
	printf("               Expected: %lu\n",it->n);
      }      
      return PEN_DUMP_ELEMENT_NUMBER_NOT_MATCH;
    }

    //Check data type and read it
    if(it->getType() == intTypes::LONG){
      readUInt(pin,pos,it->pl,nElements,elementMem);
    }
    else if(it->getType() == intTypes::LLONG){
      readUInt(pin,pos,it->pll,nElements,elementMem);
    }
    else if(it->getType() == intTypes::SHORT){
      readUInt(pin,pos,it->ps,nElements,elementMem);      
    }
    else if(it->getType() == intTypes::INT){
      readUInt(pin,pos,it->pi,nElements,elementMem);
    }
    else{
      return PEN_DUMP_INVALID_TYPE;
    }
  }
  return PEN_DUMP_SUCCESS;
}

int pen_dump::readChar(const unsigned char* const pin,
		       size_t& pos,
		       const unsigned verbose){

  uint32_t nArrays;
  uint16_t elementMem;
  memcpy(&nArrays,&pin[pos],sizeof(uint32_t));
  pos += sizeof(uint32_t);
  memcpy(&elementMem,&pin[pos],sizeof(uint16_t));
  pos += sizeof(uint16_t);

  //Check array number
  if(nArrays != puchar.size()){
    if(verbose > 0){
      printf("pen_dump:readChar: Error: Number of arrays mismatch.\n");
      printf("                   Read: %u\n",nArrays);
      printf("               Expected: %lu\n",puchar.size());
    }    
    return PEN_DUMP_NARRAY_NOT_MATCH;
  }

  //Read all arrays
  for(size_t i = 0; i < nArrays; i++){
    
    unsigned char* p = puchar[i].p;
    //Read number of elements in array
    uint32_t nElements;
    memcpy(&nElements,&pin[pos],sizeof(uint32_t));
    pos += sizeof(uint32_t);

    //Recalculate the number of elements to read according the read element size
    //and the byte size of current machine
    nElements = nElements*elementMem/charBits;
    
    if(nElements != puchar[i].n){
      if(verbose > 0){
	printf("pen_dump:readChar: Error: Number of elements in array %lu mismatch.\n",i);
	printf("                   Read: %u\n",nElements);
	printf("               Expected: %lu\n",puchar[i].n);
      }
      return PEN_DUMP_ELEMENT_NUMBER_NOT_MATCH;
    }

    //Extract elements
    memcpy(p,&pin[pos],nElements);
    pos += nElements;
  }
  return PEN_DUMP_SUCCESS;
}

int pen_dump::readSubDumps(const unsigned char* const pin,
			   size_t& pos,
			   const unsigned verbose){

  uint32_t nSubDumps;
  memcpy(&nSubDumps,&pin[pos],sizeof(uint32_t));
  pos += sizeof(uint32_t);

  //Check sub dumps number
  if(nSubDumps != subDumps.size()){
    if(verbose > 0){
      printf("pen_dump:readSubDumps: Error: Number of subdumps mismatch.\n");
      printf("                   Read: %u\n",nSubDumps);
      printf("               Expected: %lu\n",subDumps.size());
    }
    return PEN_DUMP_NARRAY_NOT_MATCH;
  }

  //Read all sub dumps
  for(pen_dump* p : subDumps){
    int err = p->read(pin,pos,verbose);
    if(err != PEN_DUMP_SUCCESS){
      return err;
    }
  }
  
  return PEN_DUMP_SUCCESS;
}


int pen_dump::read(const unsigned char* const pin,
		   size_t& pos,
		   const unsigned verbose){

  if(pin == nullptr){
    if(verbose > 0){
      printf("pen_dump:read:Error: Input buffer is a Null pointer.\n");
    }
    return PEN_DUMP_NULL_POINTER;
  }
  
  int err;

  //Read doubles
  err = readDouble(pin,pos,verbose);
  if(err != PEN_DUMP_SUCCESS){
    if(verbose > 0){
      printf("pen_dump:read:Error: unable to read double arrays.\n");
      printf("                     Error code: %d\n",err);
    }
    return PEN_DUMP_UNABLE_TO_READ_DOUBLE_ARRAYS;
  }

  //Read integers
  err = readInt(pin,pos,verbose);
  if(err != PEN_DUMP_SUCCESS){
    if(verbose > 0){
      printf("pen_dump:read:Error: unable to read int arrays.\n");
      printf("                     Error code: %d\n",err);
    }
    return PEN_DUMP_UNABLE_TO_READ_INT_ARRAYS;
  }

  //Read unsigned integers
  err = readUnsigned(pin,pos,verbose);
  if(err != PEN_DUMP_SUCCESS){
    if(verbose > 0){
      printf("pen_dump:read:Error: unable to read unsigned int arrays.\n");
      printf("                     Error code: %d\n",err);
    }
    return PEN_DUMP_UNABLE_TO_READ_UNSIGNED_ARRAYS;
  }

  //Read chars
  err = readChar(pin,pos,verbose);
  if(err != PEN_DUMP_SUCCESS){
    if(verbose > 0){
      printf("pen_dump:read:Error: unable to read char arrays.\n");
      printf("                     Error code: %d\n",err);
    }
    return PEN_DUMP_UNABLE_TO_READ_CHAR_ARRAYS;
  }

  //Read sub dumps
  err = readSubDumps(pin,pos,verbose);
  if(err != PEN_DUMP_SUCCESS){
    if(verbose > 0){
      printf("pen_dump:read:Error: unable to read sub dumps.\n");
      printf("                     Error code: %d\n",err);
    }
    return PEN_DUMP_UNABLE_TO_READ_CHAR_ARRAYS;
  }
  
  return PEN_DUMP_SUCCESS;
}

pen_dump::~pen_dump()
{}
