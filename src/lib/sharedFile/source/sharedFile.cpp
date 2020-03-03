
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

 
#include "sharedFile.hh"

int pen_sharedFile::createReader(const unsigned ID){

  //Lock file
  std::lock_guard<std::mutex> guard(fileLock);

  //Check if this ID is being used
  if(trusted_findID(ID) >= 0){
    //Used ID
    return SHARED_FILE_USED_READER_ID;
  }

  //Store ID and set position to origin
  readers.push_back(std::make_pair(ID,origin));

  return SHARED_FILE_SUCCESS;
}

bool pen_sharedFile::isBinary(){
  //Lock file
  std::lock_guard<std::mutex> guard(fileLock);
  return binary;
}

size_t pen_sharedFile::size(){
  //Lock file
  std::lock_guard<std::mutex> guard(fileLock);
  return fsize;
}

int pen_sharedFile::tell(const unsigned ID,
			 long int &tellPos){

  //Lock file
  std::lock_guard<std::mutex> guard(fileLock);
  
  if(pfile == nullptr){
    return SHARED_FILE_NO_FILE_OPENNED;
  }
  
  //Check if this ID is being used
  int pos = trusted_findID(ID);
  if(pos < 0)
    return SHARED_FILE_UNUSED_READER_ID;

  //Set file to readers position
  if(fsetpos(pfile,&(readers[pos].second)) != 0)
    return SHARED_FILE_FSETPOS_ERROR;
  
  //Call ftell
  tellPos = ftell(pfile);

  return SHARED_FILE_SUCCESS;
}

int pen_sharedFile::seek(const unsigned ID,
			 long int offset,
			 int seekOrigin,
			 int* fseekerr){

  //Lock file
  std::lock_guard<std::mutex> guard(fileLock);
  
  if(pfile == nullptr){
    return SHARED_FILE_NO_FILE_OPENNED;
  }
  
  //Check if this ID is being used
  int pos = trusted_findID(ID);
  if(pos < 0)
    return SHARED_FILE_UNUSED_READER_ID;

  //Set cursor to current file position to handle
  //fseeks at cursor position
  if(fsetpos(pfile,&(readers[pos].second)) != 0)
    return SHARED_FILE_FSETPOS_ERROR;
  
  //Set file to specified position
  int err = fseek(pfile,offset,seekOrigin);

  //Store fseek return value
  if(fseekerr != nullptr)
    *fseekerr = err;

  //Check if fseek has been executed successfuly
  if(err != 0)
    return SHARED_FILE_FSEEK_ERROR;
    
  //Store new position
  if(fgetpos(pfile,&(readers[pos].second)) != 0){
    return SHARED_FILE_FGETPOS_ERROR;
  }
  
  return SHARED_FILE_SUCCESS;
}

int pen_sharedFile::getPos(const unsigned ID,
			   fpos_t& fpos){

  //Lock file
  std::lock_guard<std::mutex> guard(fileLock);
  
  if(pfile == nullptr){
    return SHARED_FILE_NO_FILE_OPENNED;
  }
  
  //Check if this ID is being used
  int pos = trusted_findID(ID);
  if(pos < 0)
    return SHARED_FILE_UNUSED_READER_ID;
  
  //Fill position
  fpos = readers[pos].second;

  return SHARED_FILE_SUCCESS;
}

int pen_sharedFile::setPos(const unsigned ID,
			   const fpos_t& fpos){

  //Lock file
  std::lock_guard<std::mutex> guard(fileLock);
  
  if(pfile == nullptr){
    return SHARED_FILE_NO_FILE_OPENNED;
  }
  
  //Check if this ID is being used
  int pos = trusted_findID(ID);
  if(pos < 0)
    return SHARED_FILE_UNUSED_READER_ID;
  
  //Set file to specified position
  if(fsetpos(pfile,&fpos) != 0)
    return SHARED_FILE_FSETPOS_ERROR;

  //Store input position for this reader
  readers[pos].second = fpos;

  return SHARED_FILE_SUCCESS;
}

int pen_sharedFile::open(const char* filename,
			 const bool binaryFile){

  //Lock file
  std::lock_guard<std::mutex> guard(fileLock);
  
  //Close previous file
  trusted_close();

  //Open specified file in binary mode to get the total size
  pfile = fopen(filename,"rb");
  //Check openned file
  if(pfile == nullptr)
    return SHARED_FILE_UNABLE_TO_OPEN_FILE;

  unsigned char buffer[5000];
  fsize = 0;
  size_t nread = 0;
  while((nread = fread(buffer,sizeof(unsigned char),5000,pfile)) > 0){
    fsize += nread;
  }
  fclose(pfile);
  
  //Open file
  if(binaryFile){
    pfile = fopen(filename,"rb");
  }
  else{
    pfile = fopen(filename,"r");
  }

  //Check openned file
  if(pfile == nullptr)
    return SHARED_FILE_UNABLE_TO_OPEN_FILE;

  //Get origin position
  if(fgetpos(pfile,&origin) != 0){
    trusted_close();
    return SHARED_FILE_FGETPOS_ERROR;
  }

  //Set position to origin for all readers
  trusted_rewind();

  //Save format information
  binary = binaryFile;
  
  return SHARED_FILE_SUCCESS;
}

int pen_sharedFile::read(const unsigned ID,
			 char* data,
			 const size_t dataSize,
			 size_t &nread){

  //Lock file
  std::lock_guard<std::mutex> guard(fileLock);

  if(pfile == nullptr){
    return SHARED_FILE_NO_FILE_OPENNED;
  }
  
  //Check if this ID is being used
  int pos = trusted_findID(ID);  
  if(pos < 0)
    return SHARED_FILE_UNUSED_READER_ID;

  //Set file to readers position
  if(fsetpos(pfile,&(readers[pos].second)) != 0)
    return SHARED_FILE_FSETPOS_ERROR;
  
  //Read data
  if(binary)
    nread = fread(data,sizeof(char),dataSize,pfile);
  else{
    fgets(data,dataSize,pfile);
    nread = strlen(data);
  }
  
  //Store new position
  if(fgetpos(pfile,&(readers[pos].second)) != 0){
    return SHARED_FILE_FGETPOS_ERROR;
  }
  
  return SHARED_FILE_SUCCESS;
}
int pen_sharedFile::read(const unsigned ID,
			 unsigned char* data,
			 const size_t dataSize,
			 size_t &nread){

  //Lock file
  std::lock_guard<std::mutex> guard(fileLock);

  if(pfile == nullptr){
    return SHARED_FILE_NO_FILE_OPENNED;
  }
  
  //Check if this ID is being used
  int pos = trusted_findID(ID);  
  if(pos < 0)
    return SHARED_FILE_UNUSED_READER_ID;


  //Set file to readers position
  if(fsetpos(pfile,&(readers[pos].second)) != 0)
    return SHARED_FILE_FSETPOS_ERROR;
  
  //Read data
  if(binary)
    nread = fread(data,sizeof(char),dataSize,pfile);
  else{
    fgets((char*)data,dataSize,pfile);
    nread = strlen((char*)data);
  }
  
  //Store new position
  if(fgetpos(pfile,&(readers[pos].second)) != 0){
    return SHARED_FILE_FGETPOS_ERROR;
  }
  
  return SHARED_FILE_SUCCESS;
}



void pen_sharedFile::clear(){

  //Lock file
  std::lock_guard<std::mutex> guard(fileLock);

  //Close file
  trusted_close();
  
  //Clear vector
  readers.clear();
  
}

pen_sharedFile::~pen_sharedFile(){
  //Clear 
  clear();
}
