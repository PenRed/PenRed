
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

 
#ifndef __PEN_SHARED_FILE__
#define __PEN_SHARED_FILE__

#include <cstdio>
#include <vector>
#include <mutex>
#include <algorithm>
#include <cstring>

enum sharedFile_states{
		       SHARED_FILE_SUCCESS = 0,
		       SHARED_FILE_USED_READER_ID,
		       SHARED_FILE_UNUSED_READER_ID,
		       SHARED_FILE_UNABLE_TO_OPEN_FILE,
		       SHARED_FILE_NO_FILE_OPENNED,
		       SHARED_FILE_FSETPOS_ERROR,
		       SHARED_FILE_FGETPOS_ERROR,
		       SHARED_FILE_FSEEK_ERROR
};

class pen_sharedFile{

private:
  //Files vector
  std::vector<std::pair<unsigned,fpos_t>> readers;
  //Shared file pointer
  FILE* pfile;
  //Locks for read/write operations
  fpos_t origin;
  std::mutex fileLock;
  //Save if format is binary or ASCII
  bool binary;

  //Save file size
  size_t fsize;
  
  inline void trusted_rewind(){
    //Set position to origin for all readers
    for(unsigned i = 0; i < readers.size(); i++){
      readers[i].second = origin;
    }    
  }

  inline void trusted_close(){
    
    if(pfile != nullptr){
      fclose(pfile);
      pfile = nullptr;
    }
    fsize = 0;
  }

  inline int trusted_findID(const unsigned ID){
    //Check if this ID is being used
    int pos = -1;
    for(size_t i = 0; i < readers.size(); i++){
      if(readers[i].first == ID){
	pos = i;
	break;
      }
    }
    return pos;
  }
  
public:
  
  pen_sharedFile() : pfile(nullptr),fsize(0)
  {}

  bool isBinary();

  size_t size();
  
  int createReader(const unsigned ID);

  int tell(const unsigned ID,
	    long int &pos);

  int seek(const unsigned ID,
	    long int offset,
	    int origin,
	    int* fseekerr = nullptr);
  
  int getPos(const unsigned ID,
	     fpos_t& fpos);

  int setPos(const unsigned ID,
	     const fpos_t& fpos);
  
  int open(const char* filename,
	   const bool binaryFile);

  inline void close(){

    //Lock file
    std::lock_guard<std::mutex> guard(fileLock);
    trusted_close();
    
  }
  int read(const unsigned ID,
	   char* data,
	   const size_t dataSize,
	   size_t &nread);
  int read(const unsigned ID,
	   unsigned char* data,
	   const size_t dataSize,
	   size_t &nread);
  
  void clear();
  ~pen_sharedFile();
};


#endif
