
//
//
//    Copyright (C) 2019-2021 Universitat de València - UV
//    Copyright (C) 2019-2021 Universitat Politècnica de València - UPV
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

 
#ifndef __PEN_SPLITTED_FILE__
#define __PEN_SPLITTED_FILE__

#include <cstdio>
#include <vector>
#include <mutex>
#include <algorithm>
#include <cstring>
#include <string>

enum splittedFile_states{
			 SPLITTED_FILE_SUCCESS = 0,
			 SPLITTED_FILE_USED_PARTITION_ID,
			 SPLITTED_FILE_UNABLE_TO_OPEN,
			 SPLITTED_FILE_UNUSED_ID
};

class pen_splittedFile{

private:
  //Files vector
  std::vector<std::pair<unsigned,FILE*>> partitions;
  //Locks for write operations
  static const int nPL = 15;
  std::mutex partLocks[nPL];

  char formatW[5];
  char formatA[5];
  char formatR[5];

  void lockPartitions();
  void unlockPartitions();

  void trustedFlush();
  
public:

  const std::string prefix;
  const bool binary;
  
  std::string separator;
  
  pen_splittedFile(const char* prefixIn, bool binaryIn) : prefix(prefixIn),
							  binary(binaryIn)
  {
    formatW[0] = 'w';
    formatA[0] = 'a';
    formatR[0] = 'r';
    if(binary){

      formatW[1] = 'b';
      formatA[1] = 'b';
      formatR[1] = 'b';

      formatW[2] = '\0';
      formatA[2] = '\0';
      formatR[2] = '\0';
      
    }
    else{
      formatW[1] = '\0';
      formatA[1] = '\0';
      formatR[1] = '\0';
    }
  }
  
  int createPartition(unsigned ID);
  int write(const unsigned ID,
	    const char* data,
	    const size_t dataSize);
  int write(const unsigned ID,
	    const unsigned char* data,
	    const size_t dataSize);
  
  int flush(const unsigned ID);
  int merge();
  void clear();
  ~pen_splittedFile();
};


#endif
