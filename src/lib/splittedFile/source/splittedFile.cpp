
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

 
#include "splittedFile.hh"

int pen_splittedFile::createPartition(unsigned ID){

  //Lock partitions
  lockPartitions();

  //Check if this ID is being used
  for(size_t i = 0; i < partitions.size(); i++){
    if(partitions[i].first == ID){
      //Unlock partitions
      unlockPartitions();
      return SPLITTED_FILE_USED_PARTITION_ID;
    }
  }

  //Try to open a file for the new partition
  FILE* fpart = nullptr;
  std::string filename = prefix + std::string("-") + std::to_string(ID);
  fpart = fopen(filename.c_str(),formatW);

  if(fpart == nullptr){
    //Unlock partitions
    unlockPartitions();
    return SPLITTED_FILE_UNABLE_TO_OPEN;
  }
    
  //Store the file pointer and its ID
  partitions.push_back(std::make_pair(ID,fpart));

  //Unlock partitions
  unlockPartitions();

  return SPLITTED_FILE_SUCCESS;
}

int pen_splittedFile::write(const unsigned ID,
			    const char* data,
			    const size_t dataSize){

  //Lock the corresponding write lock
  unsigned lockID = ID % nPL;
  std::lock_guard<std::mutex> guard(partLocks[lockID]);

  //Check if this ID is being used
  int pos = -1;
  for(size_t i = 0; i < partitions.size(); i++){
    if(partitions[i].first == ID){
      pos = i;
      break;
    }
  }
  
  if(pos < 0)
    return SPLITTED_FILE_UNUSED_ID;

  //Write data
  if(binary)
    fwrite(data,sizeof(char),dataSize,partitions[pos].second);
  else
    fprintf(partitions[pos].second,"%s",data);
  

  return SPLITTED_FILE_SUCCESS;
}
int pen_splittedFile::write(const unsigned ID,
			    const unsigned char* data,
			    const size_t dataSize){

  //Lock the corresponding write lock
  unsigned lockID = ID % nPL;
  std::lock_guard<std::mutex> guard(partLocks[lockID]);

  //Check if this ID is being used
  int pos = -1;
  for(size_t i = 0; i < partitions.size(); i++){
    if(partitions[i].first == ID){
      pos = i;
      break;
    }
  }
  
  if(pos < 0)
    return SPLITTED_FILE_UNUSED_ID;

  //Write data
  if(binary)
    fwrite(data,sizeof(char),dataSize,partitions[pos].second);
  else
    fprintf(partitions[pos].second,"%s",data);
  

  return SPLITTED_FILE_SUCCESS;
}

int pen_splittedFile::flush(const unsigned ID){

  //Lock the corresponding write lock
  unsigned lockID = ID % nPL;
  std::lock_guard<std::mutex> guard(partLocks[lockID]);

  //Check if this ID is being used
  int pos = -1;
  for(size_t i = 0; i < partitions.size(); i++){
    if(partitions[i].first == ID){
      pos = i;
      break;
    }
  }

  if(pos < 0)
    return SPLITTED_FILE_UNUSED_ID;

  //Flush data
  fflush(partitions[pos].second);
  
  return SPLITTED_FILE_SUCCESS;
}

void pen_splittedFile::trustedFlush(){

  //Flush all partition files
  
  for(unsigned i = 0; i < partitions.size(); i++){
    //Flush data
    fflush(partitions[i].second);    
  }
}

int pen_splittedFile::merge(){

  //Lock partitions
  lockPartitions();

  //Flush partitions
  trustedFlush();
  
  //Create merge file
  FILE* fmerge = nullptr;
  std::string mergeFName = prefix + std::string("-merged.dat");
  fmerge = fopen(mergeFName.c_str(),formatW);
  if(fmerge == nullptr){
    //Unlock partitions
    unlockPartitions();
    return SPLITTED_FILE_UNABLE_TO_OPEN;
  }
  //Sort partitions by ID
  sort(partitions.begin(), partitions.end());
  
  //Merge sorted partitions
  for(unsigned i = 0; i < partitions.size(); i++){

    //Close file to read it
    fclose(partitions[i].second);

    //Open file to merge
    FILE* fin = nullptr;
    std::string inFName = prefix +
      std::string("-") +
      std::to_string(partitions[i].first);

    fin = fopen(inFName.c_str(),formatR);

    //Read-write
    const unsigned bufferSize = 50000;
    char buffer[bufferSize];
    
    if(binary){
      size_t nread = fread(buffer,sizeof(char),bufferSize,fin);
      while(nread > 0){
	fwrite(buffer,sizeof(char),nread,fmerge);
	nread = fread(buffer,sizeof(char),bufferSize,fin);
      }
    }
    else{
      while(fgets(buffer,bufferSize,fin) != nullptr){
	fprintf(fmerge,"%s",buffer);
      }
    }

    //Close merged file
    fclose(fin);

    //Reoppen file and set cursor
    partitions[i].second = nullptr;
    partitions[i].second = fopen(inFName.c_str(),formatA);
    if(partitions[i].second == nullptr){
      fclose(fmerge);
      //Unlock partitions
      unlockPartitions();
      return SPLITTED_FILE_UNABLE_TO_OPEN;
    }
  }

  //flush merged file
  fflush(fmerge);
  
  //Close merged file
  fclose(fmerge);
  
  //Unlock partitions
  unlockPartitions();

  return SPLITTED_FILE_SUCCESS;
}

void pen_splittedFile::clear(){

  //Lock partitions
  lockPartitions();
  
  size_t nParts = partitions.size();
  for(size_t i = 0; i < nParts; i++){
    fclose(partitions[i].second);
    std::string filename = prefix +
      std::string("-") +
      std::to_string(partitions[i].first);
    remove(filename.c_str());
  }
  partitions.clear();
  
  //Unlock partitions
  unlockPartitions();
}

void pen_splittedFile::lockPartitions(){
  for(unsigned i = 0; i < nPL; i++)
    partLocks[i].lock();
}
void pen_splittedFile::unlockPartitions(){
  for(unsigned i = 0; i < nPL; i++)
    partLocks[i].unlock();  
}

pen_splittedFile::~pen_splittedFile(){
  //Merge files on destroy
  merge();
  //Clear partitions
  clear();
}
