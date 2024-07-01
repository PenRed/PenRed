
//
//
//    Copyright (C) 2024 Universitat de València - UV
//    Copyright (C) 2024 Universitat Politècnica de València - UPV
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
//        vicent.gimenez.alventosa@gmail.com (Vicent Giménez Alventosa)
//
//

#include "dataBasesCommon.hh"

namespace penred{

  namespace dataBases{

    int createStringLiteralFiles(std::istream& in,
				 const std::string& filename,
				 unsigned& subFiles,
				 const size_t maxBytes){

      //Store number of subfiles and stored bytes in current file
      unsigned nSubFile = 0;
      size_t bytesInSubFile = 0;
      
      //Open the output file 0
      std::string pathToOut(filename);
      pathToOut += "_" + std::to_string(nSubFile);
      std::ofstream fout(pathToOut);
      if(!fout.is_open()){
	printf("Unable to open output file '%s'\n", pathToOut.c_str());
	return errors::UNABLE_TO_OPEN_FILE;
      }

      //File opened, start processing buffer
      fout << "R\"***(";

      std::string line;
      while(std::getline(in,line)){

	//Check line size
	if(line.length() + bytesInSubFile > maxBytes){
	  if(bytesInSubFile == 0){
	    //Unable to store a single line
	    return errors::LINE_BIGER_THAN_CHUNK_SIZE;
	  }

	  //Maximum file size reached

	  //Restart used bytes
	  bytesInSubFile = 0;

	  //Close current file and open a new subfile
	  fout << ")***\"" << std::endl;
	  fout.close();
	  fout.clear();

	  //Increase sub file number
	  ++nSubFile;

	  //Create file path
	  pathToOut = filename + "_" + std::to_string(nSubFile);
	  fout.open(pathToOut);
	  if(!fout.is_open()){
	    printf("Unable to open output file '%s'\n", pathToOut.c_str());
	    return errors::UNABLE_TO_OPEN_FILE;
	  }

	  //Write string literal header
	  fout << "R\"***(";	  
	}

	//Add bytes in this line
	bytesInSubFile += line.length();

	//Write line
	fout << line << std::endl;
      }
      fout << ")***\"" << std::endl;

      //Close file
      fout.close();

      //Save number of subfiles
      subFiles = nSubFile+1;

      return errors::SUCCESS;
    }

    
  } //namespace dataBases
} //namespace penred
