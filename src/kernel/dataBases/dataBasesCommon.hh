
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

#ifndef _PENRED_COMMONS_DATABASES_
#define _PENRED_COMMONS_DATABASES_

#ifndef PEN_MAXDB_CHUNK_SIZE
#define PEN_MAXDB_CHUNK_SIZE 60000
#endif

#include <sstream>
#include <string>
#include <array>
#include <vector>
#include <functional>
#include <limits>
#include <fstream>

namespace penred{

  struct massFraction{
    unsigned Z;
    double fraction;

    massFraction() = default;
    constexpr massFraction(const unsigned Zin) : Z(Zin), fraction(1.0) {}
    constexpr massFraction(const unsigned Zin,
			   const double fractionIn) : Z(Zin),
						   fraction(fractionIn){}
  };

  namespace dataBases{

    constexpr size_t maxChunkBytes = PEN_MAXDB_CHUNK_SIZE;

    namespace errors{
      enum errors{
	SUCCESS = 0,
	UNABLE_TO_OPEN_FILE,
	LINE_BIGER_THAN_CHUNK_SIZE,
      };
    };

    struct materials{

      virtual std::vector<std::string> matList() const = 0;
      virtual unsigned getIndex(const std::string&) const = 0;
      virtual std::vector<massFraction> getElements(const unsigned) const = 0;
      virtual double getDensity(const unsigned) const = 0; //in g/cm**3

      inline std::vector<massFraction> getElements(const std::string& name) const{
	return getElements(getIndex(name));
      }
      inline double getDensity(const std::string& name) const{
	return getDensity(getIndex(name));
      }
    };

    constexpr const char* errorMessage(unsigned ierr){
      
      switch(ierr){
      case(errors::SUCCESS): return "success";
      case(errors::UNABLE_TO_OPEN_FILE): return "unable to open file";
      case(errors::LINE_BIGER_THAN_CHUNK_SIZE): return "line biger than maximum chunk size";
      default: return "unknown error";
      };      
    }

    int createStringLiteralFiles(std::istream& in,
				 const std::string& filename,
				 unsigned& subFiles,
				 const size_t maxBytes = maxChunkBytes);

    inline std::string toVariableName( const std::string& nameIn ){
      //Create the variable name with no points
      std::string varName = nameIn;
      std::string::size_type pointPos;
      pointPos = varName.find(".") ;
      while(pointPos != std::string::npos){
	varName.replace(pointPos,1,"_");
	pointPos = varName.find(".") ;
      }
      return varName;
    }

    class literalArrayStream{

      //This class handles a input stream using several literal strings to read from, via a
      //stringstream class. When the first one ends, the second literal string is used, and so on.
      //It is supposed that each literal string ends with a complete line, i.e. lines are not
      //break between successive literal strings.

      //To get the appropiate literal string, a function must be provided to the constructor,
      //which takes a index to retrieve the appropiate data chunk. This allows to split
      //huge data base files in many chunks. This function can be changed with the "str" method.
      
    private:
      std::stringstream stream;
      size_t nextString;
      std::function<const char* (unsigned)> getData;

      inline bool handleEof(){
	//If eof is reached, change to the next data chunk.
	//Rretuns true on data change, false otherwise
	if(stream.eof()){
	  const char* const nextData = getData(nextString);
	  if(nextData != nullptr){
	    stream.clear();
	    stream.str(nextData);
	    nextString += 1;
	    return true;
	  }
	}
	return false;
      }
    public:
      
      inline literalArrayStream(const std::function<const char* (unsigned)>& getDataIn) :
	nextString(1),
	getData(getDataIn){
	const char* const firstData = getData(0);
	if(firstData == nullptr)
	  stream.str(std::string());
	else
	  stream.str(firstData);
      }

      inline unsigned partition() const {
	if(nextString > 0)
	  return nextString-1;
	return 0;
      }
      
      template<class T>
      inline literalArrayStream& operator>>(T& obj){

	stream >> obj;
	if(handleEof()){
	  stream >> obj;
	}
	
	return *this;
      }

      inline literalArrayStream& ignore(std::streamsize count = 1){
	stream.ignore(count);
	//No eof check, the whole ignore should be in the same line
	return *this;
      }
      inline literalArrayStream& ignoreLine(){
	stream.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
	if(handleEof()){
	  //The line has not been ignored, skip it
	  stream.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
	}
	return *this;
      }

      //Define getline for literalArrayStream
      inline literalArrayStream& getline(std::string& str){
	std::getline(stream, str);
	if(handleEof()){
	  std::getline(stream, str);
	}
	return *this;
      }

      inline void str(const std::function<const char* (unsigned)>& getDataIn){
	getData = getDataIn;
	nextString = 1;
	stream.clear();
	const char* const firstData = getData(0);
	if(firstData == nullptr)
	  stream.str(std::string());
	else
	  stream.str(firstData);	
      }
      
      inline bool operator!() const{
	return !stream;
      }

      inline explicit operator bool() const{
	return bool(stream);
      }

      inline bool eof(){
	if(stream.eof()){
	  if(handleEof()){
	    return false;
	  }
	  return true;
	}
	return false;
      }
    };
    
  } // namespace dataBases
  
} // namespace penred

#endif
