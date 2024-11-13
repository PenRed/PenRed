
//
//
//    Copyright (C) 2019-2024 Universitat de València - UV
//    Copyright (C) 2019-2024 Universitat Politècnica de València - UPV
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

 
#ifndef __PENRED_INTERNAL_DATA_CLASSES__
#define __PENRED_INTERNAL_DATA_CLASSES__

#include "pen_parser.hh"
#include "pen_reader.hh"
#include <sstream>
#include <cstdarg>

namespace penred{
  namespace logs{

    enum logType{
      CONFIGURATION = 0,
      SIMULATION
    };

    class logger;
    class loggerStream;

    class loggerFile{
      friend class logger;

      FILE* flog;

      inline void close(){
	if(flog != nullptr){
	  fclose(flog);
	  flog = nullptr;
	}
      }

      inline void open(const char* filename,
		       bool append = false) {
	close();
	flog = fopen(filename, append ? "a" : "w");
      }

      inline bool isOpen() const { return flog != nullptr; }

      inline FILE* getFile() { return flog; }
      inline const FILE* readFile() const { return flog; }

    public:
      loggerFile() : flog(nullptr) {}
      
      ~loggerFile(){
	close();
      }
    };
    
    class logger{

      class loggerStream{
	friend class logger;
      private:
	const logger& l;
      public:

	constexpr loggerStream(logger& lin) : l(lin){}
      
	template <class T>
	const loggerStream& operator<<(const T& o) const {
	  //Save the object in a stream
	  std::stringstream stream;
	  stream << o;
	  //Print it
	  l.printf("%s", stream.str().c_str());
	  return *this;
	}

	//Handle endl
	const loggerStream& operator<<(std::ostream&(*)(std::ostream&)) const {
	  l.printf("\n");
	  flush();
	  return *this;
	}

	void flush() const{
	  l.flush();
	}
      };
      
    public:
      static constexpr const size_t nlogs = 10;
      loggerStream cout;
    private:
  
      size_t defaultLog;
      static std::vector<loggerFile> _flog;
  
    public:

      logger() : cout(*this), defaultLog(CONFIGURATION) {}
      logger(size_t ilog) : cout(*this), defaultLog(ilog) {}

      static inline int printf(const size_t ilog, const char* format, ...){
	int result;
	va_list args;
	va_start(args, format);
	if(_flog[ilog].isOpen()){
	  result = vfprintf(_flog[ilog].getFile(), format, args);      
	  fflush(_flog[ilog].getFile());
	}else{
	  result = vprintf(format, args);
	  fflush(stdout);
	}
	va_end(args);
	return result;
      }
  
      inline int printf(const char* format, ...) const {
	int result;
	va_list args;
	va_start(args, format);
	if(_flog[defaultLog].isOpen()){
	  result = vfprintf(_flog[defaultLog].getFile(), format, args);      
	  fflush(_flog[defaultLog].getFile());
	}else{
	  result = vprintf(format, args);
	  fflush(stdout);
	}
	va_end(args);
	return result;
      }

      static inline void setLogFile(const size_t ilog,
				    const char* filename,
				    const bool append = false) {
	if(!_flog[ilog].isOpen()){
	  _flog[ilog].open(filename, append ? "a" : "w");
	}
      }

      static inline void setConfigurationLogFile(const char* filename,
						 const bool append = false) {
	setLogFile(CONFIGURATION,filename,append);
      }
      static inline void setSimulationLogFile(const char* filename,
					      const bool append = false) {
	setLogFile(SIMULATION,filename,append);
      }
      
      static inline void closeLogFile(const size_t ilog) {
	_flog[ilog].close();
      }
      static inline void closeAndSetLogFile(const size_t ilog,
					    const char* filename,
					    const bool append = false) {
	_flog[ilog].close();
	setLogFile(ilog, filename, append);
      }

      inline void setDefaultLog(const size_t ilog){
	if(ilog < nlogs)
	  defaultLog = ilog;
	else
	  defaultLog = 0;
      }

      static void flush(const size_t ilog){
	if(_flog[ilog].isOpen()){
	  fflush(_flog[ilog].getFile());
	}else{
	  fflush(stdout);
	}
      }

      void flush() const{
	flush(defaultLog);
      }

      ~logger(){
    
      }
    };
    
  } //namespace logs
} //namespace penred
  
#endif
