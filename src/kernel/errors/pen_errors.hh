
//
//
//    Copyright (C) 2019 Universitat de València - UV
//    Copyright (C) 2019 Universitat Politècnica de València - UPV
//    Copyright (C) 2026 Vicent Giménez Alventosa
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



#ifndef __PEN_ERRORS_
#define __PEN_ERRORS_

#include <memory>
#include <functional>
#include <string>
#include "eenum.h"  //Enumeration of errors

void penError(pen_errCode err);
pen_errCode penGetError();
void penClearError();
const char* penErrorString(const pen_errCode err);

namespace penred{

  namespace errors{

    //Define template error messages. This function must be overloaded
    //to support type-specific error message reading
    template<class Type>
    constexpr const char* errorMessage(const int) noexcept {
      return "Unknown type and error";
    }
    
    struct Error{

    private:
      std::shared_ptr<Error> trace;

    protected:
      std::function<const char*(int)> errorMessageFunc;
      const char* specificErrorMessage() const{
	return errorMessageFunc ? errorMessageFunc(code) : "Error message function not set";
      }
      
    public:
      int code;
      std::string description;

      Error() : code(0) {}
      
      explicit Error(const std::function<const char*(int)>& errMsgFunc) : errorMessageFunc(errMsgFunc),
									  code(0)
      {}

      virtual inline std::shared_ptr<Error> clone() const noexcept{
	return std::make_shared<Error>(*this);
      }
      
      inline const char* codeMessage() const {
	return specificErrorMessage();
      }

      inline void setTrace(const Error& err) noexcept{
	trace = err.clone();
      }
      inline std::string stringify(const std::string& indent = "") const noexcept{
	std::string result = indent + std::string("- Error message: ") + codeMessage();
	
	if(!description.empty()){
	  result += ":\n" + indent + "    Description:\n" + indent + "        " + description;
	}
	result += "\n";

	if(trace)
	  result += trace->stringify(indent + "  ");
	return result;
      }

      inline explicit operator bool() const noexcept {
        return code != 0;  // Error exists if code is non-zero
      }

      virtual ~Error() = default;
    };

    //Error comparison operators
    constexpr bool operator==(const Error& lhs, const int rhs) noexcept
    {
      return lhs.code == rhs;
    }
 
    constexpr bool operator!=(const Error& lhs, const int rhs) noexcept
    {
      return !(lhs.code == rhs);
    }    

    template<class Type>
    struct SpecificError : Error{
      
    public:
      SpecificError() : Error(errorMessage<Type>) {}
      
      inline std::shared_ptr<Error> clone() const noexcept{
	return std::make_shared<SpecificError<Type>>(*this);
      }
    };
    
  } //namespace errors
  
} //namespace penred

#endif
