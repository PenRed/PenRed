//
//
//    Copyright (C) 2024 Universitat de València - UV
//    Copyright (C) 2024 Universitat Politècnica de València - UPV
//    Copyright (C) 2025-2026 Vicent Giménez Alventosa
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
//        sanolgi@upvnet.upv.es (Sandra Oliver Gil)
//    
//

#include "functions.hh"

namespace py = pybind11;

std::string dict2SectionStringWithPrefix(const py::dict& dict, const std::string& prefixIn){

  //Convert a python dictionary to a string compatible with pen_parserSection
    
  std::string result;
  std::string prefix = prefixIn;
  if(!prefix.empty()){
    if(prefix.back() != '/'){
      //Append a slash
      prefix.append(1,'/');
    }
  }
    
  for(auto it : dict){
    //Check data type

    const std::string key = py::str(it.first).cast<std::string>();

    if(py::isinstance<py::bool_>(it.second)){
      //Boolean case
      if(it.second.cast<bool>())
	result += prefix + key + " true\n";
      else
	result += prefix + key + " false\n";
    }    
    else if(py::isinstance<py::int_>(it.second)){
      //Number case
      result += prefix + key + " " + std::to_string(it.second.cast<int>()) + "\n";
    }
    else if(py::isinstance<py::float_>(it.second)){
      //Float  case
      char aux[20];
      snprintf(aux, 20, " %15.5E\n", it.second.cast<double>());
      result += prefix + key + aux;
    }
    else if(py::isinstance<py::tuple>(it.second) ||
	    py::isinstance<py::list>(it.second)){
      //Array case
      result += prefix + key + " [";
      //Print each array element
      bool first = true;
      for(auto e : it.second){
	//Ensure the element is a number
	if(!py::isinstance<py::int_>(e) &&
	   !py::isinstance<py::float_>(e)){
	  printf("dict2SectionString: Error: The array at '%s' contains "
		 "a non numeric element. Will be skipped.\n",
		 (prefix + key).c_str());
	  continue;
	}
	if(py::isinstance<py::bool_>(e)){
	  if(it.second.cast<bool>())
	    result += (first ? " true" : ", true");
	  else
	    result += (first ? " false" : ", false");
	}
	else if(py::isinstance<py::int_>(e))
	  result += (first ? " " : ", ") + std::to_string(e.cast<int>());
	else{
	  char aux[25];
	  snprintf(aux, 25, "%s%15.5E", (first ? " " : ", "), e.cast<double>());	  
	  result += aux;
	}
	first = false;
      }
      result += " ]\n";
    }
    else if(py::isinstance<py::str>(it.second)){      
      //String case
      result += prefix + key + " \"" + it.second.cast<std::string>() + "\"\n";
    }
    else if(py::isinstance<py::dict>(it.second)){
      //Dictionary case. Add the prefix and parse it
      result += dict2SectionStringWithPrefix(it.second.cast<py::dict>(), prefix + key);
    }
    else{
      printf("dict2SectionString: Error: incompatible element at '%s'. "
	     "Only numbers, boolean, string, list, tuples and "
	     "dictionaries are allowed. Will be skipped.\n",
	     (prefix + key).c_str());
    }
  }
  return result;
}

std::string dict2SectionString(const py::dict& dict){
  return dict2SectionStringWithPrefix(dict, "");
}

size_t assertShapes(const py::array_t<double>& arr1,
		    const py::array_t<double>& arr2,
		    std::vector<unsigned long>& dimsSizes) {
  
  // Get shape information
  py::buffer_info buf1 = arr1.request();
  py::buffer_info buf2 = arr2.request();
    
  // Compare number of dimensions
  if (buf1.ndim != buf2.ndim) {
    std::string errorMsg("Different number of dimensions: ");
    errorMsg += std::to_string(buf1.ndim);
    errorMsg += " vs ";
    errorMsg += std::to_string(buf2.ndim);
    throw py::buffer_error(errorMsg);
  }
    
  // Compare shape
  for (long int i = 0; i < buf1.ndim; ++i) {
    if (buf1.shape[i] != buf2.shape[i] || buf1.shape[i] <= 0 ) {
      std::string errorMsg("Dimension ");
      errorMsg += std::to_string(i);
      errorMsg += "mismatch: ";
      errorMsg += std::to_string(buf1.shape[i]);
      errorMsg += " vs ";
      errorMsg += std::to_string(buf2.shape[i]);
      throw py::buffer_error(errorMsg);	  
      std::cout << "Dimension " << i << " mismatch: " 
		<< buf1.shape[i] << " vs " << buf2.shape[i] << std::endl;
    }
    dimsSizes.push_back(static_cast<unsigned long>(buf1.shape[i]));
  }

  return buf1.ndim;
}
