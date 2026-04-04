//
//
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

#include <fstream>
#include "functions.hh"

namespace py = pybind11;

PYBIND11_MODULE(data,m){

  m.doc() = "penred data module";

  m.def("dict2SectionString", &dict2SectionString,
	py::arg("conf"),	
	R"(

Converts a dictionary to a compatible penRed configuration section string.

Args:
    conf (dict) : configuration dictionary to be converted.
    
Returns:
    String containing the converted dictionary

    )");

  m.def("configFile2YAML",
	[](const std::string& filename) -> py::str{
	  if(filename.empty())
	    return std::string("");
	  else{

	    //Parse configuration file
	    pen_parserSection config;
	    std::string errorLine;
	    unsigned long errorLineNum;
	    int err = parseFile(filename.c_str(),config,errorLine,errorLineNum);
	    
	    if(err != INTDATA_SUCCESS){
	      printf("Error parsing configuration.\n");
	      printf("Error code: %d\n",err);
	      printf("Error message: %s\n",pen_parserError(err));
	      printf("Error located at line %lu, at text: %s\n",
		     errorLineNum,errorLine.c_str());
	      return std::string("");
	    }

	    //Create YAML string
	    return config.stringifyYAML();
	  }
	},
	R"(

Reads a configuration file and returns a YAML string with the read information.

Args:
    filename (str) : File to be read.
    
Returns:
    String containing the information in YAML format

    )");

  m.def("readResults", [](const std::string& filename, const bool extractInfo) -> py::tuple{

  //Create a results structure with maximum dimensions
  penred::measurements::results<double, penred::measurements::maxDims> reader;
  
  //Read data
  std::ifstream fin(filename, std::ifstream::in);
  if(!fin){
    std::string errorMsg("Unable to open file ");
    errorMsg += filename;
    throw py::value_error(errorMsg.c_str());
  }
  
  int err = reader.read(fin);
  fin.close();
  if(err != 0){
    std::string errorMsg("Error reading data file. ");
    errorMsg += penred::measurements::errorToString(err);
    throw py::value_error(errorMsg.c_str());
  }

  //Convert results to numpy arrays and return them
  return result2numpy(reader, extractInfo, true);

  },
    py::arg("filename"),
    py::arg("extract_info") = false,
	R"(

Reads a penRed standard results file and convert it to numpy arrays.

Args:
    filename (str) : Results filename.
    extract_info (bool) : If enabled, the information of each dimension will be read as well. Defaults to false.
    
Returns:
    A tuple containing the data, uncertainties and, if enabled, information of each dimension (lower grid limit, upper grid limits, dimension header), value header, file description

    )");

  m.def("printResults", [](const std::string& filename,
			   const py::array_t<double>& values,
			   const py::array_t<double>& sigma,
			   const py::list& info,
			   const std::string& valueHeader,
			   const unsigned nSigma,
			   const bool printCoordinates,
			   const bool printBinNumber) -> void{
    
    //Check dimensions
    std::vector<unsigned long> dimBins;
    const size_t nDim = assertShapes(values, sigma, dimBins);

    //Ensure the information is a list
    if (!py::isinstance<py::list>(info) && !py::isinstance<py::tuple>(info)) {
      throw py::type_error("Argument 'info' must be a list");
    }

    //Check information length
    if(info.size() != nDim){
      throw py::value_error("Data and information dimensions mismatch");      
    }    

    //Check and extract information
    std::vector<std::pair<double, double>> limits;
    std::vector<std::string> headers;
    for(size_t i = 0; i < info.size(); ++i){
      //Check if the element is a tuple
      if (!py::isinstance<py::tuple>(info[i])) {
	throw py::type_error("Information element " + std::to_string(i) + " is not a tuple");
      }

      py::tuple t = info[i];
      if(t.size() != 3){
	throw py::value_error("Tuple " + std::to_string(i) + " has " + 
			      std::to_string(t.size()) + " elements, expected 3");
      }

      // Check first element is numeric
      if (!py::isinstance<py::int_>(t[0]) && !py::isinstance<py::float_>(t[0])) {
	throw py::type_error("Tuple " + std::to_string(i) + 
			     ", first element must be numeric (int or float)");
      }
        
      // Check second element is numeric
      if (!py::isinstance<py::int_>(t[1]) && !py::isinstance<py::float_>(t[1])) {
	throw py::type_error("Tuple " + std::to_string(i) + 
			     ", second element must be numeric (int or float)");
      }

      // Check third element is string
      if (!py::isinstance<py::str>(t[2])) {
	throw py::type_error("Tuple " + std::to_string(i) + 
			     ", third element must be a string");
      }

      // Extract values
      double min = t[0].cast<double>();
      double max = t[1].cast<double>();
      std::string header = t[2].cast<std::string>();

      limits.emplace_back(min,max);
      headers.push_back(std::move(header));
    }
    
    //Create a results structure with maximum dimensions
    penred::measurements::results<double, penred::measurements::maxDims> results;

    //Reverse bins per dimension and limits to correct numpy ordering
    std::reverse(dimBins.begin(), dimBins.end());
    std::reverse(limits.begin(), limits.end());
    std::reverse(headers.begin(), headers.end());

    //Init it
    int err = results.init(dimBins, limits, flattenArray(values), flattenArray(sigma));
    if(err != penred::measurements::errors::SUCCESS){
      throw py::value_error(penred::measurements::errorToString(err));
    }

    //Set headers
    for(size_t i = 0; i < headers.size(); ++i){
      results.setDimHeader(i,headers[i]);
    }
    results.setValueHeader(valueHeader);

    //Open output file
    FILE* fout = fopen(filename.c_str(),"w");
    if(fout == nullptr){
      std::string errorMsg("Unable to open file ");
      errorMsg += filename;
      throw std::runtime_error(errorMsg.c_str());
    }

    //Print data
    results.print(fout, nSigma, printCoordinates, printBinNumber, true);
    fclose(fout);
  },
	py::arg("filename"),
	py::arg("values"),
	py::arg("sigma"),
	py::arg("info"),
	py::arg("value_header") = "Value",
	py::arg("print_sigmas") = 2,
	py::arg("print_coordinates") = true,
	py::arg("print_bins") = true,
	R"(

Reads a penRed standard results file and convert it to numpy arrays.

Args:
    filename (str) : Results filename.
    values (vector) : Numpy vector with the values to be saved
    sigma (vector) : Numpy vector with the uncertainty of each value corresponding to one sigma
    info (list): List of tuples where each element contains the information of the corresponding dimension. Each element must use the following format: (min, max, header), where *min* and *max* are the minimum and maximum grid value for this dimension, respectively, and *header* the text header for dimension's column.
    print_sigmas (int): Specify the number of printed sigmas, i.e. the uncertainties column will be transformed as print_sigmas*sigma
    print_coordinates (bool): If enabled, the coordinates of each dimension will be printed
    print_bins (bool): If enabled, the bin number of each dimension will be printed

Returns:
    None

Raises:
    TypeError: Incompatible types have been provided.
    ValueError: Dimensions mismatch or incompatible values have been provided.
    RuntimeError: Unable to open output file
    )");  
}
