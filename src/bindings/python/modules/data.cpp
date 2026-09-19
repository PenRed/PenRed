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

#include "functions.hh"

namespace py = pybind11;

// ** Template helpers to register results-based classes for different dimensions

// * Dimension-specific methods adders

// No special methods
template <typename T, std::size_t D, typename PyClass>
struct specialMethods{
  static void add(PyClass&) {}
};


// 1D case
template <typename T, typename PyClass>
struct specialMethods<T, 1, PyClass>{
  static void add(PyClass& cls){
    cls.def("extractValue1D",
            [](const penred::measurements::results<T, 1>& obj,
               const double x) -> py::tuple{
              
              T value;
              double uncertainty;

              int err = obj.extractValue(x, value, uncertainty);
              if(err != penred::measurements::errors::SUCCESS) {
                std::string errorMsg("Error extracting interpolated value from results. ");
                errorMsg += penred::measurements::errorToString(err);
                throw pybind11::value_error(errorMsg.c_str());
              }

              pybind11::tuple ret(2);
              ret[0] = value;
              ret[1] = uncertainty;
              return ret;              
            },
            py::arg("x"),
            R"(
Interpolate the value within the 1D grid at the specified position

Args:
    x (float) : Positions within the grid.

Returns:
    A tuple containing the interpolated value and uncertainty.

Raises:
    ValueError: If the interpolation fails.
)");
  }
};

// Partial specialization for 2D.
template <typename T, typename PyClass>
struct specialMethods<T, 2, PyClass>{
  static void add(PyClass& cls){
    cls.def("extractValue2D",
            [](const penred::measurements::results<T, 2>& obj,
               const double x,
               const double y) -> py::tuple{
              
              T value;
              double uncertainty;

              int err = obj.extractValue(x, y, value, uncertainty);
              if(err != penred::measurements::errors::SUCCESS) {
                std::string errorMsg("Error extracting interpolated value from results. ");
                errorMsg += penred::measurements::errorToString(err);
                throw pybind11::value_error(errorMsg.c_str());
              }

              pybind11::tuple ret(2);
              ret[0] = value;
              ret[1] = uncertainty;
              return ret;              
            },
            py::arg("x"),
            py::arg("y"),
            R"(
Interpolate the value within the 2D grid at the specified position

Args:
    x (float) : Positions within the grid for the first dimension.
    y (float) : Positions within the grid for the second dimension.

Returns:
    A tuple containing the interpolated value and uncertainty.

Raises:
    ValueError: If the interpolation fails.
)");
  }
};

// Partial specialization for 3D.
template <typename T, typename PyClass>
struct specialMethods<T, 3, PyClass>{
  static void add(PyClass& cls){
    cls.def("extractValue3D",
            [](const penred::measurements::results<T, 3>& obj,
               const double x,
               const double y,
               const double z) -> py::tuple{
              
              T value;
              double uncertainty;

              int err = obj.extractValue(x, y, z, value, uncertainty);
              if(err != penred::measurements::errors::SUCCESS) {
                std::string errorMsg("Error extracting interpolated value from results. ");
                errorMsg += penred::measurements::errorToString(err);
                throw pybind11::value_error(errorMsg.c_str());
              }

              pybind11::tuple ret(2);
              ret[0] = value;
              ret[1] = uncertainty;
              return ret;              
            },
            py::arg("x"),
            py::arg("y"),
            py::arg("z"),
            R"(
Interpolate the value within the 3D grid at the specified position

Args:
    x (float) : Positions within the grid for the first dimension.
    y (float) : Positions within the grid for the second dimension.
    z (float) : Positions within the grid for the third dimension.

Returns:
    A tuple containing the interpolated value and uncertainty.

Raises:
    ValueError: If the interpolation fails.
)");
  }
};

// * Register templates

template <typename T, std::size_t D>
void registerResults(py::module_& m, const std::string& name){

  using R = penred::measurements::results<T, D>;

  const std::string classDoc =
    "Results container for a maximum of " + std::to_string(D) + "D data.\n"
    "\n"
    "Stores values and their uncertainties on a " +
    std::to_string(D) + "-dimensional grid.\n";
  
  py::class_<R> cls(m, name.c_str(), classDoc.c_str());

  cls.def(py::init<>())
    .def("title", [](const R& obj) -> std::string{
      return obj.title;
    },
         R"(
Reads the results' title

Args:
    None

Returns:
    A string containing the results title

Raises:
    None
)")
    .def("description", [](const R& obj) -> std::string{
      return obj.description;
    },
         R"(
Reads the results' description

Args:
    None

Returns:
    A string containing the results' description

Raises:
    None
)")
    .def("maxDim", [](const R& /*obj*/) -> int{
      return D;
    },
         R"(

Args:
    None

Returns:
    The maximum number of dimensions

Raises:
    None
)")
    .def("load", &resultsLoad<T, D>,
         py::arg("filename"),
         R"(
Load results from a text file

Args:
    filename (str) : Data file.

Returns:
    None

Raises:
    ValueError: If file parsing fails or the fail is unrecheable.
)")
    .def("save", [](const R& obj,
                    const std::string& filename,
                    const unsigned nSigma,
                    const bool printCoordinates,
                    const bool printBinNumber) -> void{
      //Open output file
      FILE* fout = fopen(filename.c_str(),"w");
      if(fout == nullptr){
        std::string errorMsg("Unable to open file ");
        errorMsg += filename;
        throw std::runtime_error(errorMsg.c_str());
      }

      //Print data
      obj.print(fout, nSigma, printCoordinates, printBinNumber, true);
      fclose(fout);
    },
         py::arg("filename"),
         py::arg("print_sigmas") = 2,
         py::arg("print_coordinates") = true,
         py::arg("print_bins") = true,
         R"(
Writes the stored data as standard penRed's results file.

Args:
    filename (str) : Results filename.
    print_sigmas (int): Specify the number of printed sigmas, i.e. the uncertainties column will be transformed as print_sigmas*sigma
    print_coordinates (bool): If enabled, the coordinates of each dimension will be printed
    print_bins (bool): If enabled, the bin number of each dimension will be printed
Returns:
    None

Raises:
    ValueError: If file parsing fails or the fail is unrecheable.
)")
    
    .def("bins",
         [](const R& obj, const unsigned dim) -> int {
           if(dim >= D)
             return obj.getNBins();
           else
             return obj.getNBins(dim);
         },
         py::arg("dimension") = D+1,
         R"(

Args:
    dimension (int): The dimension for which the number of bins is requested. If this value exceeds the maximum number of dimensions, or it is not provided, the total number of bins across all dimensions is returned instead.

Returns:
    int: The number of bins in the specified dimension, or the total number of bins if the dimension is out of bounds.

Raises:
    None
)")
    
    .def("data",
         [](const R& obj) -> py::tuple {
           return result2numpy(obj, true, false);
         },
         R"(
Extract result's data as numpy arrays

Args:
    None

Returns:
    A tuple containing three elements. First, two numpy arrays with the values and uncertainties respectivelly, and then a list with the dimension's information.

Raises:
    ValueError: If transformation fails. If that happens, please, report the error.
)")
    .def("fill", &numpy2result<T, D>,
         py::arg("values"),
         py::arg("sigma"),
         py::arg("info"),
         py::arg("value_header") = "Value",
         R"(
Fills the results object with the provided data.

Args:
    values (vector) : Numpy vector with the values to be saved
    sigma (vector) : Numpy vector with the uncertainty of each value corresponding to one standard deviation
    info (list): List of tuples where each element contains the information of the corresponding dimension. Each element must use the following format: (min, max, header), where *min* and *max* are the minimum and maximum grid value for this dimension, respectively, and *header* the text header for dimension's column.
    value_header (str): Header assigned to *values*. An example for a energetic spectrum could be 'E (eV)'.

Returns:
    None

Raises:
    TypeError: Incompatible types have been provided.
    ValueError: Dimensions mismatch or incompatible values have been provided.
)")
    .def("extractValue", &extractValue<T,D>,
         py::arg("position"),
         R"(
Interpolate the value within the grid at the specified position in each dimension

Args:
    position (list) : List of positions for each dimension within the grid.

Returns:
    A tuple containing the interpolated value and uncertainty.

Raises:
    ValueError: If the interpolation fails.
)")

    .def("extractSpectrum1D", &extractSpectrum1D<T,D>,
         py::arg("spectrum_dim"),
         py::arg("positions"),
         R"(
Creates a 1D spectrum for the specified dimension and position within the grid interpolation values and uncertainties in each dimension

Args:
    spectrum_dim (int) : Dimension to extract the spectrum from.
    positions (list): Positions for the remaining dimensions. The spectrum dimension is not included. Therefore, the position list must have nDim-1 elements, where nDim is the maximum dimensions for the current results object.

Returns:
    A tuple containing the interpolated spectrum values and uncertainties.

Raises:
    ValueError: If the interpolation fails.
)");

  // Add the dimension-specific methods.
  specialMethods<T, D, py::class_<R>>::add(cls);
}

// Forward declaration of the recursion.
template <typename T, std::size_t D, std::size_t MaxD>
struct registerResultsRange;

// Base case: D == MaxD (nothing left to register).
template <typename T, std::size_t MaxD>
struct registerResultsRange<T, MaxD, MaxD>{
  static void apply(py::module_&, const std::string&) {}
};

// Recursive case: register D, then recurse on D+1.
template <typename T, std::size_t D, std::size_t MaxD>
struct registerResultsRange{
  static void apply(py::module_& m, const std::string& base){
    registerResults<T, D>(m, base + std::to_string(D) + "D");
    registerResultsRange<T, D + 1, MaxD>::apply(m, base);
  }
};

// Convenience entry point.
template <typename T, std::size_t MinD, std::size_t MaxD>
void registerResultsRangeAll(py::module_& m, const std::string& base)
{
    registerResultsRange<T, MinD, MaxD>::apply(m, base);
}

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
    
    penred::measurements::results<double, penred::measurements::maxDims> results;
    numpy2result(results, values,sigma, info, valueHeader);

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

Print the provided data as standard penRed's results file.

Args:
    filename (str) : Results filename.
    values (vector) : Numpy vector with the values to be saved
    sigma (vector) : Numpy vector with the uncertainty of each value corresponding to one sigma
    info (list): List of tuples where each element contains the information of the corresponding dimension. Each element must use the following format: (min, max, header), where *min* and *max* are the minimum and maximum grid value for this dimension, respectively, and *header* the text header for dimension's column.
    value_header (str): Header assigned to *values*. An example for a energetic spectrum could be 'E (eV)'.
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

  // Registers results1D, results2D, ..., results{maxDims}D
  registerResultsRangeAll<double, /*Min=*/1, /*Max=*/11>(m, "results");
  registerResults<double, penred::measurements::maxDims>(m, "results");
  
}
