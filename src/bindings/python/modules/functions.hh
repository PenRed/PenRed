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

#ifndef __PYPENRED_FUNCTIONS__
#define __PYPENRED_FUNCTIONS__

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <algorithm>
#include <fstream>
#include "pen_data.hh"
#include "math_classes.hh"

std::string dict2SectionStringWithPrefix(const pybind11::dict& dict, const std::string& prefixIn);

std::string dict2SectionString(const pybind11::dict& dict);

inline int dict2section(const pybind11::dict& dict,
			pen_parserSection& result,
			std::string& errorString){
  
  std::string text = dict2SectionStringWithPrefix(dict, "");

  unsigned long errorLine;
  return parseString(text, result, errorString, errorLine);
}

// Numpy utilities

template<typename T>
size_t assertShapes(const pybind11::array_t<T>& arr1,
                    const pybind11::array_t<T>& arr2,
                    std::vector<unsigned long>& dimsSizes) {
  
  // Get shape information
  pybind11::buffer_info buf1 = arr1.request();
  pybind11::buffer_info buf2 = arr2.request();
    
  // Compare number of dimensions
  if (buf1.ndim != buf2.ndim) {
    std::string errorMsg("Different number of dimensions: ");
    errorMsg += std::to_string(buf1.ndim);
    errorMsg += " vs ";
    errorMsg += std::to_string(buf2.ndim);
    throw pybind11::buffer_error(errorMsg);
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
      throw pybind11::buffer_error(errorMsg);	  
      std::cout << "Dimension " << i << " mismatch: " 
                << buf1.shape[i] << " vs " << buf2.shape[i] << std::endl;
    }
    dimsSizes.push_back(static_cast<unsigned long>(buf1.shape[i]));
  }

  return buf1.ndim;
}

template<typename T>
std::vector<T> flattenArray(const pybind11::array_t<T>& arr) {
  pybind11::buffer_info buf = arr.request();
    
  // Get pointer to data
  T* ptr = static_cast<T*>(buf.ptr);
    
  // Create vector from data
  std::vector<T> result(ptr, ptr + buf.size);
    
  return result;
}

// + Results value extraction

template<typename T, size_t dim>
pybind11::tuple result2numpy(const penred::measurements::results<T, dim>& results, const bool extractInfo, const bool onlyEffective){

  //Get bins in each dimension
  std::array<unsigned long, dim> nBins = results.readDimBins();
  //Reverse bins to fit the numpy ordering for dimensions
  std::array<unsigned long, dim> nBinsReverse = nBins;
  std::reverse(nBinsReverse.begin(), nBinsReverse.end());

  //Get bins for effective dimensions (whith nbins > 1)
  std::vector<unsigned long> nEffectiveBins;
  for(unsigned long dimBins : nBinsReverse){
    if(dimBins > 1){
      nEffectiveBins.push_back(dimBins);
    }
  }

  //Calculate the final tuple size
  size_t resTupleSize = 2;
  if(extractInfo){
    if(onlyEffective)
      resTupleSize += nEffectiveBins.size() + 2;
    else
      resTupleSize += dim + 2;
  }

  //Create the results tuple
  pybind11::tuple pyRes(resTupleSize);  

  //Save the data (values and uncertainties) along with the dimensions
  if(onlyEffective){
    pyRes[0] = pybind11::array_t<T>(nEffectiveBins, results.data.data());
    pyRes[1] = pybind11::array_t<T>(nEffectiveBins, results.sigma.data());
  }
  else{
    pyRes[0] = pybind11::array_t<T>(nBinsReverse, results.data.data());
    pyRes[1] = pybind11::array_t<T>(nBinsReverse, results.sigma.data());
  }
  
  if(extractInfo){
    //Get interval and description information for each dimension
    const std::array<std::pair<double, double>, dim> limits = results.readLimits();

    unsigned ituple = 2;
    for(int i = static_cast<int>(dim)-1; i >= 0; --i){

      if(onlyEffective && nBins[i] <= 1){
        //Skip "empty" dimensions
        //printf("Skipping: %s (%lu)\n", results.readDimHeader(i).c_str(), nBins[i]);
        continue;
      }
      
      pybind11::tuple dimInfo(3);
      dimInfo[0] = limits[i].first;
      dimInfo[1] = limits[i].second;
      const std::string header = results.readDimHeader(i);
      dimInfo[2] = header;

      //Append this dimension to returned results
      pyRes[ituple++] = dimInfo;
    }

    //Append value and description info
    pyRes[ituple++] = results.readValueHeader();
    pyRes[ituple++] = results.description;
  }

  return pyRes;
}

template<typename T>
pybind11::array_t<T> result2numpy(const std::vector<T>& results, const bool, const bool){

  pybind11::tuple pyRes(1);  
  pyRes[0] = pybind11::array_t<T>(results.size(), results.data());

  return pyRes;
}

// + Results constructor extraction

template<typename T, size_t dim>
void numpy2result(penred::measurements::results<T, dim>& obj,
                  const pybind11::array_t<T>& values,
                  const pybind11::array_t<T>& sigma,
                  const pybind11::list& info,
                  const std::string& valueHeader){

  //Check dimensions
  std::vector<unsigned long> dimBins;
  const size_t nDim = assertShapes(values, sigma, dimBins);

  if(nDim > dim){
    throw pybind11::type_error("Number of data dimensions is larger than results dimensions");
  }

  //Ensure the information is a list
  if (!pybind11::isinstance<pybind11::list>(info) && !pybind11::isinstance<pybind11::tuple>(info)) {
    throw pybind11::type_error("Argument 'info' must be a list");
  }

  //Check information length
  if(info.size() != nDim){
    throw pybind11::value_error("Data and information dimensions mismatch");      
  }    

  //Check and extract information
  std::vector<std::pair<double, double>> limits;
  std::vector<std::string> headers;
  for(size_t i = 0; i < info.size(); ++i){
    //Check if the element is a tuple
    if (!pybind11::isinstance<pybind11::tuple>(info[i])) {
      throw pybind11::type_error("Information element " + std::to_string(i) + " is not a tuple");
    }

    pybind11::tuple t = info[i];
    if(t.size() != 3){
      throw pybind11::value_error("Tuple " + std::to_string(i) + " has " + 
                                  std::to_string(t.size()) + " elements, expected 3");
    }

    // Check first element is numeric
    if (!pybind11::isinstance<pybind11::int_>(t[0]) && !pybind11::isinstance<pybind11::float_>(t[0])) {
      throw pybind11::type_error("Tuple " + std::to_string(i) + 
                                 ", first element must be numeric (int or float)");
    }
        
    // Check second element is numeric
    if (!pybind11::isinstance<pybind11::int_>(t[1]) && !pybind11::isinstance<pybind11::float_>(t[1])) {
      throw pybind11::type_error("Tuple " + std::to_string(i) + 
                                 ", second element must be numeric (int or float)");
    }

    // Check third element is string
    if (!pybind11::isinstance<pybind11::str>(t[2])) {
      throw pybind11::type_error("Tuple " + std::to_string(i) + 
                                 ", third element must be a string");
    }

    // Extract values
    double min = t[0].cast<double>();
    double max = t[1].cast<double>();
    std::string header = t[2].cast<std::string>();

    limits.emplace_back(min,max);
    headers.push_back(std::move(header));
  }
    
  //Reverse bins per dimension and limits to correct numpy ordering
  std::reverse(dimBins.begin(), dimBins.end());
  std::reverse(limits.begin(), limits.end());
  std::reverse(headers.begin(), headers.end());

  //Init it
  int err = obj.init(dimBins, limits, flattenArray(values), flattenArray(sigma));
  if(err != penred::measurements::errors::SUCCESS){
    throw pybind11::value_error(penred::measurements::errorToString(err));
  }

  //Set headers
  for(size_t i = 0; i < headers.size(); ++i){
    obj.setDimHeader(i,headers[i]);
  }
  obj.setValueHeader(valueHeader);
}

// + Results load
template<typename T, size_t dim>
void resultsLoad(penred::measurements::results<T, dim>& obj, const std::string& filename){

  //Read data
  std::ifstream fin(filename, std::ifstream::in);
  if(!fin){
    std::string errorMsg("Unable to open file ");
    errorMsg += filename;
    throw pybind11::value_error(errorMsg.c_str());
  }
  
  int err = obj.read(fin);
  fin.close();
  if(err != 0){
    std::string errorMsg("Error reading data file. ");
    errorMsg += penred::measurements::errorToString(err);
    throw pybind11::value_error(errorMsg.c_str());
  }
}

// + Results extract values
template<typename T, size_t dim>
pybind11::tuple extractValue(const penred::measurements::results<T, dim>& obj,
                             const std::array<double, dim>& position) {
  T value;
  double uncertainty;

  int err = obj.extractValue(position, value, uncertainty);
  if(err != penred::measurements::errors::SUCCESS) {
    std::string errorMsg("Error extracting interpolated value from results. ");
    errorMsg += penred::measurements::errorToString(err);
    throw pybind11::value_error(errorMsg.c_str());
  }

  pybind11::tuple ret(2);
  ret[0] = value;
  ret[1] = uncertainty;
  return ret;
}

template<typename T, size_t dim>
pybind11::tuple extractSpectrum1D(const penred::measurements::results<T, dim>& obj,
                                  const unsigned spectrumDim,
                                  const std::array<double, dim-1>& position){
  penred::measurements::results<T, 1> spectrum;

  int err = obj.extractSpectrum1D(spectrumDim, position, spectrum);
  if(err != penred::measurements::errors::SUCCESS) {
    std::string errorMsg("Error extracting interpolated spectrum from results. ");
    errorMsg += penred::measurements::errorToString(err);
    throw pybind11::value_error(errorMsg.c_str());
  }

  return result2numpy(spectrum, false, false);
}

#endif
