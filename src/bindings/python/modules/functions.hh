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
#include <type_traits>
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

template<typename T, size_t dim, class Grid>
pybind11::tuple result2numpy(const penred::measurements::results<T, dim, Grid>& results, const bool extractInfo, const bool onlyEffective){

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
    const auto edges = results.readGrid().getInitEdges();

    unsigned ituple = 2;
    for(int i = static_cast<int>(dim)-1; i >= 0; --i){

      if(onlyEffective && nBins[i] <= 1){
        //Skip "empty" dimensions
        //printf("Skipping: %s (%lu)\n", results.readDimHeader(i).c_str(), nBins[i]);
        continue;
      }
      
      if(std::is_same<Grid, penred::grids::RegularGrid<dim>>::value){
        pybind11::tuple dimInfo(3);
        dimInfo[0] = edges[i][0];
        dimInfo[1] = edges[i][1];
        dimInfo[2] = results.readDimHeader(i);
        pyRes[ituple++] = dimInfo;
      }
      else{
        pybind11::tuple dimInfo(2);
        dimInfo[0] = edges[i];
        dimInfo[1] = results.readDimHeader(i);
        pyRes[ituple++] = dimInfo;
      }
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

template<typename T, size_t dim, class Grid>
void numpy2result(penred::measurements::results<T, dim, Grid>& obj,
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
  std::vector<std::vector<double>> edges;
  std::vector<std::string> headers;
  for(size_t i = 0; i < info.size(); ++i){
    //Check if the element is a tuple
    if (!pybind11::isinstance<pybind11::tuple>(info[i])) {
      throw pybind11::type_error("Information element " + std::to_string(i) + " is not a tuple");
    }

    pybind11::tuple t = info[i];
    if(std::is_same<Grid, penred::grids::RegularGrid<dim>>::value){
      if(t.size() != 3){
        throw pybind11::value_error("Regular-grid information tuples must contain (min, max, header)");
      }
      if((!pybind11::isinstance<pybind11::int_>(t[0]) &&
          !pybind11::isinstance<pybind11::float_>(t[0])) ||
         (!pybind11::isinstance<pybind11::int_>(t[1]) &&
          !pybind11::isinstance<pybind11::float_>(t[1]))){
        throw pybind11::type_error("Regular-grid limits must be numeric");
      }
      if(!pybind11::isinstance<pybind11::str>(t[2])){
        throw pybind11::type_error("Dimension header must be a string");
      }
      edges.push_back({t[0].cast<double>(), t[1].cast<double>()});
      headers.push_back(t[2].cast<std::string>());
    }
    else{
      if(t.size() != 2){
        throw pybind11::value_error("Variable-grid information tuples must contain (edges, header)");
      }
      if(!pybind11::isinstance<pybind11::str>(t[1])){
        throw pybind11::type_error("Dimension header must be a string");
      }
      pybind11::sequence edgeSequence = t[0].cast<pybind11::sequence>();
      std::vector<double> dimensionEdges;
      for(const pybind11::handle& edge : edgeSequence){
        if(!pybind11::isinstance<pybind11::int_>(edge) &&
           !pybind11::isinstance<pybind11::float_>(edge)){
          throw pybind11::type_error("Variable-grid edges must be numeric");
        }
        dimensionEdges.push_back(edge.cast<double>());
      }
      edges.push_back(std::move(dimensionEdges));
      headers.push_back(t[1].cast<std::string>());
    }
  }
    
  //Reverse bins per dimension and edges to correct numpy ordering
  std::reverse(dimBins.begin(), dimBins.end());
  std::reverse(edges.begin(), edges.end());
  std::reverse(headers.begin(), headers.end());

  //Complete dimensions omitted by the input with a single default bin.
  std::array<unsigned long, dim> bins;
  std::array<std::vector<double>, dim> gridEdges;
  bins.fill(1ul);
  for(std::vector<double>& dimensionEdges : gridEdges)
    dimensionEdges = {-1.0e35, 1.0e35};
  for(size_t i = 0; i < dimBins.size(); ++i){
    bins[i] = dimBins[i];
    gridEdges[i] = edges[i];
  }

  //Init it using the common edge-based interface. RegularGrid interprets each
  //two-element edge vector as its lower and upper limits.
  int err = obj.init(bins, gridEdges);
  if(err == penred::measurements::errors::SUCCESS){
    if(obj.data.size() != static_cast<size_t>(values.size()) ||
       obj.sigma.size() != static_cast<size_t>(sigma.size())){
      err = penred::measurements::errors::DIMENSION_MISMATCH;
    }
    else{
      obj.data = flattenArray(values);
      obj.sigma = flattenArray(sigma);
    }
  }
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
template<typename T, size_t dim, class Grid>
void resultsLoad(penred::measurements::results<T, dim, Grid>& obj, const std::string& filename){

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
template<typename T, size_t dim, class Grid>
pybind11::tuple extractValue(const penred::measurements::results<T, dim, Grid>& obj,
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

template<typename T, size_t dim, class Grid>
pybind11::tuple extractSpectrum1D(const penred::measurements::results<T, dim, Grid>& obj,
                                  const unsigned spectrumDim,
                                  const std::array<double, dim-1>& position){
  penred::measurements::results<T, 1, penred::grids::rebind_t<1, Grid>> spectrum;

  int err = obj.extractSpectrum1D(spectrumDim, position, spectrum);
  if(err != penred::measurements::errors::SUCCESS) {
    std::string errorMsg("Error extracting interpolated spectrum from results. ");
    errorMsg += penred::measurements::errorToString(err);
    throw pybind11::value_error(errorMsg.c_str());
  }

  return result2numpy(spectrum, false, false);
}

#endif
