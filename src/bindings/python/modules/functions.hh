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

size_t assertShapes(const pybind11::array_t<double>& arr1,
		    const pybind11::array_t<double>& arr2,
		    std::vector<unsigned long>& dimsSizes);

template<typename T>
std::vector<T> flattenArray(const pybind11::array_t<T>& arr) {
    pybind11::buffer_info buf = arr.request();
    
    // Get pointer to data
    T* ptr = static_cast<T*>(buf.ptr);
    
    // Create vector from data
    std::vector<T> result(ptr, ptr + buf.size);
    
    return result;
}

#endif
