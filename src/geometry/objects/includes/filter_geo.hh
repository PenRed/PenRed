
//
//
//    Copyright (C) 2023 Universitat de València - UV
//    Copyright (C) 2023 Universitat Politècnica de València - UPV
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
//    
//

 
#ifndef __PEN_FILTER_GEOMETRY__
#define __PEN_FILTER_GEOMETRY__

#include "geometry_classes.hh"

struct pen_filterBody : public pen_baseBody{
  std::string name;
  double origin, limit;
  double width;

  inline bool isIn(const pen_particleState& state) const {
    return state.Z >= origin && state.Z < limit;
  }
};

class pen_filterGeo : public abc_geometry<pen_filterBody>{
  DECLARE_GEOMETRY(pen_filterGeo)

  double origin;
  
  public:

  enum errors{
    PEN_FILTER_GEO_SUCCESS = 0,
    PEN_FILTER_GEO_NO_FILTERS,
    PEN_FILTER_GEO_INVALID_NUMBER_OF_FILTERS,
    PEN_FILTER_GEO_MISSING_PARAMETER,
    PEN_FILTER_GEO_BAD_VALUE,
  };

  static constexpr const char* errorMessage(const int val) noexcept {
    switch(val){
    case PEN_FILTER_GEO_SUCCESS: return "Success";
    case PEN_FILTER_GEO_NO_FILTERS: return "No filters provided";
    case PEN_FILTER_GEO_INVALID_NUMBER_OF_FILTERS: return "Invalid number of filters";
    case PEN_FILTER_GEO_MISSING_PARAMETER: return "Missing parameter";
    case PEN_FILTER_GEO_BAD_VALUE: return "Bad parameter value";
    default: return "Unknown error";
    }
  }  
  
  pen_filterGeo() : origin(0.0) {
  }
  
  penred::errors::Error specificConfigure(const pen_parserSection& config, unsigned verbose);
  
  void locateLocal(pen_particleState&) const final override;
  void stepLocal(pen_particleState&,
                 double,
                 double &,
                 double &,
                 int &) const final override;

  inline unsigned getIBody(const char*) const {return getElements();}

  inline std::string getBodyName(const unsigned) const {return std::string("NONE");}
  
};

//Define error message function
template<>
constexpr const char* penred::errors::errorMessage<pen_filterGeo>(const int val) noexcept {
  return pen_filterGeo::errorMessage(val);
}

#endif
