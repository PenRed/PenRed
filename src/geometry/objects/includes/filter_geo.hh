
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

enum pen_filterGeoErr{
  PEN_FILTER_GEO_SUCCESS = 0,
  PEN_FILTER_GEO_NO_FILTERS,
  PEN_FILTER_GEO_INVALID_NUMBER_OF_FILTERS,
  PEN_FILTER_GEO_MISSING_POSITION,
  PEN_FILTER_GEO_INVALID_POSITION,
  PEN_FILTER_GEO_MISSING_WIDTH,
  PEN_FILTER_GEO_INVALID_WIDTH,
  PEN_FILTER_GEO_MISSING_MATERIAL,
  PEN_FILTER_GEO_INVALID_MATERIAL,
  PEN_FILTER_GEO_INVALID_DETECTOR,
  PEN_FILTER_GEO_KNOWN_PARTICLE,
  PEN_FILTER_GEO_BAD_EABS,
};

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
  pen_filterGeo() : origin(0.0) {
  }
  
  int configure(const pen_parserSection& config, unsigned verbose);
  
  void locate(pen_particleState&) const;
  void step(pen_particleState&,
	    double,
	    double &,
	    double &,
	    int &) const;

  inline unsigned getIBody(const char*) const {return getElements();}

  inline std::string getBodyName(const unsigned) const {return std::string("NONE");}
  
};

#endif
