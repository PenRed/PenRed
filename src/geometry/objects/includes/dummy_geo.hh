
//
//
//    Copyright (C) 2019 Universitat de València - UV
//    Copyright (C) 2019 Universitat Politècnica de València - UPV
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

 
#ifndef __PEN_DUMMY_GEOMETRY__
#define __PEN_DUMMY_GEOMETRY__

#include "geometry_classes.hh"

class pen_dummyGeo : public abc_geometry<pen_baseBody>{
  DECLARE_GEOMETRY(pen_dummyGeo)

  public:
  pen_dummyGeo() {
    configStatus = 0;
  }
  
  int configure(const pen_parserSection& /*config*/, unsigned /*verbose*/){
    return 0;
  }
  void locate(pen_particleState&) const;
  void step(pen_particleState&,
	    double,
	    double &,
	    double &,
	    int &) const;

  inline unsigned getIBody(const char*) const {return getElements();}
};

#endif
