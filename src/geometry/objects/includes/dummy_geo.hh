
//
//
//    Copyright (C) 2019-2021 Universitat de València - UV
//    Copyright (C) 2019-2021 Universitat Politècnica de València - UPV
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
  
  int configure(const pen_parserSection& /*config*/, unsigned verbose){
    NBODYS = 1; //The geometry consists of a single body
    bodies[0].MATER = 1; //of material 1
    bodies[0].DSMAX = 1.0e35; //with infinite eabs

    if(verbose > 1){
      printf("Dummy body:\n");
      printf("  Material: %d\n",bodies[0].MATER);
      printf("  DSMAX   : %E\n",bodies[0].DSMAX);
      printf("  KDET    : %d\n",bodies[0].KDET);
      printf("  EABS    :\n");
      for(unsigned k = 0; k < constants::nParTypes; k++)
	printf("  %20.20s: %14.5E\n",particleName(k),bodies[0].localEABS[k]);

    }
    
    return 0;
  }
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
