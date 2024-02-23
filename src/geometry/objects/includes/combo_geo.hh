
//
//
//    Copyright (C) 2023-2024 Universitat de València - UV
//    Copyright (C) 2023-2024 Universitat Politècnica de València - UPV
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
//        sanolgi@upvnet.upv.es              (Sandra Oliver Gil)
//    
//

 
#ifndef __PEN_COMBO_GEOMETRY__
#define __PEN_COMBO_GEOMETRY__

#include <string>
#include "geometry_classes.hh"

struct pen_comboBody : public pen_baseBody{

  //Saves the geometry number to which this body belongs
  unsigned geoIndex;

  //The name of each body is assigned according to the following pattern:
  //     <geoIndex>_<bodyIndex>
  //
  // where <geoIndex> is the geometry number to which the body belongs and
  // <bodyName> is its name in the corresponding geometry
  std::string name;
};

class pen_comboGeo : public abc_geometry<pen_comboBody>{
  DECLARE_GEOMETRY(pen_comboGeo)

  private:

  std::vector<wrapper_geometry*> geometries; //Vector of combined geometries
  std::vector<unsigned> firstIBody; //First body index for each geometry
  
  public:
  pen_comboGeo() : abc_geometry<pen_comboBody>() {
    configStatus = 0;
  }

  const std::vector<wrapper_geometry*>& readGeometries() const { return geometries; }
  
  int configure(const pen_parserSection& /*config*/, unsigned verbose);
  void locate(pen_particleState&) const;
  void step(pen_particleState&,
	    double,
	    double &,
	    double &,
	    int &) const;
  
  unsigned getIBody(const char*) const final override;
  std::string getBodyName(const unsigned ibody) const final override;

  //Methods for nested geometries
  inline size_t nInternalGeometries() const final override {
    size_t nInternal = 0;
    for(const wrapper_geometry* p : geometries){
      nInternal += p->nInternalGeometries();
    }
    nInternal += geometries.size();
    return nInternal;
  }
  inline const wrapper_geometry* getInternalGeo(const size_t index) const final override {

    size_t actualPos = 0;
    for(size_t igeo = 0; igeo < geometries.size(); ++igeo){

      //Check inmidiatly internal geometry
      if(actualPos == index){
	return geometries[igeo];
      }

      //Increase actual position
      ++actualPos;

      //Check internal geometries of internal geometry
      const wrapper_geometry* pIGeo = geometries[igeo]->getInternalGeo(index-actualPos);
      if(pIGeo != nullptr)
	return pIGeo;      

      //Increase actual position skipping internal geometries
      actualPos += geometries[igeo]->nInternalGeometries();
    }

    return nullptr;
  }
  
  
};

#endif
