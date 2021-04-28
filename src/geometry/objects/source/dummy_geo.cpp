
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


#include "dummy_geo.hh"

void pen_dummyGeo::locate(pen_particleState& state) const{
  state.MAT = bodies[0].MATER;
  state.IBODY = 0;  
}

void pen_dummyGeo::step(pen_particleState& state, double DS, double &DSEF, double &DSTOT, int &NCROSS) const{

  DSEF = DS;
  DSTOT = DS;
  NCROSS = 0;
  state.X += DS*state.U;
  state.Y += DS*state.V;
  state.Z += DS*state.W;
}

REGISTER_GEOMETRY(pen_dummyGeo,DUMMY)
