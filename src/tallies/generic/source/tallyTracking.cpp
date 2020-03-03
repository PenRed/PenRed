
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


#include "tallyTracking.hh"


void pen_tallyTracking::flush(){
}

void printstate(FILE* fout, const pen_particleState& state){
  fprintf(fout,"%s\n",state.stringifyBase().c_str());
  fflush(fout);
}

void pen_tallyTracking::tally_beginSim(){
  fout = fopen("tracking.dat","w");
  fprintf(fout,"%s\n\n",baseStateHeader());  
  active = true;
}

void pen_tallyTracking::tally_beginPart(const double /*nhist*/,
					const unsigned /*kdet*/,
					const pen_KPAR kpar,
					const pen_particleState& state){

  if(active){
    fprintf(fout,"#Begin particle simulation of kpar %d:\n",kpar);
    printstate(fout,state);
  }
}

void pen_tallyTracking::tally_beginHist(const double nhist,
					const unsigned /*kdet*/,
					const pen_KPAR kpar,
					const pen_particleState& state){

  if(nhist < nhists){
    fprintf(fout,"#Begin history %20.0f:\n",nhist);
    fprintf(fout,"#      kpar -> %d\n",kpar);
    fprintf(fout,"#X,Y,Z (cm) -> %12.4E, %12.4E, %12.4E \n",
	    state.X,state.Y,state.Z);
    fprintf(fout,"#    E (eV) -> %12.4E\n",state.E);
    printstate(fout,state);
  } else{
    active = false;
  }
}

void pen_tallyTracking::tally_move2geo(const double /*nhist*/,
				       const unsigned /*kdet*/,
				       const pen_KPAR /*kpar*/,
				       const pen_particleState& state,
				       const double /*dsef*/,
				       const double /*dstot*/){

  if(active){
    fprintf(fout,"# Moved to geometry.\n");
    printstate(fout,state);
  }  
}

void pen_tallyTracking::tally_endPart(const double /*nhist*/,
				    const pen_KPAR /*kpar*/,
				    const pen_particleState& state){
  if(active){
    fprintf(fout,"# Particle simulation ended.\n");
    printstate(fout,state);
  }
}

void pen_tallyTracking::tally_step(const double /*nhist*/,
				   const pen_KPAR /*kpar*/,
				   const pen_particleState& state,
				   const tally_StepData& stepData){
  if(active){
    fprintf(fout,"# Particle moved dsef = %12.4E cm, dstot = %12.4E cm, deSoft = %12.4E eV.\n",stepData.dsef,stepData.dstot,stepData.softDE);
    fprintf(fout,"# Energy deposited at (%12.4E,%12.4E,%12.4E).\n",
	    stepData.softX,stepData.softY,stepData.softZ);
    fprintf(fout,"# Previous IBody: %u, actual IBody: %u\n",
	    stepData.originIBODY,state.IBODY);
    printstate(fout,state);
  }

}

void pen_tallyTracking::tally_interfCross(const double /*nhist*/,
					  const unsigned /*kdet*/,
					  const pen_KPAR /*kpar*/,
					  const pen_particleState& /*state*/){
  if(active){
    fprintf(fout,"# Particle crossed an interface\n");
  }

}

void pen_tallyTracking::tally_jump(const double /*nhist*/,
				 const pen_KPAR /*kpar*/,
				 const pen_particleState& state,
				 const double ds){
  if(active){
    fprintf(fout,"# Particle will jump %12.4Ecm.\n",ds);
    printstate(fout,state);
  }
}

void pen_tallyTracking::tally_knock(const double /*nhist*/,
				  const pen_KPAR /*kpar*/,
				  const pen_particleState& state,
				  const int icol){
  if(active){
    fprintf(fout,"# Particle suffers interaction %d.\n",icol);
    printstate(fout,state);
  }
}

int pen_tallyTracking::configure(const wrapper_geometry& /*geometry*/,
				 const abc_material* const /*materials*/[constants::MAXMAT],
				 const pen_parserSection& config,
				 const unsigned verbose){
    
  int err;     
  err = config.read("nhists", nhists);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("Tracking:configure:unable to read 'nhists' in configuration. Integrer expected");
    }
    return -1;
  }

  active = true;

  if(verbose > 1){
    printf("Histories to track: %d\n",nhists);
  }  

  
  return 0;
}



void pen_tallyTracking::saveData(const double /*nhist*/) const{}
int pen_tallyTracking::sumTally(const pen_tallyTracking& /*tally*/){return 0;}



REGISTER_COMMON_TALLY(pen_tallyTracking, TRACK)
