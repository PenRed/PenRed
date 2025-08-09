
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
//        vicent.gimenez.alventosa@gmail.com  (Vicent Giménez Alventosa)
//        vicente.gimenez@uv.es (Vicente Giménez Gómez)
//        sanolgi@upvnet.upv.es (Sandra Oliver Gil)
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

  if(debug){
    fout = fopen("tracking-debug.dat","w");
    fprintf(fout,"%s\n\n",baseStateHeader());
  }

  for(size_t i = 0; i < constants::nParTypes; ++i){
    std::string filename = std::string(particleName(i)) + "-tracks.dat";
    trackFiles[i] = fopen(filename.c_str(), "w");
    fprintf(trackFiles[i],"#%s\n",baseStateHeader());
  }
  
  active = true;
  sampledToPrint = false;
}

void pen_tallyTracking::tally_endSim(const unsigned long long){

  if(debug){
    if(fout != nullptr){
      fclose(fout);
      fout = nullptr;
    }
  }
  
  for(size_t i = 0; i < constants::nParTypes; ++i){
    if(trackFiles[i] != nullptr){
      fclose(trackFiles[i]);
      trackFiles[i] = nullptr;
    }
  }

  active = false;
}

void pen_tallyTracking::tally_beginPart(const unsigned long long /*nhist*/,
					const unsigned /*kdet*/,
					const pen_KPAR kpar,
					const pen_particleState& state){

  if(active){
    
    if(debug){
      fprintf(fout,"#Begin particle simulation of kpar %d:\n",kpar);
      printstate(fout,state);
    }

    //Check if a previous sampled state must be printed
    if(sampledToPrint){
      printstate(trackFiles[kpar],auxState);
      sampledToPrint = false;
    }
    
    //Print particle start position
    printstate(trackFiles[kpar],state);    
    
  }
}

void pen_tallyTracking::tally_sampledPart(const unsigned long long nhist,
					  const unsigned long long dhist,
					  const unsigned /*kdet*/,
					  const pen_KPAR kpar,
					  const pen_particleState& state){

  if(nhist < nhists){

    if(debug){
      if(dhist > 0){
	fprintf(fout,"#Begin history %llu:\n",nhist);
	fprintf(fout,"#      kpar -> %d\n",kpar);
	fprintf(fout,"#X,Y,Z (cm) -> %12.4E, %12.4E, %12.4E \n",
		state.X,state.Y,state.Z);
	fprintf(fout,"#    E (eV) -> %12.4E\n",state.E);
	printstate(fout,state);
      } else{

	fprintf(fout,"#Sampled particle in the same history %llu:\n",nhist);
	fprintf(fout,"#      kpar -> %d\n",kpar);
	fprintf(fout,"#X,Y,Z (cm) -> %12.4E, %12.4E, %12.4E \n",
		state.X,state.Y,state.Z);
	fprintf(fout,"#    E (eV) -> %12.4E\n",state.E);
	printstate(fout,state);      
      
      }
    }

    //Print sample position if it is in void
    if(state.MAT == 0){
      //Save state to print it if reach the geometry system
      auxState = state;
      sampledToPrint = true;
    }
    
  } else{
    active = false;
  }
}

void pen_tallyTracking::tally_move2geo(const unsigned long long /*nhist*/,
				       const unsigned /*kdet*/,
				       const pen_KPAR /*kpar*/,
				       const pen_particleState& state,
				       const double /*dsef*/,
				       const double /*dstot*/){

  if(active && debug){
    fprintf(fout,"# Moved to geometry.\n");
    printstate(fout,state);
  }  
}

void pen_tallyTracking::tally_endPart(const unsigned long long /*nhist*/,
				    const pen_KPAR kpar,
				    const pen_particleState& state){
  if(active){

    if(debug){
      fprintf(fout,"# Particle simulation ended.\n");
      printstate(fout,state);
    }

    //Let empty lines to separate particles
    if(!sampledToPrint){
      fprintf(trackFiles[kpar], "\n\n");
    }
    else{
      sampledToPrint = false;
    }
  }
}

void pen_tallyTracking::tally_step(const unsigned long long /*nhist*/,
				   const pen_KPAR kpar,
				   const pen_particleState& state,
				   const tally_StepData& stepData){
  if(active){
    if(debug){
      fprintf(fout,"# Particle moved dsef = %12.4E cm, dstot = %12.4E cm, deSoft = %12.4E eV.\n",
	      stepData.dsef,stepData.dstot,stepData.softDE);
      fprintf(fout,"# Energy deposited at (%12.4E,%12.4E,%12.4E).\n",
	      stepData.softX,stepData.softY,stepData.softZ);
      fprintf(fout,"# Previous IBody: %u, actual IBody: %u\n",
	      stepData.originIBODY,state.IBODY);
      printstate(fout,state);
    }

    //Check if the particle scapes the geometry system
    if(state.MAT == 0){
      //Create a new state moving the particle to the geometry boundary
      auxState = state;
      auxState.X = xlast + auxState.U*stepData.dsef;
      auxState.Y = ylast + auxState.V*stepData.dsef;
      auxState.Z = zlast + auxState.W*stepData.dsef;
      printstate(trackFiles[kpar],auxState);
    }else{
      //Print state after moving
      printstate(trackFiles[kpar],state);
    }
  }

}

void pen_tallyTracking::tally_interfCross(const unsigned long long /*nhist*/,
					  const unsigned /*kdet*/,
					  const pen_KPAR /*kpar*/,
					  const pen_particleState& /*state*/){
  if(active && debug){
    fprintf(fout,"# Particle crossed an interface\n");
  }

}

void pen_tallyTracking::tally_jump(const unsigned long long /*nhist*/,
				 const pen_KPAR /*kpar*/,
				 const pen_particleState& state,
				 const double ds){
  if(active){

    if(debug){
      fprintf(fout,"# Particle will jump %12.4Ecm.\n",ds);
      printstate(fout,state);
    }

    //Save particle last position
    xlast = state.X;
    ylast = state.Y;
    zlast = state.Z;
  }
}

void pen_tallyTracking::tally_knock(const unsigned long long /*nhist*/,
				  const pen_KPAR kpar,
				  const pen_particleState& state,
				  const int icol){
  if(active){
    if(debug){
      fprintf(fout,"# Particle suffers interaction %d.\n",icol);
      printstate(fout,state);
    }

    //Print state after interaction
    printstate(trackFiles[kpar],state);

  }
}

int pen_tallyTracking::configure(const wrapper_geometry& /*geometry*/,
				 const abc_material* const /*materials*/[constants::MAXMAT],
				 const pen_parserSection& config,
				 const unsigned verbose){
    
  int err;
  int nhistsAux;
  err = config.read("nhists", nhistsAux);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("Tracking:configure:unable to read 'nhists' in configuration. Integrer expected");
    }
    return -1;
  }
  active = true;

  if(nhistsAux <= 0){
    if(verbose > 0){
      printf("Tracking:configure: Invalid value to 'nhists', must be greater than 0.\n");
    }
    return -2;
  }
  
  nhists = static_cast<unsigned long long>(nhistsAux);
  if(verbose > 1){
    printf("Histories to track: %llu\n",nhists);
  }

  //Check if debug file is requested
  err = config.read("debug", debug);
  if(err != INTDATA_SUCCESS){
    debug = false;
  }else{
    if(verbose > 1){
      printf("Tracking debug is %s\n", debug ? "enabled" : "disabled");
    }
  }
  
  return 0;
}



void pen_tallyTracking::saveData(const unsigned long long /*nhist*/) const{}
int pen_tallyTracking::sumTally(const pen_tallyTracking& /*tally*/){return 0;}



REGISTER_COMMON_TALLY(pen_tallyTracking)
