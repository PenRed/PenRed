
//
//
//    Copyright (C) 2025 Vicent Giménez Alventosa
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
//    
//


#include "tallyTrackingMemory.hh"


void pen_tallyTrackingMemory::flush(){
}

void pen_tallyTrackingMemory::tally_beginSim(){
  
  active = true;
  sampledToPrint = false;
  nStates = 0;
}

void pen_tallyTrackingMemory::tally_endSim(const unsigned long long){
  
  active = false;
}

void pen_tallyTrackingMemory::tally_beginPart(const unsigned long long /*nhist*/,
					      const unsigned /*kdet*/,
					      const pen_KPAR kpar,
					      const pen_particleState& state){

  if(active && enabledKpar[kpar]){
    
    //Check if a previous sampled state must be printed
    if(sampledToPrint){
      encodeState(auxState);
      sampledToPrint = false;
    }
    else{    
      //Print particle start position
      encodeState(state);    
    }
  }
}

void pen_tallyTrackingMemory::tally_sampledPart(const unsigned long long nhist,
						const unsigned long long /*dhist*/,
						const unsigned /*kdet*/,
						const pen_KPAR /*kpar*/,
						const pen_particleState& state){

  if(nhist < lastHist){
    
    //Print sample position if it is in void
    if(state.MAT == 0){
      //Save state to print it if reach the geometry system
      auxState = state;
      sampledToPrint = true;
    }
    
  } else {
    active = false;
  }
}

void pen_tallyTrackingMemory::tally_endPart(const unsigned long long /*nhist*/,
					    const pen_KPAR kpar,
					    const pen_particleState& state){
  sampledToPrint = false;
  if(state.MAT > 0 && enabledKpar[kpar]){
    encodeState(state);
    //Flag end of track
    encodeState(pen_particleState());
  }
}

void pen_tallyTrackingMemory::tally_step(const unsigned long long /*nhist*/,
					 const pen_KPAR kpar,
					 const pen_particleState& state,
					 const tally_StepData& stepData){
  if(active && enabledKpar[kpar]){

    //Check if the particle scapes the geometry system
    if(state.MAT == 0){
      //Create a new state moving the particle to the geometry boundary
      auxState = state;
      auxState.X = xlast + auxState.U*stepData.dsef;
      auxState.Y = ylast + auxState.V*stepData.dsef;
      auxState.Z = zlast + auxState.W*stepData.dsef;
      encodeState(auxState);
      //Flag end of track
      if(nStates > 0 && tracks[1 + (nStates-1)*STATESIZE] < 1.0)
	return; //The last state is already a track end      
      encodeState(pen_particleState());
    }
  }

}

void pen_tallyTrackingMemory::tally_jump(const unsigned long long /*nhist*/,
					 const pen_KPAR /*kpar*/,
					 const pen_particleState& state,
					 const double /*ds*/){
  if(active){

    //Save particle last position
    xlast = state.X;
    ylast = state.Y;
    zlast = state.Z;
  }
}

void pen_tallyTrackingMemory::tally_knock(const unsigned long long /*nhist*/,
					  const pen_KPAR kpar,
					  const pen_particleState& state,
					  const int /*icol*/){
  if(active && enabledKpar[kpar]){
    //Print state after interaction
    encodeState(state);
  }
}

int pen_tallyTrackingMemory::configure(const wrapper_geometry& /*geometry*/,
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

  for(size_t i = 0; i < constants::nParTypes; ++i){
    std::string enableKey("enable/");
    enableKey += particleName(i);
    if(config.read(enableKey.c_str(), enabledKpar[i]) != INTDATA_SUCCESS){
      enabledKpar[i] = true;
    }
    printf("Recording %s %s\n",
	   particleName(i), enabledKpar[i] ? "enabled" : "disabled");
  }
  
  return 0;
}



void pen_tallyTrackingMemory::saveData(const unsigned long long /*nhist*/) const{}
int pen_tallyTrackingMemory::sumTally(const pen_tallyTrackingMemory& /*tally*/){return 0;}



REGISTER_COMMON_TALLY(pen_tallyTrackingMemory)
