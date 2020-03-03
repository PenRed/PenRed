
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

 
#ifndef __PEN_LOOPS__
#define __PEN_LOOPS__

#include "PenRed.hh"
#include "pen_samplers.hh"
#include "pen_tallies.hh"
#include "pen_geometries.hh"

template <class particleType>

bool move2geo(const double nhist,
	      particleType& particle,
	      pen_commonTallyCluster& tallies){

  pen_particleState& state = particle.getBaseState();
  //Check if particle has been generated at vacuum
  if(state.MAT == 0){
    //Move it to geometry system
    particle.jumpVolume();
    
    //Call tallies with move2geo collect function
    unsigned kdet = particle.getDET();
    tallies.run_move2geo(nhist,kdet,particle.kpar,
			 state,particle.DSef(),particle.DStot());

    if(state.MAT == 0){
      tallies.run_endPart(nhist,particle.kpar,state);
      return false;
    }
  }
  else{
    particle.updateMat();
    particle.updateBody();
  }
  return true;
}


template <class particleType,
	  class stateType>

bool move2geo(const double nhist,
	      particleType& particle,
	      pen_commonTallyCluster& tallies,
	      pen_specificTallyCluster<stateType>& specificTallies){

  stateType& state = particle.getState();
  //Check if particle has been generated at vacuum
  if(state.MAT == 0){
    //Move it to geometry system
    particle.jumpVolume();
    
    //Call tallies with move2geo collect function
    unsigned kdet = particle.getDET();
    tallies.run_move2geo(nhist,kdet,particle.kpar,state,
			 particle.DSef(),particle.DStot());
    specificTallies.run_move2geo(nhist,kdet,particle.kpar,
				 state,particle.DSef(),particle.DStot());
    
    if(state.MAT == 0){
      tallies.run_endPart(nhist,particle.kpar,state);
      specificTallies.run_endPart(nhist,particle.kpar,state);
      return false;
    }
  }
  else{
    particle.updateMat();
    particle.updateBody();
  }  
  return true;
}


template <class particleType>

inline bool absorb(const double nhist,
		   particleType& particle,
		   pen_commonTallyCluster& tallies,
		   pen_rand& penRand){

  const pen_particleState& state = particle.readBaseState();

  //Check if particle must be absorbed
  if(state.E < particle.getEABS()){ 

    //Check if the particle has already been annihilated
    if(state.E > 1.0e-6){
      // run annihilation process
      double Eprod = particle.annihilationEDep;
      particle.annihilate(penRand);
      tallies.run_localEdep(nhist,particle.kpar,state,state.E+Eprod);
    }
    //Call tallies with end particle collect function
    tallies.run_endPart(nhist,particle.kpar,state);
    return true;
  }

  //Particle can't be absorbed yet
  return false;
}


template <class particleType,
	  class stateType>

inline bool absorb(const double nhist,
		   particleType& particle,
		   pen_commonTallyCluster& tallies,
		   pen_specificTallyCluster<stateType>& specificTallies,
		   pen_rand& penRand){

  const stateType& state = particle.readState();

  //Check if particle must be absorbed
  if(state.E < particle.getEABS()){ 

    //Check if the particle has already been annihilated
    if(state.E > 1.0e-6){
      // run annihilation process
      double Eprod = particle.annihilationEDep;
      particle.annihilate(penRand);
      tallies.run_localEdep(nhist,particle.kpar,state,state.E+Eprod);
      specificTallies.run_localEdep(nhist,particle.kpar,state,state.E+Eprod);
    }
    //Call tallies with end particle collect function
    tallies.run_endPart(nhist,particle.kpar,state);
    specificTallies.run_endPart(nhist,particle.kpar,state);      
    return true;
  }

  //Particle can't be absorbed yet
  return false;
}


template <class particleType>

void simulatePart(const double nhist,
		  particleType& particle,
		  pen_commonTallyCluster& tallies,
		  pen_rand& randoms){

  pen_particleState& state = particle.getBaseState();
  unsigned kdet = particle.getDET();
  const pen_KPAR kpar = particle.kpar;

  //Check if particle scaped from geometry
  if(state.MAT == 0){
    //Call tallies with end particle collect function
    tallies.run_endPart(nhist,particle.kpar,state);      
    return;
  }  
  
  //Check if the particle must be absorbed
  if(absorb(nhist,particle,tallies,randoms))
    return;

  particle.START();
  
  for(;;){

    double ds, de;
    int icol;

    //Create step data structure
    tally_StepData stepData;
    
    //Check for interaction forcing
    bool iforced = particle.readContext().isForcing(kpar,state.IBODY,state.WGHT);

    //Calculate path to next iteration
    if(iforced)
      particle.JUMPF(ds, randoms, particle.getDSMAX());
    else
      particle.JUMP(ds, randoms, particle.getDSMAX());

    //Call tallies with jump collect function
    tallies.run_jump(nhist,kpar,state,ds);

    //Move particle across geometry
    double lastEABS = particle.getEABS();
    particle.move(ds,stepData.softDE,
		  stepData.softX,stepData.softY,stepData.softZ,
		  randoms);    

    //Update step data
    stepData.update(particle); 
    //Call tallies with step collect function
    tallies.run_step(nhist,kpar,state,stepData);

    //Check if the particle must be absorbed at original material and body
    if(state.E < lastEABS){

      double Eprod = particle.annihilationEDep;
      if(state.MAT > 0){
	// run annihilation process
	particle.annihilate(randoms);
      }
      //Set to the particle state the origin material and body
      state.IBODY = particle.lastBody();
      state.MAT = particle.lastMat();

      //Set position to the end of dsef
      double XL,YL,ZL;
      particle.lastPos(XL,YL,ZL);
      state.X = XL+state.U*particle.DSef();
      state.Y = YL+state.V*particle.DSef();
      state.Z = ZL+state.W*particle.DSef();

      tallies.run_localEdep(nhist,kpar,state,state.E+Eprod);
      //Call tallies with end particle collect function
      tallies.run_endPart(nhist,kpar,state);
      return;
    }
      
    //Check if particle cross some interface
    if(particle.NCross() != 0){
      
      //Interface crossed, Update kdet
      kdet = particle.getDET();
      
      //Call tallies with cross interface collect function
      tallies.run_interfCross(nhist,kdet,kpar,state);      

      //Check if material has been changed
      if(particle.lastMat() != state.MAT){
	//Call changed material collect function
	tallies.run_matChange(nhist,kpar,state,particle.lastMat());
	
	if(state.MAT == 0){
	  //New material is vacuum
	  //Call tallies with end particle collect function
	  tallies.run_endPart(nhist,kpar,state);
	  return;
	}
      }
	
      // Particle has been stopped at interface,
      // restart particle state
      particle.START();

      //Check if the particle must be absorbed
      if(absorb(nhist,particle,tallies,randoms))
	return;
    }
    else{

      //No interface crossed, simulate interaction
      if(iforced)
	particle.KNOCKF(de,icol,randoms);
      else
	particle.KNOCK(de,icol,randoms);

      //Call tallies with knock collect function
      tallies.run_knock(nhist,kpar,state,icol);
      //Call tallies with eloss collect function
      tallies.run_localEdep(nhist,kpar,state,de);
	
      //Check if the particle must be absorbed
      if(absorb(nhist,particle,tallies,randoms)){
	return;
      }
    }
  }
}


/*

template <class particleType,
	  class stateType>

void simulatePart(const double nhist,
		  const wrapper_geometry& geometry,
		  particleType& particle,
		  pen_commonTallyCluster& tallies,
		  pen_specificTallyCluster<stateType>& specificTallies,
		  pen_rand& randoms){

  stateType& state = particle.getState();
  unsigned lastMat = state.MAT;
  unsigned lastIBODY = state.IBODY;
  unsigned kdet = geometry.getDET(lastIBODY);
  const pen_KPAR& kpar = particle.kpar;

  //Call tallies with begin particle collect function
  tallies.run_beginPart(nhist,kdet,kpar,state);
  specificTallies.run_beginPart(nhist,kdet,kpar,state);
  
  
  //Check if the particle must be absorbed
  if(absorb(nhist,particle,tallies,specificTallies,randoms))
    return;

  particle.START();
  
  for(;;){
    
    double ds, dsef, dstot, de;
    int ncross, icol;
    double softX, softY, softZ;

    //Check for interaction forcing
    bool iforced = particle.readContext().isForcing(kpar,state.IBODY,state.WGHT);
    
    //Calculate path to next iteration
    if(iforced)      
      particle.JUMPF(ds, randoms, geometry.getDSMAX(state.IBODY));
    else
      particle.JUMP(ds, randoms, geometry.getDSMAX(state.IBODY));    
    
    //Call tallies with jump collect function
    tallies.run_jump(nhist,kpar,state,ds);
    specificTallies.run_jump(nhist,kpar,state,ds);

    //Save origin body index
    lastIBODY = state.IBODY;
    //Move particle across geometry
    double lastEABS = particle.getEABS();
    geometry.step(state, ds, dsef, dstot, ncross);      
    //Add flight time
    if(state.LAGE){particle.dpage(dsef, dstot);}

    //Considere soft energy loss
    if(particle.reqSoftELoss() == 1){
      particle.softEloss(dsef,lastMat,softX,softY,softZ,de,randoms);
      //Call tallies with step collect function
      tallies.run_step(nhist,kpar,state,dsef,dstot,de,
		       softX,softY,softZ,lastIBODY,lastMat);
      specificTallies.run_step(nhist,kpar,state,dsef,dstot,de,
			       softX,softY,softZ,lastIBODY,lastMat);

      //Check if the particle must be absorbed
      if(state.E < lastEABS){

	double Eprod = particle.annihilationEDep;
	if(state.MAT > 0){
	  // run annihilation process
	  particle.annihilate(randoms);
	}	
	//Set to the particle state the origin material and body
	state.IBODY = lastIBODY;
	state.MAT = lastMat;
	tallies.run_localEdep(nhist,kpar,state,state.E+Eprod);
	specificTallies.run_localEdep(nhist,kpar,state,state.E+Eprod);
	//Call tallies with end particle collect function
	tallies.run_endPart(nhist,kpar,state);
	specificTallies.run_endPart(nhist,kpar,state);
	return;
      }
    }
    else{
      //Call tallies with step collect function
      tallies.run_step(nhist,kpar,state,dsef,dstot,
		       0.0,0.0,0.0,0.0,lastIBODY,lastMat);
      specificTallies.run_step(nhist,kpar,state,dsef,dstot,
			       0.0,0.0,0.0,0.0,lastIBODY,lastMat);
    }
    
    //Check if particle cross some interface
    if(ncross != 0){
      //Interface crossed, update kdet
      kdet = geometry.getDET(state.IBODY);
      //Call tallies with cross interface collect function
      tallies.run_interfCross(nhist,kdet,kpar,state);      
      specificTallies.run_interfCross(nhist,kdet,kpar,state);      

      if(lastMat != state.MAT){
	//Call changed material collect function
	tallies.run_matChange(nhist,kpar,state,lastMat);
	specificTallies.run_matChange(nhist,kpar,state,lastMat);

	lastMat = state.MAT;
	
	if(lastMat == 0){
	  //New material is vacuum
	  //Call tallies with end particle collect function
	  tallies.run_endPart(nhist,kpar,state);
	  specificTallies.run_endPart(nhist,kpar,state);

	  return;
	}
      }
	
      // Particle has been stopped at interface,
      // restart particle state
      particle.START();

      //Check if the particle must be absorbed
      if(absorb(nhist,particle,tallies,specificTallies,randoms))
	return;
    }
    else{
      //No interface crossed, simulate interaction
      if(iforced)
	particle.KNOCKF(de,icol,randoms);
      else
	particle.KNOCK(de,icol,randoms);

      //Call tallies with knock collect function
      tallies.run_knock(nhist,kpar,state,icol);
      specificTallies.run_knock(nhist,kpar,state,icol);
      //Call tallies with eloss collect function
      tallies.run_localEdep(nhist,kpar,state,de);
      specificTallies.run_localEdep(nhist,kpar,state,de);

      
      //Check if the particle must be absorbed
      if(absorb(nhist,particle,tallies,specificTallies,randoms))
	return;
    }
  }
}
*/

#endif
