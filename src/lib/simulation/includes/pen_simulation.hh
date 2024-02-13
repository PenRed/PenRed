
//
//
//    Copyright (C) 2019-2024 Universitat de València - UV
//    Copyright (C) 2019-2024 Universitat Politècnica de València - UPV
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
//        sanolgi@upvnet.upv.es
//        vicente.gimenez@uv.es
//    
//

 
#ifndef __PEN_SIMULATION__
#define __PEN_SIMULATION__

#include "PenRed.hh"
#include "pen_samplers.hh"
#include "pen_tallies.hh"
#include "pen_geometries.hh"
#include "pen_vr.hh"

namespace penred{

  namespace simulation{

    typedef std::function<void(const pen_particleState&,const pen_KPAR,const int)> simFuncType;

    enum finishTypes{
      NO_CONDITION = 0,
      DETECTOR_REACHED = 1,
    };

    //Move to geo function family:
    //
    // This function check if a sampled particle has been created in the void.
    // In that case, it tries to move the particle to reach a non void region
    //
    template <class particleType>
    bool move2geo(const unsigned long long& nhist,
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
    bool move2geo(const unsigned long long& nhist,
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
    bool move2geo(particleType& particle){

      pen_particleState& state = particle.getBaseState();
      //Check if particle has been generated at vacuum
      if(state.MAT == 0){
	//Move it to geometry system
	particle.jumpVolume();

	if(state.MAT == 0){
	  return false;
	}
      }
      else{
	particle.updateMat();
	particle.updateBody();
      }
      return true;
    }

    //Absorb function family:
    //
    // These functions check if the particle must be absrobed
    // in its current state
    //
    template <class particleType>
    inline bool absorb(const unsigned long long& nhist,
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
    inline bool absorb(const unsigned long long& nhist,
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
    inline bool absorb(particleType& particle,
		       pen_rand& penRand){

      const pen_particleState& state = particle.readBaseState();

      //Check if particle must be absorbed
      if(state.E < particle.getEABS()){ 

	//Check if the particle has already been annihilated
	if(state.E > 1.0e-6){
	  // run annihilation process
	  particle.annihilate(penRand);
	}
	return true;
      }

      //Particle can't be absorbed yet
      return false;
    }

    //Simulate particle function family
    //
    // These functions simulates a particle already situated in
    // a non void region. This shoud be satisfied if the particle
    // to be simulated is a secondary produced by a previous particle
    // simulation, or the particle to be simulated has been sampled and
    // "moved" using the move2geo function. 
    //
    template <class particleType>
    void simulatePart(const unsigned long long& nhist,
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

	    // Particle has been stopped at interface,
	    // restart particle state
	    particle.START();

	    //Check if the particle must be absorbed
	    if(absorb(nhist,particle,tallies,randoms))
	      return;

	    //VR
	    double deVR = particle.vr_matChange(nhist,randoms,2);	
	    //VR
	    deVR += particle.vr_interfCross(nhist,randoms,2);

	    //Tally the deposited energy by VR
	    if(deVR > 0.0){
	      //We must compensate the particle weight to score correctly
	      //the energy stored in the secondary stacks by the VR call
	      double originWGHT = state.WGHT;
	      state.WGHT = 1.0;
	      tallies.run_localEdep(nhist,kpar,state,deVR);
	      state.WGHT = originWGHT;
	    }

	    //Check if the particle must be absorbed because VR
	    if(absorb(nhist,particle,tallies,randoms))
	      return;
	
	  }
	  else{
	    // Particle has been stopped at interface,
	    // restart particle state
	    particle.START();

	    //Check if the particle must be absorbed
	    if(absorb(nhist,particle,tallies,randoms))
	      return;
      
	    //VR
	    double deVR = particle.vr_interfCross(nhist,randoms,2);

	    //Tally the deposited energy by VR
	    if(deVR > 0.0){
	      //We must compensate the particle weight to score correctly
	      //the energy stored in the secondary stacks by the VR call
	      double originWGHT = state.WGHT;
	      state.WGHT = 1.0;
	      tallies.run_localEdep(nhist,kpar,state,deVR);
	      state.WGHT = originWGHT;
	    }

	    //Check if the particle must be absorbed because VR
	    if(absorb(nhist,particle,tallies,randoms))
	      return;	
	  }
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

    template <class particleType>
    int simulatePartCond(const unsigned long long& nhist,
			 particleType& particle,
			 pen_rand& randoms,
			 const finishTypes& finishType,
			 const unsigned& finishValue){
  
      pen_particleState& state = particle.getBaseState();
      const pen_KPAR kpar = particle.kpar;

      //Check if particle scaped from geometry
      if(state.MAT == 0){
	return 0;
      }  
  
      //Check if the particle must be absorbed
      if(absorb(particle,randoms))
	return 0;

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

	//Move particle across geometry
	double lastEABS = particle.getEABS();
	particle.move(ds,stepData.softDE,
		      stepData.softX,stepData.softY,stepData.softZ,
		      randoms);    

	//Update step data
	stepData.update(particle);

	//Check if the particle must be absorbed at original material and body
	if(state.E < lastEABS){

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

	  return 0;
	}
      
	//Check if particle cross some interface
	if(particle.NCross() != 0){
      
	  if(finishType == 1){      
	    //Finish on detector reached condition
	    unsigned kdet = particle.getDET();
	    if(finishValue == kdet)
	      return 1;
	  }
	
	  //Check if material has been changed
	  if(particle.lastMat() != state.MAT){
	
	    if(state.MAT == 0){
	      //New material is vacuum
	      return 0;
	    }

	    // Particle has been stopped at interface,
	    // restart particle state
	    particle.START();

	    //Check if the particle must be absorbed
	    if(absorb(particle,randoms))
	      return 0;

	    //VR
	    particle.vr_matChange(nhist,randoms,2);	
	    //VR
	    particle.vr_interfCross(nhist,randoms,2);

	    //Check if the particle must be absorbed because VR
	    if(absorb(particle,randoms))
	      return 0;
	
	  }
	  else{
	    // Particle has been stopped at interface,
	    // restart particle state
	    particle.START();

	    //Check if the particle must be absorbed
	    if(absorb(particle,randoms))
	      return 0;
      
	    //VR
	    particle.vr_interfCross(nhist,randoms,2);

	    //Check if the particle must be absorbed because VR
	    if(absorb(particle,randoms))
	      return 0;
	  }
	}
	else{

	  //No interface crossed, simulate interaction
	  if(iforced)
	    particle.KNOCKF(de,icol,randoms);
	  else
	    particle.KNOCK(de,icol,randoms);
	
	  //Check if the particle must be absorbed
	  if(absorb(particle,randoms)){
	    return 0;
	  }
	}
      }
    }

    //Simulate stack function family:
    //
    // These functions simulates "n" seconday particles stored
    // in the particle main stack belonging the particle "p" 
    //
    template<class particleType>
    inline unsigned simulateStack(const unsigned long long& hist,
				  particleType& p,
				  const unsigned& n,
				  pen_commonTallyCluster& tallies,
				  pen_rand& randoms){
      //Check the number of particles to simulate
      const unsigned toSim = std::min(n,p.nStacked());
  
      for(unsigned i = 0; i < toSim; ++i){
	//Get next state from particle stack
	p.setStateFromStack();

	//Get kdet
	unsigned kdet = p.getDET();
    
	//VR
	p.vr_particleStack(hist,randoms,2);

	//Simulate particle
	tallies.run_beginPart(hist,kdet,p.getKpar(),p.readState());
	simulatePart(hist,p,tallies,randoms);    
      }

      return toSim;
    }

    template<class particleType>
    inline unsigned simulateStackCond(const unsigned long long& hist,
				      particleType& p,
				      const unsigned& n,
				      pen_rand& randoms,
				      const finishTypes& finishType,
				      const unsigned& finishValue,
				      simFuncType& f){
  
      //Check the number of particles to simulate
      const unsigned toSim = std::min(n,p.nStacked());
  
      for(unsigned i = 0; i < toSim; ++i){
	//Get next state from particle stack
	p.setStateFromStack();

	//VR
	p.vr_particleStack(hist,randoms,2);

	//Simulate particle
	int ret = simulatePartCond(hist,p,randoms,finishType,finishValue);
	f(p.readState(), p.getKpar(), ret);
      }

      return toSim;
    }

    //Simulate stacks N function family:
    //
    // These functions simulates "n" seconday particles of each
    // particle main stack belonging to each particle provided in the argument
    //
    template<class particleType>
    inline unsigned simulateStacksN(const unsigned long long& hist,
				    const unsigned& n,
				    pen_commonTallyCluster& tallies,
				    pen_rand& randoms,
				    particleType& p){
  
      //Try to simulate "n" particles in the stack
      simulateStack(hist, p, n, tallies, randoms);

      //Return the number of remaining particles in the stack
      //Notice that this is the last particle stack to process
      return p.nStacked();
    }

    template<class particleType, class... otherParticleTypes>
    inline typename std::enable_if<sizeof...(otherParticleTypes) >= 1, unsigned>::type
    simulateStacksN(const unsigned long long& hist,
		    const unsigned& n,
		    pen_commonTallyCluster& tallies,
		    pen_rand& randoms,
		    particleType& p,
		    otherParticleTypes&... others){
  
      //Try to simulate "n" particles in the stack
      simulateStack(hist, p, n, tallies, randoms);

      //Simulate other particle stacks
      unsigned remainingOthers = simulateStacksN(hist,n,tallies,randoms,others...);

      //Return the number of remaining particles in the stack
      //plus the other stacks after simulating "n" particles
      return p.nStacked() + remainingOthers;
    }

    template<class particleType>
    inline unsigned simulateStacksNCond(const unsigned long long& hist,
					const unsigned& n,
					pen_rand& randoms,
					const finishTypes& finishType,
					const unsigned& finishValue,
					simFuncType& f,
					particleType& p){
  
      //Try to simulate "n" particles in the stack
      simulateStackCond(hist, p, n, randoms, finishType, finishValue, f);

      //Return the number of remaining particles in the stack
      //Notice that this is the last particle stack to process
	return p.nStacked();      
    }

    template<class particleType, class... otherParticleTypes>
    inline typename std::enable_if<sizeof...(otherParticleTypes) >= 1, unsigned>::type
    simulateStacksNCond(const unsigned long long& hist,
			const unsigned& n,
			pen_rand& randoms,
			const finishTypes& finishType,
			const unsigned& finishValue,
			simFuncType& f,
			particleType& p,
			otherParticleTypes&... others){
  
      //Try to simulate "n" particles in the stack
      simulateStackCond(hist, p, n, randoms, finishType, finishValue, f);

      //Simulate other particle stacks
      unsigned remainingOthers = simulateStacksNCond(hist,n,randoms,
						     finishType,finishValue,
						     f,others...);
      
      //Return the number of remaining particles in the stack
      //plus the other stacks after simulating "n" particles
      return p.nStacked() + remainingOthers;
      
    }

    
    //Simulate stacks function family:
    //
    // These functions simulates all the particles in the main
    // stack belonging to each particle provided in the argument
    //
    template<class... particleTypes>
    inline void simulateStacks(const unsigned long long& hist,
			       pen_commonTallyCluster& tallies,
			       pen_rand& randoms,
			       particleTypes&... particles){

      unsigned remaining;
      do{
	//Simulate stack particles in chunks of 100 particles
	remaining = simulateStacksN(hist, 100, tallies, randoms, particles...);
      }while(remaining > 0);
  
    }

    template<class... particleTypes>
    inline void simulateStacksCond(const unsigned long long& hist,
				   pen_rand& randoms,
				   const finishTypes& finishType,
				   const unsigned& finishValue,
				   simFuncType& f,
				   particleTypes&... particles){

      unsigned remaining;
      do{
	//Simulate stack particles in chunks of 100 particles
	remaining = simulateStacksNCond(hist, 100, randoms,
					finishType, finishValue, f, particles...);
      }while(remaining > 0);
  
    }

    //Simulate shower function family
    //
    // These functions simulate a complete shower, starting from a primary
    // particle and ending with all the secondary generated particles
    //
    template<class primaryType, class... secondaryTypes>
    void simulateShower(const unsigned long long& hist,
			pen_commonTallyCluster& tallies,
			pen_rand& randoms,
			primaryType& p,
			secondaryTypes&... secondary){

      //Init Page
      p.page0();

      //Try to move the generated particle to geometry system
      if(move2geo(hist,p,tallies)){

	//Get detector ID
	unsigned kdet = p.getDET();	
	
	//Tally particle begins
	tallies.run_beginPart(hist,kdet,p.getKpar(),p.readState());
    
	//Simulate particle
	simulatePart(hist,p,tallies,randoms);

	//Simulate secondary particles
	simulateStacks(hist, tallies, randoms, p, secondary...);
      }

  
  
    }

    template<class primaryType, class... secondaryTypes>
    void simulateShowerCond(const unsigned long long& hist,
			    pen_rand& randoms,
			    const finishTypes& finishType,
			    const unsigned& finishValue,
			    simFuncType& f,
			    primaryType& p,
			    secondaryTypes&... secondary){

      //Init Page
      p.page0();

      //Try to move the generated particle to geometry system
      if(move2geo(p)){

	//Simulate particle
	int ret = simulatePartCond(hist,p,randoms,finishType,finishValue);
	//Call tally function
	f(p.readState(), p.getKpar(), ret);

	//Simulate secondary particles
	simulateStacksCond(hist, randoms, finishType, finishValue, f, p, secondary...);
      }

  
  
    }

  };
};
#endif
