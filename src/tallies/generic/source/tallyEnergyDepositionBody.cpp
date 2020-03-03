
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


#include "tallyEnergyDepositionBody.hh"

void pen_EdepBody::flush()
{
    for(int i = 0; i < nBody; i++) //For all existing bodies
    {
      if(edptmp[i] == 0.0){continue;}   //Skip void counters
      edep[i] += edptmp[i];    //Transfer temporary counter to mean and variance
      edep2[i] +=  edptmp[i]*edptmp[i];
      edptmp[i] = 0.0;                  //Clear counter to start a fresh history
    }
}

 
void pen_EdepBody::tally_localEdep(const double /*nhist*/,
				  const pen_KPAR /*kpar*/,
				  const pen_particleState& state,
				  const double dE){

  //Energy deposited at body
  edptmp[state.IBODY] += dE*state.WGHT; 
}

void pen_EdepBody::tally_beginPart(const double /*nhist*/,
				  const unsigned /*kdet*/,
				  const pen_KPAR /*kpar*/,
				  const pen_particleState& state){

  //Extract energy from body to create new particle
  edptmp[state.IBODY] -= state.E*state.WGHT;

}

void pen_EdepBody::tally_beginHist(const double /*nhist*/,
				  const unsigned /*kdet*/,
				  const pen_KPAR /*kpar*/,
				  const pen_particleState& state){

  //Add energy to body to compensate substracted one when
  //beginPartwill be called.
  if(state.MAT > 0){
    //Ensure that primary particle has not been created at void volume.
    //If this happens, energy will be added at move2geo call.
    edptmp[state.IBODY] += state.E*state.WGHT;
  }
}

void pen_EdepBody::tally_step(const double /*nhist*/,
			      const pen_KPAR /*kpar*/,
			      const pen_particleState& state,
			      const tally_StepData& stepData){

  if(stepData.softDE == 0.0){return;}  //Nothing to deposit.
  //Energy deposited at body
  edptmp[stepData.originIBODY] += stepData.softDE*state.WGHT;
}

void pen_EdepBody::tally_move2geo(const double /*nhist*/,
				 const unsigned /*kdet*/,
				 const pen_KPAR /*kpar*/,
				 const pen_particleState& state,
				 const double /*dsef*/,
				 const double /*dstot*/){
  //Primary particle has been created at void volume. Check if
  //after step call particle reached non void volume.
  if(state.MAT > 0){
    //Non void volume reached. Add particle energy to compensate
    //substracted one when beginPart will be called.
    edptmp[state.IBODY] += state.E*state.WGHT;
  }  
}



void pen_EdepBody::tally_endHist(const double /*nhist*/){

    flush();
}

int pen_EdepBody::configure(const wrapper_geometry& /*geometry*/,
			   const abc_material* const /*materials*/[pen_geoconst::NB],
			   const pen_parserSection& config,
			   const unsigned verbose){
    
  int err;     
  err = config.read("nBody", nBody);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("EdepBody:configure:unable to read 'nBody' in configuration. Integrer expected");
    }
    return -1;
  }
    
  //Clear counters:
  for(unsigned int j = 0; j < pen_geoconst::NB; j++){
    edptmp[j] = 0.0;
    edep[j] = 0.0;
    edep2[j] = 0.0;
  }

  if(verbose > 1){
    printf("Detection material: %d\n",nBody);
  }  

  //Register data to dump
  dump.toDump(edptmp,nBody);
  dump.toDump(edep,nBody);
  dump.toDump(edep2,nBody);
  dump.toDump(&nBody,1);
  
  return 0;
}

void pen_EdepBody::saveData(const double nhist) const{
    
  FILE* out;
  double q, q2, sigma, invn;
    
  out = fopen("bodyEnergyDeposition.dat", "w");
  if(out == NULL){
    
    printf(" *********************************************\n");
    printf(" EdepBody:saveData:ERROR: cannot open output data file\n");
    printf(" *********************************************\n");
    return;
  }
    
  fprintf(out, "#------------------------------------------------------------\n");
  fprintf(out, "# PenRed: Body energy deposition report\n");
  fprintf(out, "# Units are eV per history\n");
  fprintf(out, "#\n");
  fprintf(out, "# Body : Energy (eV/hist) : +-2sigma\n");  
  
  invn = 1.0/nhist;
  for(int i = 0; i < nBody; i++)
    {
      q  = edep[i]*invn;
      q2 = edep2[i]*invn;
      sigma = (q2-(q*q))*invn;
      if(sigma > 0.0){ sigma = sqrt(sigma);}
      else{sigma = 0.0;}

      fprintf(out, " %3d      %12.5E     %8.1E\n", i+1,q,2.0*sigma);
    }
  fclose(out);
}

int pen_EdepBody::sumTally(const pen_EdepBody& tally){

  if(nBody != tally.nBody)
    return -1;
  
  for(int i = 0; i < nBody; ++i){
    edep[i] += tally.edep[i];
  }
  for(int i = 0; i < nBody; ++i){
    edep2[i] += tally.edep2[i];
  }

  return 0;
}

REGISTER_COMMON_TALLY(pen_EdepBody, EDEP_BODY)
