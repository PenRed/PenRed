
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


#include "tallySecondaryGen.hh"

void pen_tallySecondary::flush(){
    
      // Primary particles

  for(unsigned i = 0; i < constants::nParTypes; ++i){
    if(absPrimTmp[i] > 0.0){
      absPrim[i] += absPrimTmp[i];
      absPrim2[i] += absPrimTmp[i]*absPrimTmp[i];
      absPrimTmp[i] = 0.0;
    }
  }

  for(unsigned i = 0; i < constants::nParTypes; ++i){
    if(upBoundPrimTmp[i] > 0.0){
      upBoundPrim[i] += upBoundPrimTmp[i];
      upBoundPrim2[i] += upBoundPrimTmp[i]*upBoundPrimTmp[i];
      upBoundPrimTmp[i] = 0.0;
    }
  }

  for(unsigned i = 0; i < constants::nParTypes; ++i){
    if(downBoundPrimTmp[i] > 0.0){
      downBoundPrim[i] += downBoundPrimTmp[i];
      downBoundPrim2[i] += downBoundPrimTmp[i]*downBoundPrimTmp[i];
      downBoundPrimTmp[i] = 0.0;
    }
  }  

  // Secondary particles
  
  for(unsigned i = 0; i < constants::nParTypes; ++i){
    if(absSecTmp[i] > 0.0){
      absSec[i] += absSecTmp[i];
      absSec2[i] += absSecTmp[i]*absSecTmp[i];
      absSecTmp[i] = 0.0;
    }
  }

  for(unsigned i = 0; i < constants::nParTypes; ++i){
    if(upBoundSecTmp[i] > 0.0){
      upBoundSec[i] += upBoundSecTmp[i];
      upBoundSec2[i] += upBoundSecTmp[i]*upBoundSecTmp[i];
      upBoundSecTmp[i] = 0.0;
    }
  }

  for(unsigned i = 0; i < constants::nParTypes; ++i){
    if(downBoundSecTmp[i] > 0.0){
      downBoundSec[i] += downBoundSecTmp[i];
      downBoundSec2[i] += downBoundSecTmp[i]*downBoundSecTmp[i];
      downBoundSecTmp[i] = 0.0;
    }
  }
}


void pen_tallySecondary::tally_endPart(const unsigned long long /*nhist*/,
				       const pen_KPAR kpar,
				       const pen_particleState& state){

  //Check if is a secondary particle
  if(state.ILB[0] > 1){
    //Is a secondary particle
    //Check if has scaped from geometry system
    if(state.MAT == 0){
      //Check if the particle scape upbound or downbound
      if(state.W > 0.0){
	//Upbound
	upBoundSecTmp[kpar] += state.WGHT;
      }
      else{
	//Downbound
	downBoundSecTmp[kpar] += state.WGHT;	
      }
    }
    //So, has been absorbed
    else{
      absSecTmp[kpar] += state.WGHT;
    }
  }
  else{
    //Is a primary particle
    //Check if has scaped from geometry system
    if(state.MAT == 0){
      //Check if the particle scape upbound or downbound
      if(state.W > 0.0){
	//Upbound
	upBoundPrimTmp[kpar] += state.WGHT;
      }
      else{
	//Downbound
	downBoundPrimTmp[kpar] += state.WGHT;	
      }
    }
    //So, has been absorbed
    else{
      absPrimTmp[kpar] += state.WGHT;
    }
  }
}


void pen_tallySecondary::tally_endHist(const unsigned long long /*nhist*/){
    
    flush();

}

int pen_tallySecondary::configure(const wrapper_geometry& /*geometry*/,
				  const abc_material* const /*materials*/[constants::MAXMAT],
				  const pen_parserSection& /*config*/,
				  const unsigned verbose){
  if(verbose > 1)
    printf("Primary and secondary particles register enabled.\n");

  //Primary
  
  for(unsigned i = 0; i < constants::nParTypes; i++){
    absPrim[i] = 0.0;
    absPrim2[i] = 0.0;
    absPrimTmp[i] = 0.0;
  }

  for(unsigned i = 0; i < constants::nParTypes; i++){
    upBoundPrim[i] = 0.0;
    upBoundPrim2[i] = 0.0;
    upBoundPrimTmp[i] = 0.0;
  }

  for(unsigned i = 0; i < constants::nParTypes; i++){
    downBoundPrim[i] = 0.0;
    downBoundPrim2[i] = 0.0;
    downBoundPrimTmp[i] = 0.0;
  }

  //Secondary
  
  for(unsigned i = 0; i < constants::nParTypes; i++){
    absSec[i] = 0.0;
    absSec2[i] = 0.0;
    absSecTmp[i] = 0.0;
  }

  for(unsigned i = 0; i < constants::nParTypes; i++){
    upBoundSec[i] = 0.0;
    upBoundSec2[i] = 0.0;
    upBoundSecTmp[i] = 0.0;
  }

  for(unsigned i = 0; i < constants::nParTypes; i++){
    downBoundSec[i] = 0.0;
    downBoundSec2[i] = 0.0;
    downBoundSecTmp[i] = 0.0;
  }
    
    //Register data to dump primary particles
    dump.toDump(absPrim,constants::nParTypes);
    dump.toDump(absPrim2,constants::nParTypes);
    dump.toDump(upBoundPrim,constants::nParTypes);
    dump.toDump(upBoundPrim2,constants::nParTypes);
    dump.toDump(downBoundPrim,constants::nParTypes);
    dump.toDump(downBoundPrim2,constants::nParTypes);
    
    
    //Register data to dump seconsary particles
    dump.toDump(absSec,constants::nParTypes);
    dump.toDump(absSec2,constants::nParTypes);
    dump.toDump(upBoundSec,constants::nParTypes);
    dump.toDump(upBoundSec2,constants::nParTypes);
    dump.toDump(downBoundSec,constants::nParTypes);
    dump.toDump(downBoundSec2,constants::nParTypes);
    

  return 0;
}

void pen_tallySecondary::saveData(const unsigned long long nhist) const{

  FILE* fout = nullptr;
  fout = fopen("particleGeneration.dat","w");
  if(fout == nullptr){
    printf("pen_tallySecondary:saveData: Unable to open results file.\n");
    return;
  }

  double invn = 1.0/static_cast<double>(nhist);

  fprintf(fout,"#----------------------------------------------------------\n");
  fprintf(fout,"# PenRed: Particle generation report\n");
  fprintf(fout,"# \n");
  fprintf(fout,"# Primary particles\n");
  fprintf(fout,"# \n");
  fprintf(fout,"                           _______________________________________________ \n");
  fprintf(fout,"                          |    upbound    |   downbound   |    absorbed   |\n");
  fprintf(fout,"__________________________|_______________|_______________|_______________|\n");
  for(unsigned i = 0; i < constants::nParTypes; i++){

    //Check if this particle type has record some value
    if(upBoundPrim[i] == 0.0 && downBoundPrim[i] == 0.0 && absPrim[i] == 0.0)
      continue;
          
    fprintf(fout,"%25.25s | %12.5E  | %12.5E  | %12.5E  |\n",
	    particleName(i),upBoundPrim[i],downBoundPrim[i],absPrim[i]);    
    fprintf(fout,"__________________________|_______________|_______________|_______________|\n");
  }
  fprintf(fout,"\n");

  fprintf(fout,"# Primary particles probabilities\n");
  fprintf(fout,"# \n");
  fprintf(fout,"                           _______________________________________________ \n");
  fprintf(fout,"                          |    upbound    |   downbound   |    absorbed   |\n");
  fprintf(fout,"__________________________|_______________|_______________|_______________|\n");
  for(unsigned i = 0; i < constants::nParTypes; i++){

    //Check if this particle type has record some value
    if(upBoundPrim[i] == 0.0 && downBoundPrim[i] == 0.0 && absPrim[i] == 0.0)
      continue;

    double sigmaUp   = invn*sqrt(fabs(upBoundPrim2[i])-pow(upBoundPrim[i],2)*invn);
    double sigmaDown = invn*sqrt(fabs(downBoundPrim2[i])-pow(downBoundPrim[i],2)*invn);;
    double sigmaAbs  = invn*sqrt(fabs(absPrim2[i])-pow(absPrim[i],2)*invn);;
    
    fprintf(fout,"%25.25s | %12.5E  | %12.5E  | %12.5E  |\n",
	    particleName(i),upBoundPrim[i]*invn,downBoundPrim[i]*invn,absPrim[i]*invn);    
    fprintf(fout,"%25.25s |  +- %8.2E  |  +- %8.2E  |  +- %8.2E  |\n",
	    " ",2.0*sigmaUp,2.0*sigmaDown,2.0*sigmaAbs);
    fprintf(fout,"__________________________|_______________|_______________|_______________|\n");
  }
  fprintf(fout,"\n");
  
  
  fprintf(fout,"# Secondary particles probabilities\n");
  fprintf(fout,"# \n");    
  fprintf(fout,"                           _______________________________________________ \n");
  fprintf(fout,"                          |    upbound    |   downbound   |    absorbed   |\n");
  fprintf(fout,"__________________________|_______________|_______________|_______________|\n");
  for(unsigned i = 0; i < constants::nParTypes; i++){

    //Check if this particle type has record some value
    if(upBoundSec[i] == 0.0 && downBoundSec[i] == 0.0 && absSec[i] == 0.0)
      continue;
    
    //Calculate final values
    double qUp = upBoundSec[i]*invn;
    double sigmaUp = (upBoundSec2[i]*invn - qUp*qUp)*invn;
    sigmaUp = sqrt(sigmaUp > 0.0 ? sigmaUp : 0.0);

    double qDown = downBoundSec[i]*invn;
    double sigmaDown = (downBoundSec2[i]*invn - qDown*qDown)*invn;
    sigmaDown = sqrt(sigmaDown > 0.0 ? sigmaDown : 0.0);

    double qAbs = absSec[i]*invn;
    double sigmaAbs = (absSec2[i]*invn - qAbs*qAbs)*invn;
    sigmaAbs = sqrt(sigmaAbs > 0.0 ? sigmaAbs : 0.0);
    
    fprintf(fout,"%25.25s | %12.5E  | %12.5E  | %12.5E  |\n",particleName(i),qUp,qDown,qAbs);
    fprintf(fout,"%25.25s |  +- %8.2E  |  +- %8.2E  |  +- %8.2E  |\n",
	    " ",2.0*sigmaUp,2.0*sigmaDown,2.0*sigmaAbs);
    
    fprintf(fout,"__________________________|_______________|_______________|_______________|\n");
  }    
  fprintf(fout,"\n");

  fclose(fout);
}


int pen_tallySecondary::sumTally(const pen_tallySecondary& tally){
    
    //Primary particles
    //******************
    
    for(unsigned k = 0; k < constants::nParTypes; ++k)
    {
        absPrim[k] += tally.absPrim[k];
    
    }
    
    
    for(unsigned k = 0; k < constants::nParTypes; ++k)
    {
        absPrim2[k] += tally.absPrim2[k];
    
    }
    
    
    for(unsigned k = 0; k < constants::nParTypes; ++k)
    {
        upBoundPrim[k] += tally.upBoundPrim[k];
    
    }
    
        
    for(unsigned k = 0; k < constants::nParTypes; ++k)
    {
        upBoundPrim2[k] += tally.upBoundPrim2[k];
    
    }
    
    
    
    for(unsigned k = 0; k < constants::nParTypes; ++k)
    {
        downBoundPrim[k] += tally.downBoundPrim[k];
    
    }
    
        
    for(unsigned k = 0; k < constants::nParTypes; ++k)
    {
        downBoundPrim2[k] += tally.downBoundPrim2[k];
    
    }
    
    
    
    
    //Secondary particles
    //******************
    
    for(unsigned k = 0; k < constants::nParTypes; ++k)
    {
        absSec[k] += tally.absSec[k];
    
    }
    
    
    for(unsigned k = 0; k < constants::nParTypes; ++k)
    {
        absSec2[k] += tally.absSec2[k];
    
    }
    
    
    for(unsigned k = 0; k < constants::nParTypes; ++k)
    {
        upBoundSec[k] += tally.upBoundSec[k];
    
    }
    
        
    for(unsigned k = 0; k < constants::nParTypes; ++k)
    {
        upBoundSec2[k] += tally.upBoundSec2[k];
    
    }
    
    
    
    for(unsigned k = 0; k < constants::nParTypes; ++k)
    {
        downBoundSec[k] += tally.downBoundSec[k];
    
    }
    
        
    for(unsigned k = 0; k < constants::nParTypes; ++k)
    {
        downBoundSec2[k] += tally.downBoundSec2[k];
    
    }
    
    
    return 0;
}




REGISTER_COMMON_TALLY(pen_tallySecondary)
