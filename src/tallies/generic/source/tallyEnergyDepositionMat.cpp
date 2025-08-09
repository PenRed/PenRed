
//
//
//    Copyright (C) 2019-2021 Universitat de València - UV
//    Copyright (C) 2019-2021 Universitat Politècnica de València - UPV
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
//        vicente.gimenez@uv.es (Vicente Giménez Gómez)
//        sanolgi@upvnet.upv.es (Sandra Oliver Gil)
//    
//


#include "tallyEnergyDepositionMat.hh"

pen_EdepMat::pen_EdepMat() : pen_genericTally( USE_LOCALEDEP |
					       USE_BEGINPART |
					       USE_SAMPLEDPART |
					       USE_STEP |
					       USE_ENDHIST |
					       USE_MOVE2GEO),
			     nmat(0)
{
  setResultsGenerator<0>
    ([this](const unsigned long long nhists) -> penred::measurements::results<double, 1>{

      if(nmat <= 0){
	penred::measurements::results<double, 1> results;
	results.description = "Error: No material provided. "
	  "Configure the tally before getting results.";
	return results;
      }

      double invn = 1.0/static_cast<double>(nhists);
	
      //Create results
      penred::measurements::results<double, 1> results;
      results.initFromLists({static_cast<unsigned long>(nmat)},
			    {penred::measurements::limitsType(0.0, static_cast<double>(nmat))});
	  
      results.description = "PenRed: Material energy deposition report.\n\n";
  
      results.setDimHeader(0, "Material");
      results.setValueHeader("Energy (eV/hist)");

      for(int i = 0; i < nmat; i++)
	{
	  double q  = edep[i]*invn;
	  double q2 = edep2[i]*invn;
	  double sigma = (q2-(q*q))*invn;
	  if(sigma > 0.0){ sigma = sqrt(sigma);}
	  else{sigma = 0.0;}

	  results.data[i] = q;
	  results.sigma[i] = sigma;
	}

      return results;
    });      
}

void pen_EdepMat::flush()
{
    for(int i = 0; i < nmat; i++) //For all existing materials
    {
      if(edptmp[i] == 0.0){continue;}   //Skip void counters
      edep[i] += edptmp[i];    //Transfer temporary counter to mean and variance
      edep2[i] +=  edptmp[i]*edptmp[i];
      edptmp[i] = 0.0;                  //Clear counter to start a fresh history
    }
}

 
void pen_EdepMat::tally_localEdep(const unsigned long long /*nhist*/,
				  const pen_KPAR /*kpar*/,
				  const pen_particleState& state,
				  const double dE){

  //Energy deposited at material
  edptmp[state.MAT-1] += dE*state.WGHT; 
}

void pen_EdepMat::tally_beginPart(const unsigned long long /*nhist*/,
				  const unsigned /*kdet*/,
				  const pen_KPAR /*kpar*/,
				  const pen_particleState& state){

  //Extract energy from material to create new particle
  edptmp[state.MAT-1] -= state.E*state.WGHT;

}

void pen_EdepMat::tally_sampledPart(const unsigned long long /*nhist*/,
				    const unsigned long long /*dhist*/,
				    const unsigned /*kdet*/,
				    const pen_KPAR /*kpar*/,
				    const pen_particleState& state){

  //Add energy to material to compensate substracted one when
  //beginPartwill be called.
  if(state.MAT > 0){
    //Ensure that primary particle has not been created at void volume.
    //If this happens, energy will be added at move2geo call.
    edptmp[state.MAT-1] += state.E*state.WGHT;
  }
}

void pen_EdepMat::tally_step(const unsigned long long /*nhist*/,
			     const pen_KPAR /*kpar*/,
			     const pen_particleState& state,
			     const tally_StepData& stepData){

  if(stepData.softDE == 0.0){return;}  //Nothing to deposit.
  //Energy deposited at material
  edptmp[stepData.originMAT-1] += stepData.softDE*state.WGHT;
}

void pen_EdepMat::tally_move2geo(const unsigned long long /*nhist*/,
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
    edptmp[state.MAT-1] += state.E*state.WGHT;
  }  
}



void pen_EdepMat::tally_endHist(const unsigned long long /*nhist*/){

    flush();
}

int pen_EdepMat::configure(const wrapper_geometry& /*geometry*/,
			   const abc_material* const /*materials*/[constants::MAXMAT],
			   const pen_parserSection& /*config*/,
			   const unsigned verbose){

  nmat = constants::MAXMAT;
    
  //Clear counters:
  for(unsigned int j = 0; j < constants::MAXMAT; j++){
    edptmp[j] = 0.0;
    edep[j] = 0.0;
    edep2[j] = 0.0;
  }

  if(verbose > 1){
    printf("Detection material: %d\n",nmat);
  }  

  //Register data to dump
  dump.toDump(edptmp,nmat);
  dump.toDump(edep,nmat);
  dump.toDump(edep2,nmat);
  dump.toDump(&nmat,1);
  
  return 0;
}

void pen_EdepMat::saveData(const unsigned long long nhist) const{
    
  FILE* out;
  int i;
  double q, q2, sigma, invn;
    
  out = fopen("materialEnergyDeposition.dat", "w");
  if(out == NULL){
    
    printf(" *********************************************\n");
    printf(" EdepMat:saveData:ERROR: cannot open output data file\n");
    printf(" *********************************************\n");
    return;
  }
    
  fprintf(out, "#------------------------------------------------------------\n");
  fprintf(out, "# PenRed: Material energy deposition\n");
  fprintf(out, "# Units are eV per history\n");
  fprintf(out, "#\n");
  fprintf(out, "# Material : Energy (eV/hist) : +-2sigma\n");  
  
  invn = 1.0/static_cast<double>(nhist);
  for(i = 0; i < nmat; i++)
    {
      q  = edep[i]*invn;
      q2 = edep2[i]*invn;
      sigma = (q2-(q*q))*invn;
      if(sigma > 0.0){ sigma = sqrt(sigma);}
      else{sigma = 0.0;}

      fprintf(out, "    %3d       %12.5E      %8.1E\n", i+1,q,2.0*sigma);
    }
  fclose(out);
}

int pen_EdepMat::sumTally(const pen_EdepMat& tally){

  if(nmat != tally.nmat)
    return -1;
  
  for(int i = 0; i < nmat; ++i){
    edep[i] += tally.edep[i];
  }
  for(int i = 0; i < nmat; ++i){
    edep2[i] += tally.edep2[i];
  }

  return 0;
}

REGISTER_COMMON_TALLY(pen_EdepMat)
