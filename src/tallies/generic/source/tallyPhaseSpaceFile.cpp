
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


#include "tallyPhaseSpaceFile.hh"  

std::vector<
  std::pair<std::string,std::shared_ptr<pen_splittedFile>>
  > pen_tallyPhaseSpaceFile::splittedFiles;

std::mutex pen_tallyPhaseSpaceFile::SFlock;


void pen_tallyPhaseSpaceFile::flush(){
}

void pen_tallyPhaseSpaceFile::tally_lastHist(const double lastHist){
  //Set last hist to phase space file
  psf.setLast(lastHist);
}

void pen_tallyPhaseSpaceFile::tally_interfCross(const double nhist,
						const unsigned kdet,
						const pen_KPAR kpar,
						const pen_particleState& state){

  if((int)kdet == detector){
    if(!inside && state.E >= emin && state.E <= emax)
      store(nhist,kpar,state);

    inside = true;
  }
  else{
    inside = false;
  }
}

void pen_tallyPhaseSpaceFile::tally_matChange(const double nhist,
					      const pen_KPAR kpar,
					      const pen_particleState& state,
					      const unsigned /*prevMat*/){
  if((int)state.MAT == material){
    if(state.E >= emin && state.E <= emax)
      store(nhist,kpar,state);
  }
}

void pen_tallyPhaseSpaceFile::tally_move2geo(const double nhist,
					     const unsigned kdet,
					     const pen_KPAR kpar,
					     const pen_particleState& state,
					     const double /*dsef*/,
					     const double /*dstot*/){

  if((int)kdet == detector || (int)state.MAT == material){
    inside = true;
    if(state.E >= emin && state.E <= emax)
      store(nhist,kpar,state);
  }else{
    inside = false;
  }
  
}

void pen_tallyPhaseSpaceFile::tally_beginPart(const double /*nhist*/,
					      const unsigned kdet,
					      const pen_KPAR /*kpar*/,
					      const pen_particleState& /*state*/){
  if((int)kdet == detector){
    inside = true;
  }else{
    inside = false;
  }
}

void pen_tallyPhaseSpaceFile::tally_endSim(const double /*nhist*/){
  //Dump residual states
  dump(true);
  pSF.get()->flush(getThread());
}

int pen_tallyPhaseSpaceFile::configure(const wrapper_geometry& /*geometry*/,
				 const abc_material* const /*materials*/[constants::MAXMAT],
				 const pen_parserSection& config,
				 const unsigned verbose){
    
  int err;

  //Clear phase space file
  psf.clear();
  
  //Try to create a splitted file with this tally name

  //First, lock splitted files vector
  std::lock_guard<std::mutex> guard(SFlock);

  //Next, check if a splitted file with our name exists
  bool found = false;
  size_t nSF = splittedFiles.size();
  for(unsigned i = 0; i < nSF; i++){
    if(splittedFiles[i].first.compare(readName()) == 0){
      //This splitted file has already created. Get its pointer
      pSF = splittedFiles[i].second;
      found = true;
      break;
    }
  }

  if(!found){
    //Create the required splitted file

    std::string filename = readName();
  // ******************************* MPI ************************************ //
#ifdef _PEN_USE_MPI_

    //Get process rank
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    //Add MPI rank to filename 
    filename = std::string("MPI") + std::to_string(rank) +
      std::string("-") + filename;
#endif
  // ***************************** MPI END ********************************** //
    
    splittedFiles.push_back(std::make_pair(readName(),std::make_shared<
					   pen_splittedFile>(filename.c_str(),
							     true
							     )
					   )
			    );

    //Store a shared pointer
    pSF = splittedFiles[nSF].second;
  }
  
  //try to create a partition for our thread
  err = pSF.get()->createPartition(getThread());
  if(err != SPLITTED_FILE_SUCCESS){
    if(verbose > 0){
      printf("PhaseSpaceFile:configure: Error: Unable to create a partition for tally %s thread %u.\n",readName().c_str(),getThread());
      printf("                 Error code: %d\n",err);
    }
    return -1;
  }

  
 // Detector
 //*************** 
  int auxkDet;
  err = config.read("detector", auxkDet);
  if(err != INTDATA_SUCCESS){
    if(verbose > 1){
      printf("No detector specified. Integrer expected\n");
    }
  }
  else if(auxkDet <= 0){
    if(verbose > 0){
      printf("PhaseSpaceFile:configure: Error: Detector index must be greater than 0.\n");
    }
    return -2;    
  }
  else{
    detector = auxkDet;
    if(verbose > 1){
      printf("Detector:\n");
      printf(" %u \n\n",detector);
    }    
  }


 // Detector material
 //******************* 
  int auxMat;
  err = config.read("material", auxMat);
  if(err != INTDATA_SUCCESS){
    if(verbose > 1){
      printf("No detector material specified. Integrer expected\n");
    }
  }
  else if(auxMat <= 0){
    if(verbose > 0){
      printf("PhaseSpaceFile:configure: Error: Detector material index must be greater than 0.\n");
    }
    return -3;    
  }
  else{
    material = auxMat;
    if(verbose > 1){
      printf("Detector material:\n");
      printf(" %u \n\n",auxMat);
    }
  }

  if(material <= 0 && detector <= 0){
    if(verbose > 0){
      printf("PhaseSpaceFile:configure: Error: No detector nor material specified.\n");      
    }
    return -4;
  }

  if(material > 0 && detector > 0){
    if(verbose > 0){
      printf("PhaseSpaceFile:configure: Error: Detector and material specified. Only one can be set as SPF detector.\n");      
    }
    return -4;
  }

 //Energies
 //****************

  // Minimum energy
  //***************    
  err = config.read("emin", emin);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("PhaseSpaceFile:configure: Error: Unable to read 'emin' in configuration. Double expected\n");
    }
    return -5;
  }

  // Maximum energy
  //***************    
    
  err = config.read("emax", emax);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("PhaseSpaceFile:configure: Error: Unable to read 'emax' in configuration. Double expected\n");
    }
    return -6;
  }

  if(emin >= emax){
    if(verbose > 0){
      printf("PhaseSpaceFile:configure: Error: Minimum energy must be lower than maximum one.\n");
      printf("            Emin: %14.5E\n",emin);
      printf("            Emax: %14.5E\n",emax);
    }
    return -7;
  }
  
  if(verbose > 1){
    printf("PSF limits [Emin,Emax] (eV):\n");
    printf(" %12.5E %12.5E \n\n",emin,emax);
  }
  
  return 0;
}


void pen_tallyPhaseSpaceFile::saveData(const double /*nhist*/) const{
  
}
int pen_tallyPhaseSpaceFile::sumTally(const pen_tallyPhaseSpaceFile& /*tally*/){return 0;}


REGISTER_COMMON_TALLY(pen_tallyPhaseSpaceFile, PSF)
