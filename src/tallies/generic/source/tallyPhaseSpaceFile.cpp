
//
//
//    Copyright (C) 2019-2025 Universitat de València - UV
//    Copyright (C) 2019-2025 Universitat Politècnica de València - UPV
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
//        vicent.gimenez.alventosa@gmail.com
//        sanolgi@upvnet.upv.es
//        vicente.gimenez@uv.es
//    
//


#include "tallyPhaseSpaceFile.hh"  

void pen_tallyPhaseSpaceFile::flush(){
}

void pen_tallyPhaseSpaceFile::tally_lastHist(const unsigned long long lastHist){
  //Set last hist to phase space file
  psf.setLast(lastHist);
}

void pen_tallyPhaseSpaceFile::tally_interfCross(const unsigned long long nhist,
						const unsigned kdet,
						const pen_KPAR kpar,
						const pen_particleState& state){

  if((int)kdet == detector){
    if(!inside && enabledKpars[kpar] && state.E >= emin && state.E <= emax)
      store(nhist,kpar,state);

    inside = true;
  }
  else{
    inside = false;
  }
}

void pen_tallyPhaseSpaceFile::tally_matChange(const unsigned long long nhist,
					      const pen_KPAR kpar,
					      const pen_particleState& state,
					      const unsigned /*prevMat*/){
  if((int)state.MAT == material){
    if(enabledKpars[kpar] && state.E >= emin && state.E <= emax)
      store(nhist,kpar,state);
  }
}

void pen_tallyPhaseSpaceFile::tally_move2geo(const unsigned long long nhist,
					     const unsigned kdet,
					     const pen_KPAR kpar,
					     const pen_particleState& state,
					     const double /*dsef*/,
					     const double /*dstot*/){

  if((int)kdet == detector || (int)state.MAT == material){
    inside = true;
    if(enabledKpars[kpar] && state.E >= emin && state.E <= emax)
      store(nhist,kpar,state);
  }else{
    inside = false;
  }
  
}

void pen_tallyPhaseSpaceFile::tally_beginPart(const unsigned long long /*nhist*/,
					      const unsigned kdet,
					      const pen_KPAR /*kpar*/,
					      const pen_particleState& /*state*/){
  if((int)kdet == detector){
    inside = true;
  }else{
    inside = false;
  }
}

void pen_tallyPhaseSpaceFile::tally_endSim(const unsigned long long /*nhist*/){
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
  
  //Create the splitted file in thread 0
  if(getThread() == 0){
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

    pSF = std::make_shared<pen_splittedFile>(filename.c_str(),true);
    
    //try to create a partition for thread 0
    err = pSF.get()->createPartition(getThread());
    if(err != SPLITTED_FILE_SUCCESS){
      if(verbose > 0){
	printf("PhaseSpaceFile:configure: Error: Unable to create a partition for "
	       "tally %s thread 0.\n",readName().c_str());
	printf("                 Error code: %d\n",err);
      }
      return -1;
    }
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

  //Check if particles must be enabled or disabled by default
  bool defaultEnabled;
  if(config.read("particles/default", defaultEnabled) == INTDATA_SUCCESS){
    if(verbose > 1){
      printf(" Default behaviour for particle type recording has been "
	     "set to: %s\n", defaultEnabled ? "enabled" : "disabled");
    }
    std::fill(enabledKpars.begin(), enabledKpars.end(), defaultEnabled);
  }else{
    std::fill(enabledKpars.begin(), enabledKpars.end(), true);    
  }


  //Check if someone must be disabled
  std::vector<std::string> enabledPartNames;
  config.ls("particles",enabledPartNames);

  for(const std::string& partName : enabledPartNames){

    if(partName.compare("default") == 0)
      continue;
    
    //Extract kpar
    unsigned enableKpar = particleID(partName.c_str());
    
    if(enableKpar == pen_KPAR::ALWAYS_AT_END){
      //Unknown particle name
      if(verbose > 1){
	printf("tallyPhaseSpaceFile:configure: Warning: "
	       "Unknow particle name '%s'\n",partName.c_str());
      }
    }else{
      //Enable or disable this particle type
      std::string enablePath("particles/");
      enablePath += partName;
      bool enabledkpar;
      if(config.read(enablePath,enabledkpar) == INTDATA_SUCCESS){
	enabledKpars[enableKpar] = enabledkpar;
	
	if(verbose > 1){
	  printf(" Particle '%s' recording: %s",
		 partName.c_str(),
		 enabledkpar ? "enabled" : "disabled");
	}
      }else{
	if(verbose > 0){
	  printf("PhaseSpaceFile:configure: Error: Unable to read "
		 "'%s' in configuration. Boolean expected\n",
		 enablePath.c_str());
	}
	return -8;
      }
    }
  }
  
  return 0;
}


void pen_tallyPhaseSpaceFile::saveData(const unsigned long long /*nhist*/) const{
  
}
int pen_tallyPhaseSpaceFile::sumTally(const pen_tallyPhaseSpaceFile& /*tally*/){return 0;}


REGISTER_COMMON_TALLY(pen_tallyPhaseSpaceFile, PSF)
