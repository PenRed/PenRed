
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
//        vicente.gimenez@uv.es
//    
//

 
#include "pen_context.hh"

// Material reader
int pen_contextReaderMat::beginSectionFamily(const std::string& pathInSection,
					  const size_t,
					  const unsigned){

  //First level sections
  if(family == -1){
    if(pathInSection.compare("material-ranges") == 0){
      family = 0;
      return SUCCESS;
    }
    else if(pathInSection.compare("material-eabs") == 0){
      family = 1;
      return SUCCESS;
    }
    else if(pathInSection.compare("materials") == 0){
      family = 2;
      return SUCCESS;
    }
  }
  //Sections in "materials"
  else if(family == 2){
    if(pathInSection.compare("elements") == 0){
      family = 3;
      return SUCCESS;
    }
    else if(pathInSection.compare("range") == 0){
      family = 4;
      return SUCCESS;
    }
    else if(pathInSection.compare("eabs") == 0){
      family = 5;
      return SUCCESS;
    }    
  }
  
  return UNHANDLED;
}

int pen_contextReaderMat::endSectionFamily(const unsigned){

  if(family == 0 || family == 1 || family == 2){
    //End material-ranges, material-eabs, materials section
    family = -1;
    return SUCCESS;
  }
  //End materials/${subsection} subsections. Return to "material" subsection
  else if(family == 3 || family == 4 || family == 5){
    //End elements, range or eabs subsections
    family = 2;
    return SUCCESS;
  }
  
  return UNHANDLED;
}

int pen_contextReaderMat::beginSection(const std::string& name,
				       const unsigned verbose){

  if(family == 0 || family == 1 || family == 4 || family == 5){
    //Subsection in material-ranges, material-eabs, material range or eabs
    
    //Get particle ID from particle "name"
    actualPartID = particleID(name.c_str());
    if(actualPartID == ALWAYS_AT_END){
      if(verbose > 0){
	printf("Error: Unknown particle type %s\n"
	       " Known particles:\n", name.c_str());
	for(size_t i = 0; i < constants::nParTypes; ++i){
	  printf("    %s\n", particleName(i));
	}
	
      }
      return UNKNOWN_PARTICLE;
    }
    return SUCCESS;
  }
  else if(family == 2){
    //Subsection in materials
    if(name.empty()){
      if(verbose > 0){
	printf("Error: Empty material name");
      }
      return INVALID_NAME;
    }
    mats.emplace_back(name);
    return SUCCESS;
  }
  else if(family == 3){
    //Subsection in materials/${subsection}/elements

    //Try to convert the section name to an integer
    unsigned Z;
    if(sscanf(name.c_str(), "%u", &Z) != 1){
      if(verbose > 0){
	printf("Error: Unknown atomic number: %s\n", name.c_str());
      }
      return INVALID_ATOMIC_NUMBER;      
    }
    
    mats.back().composition.emplace_back(Z);
    return SUCCESS;
  }
  
  return UNHANDLED;
}

int pen_contextReaderMat::endSection(const unsigned){ return SUCCESS; }  

int pen_contextReaderMat::storeElement(const std::string& pathInSection,
				       const pen_parserData& element,
				       const unsigned verbose){

  if(family == 0){
    //Subsection in material-ranges
    
    if(pathInSection.empty()){
      if(actualPartID >= ALWAYS_AT_END){
	if(verbose > 0){
	  printf("Error: Unexpected error, unknown "
		 "particle id %u\n", actualPartID);
	}
	return UNKNOWN_PARTICLE;
      }else{
	maxRanges[actualPartID] = element;
	return SUCCESS;
      }
    }
  }
  else if(family == 1){
    //Subsection in material-eabs
    
    if(pathInSection.empty()){
      if(actualPartID >= ALWAYS_AT_END){
	if(verbose > 0){
	  printf("Error: Unexpected error, unknown "
		 "particle id %u\n", actualPartID);
	}
	return UNKNOWN_PARTICLE;
      }else{
	defaultEabs[actualPartID] = element;
	return SUCCESS;
      }
    }
  }
  else if(family == 2){
    //Subsection in materials
    
    if(pathInSection.compare("number") == 0){
      mats.back().index = element;
      return SUCCESS;
    }
    else if(pathInSection.compare("force-creation") == 0){
      mats.back().forceCreation = element;
      return SUCCESS;
    }
    else if(pathInSection.compare("density") == 0){
      mats.back().density = element;
      return SUCCESS;
    }
    else if(pathInSection.compare("C1") == 0){
      mats.back().C1 = element;
      return SUCCESS;
    }
    else if(pathInSection.compare("C2") == 0){
      mats.back().C2 = element;
      return SUCCESS;
    }
    else if(pathInSection.compare("WCC") == 0){
      mats.back().WCC = element;
      return SUCCESS;
    }
    else if(pathInSection.compare("WCR") == 0){
      mats.back().WCR = element;
      return SUCCESS;
    }

  }
  else if(family == 3){
    //Subsection in materials/${subsection}/elements
    if(pathInSection.empty()){
      mats.back().composition.back().fraction = element;
      return SUCCESS;
    }
  }
  else if(family == 4){
    //Subsection in materials/${subsection}/range
    if(pathInSection.empty()){
      mats.back().maxRanges[actualPartID] = element;
      return SUCCESS;
    }
  }
  else if(family == 5){
    //Subsection in materials/${subsection}/eabs
    if(pathInSection.empty()){
      mats.back().eabs[actualPartID] = element;
      return SUCCESS;
    }
  }
  
  return UNHANDLED;
}

int pen_contextReaderMat::storeString(const std::string& pathInSection,
				      const std::string& element,
				      const unsigned){
  if(family == 2){
    if(pathInSection.compare("filename") == 0){
      if(element.compare("-") != 0)
	mats.back().filename = element;
      return SUCCESS;
    }
  }else if(family == -1){
    if(pathInSection.compare("context-log") == 0){
      contextlogfile = element;
      return SUCCESS;
    }
  }
  return UNHANDLED;  
}

// VR reader
int pen_contextReaderVR::beginSectionFamily(const std::string& pathInSection,
					    const size_t,
					    const unsigned){

  //First level sections
  if(family == -1){
    if(pathInSection.compare("VR/IForcing") == 0){
      family = 0;
      return SUCCESS;
    }
    else if(pathInSection.compare("VR/bremss") == 0){
      family = 3;
      return SUCCESS;
    }
  }
  //Sections in "VR/IForcing"
  else if(family == 0){
    if(pathInSection.compare("bodies") == 0){
      family = 1;
      return SUCCESS;
    }
    else if(pathInSection.compare("materials") == 0){
      family = 2;
      return SUCCESS;
    }
  }
  //Sections in "VR/bremss"
  else if(family == 3){
    if(pathInSection.compare("bodies") == 0){
      family = 4;
      return SUCCESS;
    }
    else if(pathInSection.compare("materials") == 0){
      family = 5;
      return SUCCESS;
    }
  }
  return UNHANDLED;
}

int pen_contextReaderVR::endSectionFamily(const unsigned){

  if(family == 0 || family == 3){
    //VR/IForcing or VR/bremss sections
    family = -1;
    return SUCCESS;
  }
  //End VR/IForcing/${subsection} subsections. Return to "VR/IForcing" subsection
  else if(family == 1 || family == 2){
    family = 0;
    return SUCCESS;
  }
  //End VR/bremss/${subsection} subsections. Return to "VR/bremss" subsection
  else if(family == 4 || family == 5){
    family = 3;
    return SUCCESS;
  }
  
  return UNHANDLED;
}

int pen_contextReaderVR::beginSection(const std::string& name,
				      const unsigned verbose){

  //Interaction forcing sections
  if(family == 0){
    //VR/IForcing begins
    iforcing.emplace_back(name);
    return SUCCESS;
  }
  else if(family == 1){
    //VR/IForcing/${subsection}/bodies begins
    iforcing.back().enabledBodies.emplace_back(name);
    return SUCCESS;
  }
  else if(family == 2){
    //VR/IForcing/${subsection}/materials begins

    //Convert name to material index
    unsigned imat;
    if(sscanf(name.c_str(), "%u", &imat) != 1){
      if(verbose > 0){
	printf("Error: Invalid material index: %s\n", name.c_str());
      }
      return INVALID_MATERIAL_INDEX;      
    }
    if(imat <= 0 || imat >= constants::MAXMAT){
      if(verbose > 0){
	printf("Error: Material index must be in the interval (0,%u)\n",
	       constants::MAXMAT);
      }
      return INVALID_MATERIAL_INDEX;      
    }
    iforcing.back().enabledMats.emplace_back(imat);
    return SUCCESS;
  }
  //Bremsstrahlung splitting sections
  else if(family == 3){
    //VR/bremss begins
    bremss.emplace_back(name);
    return SUCCESS;
  }
  else if(family == 4){
    //VR/bremss/${subsection}/bodies begins
    bremss.back().enabledBodies.emplace_back(name);
    return SUCCESS;
  }
  else if(family == 5){
    //VR/bremss/${subsection}/materials begins

    //Convert name to material index
    unsigned imat;
    if(sscanf(name.c_str(), "%u", &imat) != 1){
      if(verbose > 0){
	printf("Error: Invalid material index: %s\n", name.c_str());
      }
      return INVALID_MATERIAL_INDEX;      
    }
    if(imat <= 0 || imat >= constants::MAXMAT){
      if(verbose > 0){
	printf("Error: Material index must be in the interval (0,%u)\n",
	       constants::MAXMAT);
      }
      return INVALID_MATERIAL_INDEX;      
    }
    bremss.back().enabledMats.emplace_back(imat);
    return SUCCESS;
  }
  
  return UNHANDLED;
}

int pen_contextReaderVR::endSection(const unsigned){ return SUCCESS; }  

int pen_contextReaderVR::storeElement(const std::string& pathInSection,
				       const pen_parserData& element,
				       const unsigned verbose){

  //Interaction forcing sections
  if(family == 0){
    if(pathInSection.compare("interaction") == 0){
      iforcing.back().interaction = element;
      if(iforcing.back().interaction >= constants::MAXINTERACTIONS){
	if(verbose > 0)
	  printf("Error: Interaction index (%u) out of range.\n",
		 iforcing.back().interaction);
	return INVALID_INTERACTION_INDEX;
      }
      return SUCCESS;
    }
    else if(pathInSection.compare("factor") == 0){
      iforcing.back().factor = element;
      return SUCCESS;
    }
    else if(pathInSection.compare("min-weight") == 0){
      iforcing.back().minW = element;
      return SUCCESS;
    }    
    else if(pathInSection.compare("max-weight") == 0){
      iforcing.back().maxW = element;
      return SUCCESS;
    }
  }
  else if(family == 1){
    //VR/IForcing/${subsection}/bodies
    if(pathInSection.empty()){
      if(!element){
	iforcing.back().enabledBodies.pop_back();
      }
      return SUCCESS;
    }
  }
  else if(family == 2){
    //VR/IForcing/${subsection}/materials
    if(pathInSection.empty()){
      if(!element){
	iforcing.back().enabledMats.pop_back();
      }
      return SUCCESS;
    }
  }
  //Bremsstrahlung splitting sections
  else if(family == 3){
    if(pathInSection.compare("splitting") == 0){
      bremss.back().factor = element;
      return SUCCESS;
    }
  }
  else if(family == 4){
    //VR/bremss/${subsection}/bodies
    if(pathInSection.empty()){
      if(!element){
	bremss.back().enabledBodies.pop_back();
      }
      return SUCCESS;
    }
  }
  else if(family == 5){
    //VR/bremss/${subsection}/materials
    if(pathInSection.empty()){
      if(!element){
	bremss.back().enabledMats.pop_back();
      }
      return SUCCESS;
    }
  }
  
  return UNHANDLED;
}

int pen_contextReaderVR::storeString(const std::string& pathInSection,
				     const std::string& element,
				     const unsigned verbose){
  if(family == 0){
    if(pathInSection.compare("particle") == 0){
      unsigned kpar = particleID(element.c_str());
      if(kpar == ALWAYS_AT_END){
	if(verbose > 0){
	  printf("Error: Unknown particle id %u\n",
		 kpar);
	}
	return UNKNOWN_PARTICLE;	
      }
      iforcing.back().particleType = static_cast<pen_KPAR>(kpar);
      return SUCCESS;
    }
  }
  return UNHANDLED;  
}

int pen_context::init(double EMAX, FILE *IWR, int INFO, std::string PMFILE[constants::MAXMAT])
{
  
  //  Input of material data and initialisation of simulation routines.
  //
  //  Each material is defined through an input file, which is created by
  //  the program MATERIAL using information contained in the database.
  //  This file can be modified by the user if more accurate interaction
  //  data are available.
  //
  //  Input arguments:
  //    EMAX ... maximum particle energy (kinetic energy for electrons and
  //             positrons) used in the simulation. Note: Positrons with
  //             energy E may produce photons with energy E+1.022E6.
  //    NMATER .... number of materials in the geometry.
  //    PMFILE .... array of MAXMAT strings. The first NMATER
  //             elements are the filenames of the material data files.
  //             The file PMFILE(M) contains radiation interaction data
  //             for material M (that is, the order is important!).
  //    IWR .... output unit.
  //    INFO ... determines the amount of information that is written on
  //             the output file,
  //               INFO=1 (or less), minimal (composition data only).
  //               INFO=2, medium (same information as in the material
  //                 definition data file, useful to check that the struc-
  //                 ture of the latter is correct).
  //               INFO=3 or larger, full information, including tables of
  //                 interaction properties used in the simulation.
  //
  //  For the preliminary computations, PEINIT needs to know the absorption
  //  energies EABS(KPAR,M) and the simulation parameters C1(M), C2(M),
  //  WCC(M) and WCR(M). This information is introduced through the module
  //  PENELOPE_mod, which has to be loaded before calling PEINIT.
  //
  
  char LIT[3];
  //  ****  Simulation parameters.
  double EABS0[constants::nParTypes][constants::MAXMAT];
  
  if (IWR != nullptr) fprintf(IWR  ,"\n **********************************\n **   PENELOPE  (version 2018)   **");
  if (IWR != nullptr) fprintf(IWR  ,"\n **********************************\n");

  //Check material absortion energies
  for(unsigned M = 0; M < getNMats(); M++)
    {
      pen_material& mat = getBaseMaterial(M);
      
      if(mat.EABS[PEN_ELECTRON] < 49.999)
	{
	  if (IWR != nullptr) fprintf(IWR , "\nMaterial %3d, EABS(%d,%2d) = %11.4E eV\n *** ERROR: electron absorption energy cannot be less than %.5E eV\n",
		  PEN_ELECTRON,M, M, mat.EABS[PEN_ELECTRON], double(constants::MINEABSE));
	  penError(ERR_PEINIT_ELECTRON_ENERGY);
	  return ERR_PEINIT_ELECTRON_ENERGY;
	}
      EABS0[PEN_ELECTRON][M] = mat.EABS[PEN_ELECTRON];
      if(mat.EABS[PEN_PHOTON] < 49.999)
	{
	  if (IWR != nullptr) fprintf(IWR , "\nMaterial %3d, EABS(%d,%2d) = %11.4E eV\n *** ERROR: photon absorption energy cannot be less than %.5E eV\n",
		  PEN_PHOTON,M, M, mat.EABS[PEN_PHOTON], double(constants::MINEABSPh));
	  penError(ERR_PEINIT_GAMMA_ENERGY);
	  return ERR_PEINIT_GAMMA_ENERGY;
	}
      EABS0[PEN_PHOTON][M] = mat.EABS[PEN_PHOTON]; 
      if(mat.EABS[PEN_POSITRON] < 49.999)
	{
	  if (IWR != nullptr) fprintf(IWR , "\nMaterial %3d, EABS(%d,%2d) = %11.4E eV\n *** ERROR: positron absorption energy cannot be less than %.5E eV\n",
		  PEN_POSITRON,M, M, mat.EABS[PEN_POSITRON], double(constants::MINEABSPo));
	  penError(ERR_PEINIT_POSITRON_ENERGY);
	  return ERR_PEINIT_POSITRON_ENERGY;
	}
      EABS0[PEN_POSITRON][M] = mat.EABS[PEN_POSITRON];
    }

  //  ****  Lower limit of the energy grid.

  double EMIN = 1.0E35;
  for(unsigned M = 0; M < getNMats(); M++)
    {
      pen_material& mat = getBaseMaterial(M);

      if(mat.EABS[PEN_ELECTRON] >= EMAX)
	{
	  mat.EABS[PEN_ELECTRON] = EMAX*(double)0.9999;
	  if (IWR != nullptr) fprintf(IWR , "\n WARNING: EABS(%d,%2d) has been modified because is greater than EMAX\n", PEN_ELECTRON,M);
	}
      if(mat.EABS[PEN_PHOTON] >= EMAX)
	{
	  mat.EABS[PEN_PHOTON] = EMAX*(double)0.9999;
	  if (IWR != nullptr) fprintf(IWR , "\n WARNING: EABS(%d,%2d) has been modified because is greater than EMAX\n", PEN_PHOTON,M);
	}
      if(mat.EABS[PEN_POSITRON] >= EMAX)
	{
	  mat.EABS[PEN_POSITRON] = EMAX*(double)0.9999;
	  if (IWR != nullptr) fprintf(IWR , "\n WARNING: EABS(%d,%2d) has been modified because is greater than EMAX\n", PEN_POSITRON,M);
	}
      if(mat.EABS[PEN_ELECTRON] < EMIN)
	{
	  EMIN = mat.EABS[PEN_ELECTRON];
	}
      if(mat.EABS[PEN_PHOTON] < EMIN)
	{
	  EMIN = mat.EABS[PEN_PHOTON];
	}
      if(mat.EABS[PEN_POSITRON] < EMIN)
	{
	  EMIN = mat.EABS[PEN_POSITRON];
	}
    }
      
  if(EMIN < (double)constants::MINEGRID){ EMIN = (double)constants::MINEGRID;}

  if (IWR != nullptr) fprintf(IWR , "\n EMIN =%11.4E eV,  EMAX =%11.4E eV\n", EMIN, EMAX);
  
  if(EMAX < EMIN+(double)constants::MINDIMGRID)
    {
      if (IWR != nullptr) fprintf(IWR , "\n ERROR: The energy interval between EMIN and EMAX is less than %.5E eV\n", constants::MINDIMGRID);
      penError(ERR_PEINIT_ENERGY_INTERVAL);
      return ERR_PEINIT_ENERGY_INTERVAL;
    }
      
  if(INFO > 2){ if (IWR != nullptr) fprintf(IWR , "\n NOTE: 1 mtu = 1 g/cm**2\n");}
		     
  grid.init(EMIN,EMAX);  // Defines the simulation energy grid.
  elements.GPHa0();  // Initialises photoelectric routines.
  elements.RELAX0();  // Initialises atomic relaxation routines.
  RNDG30(rndg3);  // Initialises the Gaussian sampling routine.
  if(penGetError() != PEN_SUCCESS){ return penGetError();}

  for(unsigned M = 0; M < getNMats(); M++)
    {
      pen_material& mat = getBaseMaterial(M);
    
      if(M == 0){ LIT[0]='s'; LIT[1]='t';}
      if(M == 1){ LIT[0]='n'; LIT[1]='d';}
      if(M == 2){ LIT[0]='r'; LIT[1]='d';}
      if(M > 2){ LIT[0]='t'; LIT[1]='h';}
      LIT[2]='\0';
      if (IWR != nullptr) fprintf(IWR , "\n\n **********************\n **  %2d%s material   **\n **********************\n", M+1, LIT);
	  
      if (IWR != nullptr) fprintf(IWR , "\n Material data file: %-20s\n", PMFILE[M].c_str());
      FILE* IRD = nullptr;
      IRD = fopen(PMFILE[M].c_str(),"r");
      if(IRD == nullptr)
	{
	  if (IWR != nullptr) fprintf(IWR , "Error: Material file %s could not be opened\n", PMFILE[M].c_str());
	  penError(ERR_PEINIT_MATERIAL_FILE);
	  return ERR_PEINIT_MATERIAL_FILE;
	}
      
      //Save material filename
      mat.setDataPath(PMFILE[M]);
      
      //  ****  Energy limits and thresholds.
	  
      if (IWR != nullptr) fprintf(IWR , "\n *** Simulation parameters:\n");
      if (IWR != nullptr) fprintf(IWR , "     Electron absorption energy =%11.4E eV\n", mat.EABS[PEN_ELECTRON]);
      if (IWR != nullptr) fprintf(IWR , "       Photon absorption energy =%11.4E eV\n", mat.EABS[PEN_PHOTON]);
      if (IWR != nullptr) fprintf(IWR , "     Positron absorption energy =%11.4E eV\n", mat.EABS[PEN_POSITRON]);

      if(fabs(mat.C1) < 0.2){mat.C1 = fabs(mat.C1);}
      else{mat.C1 = 0.2;}
	  
      if(fabs(mat.C2) < 0.2){mat.C2 = fabs(mat.C2);}
      else{mat.C2 = 0.2;}
	  
      if(fabs(mat.WCC) < EMAX){mat.WCC = fabs(mat.WCC);}
      else{mat.WCC = EMAX;}

//      if(mat.WCC > mat.EABS[PEN_ELECTRON]){ mat.WCC = mat.EABS[PEN_ELECTRON];}
//      if(mat.WCR > mat.EABS[PEN_PHOTON]){ mat.WCR = mat.EABS[PEN_PHOTON];}
      if(mat.WCR < 0.0){
	if (IWR != nullptr) fprintf(IWR , "*** Warning: soft radiative losses are switched off in material number %3d\n", M);
      }
	  
      if (IWR != nullptr) fprintf(IWR , "      C1 =%11.4E,       C2 =%11.4E\n     WCC =%11.4E eV,   WCR =%11.4E eV\n\n",
	      mat.C1, mat.C2, mat.WCC, (mat.WCR>10.0 ? mat.WCR:10.0));

      initStructs initStore;
      mat.load(IRD,IWR,initStore,elements,grid,INFO);
      fclose(IRD);

      if(penGetError() != PEN_SUCCESS){
	return penGetError();
      }
    }
      
  //  ****  Restore the user values of EABS(.), to avoid inconsistencies
  //        in the main program.

  for(unsigned M=0; M < getNMats(); M++)
    {
      pen_material& mat = getBaseMaterial(M);
      
      mat.EABS[PEN_ELECTRON]=EABS0[PEN_ELECTRON][M];
      mat.EABS[PEN_PHOTON]=EABS0[PEN_PHOTON][M];
      mat.EABS[PEN_POSITRON]=EABS0[PEN_POSITRON][M];
    }
  return PEN_SUCCESS;
}

int pen_context::configure(const double EMAX,
			   const pen_parserSection& config,
			   pen_parserSection& matInfo,
			   const unsigned verbose){

  if(verbose > 1){
    printf(" ** Context configuration:\n\n");
  }
  
  //Read material information from config section
  pen_contextReaderMat matReader;
  int err = matReader.read(config,verbose);
  if(err != pen_contextReaderMat::SUCCESS){
    return err;
  }

  //Process material data
  if(verbose > 1){
    printf("   + Global maximum ranges:\n\n");
    unsigned count = 0;
    for(unsigned i = 0; i < constants::nParTypes; ++i){
      if(matReader.maxRanges[i] > 0.0){
	printf("     - %s : %E cm\n",
	       particleName(i), matReader.maxRanges[i]);
	++count;
      }
    }
    if(count == 0){
      printf("     - No range specified\n");
    }
    printf("\n");

    printf("   + Default global absorption energies:\n\n");
    for(unsigned i = 0; i < constants::nParTypes; ++i){
      printf("     - %s : %E eV\n",
	     particleName(i), matReader.defaultEabs[i]);
      
    }
    printf("\n");

    printf("  + Number of materials: %u\n",
	   static_cast<unsigned>(matReader.mats.size()));    
  }

  // ** Check number of materials
  if(matReader.mats.size() == 0){
    if(verbose > 0){
      printf("Error: No material found at configuration. "
	     "Simulation requires, at least, one material.\n");
    }
    return NO_MATERIAL_PROVIDED;
  }
  
  if(matReader.mats.size() > constants::MAXMAT){
    if(verbose > 0){
      printf("Error: %d materials found. The maximum "
	     "number of materials is %d.\n",
	     static_cast<unsigned>(matReader.mats.size()),
	     constants::MAXMAT);
    }
    return TOO_MUCH_MATERIALS;
  }

  // ** Set materials in the context
  err = setMats<pen_material>(matReader.mats.size());
  if(err != 0){
    printf("Error creating materials: %d.\n",err);
    return UNABLE_TO_CREATE_MATERIALS;
  }

  //Create a vector with materials filenames
  std::string matFilenames[constants::MAXMAT];

  //Iterate over each material to configure it
  for(size_t i = 0; i < matReader.mats.size(); ++i){

    //Get material reader reference
    const pen_contextReaderMat::materialData& matData =
      matReader.mats[i];
    
    if(verbose > 1){
      printf("\n\n------------------------------------\n\n");
      printf(" **** Material '%s'\n\n",matData.name.c_str());
    }

    //Ensure material index is inside bounds
    if(getNMats() <= matData.index-1){
      if(verbose > 0){
	printf("Error: Material index must be lesser or equal"
	       " to the total number of materials.\n"
	       "    Error material     : %s\n"
	       "    Material index     : %u\n"
	       "    Number of materials: %u\n",
	       matData.name.c_str(), matData.index,
	       getNMats());
      }
      return MATERIAL_OUT_OF_RANGE;
    }

    //Check if this material has already been set
    if(!matFilenames[matData.index-1].empty()){
      if(verbose > 0){
	printf("Error: Material with index %u has been "
	       "defined multiple times\n", matData.index);
      }
      return MATERIAL_DEFINED_TWICE;
    }
    
    //Save material filename
    matFilenames[matData.index-1] = matData.filename;

    //Check if the file exists or if the
    //material construction is forced
    FILE* file = nullptr;
    file = fopen(matData.filename.c_str(), "r");
    bool fileAccessible = false;
    if(file != nullptr){
      fileAccessible = true;
      fclose(file);
    }

#ifdef _PEN_EMBEDDED_DATA_BASE_
    if(!fileAccessible || matData.forceCreation){

      //Create the material
      if(matData.composition.size() == 0){
	if(verbose > 0){
	  printf("Error: Unable to create material %s."
		 " No composition provided.\n",
		 matData.name.c_str());
	}
	return MISSING_MATERIAL_COMPOSITION;
      }

      //******************************************************
      //Fix the minimum number of exponent digits in old MVS versions to 2 
#ifdef _MSC_VER
  #if _MSC_VER < 1900
      unsigned int prev_exponent_format =
	_set_output_format(_TWO_DIGIT_EXPONENT);
  #endif
#endif
      //******************************************************

      //Create the material
      std::string errorString;
      err = penred::penMaterialCreator::createMat(matData.name,
						  matData.density,
						  matData.composition,
						  errorString,
						  matData.filename);
      if(err != 0){
	if(verbose > 0){
	  printf("Unable to create material %s.\n"
		 "%s\n"
		 "IRETRN =%d\n",
		 matData.name.c_str(),
		 errorString.c_str(),
		 err);
	}
	return MATERIAL_CREATION_FAILED;
      }
    }
#else
    if(!fileAccessible){
      if(verbose > 0){
	printf("Error: Unable to access material file '%s'\n",
	       matData.filename.c_str());
      }
      return UNABLE_TO_ACCESS_MATERIAL_FILE;
    }
#endif

    
    //Get context material
    pen_material& mat = getBaseMaterial(matData.index-1);

    //Set initial parameters for class II transport
    mat.C1 = matData.C1;
    mat.C2 = matData.C2;
    mat.WCC = matData.WCC;
    mat.WCR = matData.WCR;

    //Set particle absorption energies
    pen_context* rangeContext = nullptr;
    for(unsigned ip = 0; ip < constants::nParTypes; ++ip){
      double partRange = matReader.maxRanges[ip];
      if(matData.maxRanges[ip] > 0.0)
	partRange = matData.maxRanges[ip];
      double eabs = matReader.defaultEabs[ip];
      if(matData.eabs[ip] > 0.0){
	eabs = matData.eabs[ip];
      }

      //Assign absortion energy by range, explicit value or default value
      if(partRange > 0.0 && matData.eabs[ip] < 0.0){

	//Create range context if has not previously created
	if(rangeContext == nullptr)
	  rangeContext = createAuxContext(EMAX,
					  matData.filename.c_str(),
					  verbose);
	if(rangeContext == nullptr){
	  if(verbose > 0){
	    printf("Error: Unable to create auxiliary "
		   "context to calculate range limits.\n");
	  }
	  return RANGE_CONTEXT_CREATION_FAILED;
	}
	  
	// + Range
	double topE = EMAX;
	double lowE = 50;
	double objectiveRange = partRange;
	unsigned nTries = 0;
	do{
	  double midE = (topE+lowE)/2.0;
	  double range = rangeContext->range(midE,static_cast<pen_KPAR>(ip),0);
	  if(range == objectiveRange){
	    topE = midE;
	    lowE = midE;
	  }
	  else if(range > objectiveRange){
	    topE = midE;
	  }else{
	    lowE = midE;
	  }
	  ++nTries;
	}while(nTries < 1000000 && topE/lowE > 1.001);
	if(topE == EMAX)
	  mat.EABS[ip] = 1.0e35;
	else
	  mat.EABS[ip] = lowE;
	
      }
      else{
	// + Set absorption energy
	mat.EABS[ip] = eabs;
      }
    }    

    //Free range context if has been created
    if(rangeContext != nullptr){
      delete rangeContext;
    }

    //Check WCC and WCR values
    if(mat.WCC <= 0.0){
      mat.WCC = std::min(5e3,mat.EABS[PEN_ELECTRON]/100.0);
    }
    if(mat.WCR <= 0.0){
      mat.WCR = std::min(5e3,mat.EABS[PEN_PHOTON]/100.0);
    }

    if(verbose > 1){
      printf("      C1 =%11.4E       C2 =%11.4E\n",mat.C1,mat.C2);
      printf("     WCC =%11.4E eV,   WCR =%11.4E eV\n",
	     mat.WCC,(mat.WCR > 10.0E0 ? mat.WCR : 10.0E0));
      printf("  electron EABS: %11.4E eV\n",mat.EABS[PEN_ELECTRON]);
      printf("     gamma EABS: %11.4E eV\n",mat.EABS[PEN_PHOTON]);
      printf("  positron EABS: %11.4E eV\n",mat.EABS[PEN_POSITRON]);
      
      printf("\n Material filename: '%s'.\n",matData.filename.c_str());
      printf("\n Material number: %d.\n",matData.index);
      
      printf("\n\n------------------------------------\n\n");    
    }    
  }

  //Once materials have been configured, init context
  FILE* fcontext = nullptr;
  int contextVerbose = verbose;
  if(matReader.contextlogfile.length() > 0){
    fcontext = fopen(matReader.contextlogfile.c_str(),"w");
    if(fcontext == nullptr){
      printf("Error: unable to open context log file '%s'\n",
	     matReader.contextlogfile.c_str());
    }
  }else{
    fcontext = stdout;
    contextVerbose = 1;
  }
  if(init(EMAX,fcontext,contextVerbose,matFilenames) != PEN_SUCCESS){
    if(fcontext != nullptr) fclose(fcontext);
    printf("Error at context initialization.\n");
    return ERROR_AT_CONTEXT_INIT;
  }
  if(matReader.contextlogfile.length() > 0 && fcontext != nullptr) fclose(fcontext);

  //Construct material information section
  matInfo.clear();
  for(unsigned imat = 0; imat < getNMats(); imat++){
    std::string matPrefix = matReader.mats[imat].name;
    matInfo.set((matPrefix + "/ID").c_str(),
		static_cast<int>(matReader.mats[imat].index));
    matInfo.set((matPrefix + "/density").c_str(),
		readBaseMaterial(matReader.mats[imat].index-1).readDens());
  }

  return SUCCESS;
}

int pen_context::configureWithGeo(const pen_parserSection& config,
				  const unsigned verbose){

  if(readGeometry() == nullptr)
    return GEOMETRY_NOT_SET;

  //Read VR information from config section
  pen_contextReaderVR VRReader;
  int err = VRReader.read(config,verbose);
  if(err != pen_contextReaderVR::SUCCESS){
    return err;
  }

  //Get materials used by the current geometry
  bool usedMat[constants::MAXMAT+1];
  readGeometry()->usedMat(usedMat);
  
  if(VRReader.iforcing.size() == 0 && VRReader.bremss.size() == 0)
    return SUCCESS;
  
  if(verbose > 1){
    printf("  * Context variance reduction configuration:\n\n");
  }
  

  //Iterate over interaction forcing data
  for(const pen_contextReaderVR::IFdata& ifdata : VRReader.iforcing){
    
    if(verbose > 1){
      printf("    + Interaction forcing %s:\n"
	     "  Particle    | Interaction |  IF factor  | Weight Range\n"
	     " %12s       %2u           %+4f       (%f,%f)\n\n",
	     ifdata.name.c_str(),
	     particleName(ifdata.particleType),
	     ifdata.interaction,
	     ifdata.factor,
	     ifdata.minW, ifdata.maxW);
    }

    //Iterate over material indexes
    if(ifdata.enabledMats.size() > 0){
      if(verbose > 1){
	printf("      - Applied on materials: \n");
      }
      for(const unsigned& imat : ifdata.enabledMats){
	//Check if the material is used in the geometry
	if(!usedMat[imat]){
	  if(verbose > 0){
	    printf("Error: Specified material index (%u) for "
		   "interaction forcing '%s' is not used at current geometry.\n",
		   imat, ifdata.name.c_str());
	  }
	  return UNUSED_MATERIAL;
	}

	if(verbose > 1){
	  printf("        + %u \n", imat);
	}

	//Set parameters for bodies with this material index
	for(unsigned ibody = 0; ibody < readGeometry()->getBodies(); ibody++){
	
	  if(readGeometry()->getMat(ibody) != imat) continue;
	  
	  if(ibody >= NBV){
	    if(verbose > 0){
	      printf("Error: Maximum body index for IF reached (%u)\n",NBV);
	    }
	    return MAXIMUM_BODY_IF_REACHED;
	  }

	  //Set interaction forcing in context
	  setForcing(ifdata.factor, ifdata.particleType,
			     ifdata.interaction, ibody,
			     ifdata.minW, ifdata.maxW);
	}
      }

      if(verbose > 1){
	printf("\n");
      }
      
    }

    //Iterate over enabled bodies
    if(ifdata.enabledBodies.size() > 0){
      if(verbose > 1){
	printf("      - Applied on bodies: \n");
      }
    }
    for(const std::string& bodyAlias : ifdata.enabledBodies){

      if(verbose > 1){
	printf("        + %s \n", bodyAlias.c_str());
      }
      
      unsigned ibody = readGeometry()->getIBody(bodyAlias.c_str());
      if(ibody >= readGeometry()->getBodies()){
	if(verbose > 0){
	  printf("Error: Body '%s' doesn't exists in loaded geometry.\n",
		 bodyAlias.c_str());
	}
	return UNKNOWN_BODY;
      }

      if(ibody >= NBV){
	if(verbose > 0){
	  printf("Error: Maximum body index for IF is (%u)\n",NBV);
	  printf("                  specified index: %u\n",ibody);
	}
	return MAXIMUM_BODY_IF_REACHED;	  
      }

      //Set interaction forcing in context
      setForcing(ifdata.factor, ifdata.particleType,
		 ifdata.interaction, ibody,
		 ifdata.minW, ifdata.maxW);
    } 

    if(verbose > 1){
      printf("\n");
    }
    
  }

  //Iterate over bremss data
  if(VRReader.bremss.size() > 0){
    
    for(const pen_contextReaderVR::bremssData& bremss : VRReader.bremss){
      if(verbose > 1){
	printf("    + Bremsstrahlung splitting %s:\n"
	       "  Splitting \n"
	       "    %4u \n\n",
	       bremss.name.c_str(),
	       bremss.factor);
      }

      //Iterate over material indexes
      if(bremss.enabledMats.size() > 0){
	if(verbose > 1){
	  printf("      - Applied on materials: \n");
	}
      }
      for(const unsigned& imat : bremss.enabledMats){
	//Check if the material is used in the geometry
	if(!usedMat[imat]){
	  if(verbose > 0){
	    printf("Error: Specified material index (%u) for "
		   "Bremsstrahlung splitting '%s' is not used at current geometry.\n",
		   imat, bremss.name.c_str());
	  }
	  return UNUSED_MATERIAL;
	}

	if(verbose > 1){
	  printf("        + %u \n", imat);
	}

	//Set parameters for bodies with this material index
	for(unsigned ibody = 0; ibody < readGeometry()->getBodies(); ibody++){
	
	  if(readGeometry()->getMat(ibody) != imat) continue;
	  
	  if(ibody >= NBV){
	    if(verbose > 0){
	      printf("Error: Maximum body index for Bremssthralung"
		     " splitting reached (%u)\n",
		     NBV);
	    }
	    return MAXIMUM_BODY_IF_REACHED;
	  }

	  //Set bremss splitting in context
	  if(LFORCE[ibody][PEN_ELECTRON] || LFORCE[ibody][PEN_POSITRON]){
	    IBRSPL[ibody] = bremss.factor;	  
	  }else{
	    if(verbose > 1){
	      printf("Warning: Bremss factor specified with no "
		     "interaction forcing in body %s."
		     " Value will be ignored.\n",
		     readGeometry()->getBodyName(ibody).c_str());
	    }
	  }
	}
      }

      if(verbose > 1){
	printf("\n");
      }
      
      //Iterate over enabled bodies
      if(bremss.enabledBodies.size() > 0){
	if(verbose > 1){
	  printf("      - Applied on bodies: \n");
	}
      }
      for(const std::string& bodyAlias : bremss.enabledBodies){

	if(verbose > 1){
	  printf("        + %s \n", bodyAlias.c_str());
	}
      
	unsigned ibody = readGeometry()->getIBody(bodyAlias.c_str());
	if(ibody >= readGeometry()->getBodies()){
	  if(verbose > 0){
	    printf("Error: Body '%s' doesn't exists in loaded geometry.\n",
		   bodyAlias.c_str());
	  }
	  return UNKNOWN_BODY;
	}

	if(ibody >= NBV){
	  if(verbose > 0){
	    printf("Error: Maximum body index for Bremssthralung"
		   " splitting reached (%u)\n",
		   NBV);
	  }
	  return MAXIMUM_BODY_IF_REACHED;
	}

	//Set bremss splitting in context
	if(LFORCE[ibody][PEN_ELECTRON] || LFORCE[ibody][PEN_POSITRON]){
	  IBRSPL[ibody] = bremss.factor;	  
	}else{
	  if(verbose > 1){
	    printf("Warning: Bremss factor specified with no "
		   "interaction forcing in body %s."
		   " Value will be ignored.\n",
		   readGeometry()->getBodyName(ibody).c_str());
	  }
	}
	
      }

      if(verbose > 1){
	printf("\n");
      }

      
    }
    
  }

  return SUCCESS;
}


double pen_context::range(const double E, const pen_KPAR kpar, const unsigned M) const {

  //  This function computes the range (in cm) of particles of type KPAR
  //  and energy E in material M. For electrons and positrons, the output
  //  is the CSDA range. For photons, the delivered value is the mean free
  //  path (=inverse attenuation coefficient).

  const pen_material& mat = readBaseMaterial(M);

  //Check if energy is in range
  double EE;
  if(E < grid.EL){
    EE = grid.EL;
  }
  else if(E > grid.EU){
    EE = grid.EU;
  }
  else{
    EE = E;
  }

  //Get energy interval
  int KE;
  double XEL, XE, XEK;
  grid.getInterval(E,KE,XEL,XE,XEK);
  if(KE >= static_cast<int>(constants::NEGP-1)){
    KE = constants::NEGP-2;
    XEK = XE - KE;
  }

  // *** Electrons and positrons
  if(kpar == PEN_ELECTRON || kpar == PEN_POSITRON){
    return exp(mat.RANGEL[kpar][KE] +
	       (mat.RANGEL[kpar][KE+1] - mat.RANGEL[kpar][KE])*XEK);
  }

  // *** Photons

  if(kpar == PEN_PHOTON){

    // Rayleigh scattering
    
    int II = mat.IED[KE];
    int IU = mat.IEU[KE];
    do{
      int IT = (II+IU)/2;
      if(XEL > mat.ERA[IT]){
	II = IT;
      }
      else{
	IU = IT;
      }
    }while(IU-II > 1);

    double HMF1 = exp(mat.XSRA[II] + (mat.XSRA[II+1] - mat.XSRA[II])*
		      (XEL-mat.ERA[II])
		      /(mat.ERA[II+1]-mat.ERA[II]));

    // Compton scattering
    double HMF2 = exp(mat.SGCO[KE] + (mat.SGCO[KE+1] - mat.SGCO[KE])*XEK);

    // Photoelectric absorption
    double XS;
    GPHaT(EE, XS, mat, elements);
    double HMF3 = XS*mat.VMOL;

    // Pair production
    double HMF4;
    if(E > 1.022E6){
      HMF4 = exp(mat.SGPP[KE] + (mat.SGPP[KE+1] - mat.SGPP[KE])*XEK);
    }
    else{
      HMF4 = 0.0;
    }

    return 1.0/std::max(HMF1+HMF2+HMF3+HMF4,1.0e-35);
    
  }

  return 1.0e35;
}

double pen_context::avncol(const double E,
			   const pen_KPAR kpar,
			   const int icol,
			   const unsigned M) const {
  //Get material
  const pen_material& mat = readBaseMaterial(M);
  
  double averageNumberOfInteractions = 0.0; // Value returned by the function
  
  if (E < grid.EL || E > grid.EU){
    return averageNumberOfInteractions;
  }
  
  double AUX[constants::NEGP];
  
  if(kpar == PEN_ELECTRON){
    
    pen_betaE_interact icole = static_cast<pen_betaE_interact>(icol);

    /*
      Total stopping power
      mat.TSTPE[I] = log(mat.CSTPE[I]+mat.RSTPE[I]);
      where mat.CSTPE[I] and mat.RSTPE[I] are collision and radiation stopping power.
      Average number of hard collisions between energy nodes of grid
      grid.ET[i + 1] and grid.ET[i] is DS[i]/L[i], where DS[i] is electron
      path length on which electron lost energy DE = grid.ET[i + 1] - grid.ET[i], and
      L(Ei) is electron mean free path for given interaction.
      DS[i] = (-1/dE/ds)*DE, wher -dE/ds is stopping power.
      _
      \
      Total avarage number of interactions is /_((1/L(Ei)) / (-1/dE/ds)*DE
      i

    */

    if(icole == BETAe_HARD_ELASTIC){
      for (unsigned KE = 0; KE < constants::NEGP; ++KE){
	AUX[KE] = exp(mat.SEHEL[KE] - mat.TSTPE[KE]);
      }
    }
    else if(icole == BETAe_HARD_INELASTIC){
      for(unsigned KE = 0; KE < constants::NEGP; ++KE){
	AUX[KE] = (exp(mat.SEHIN[KE]) + exp(mat.SEISI[KE]))
	  /(mat.CSTPE[KE]+mat.RSTPE[KE]);
      }
    }
    else if(icole == BETAe_HARD_BREMSSTRAHLUNG){
      for (unsigned KE = 0; KE < constants::NEGP; ++KE){
	AUX[KE] = exp(mat.SEHBR[KE] - mat.TSTPE[KE]);
      }
    }
    else if(icole == BETAe_HARD_INNER_SHELL){
      for(unsigned KE = 0; KE < constants::NEGP; ++KE){
	AUX[KE] = exp(mat.SEISI[KE] - mat.TSTPE[KE]);
      }
    }
    else if (icole == BETAe_HARD_TOTAL){
      for (unsigned KE = 0; KE < constants::NEGP; ++KE){
	AUX[KE] = exp(mat.SETOT[KE] - mat.TSTPE[KE]);
      }
    }
    
    averageNumberOfInteractions =
      RMOMX(grid.ET, AUX, grid.EL, E, constants::NEGP, 0);
  }
  else if(kpar == PEN_PHOTON){
    pen_gamma_interact icolg = static_cast<pen_gamma_interact>(icol);

    int KE;
    double XEL, XE, XEK;
    grid.getInterval(E, KE, XEL, XE, XEK);

    if(KE >= static_cast<int>(constants::NEGP - 1)){
      KE = constants::NEGP - 2;
      XEK = XE - KE;
    }

    int II = mat.IED[KE];
    int IU = mat.IEU[KE];
    do{
      int IT = (II+IU)/2;
      if(XEL > mat.ERA[IT]){
	II = IT;
      }
      else{
	IU = IT;
      }
    }while(IU-II > 1);

    double HMF[4];
    
    // Rayleigh scattering
    HMF[static_cast<int>(GAMMA_RAYLEIGH)] =
      exp(mat.XSRA[II]+(mat.XSRA[II+1]-mat.XSRA[II])
	  *(XEL-mat.ERA[II])/(mat.ERA[II+1]-mat.ERA[II]));
    
    // Compton scattering
    HMF[static_cast<int>(GAMMA_COMPTON)] =
      exp(mat.SGCO[KE] + (mat.SGCO[KE+1] - mat.SGCO[KE])*XEK);

    // Photoelectric absorption
    double XS;
    GPHaT(E, XS, mat, elements);
    HMF[static_cast<int>(GAMMA_PHOTOELECTRIC)] = XS*mat.VMOL;

    // Pair production
    if(E > 1.022E6){
      HMF[static_cast<int>(GAMMA_PAIR_PRODUCTION)] =
	exp(mat.SGPP[KE] + (mat.SGPP[KE+1] - mat.SGPP[KE])*XEK);
    }
    else{
      HMF[static_cast<int>(GAMMA_PAIR_PRODUCTION)] = 0.0;
    }

    averageNumberOfInteractions = HMF[icolg]
      / std::max(HMF[static_cast<int>(GAMMA_RAYLEIGH)] 
		 + HMF[static_cast<int>(GAMMA_COMPTON)] 
		 + HMF[static_cast<int>(GAMMA_PHOTOELECTRIC)] 
		 + HMF[static_cast<int>(GAMMA_PAIR_PRODUCTION)], 1.0e-35);
  }
  else if(kpar == PEN_POSITRON){
    pen_betaP_interact icolp = static_cast<pen_betaP_interact>(icol);

    if(icolp == BETAp_HARD_ELASTIC){
      for(unsigned KE = 0; KE < constants::NEGP; ++KE){
	AUX[KE] = exp(mat.SPHEL[KE] - mat.TSTPP[KE]);
      }
    }
    else if(icolp == BETAp_HARD_INELASTIC){
      for(unsigned KE = 0; KE < constants::NEGP; ++KE){
	AUX[KE] = (exp(mat.SPHIN[KE]) + exp(mat.SPISI[KE]))
	  /(mat.CSTPP[KE]+mat.RSTPP[KE]);
      }
    }
    else if(icolp == BETAp_HARD_BREMSSTRAHLUNG){
      for(unsigned KE = 0; KE < constants::NEGP; ++KE){
	AUX[KE] = exp(mat.SPHBR[KE] - mat.TSTPP[KE]);
      }
    }
    else if(icolp == BETAp_HARD_INNER_SHELL){
      for(unsigned KE = 0; KE < constants::NEGP; ++KE){
	AUX[KE] = exp(mat.SPISI[KE] - mat.TSTPP[KE]);
      }
    }
    else if(icolp == BETAp_ANNIHILATION){
      for(unsigned KE = 0; KE < constants::NEGP; ++KE){
	AUX[KE] = exp(mat.SPAN[KE] - mat.TSTPP[KE]);
      }
    }
    else if(icolp == BETAp_HARD_TOTAL){
      for(unsigned KE = 0; KE < constants::NEGP; ++KE){
	AUX[KE] = exp(mat.SPTOT[KE] - mat.TSTPP[KE]);
      }
    }

    averageNumberOfInteractions =
      RMOMX(grid.ET, AUX, grid.EL, E, constants::NEGP, 0);
  }

  return averageNumberOfInteractions;

}

double pen_context::IMFP(const double E,
			 const pen_KPAR kpar,
			 const int icol,
			 const unsigned M) const {
    //Get material
    const pen_material& mat = readBaseMaterial(M);
    double inverse_mean_free_path = 1.0e-35; // Value returned by the function

    if(E < grid.EL || E > grid.EU){
      return inverse_mean_free_path;
    }

    int KE;
    double XEL, XE, XEK;
    grid.getInterval(E, KE, XEL, XE, XEK);
    
    
    //Compute inverse mean free path
    if(kpar == PEN_ELECTRON){
      pen_betaE_interact icole = static_cast<pen_betaE_interact>(icol);
        
      if (icole == BETAe_HARD_ELASTIC)
        {
	  inverse_mean_free_path = exp(mat.SEHEL[KE] + mat.DSEHEL[KE] * XEK);
        }
      else if (icole == BETAe_HARD_INELASTIC)
        {
	  inverse_mean_free_path = exp(mat.SEHIN[KE] + mat.DSEHIN[KE] * XEK);
        }
      else if (icole == BETAe_HARD_BREMSSTRAHLUNG)
        {
	  inverse_mean_free_path = exp(mat.SEHBR[KE] + mat.DSEHBR[KE] * XEK);
        }
      else if (icole == BETAe_HARD_INNER_SHELL)
        {
	  inverse_mean_free_path = exp(mat.SEISI[KE] + mat.DSEISI[KE] * XEK);
        }
    }
    else if (kpar == PEN_PHOTON){
      pen_gamma_interact icolg = static_cast<pen_gamma_interact>(icol);

      if (KE < 0)
        {
	  KE = 0;
	  XEK = XE - KE;
        }
      if (KE >= static_cast<int>(constants::NEGP - 1))
        {
	  KE = constants::NEGP - 2;
	  XEK = XE - KE;
        }

      if (icolg == GAMMA_RAYLEIGH){
	int II = mat.IED[KE];
	int IU = mat.IEU[KE];
	do{
	  int IT = (II+IU)/2;
	  if(XEL > mat.ERA[IT]){
	    II = IT;
	  }
	  else{
	    IU = IT;
	  }
	}while(IU-II > 1);
	
	inverse_mean_free_path =
	  exp(mat.XSRA[II]+(mat.XSRA[II+1]-mat.XSRA[II])*(XEL-mat.ERA[II])
	      /(mat.ERA[II+1]-mat.ERA[II]));
      }
      else if(icolg == GAMMA_COMPTON){
	inverse_mean_free_path = exp(mat.SGCO[KE] + mat.DSGCO[KE] * XEK);
      }
      else if(icolg == GAMMA_PHOTOELECTRIC){
	double XS;
	GPHaT(E, XS, mat, elements);
	inverse_mean_free_path = XS*mat.VMOL;
      }
      else if(icolg == GAMMA_PAIR_PRODUCTION){
	inverse_mean_free_path =
	  E < 1.023e6 ? 0.0 : exp(mat.SGPP[KE] + mat.DSGPP[KE] * XEK);
      }
    }
    else if (kpar == PEN_POSITRON){
      pen_betaP_interact icolp = static_cast<pen_betaP_interact>(icol);

      if (icolp == BETAp_HARD_ELASTIC)
        {
	  inverse_mean_free_path = exp(mat.SPHEL[KE] + mat.DSPHEL[KE] * XEK);
        }
      else if (icolp == BETAp_HARD_INELASTIC)
        {
	  inverse_mean_free_path = exp(mat.SPHIN[KE] + mat.DSPHIN[KE] * XEK);
        }
      else if (icolp == BETAp_HARD_BREMSSTRAHLUNG)
        {
	  inverse_mean_free_path = exp(mat.SPHBR[KE] + mat.DSPHBR[KE] * XEK);
        }
      else if (icolp == BETAp_HARD_INNER_SHELL)
        {
	  inverse_mean_free_path = exp(mat.SPISI[KE] + mat.DSPISI[KE] * XEK);
        }
      else if (icolp == BETAp_ANNIHILATION)
        {
	  inverse_mean_free_path = exp(mat.SPAN[KE] + mat.DSPAN[KE] * XEK);
        }
    }

    return std::max(inverse_mean_free_path,1.0e-35);

}

double pen_context::avninter(const double E,
			     const pen_KPAR kpar,
			     const int icol,
			     const unsigned M,
			     const bool calc_piecewise) const {
  double mean_number_of_interactions = 0.0;
  if (calc_piecewise == true)
    {
      mean_number_of_interactions = avncol(E, kpar, icol, M);
    }
  else
    {
      double inverse_mean_free_path = IMFP(E, kpar, icol, M);
      double plt = range(E, kpar, M);
      mean_number_of_interactions = inverse_mean_free_path * plt;
    }

  return mean_number_of_interactions;
}

double pen_context::getIF(const double forcerIn,
			  const pen_KPAR kpar,
			  const int icol,
			  const unsigned M,
			  const bool calc_piecewise) const {

/*
 forcer is the forcing factor, which must be larger than unity.
 TRICK: a negative input value of forcer, -FN, is assumed to
 mean that a particle with energy E=EPMAX should interact,
 on average, +FN times in the course of its slowing down to
 rest, for electrons and positrons, or along a mean free
 path, for photons. 
 1. forcer > 0, emfp = mfp / force; where mfp is mean free pass, emfp
    is effective mfp.  
    avncol = range / emfp(Emax) is avarage number of interactions on full range down to rest.
    This value is not known when creating configuration file that is why 
    selecting forcer is not intuitive.
 2. forcer < 0 is reinterpreted. forcer = abs(forcer).
    avncol0 = range(Emax) / mfp(Emax)
    forcer_new = forcer / avncol0
    emfp = mfp / forcer_new = mfp * avncol0/ forcer = mfp * range(Emax) / (forcer * mfp(Emax))
    avncol = range / emfp = range(Emax) * forcer * mfp(Emax) / (mfp * range(Emax)) = forcer * mfp(Emax) / mfp 
    avncol ~ forcer
 */

//#define PENEPMA_NEGATIVE_FORCER
  double forcer = forcerIn;
  if (forcer < -1.0e-6)
    {
#ifdef PENEPMA_NEGATIVE_FORCER // penepma 2018
      // Negative forcer values are re-interpreted.
      double E0 = grid.EU;
      // penepma
      double avncl = avninter(E0, kpar, icol, M, calc_piecewise);

      if (avncl > 1.0e-8)
        {
	  forcer = std::max(fabs(forcer) / avncl, 1.0);
        }
      else
        {
	  forcer = std::max(fabs(forcer), 1.0);
        }

#else        // penelope 2018
      if (kpar == PEN_PHOTON)
        {
	  double E0 = grid.EU;
	  double fp = 1.0 / IMFP(E0, kpar, icol, M);
	  double tst = std::abs(forcer);
	  if (fp > tst)
            {
	      forcer = fp / tst;
            }
	  else
            {
	      forcer = 1.0;
            }
        }
      else
        {
	  double E0 = grid.EU;
	  double avncl = avninter(E0, kpar, icol, M, calc_piecewise);

	  if (avncl > 1.0e-8)
            {
	      forcer = std::max(fabs(forcer) / avncl, 1.0);
            }
	  else
            {
	      forcer = std::max(fabs(forcer), 1.0);
            }
        }
#endif
    }

  return forcer;
}

pen_context* pen_context::createAuxContext(double EMAX,
					   const char* matFilename,
					   const unsigned verbose){
  
  //Create a context
  pen_context* context = nullptr;
  context = new pen_context();
  if(context == nullptr){
    if(verbose > 0){
      printf("createAuxContext: Error allocating auxilary context.\n");      
    }
    return nullptr;
  }

  //Set the number of materials to context (1)
  int errmat = context->setMats<pen_material>(1);
  if(errmat != 0){
    if(verbose > 0){
      printf("createAuxContext: Error at context "
	     "material creation: %d.\n",errmat);
    }
    delete context;
    return nullptr;
  }
  
  //Get the material
  pen_material& mat = context->getBaseMaterial(0);

  //Configure the material
  mat.C1=0.2;
  mat.C2=0.2;
  mat.WCC=1.0e3;
  mat.WCR=1.0e3;

  mat.EABS[PEN_ELECTRON] = 50.0E0;
  mat.EABS[PEN_PHOTON]   = 50.0E0;
  mat.EABS[PEN_POSITRON] = 50.0E0;

  FILE* fcontext = nullptr;
  
  if(verbose > 0){
    fcontext = fopen("rangeContext.rep","w");
    if(fcontext == nullptr){
      printf("createAuxContext: Error: unable to create "
	     "file 'rangeContext.rep'\n");
      delete context;
      return nullptr;
    }
  }
  
  int INFO = 1;
  std::string PMFILEstr[constants::MAXMAT];
  PMFILEstr[0].assign(matFilename);
  int err = context->init(EMAX,nullptr,INFO,PMFILEstr);
  if(err != 0){
    if(verbose > 0){
      printf("createAuxContext: Error: Unable to configure range context."
	     "More details can be found in 'rangeContext.rep' file.\n");
    }
    delete context;
    return nullptr;
  }

  return context;
  
}


void DIRECT(const double CDT, const double DF, double &U, double &V, double &W)
{
  //  This subroutine computes the new direction cosines of the particle
  //  velocity after a collision with given polar and azimuthal scattering
  //  angles.

  //  Input:  U,V,W ... initial direction cosines.
  //          CDT ..... cosine of the polar scattering angle.
  //          DF ...... azimuthal scattering angle (rad).

  //  Output: U,V,W ... new direction cosines.
  //          CDT and DF remain unchanged.

  //  ****  Ensure normalisation.
  double UV = U*U+V*V;
  double UVW = UV+W*W;
  if(fabs(UVW-1.0) > 1.0E-13)
  {
    double FNORM = 1.0/sqrt(UVW);
    U = FNORM*U;
    V = FNORM*V;
    W = FNORM*W;
    UV = U*U+V*V;
  }

  //  ****  Calculate new direction.

  double SDT;
  const double ABSCDT = fabs(CDT);
  if(1.0-ABSCDT > 1.0E-8)
  {
    SDT = sqrt(1.0-CDT*CDT);
  }
  else
  {
    SDT = sqrt(2.0*(1.0-ABSCDT));
  }

  if(SDT < 1.0E-13)
  {
    if(CDT < 0.0)
    {
      U = -U;
      V = -V;
      W = -W;
    }
  }
  else
  {
    const double SDTSDF = SDT*sin(DF);
    const double SDTCDF = SDT*cos(DF);
    if(UV > 1.0E-26)
    {
      const double SUV = sqrt(UV);
      const double UN = U/SUV;
      const double VN = V/SUV;
      const double WSDTCDF = W*SDTCDF;
  
      U = U*CDT+(UN*WSDTCDF-VN*SDTSDF);
      V = V*CDT+(VN*WSDTCDF+UN*SDTSDF);
      W = W*CDT-SUV*SDTCDF;
    }
    else
    {
      if(W > 0.0)
      {
        U = SDTCDF;
        V = SDTSDF;
        W = CDT;
      }
      else
      {
        U = -SDTCDF;
        V = -SDTSDF;
        W = -CDT;
      }
    }
  }
}


//  *********************************************************************
//                       SUBROUTINE DIRPOL
//  *********************************************************************
void DIRPOL(const double CDT,
	    double &DF,
	    const double CONS,
	    double &SP1,
	    double &SP2,
	    double &SP3,
	    double &U,
	    double &V,
	    double &W,
	    pen_rand& random)
{
  //     This subroutine computes the direction cosines _and_ the Stokes
  //  parameters of a polarised photon after scattering with a given polar
  //  angle.

  //  Input:  U,V,W ... initial direction cosines.
  //          SP1,SP2,SP3 ... initial Stokes parameters.
  //          CDT ..... cosine of the polar scattering angle.
  //          CONS .... constant in the PDF of the azimuthal angle.
  //  Output: U,V,W ... new direction cosines.
  //          SP1,SP2,SP3 ... new Stokes parameters.
  //          DF ...... azimuthal scattering angle.
  //          CDT and CONS remain unchanged.
  
  //  ****  Sampling the azimuthal scattering angle.

  double CDT2 = CDT*CDT;
  double CDT21 = CDT2+1.0;
  double PHA = CDT21+CONS;
  double PHB = 1.0-CDT2;
  double SP0MAX = PHA+PHB*sqrt(SP1*SP1+SP3*SP3+1.0E-35);
  double SDT, SDF, CDF, S2DF, C2DF, SP3P, SP0P;
  do
  {
    DF = random.rand()*constants::TWOPI;
    SDF = sin(DF);
    CDF = cos(DF);
    S2DF = 2.0*SDF*CDF;
    C2DF = CDF*CDF-SDF*SDF;
    SP3P = S2DF*SP1+C2DF*SP3;  // Stokes parameter with new zero-azimuth.
    SP0P = PHA-PHB*SP3P;
  }while(random.rand()*SP0MAX > SP0P);

  //  ****  Calculate new Stokes parameters.

  double SP1P = C2DF*SP1-S2DF*SP3;  // Stokes parameter with new zero-azimuth.
  double RSP0 = 1.0/SP0P;

  SP1 = 2.0*CDT*SP1P*RSP0;
  SP2 = (2.0+CONS)*CDT*SP2*RSP0;
  SP3 = (CDT21*SP3P-PHB)*RSP0;

  //  ****  Ensure normalisation.

  double UV = U*U+V*V;
  double UVW = UV+W*W;
  if(fabs(UVW-1.0) > 1.0E-13)
  {
    double FNORM = 1.0/sqrt(UVW);
    U = FNORM*U;
    V = FNORM*V;
    W = FNORM*W;
    UV = U*U+V*V;
  }

  //  ****  Calculate new direction.

  if(1.0-fabs(CDT) > 1.0E-8)
  {
    SDT = sqrt(PHB);
  }
  else
  {
    SDT = sqrt(2.0*(1.0-fabs(CDT)));
  }

  if(SDT < 1.0E-13)
  {
    if(CDT < 0.0)
    {
      U = -U;
      V = -V;
      W = -W;
    }
  }
  else
  {
    double SDTSDF = SDT*SDF;
    double SDTCDF = SDT*CDF;
    if(UV > 1.0E-26)
    {
      double SUV = sqrt(UV);
      double UN = U/SUV;
      double VN = V/SUV;
      U = U*CDT+(UN*W*SDTCDF-VN*SDTSDF);
      V = V*CDT+(VN*W*SDTCDF+UN*SDTSDF);
      W = W*CDT-SUV*SDTCDF;
    }
    else
    {
      if(W > 0.0)
      {
        U = SDTCDF;
        V = SDTSDF;
        W = CDT;
      }
      else
      {
        U = -SDTCDF;
        V = -SDTSDF;
        W = -CDT;
      }
    }
  }
}

//  *********************************************************************
//                       SUBROUTINE RELAX
//  *********************************************************************
void RELAX(const pen_elementDataBase& elements,
	   const pen_material& mat,
	   pen_particleState& state,
	   const int ICOL,
	   const int MODER,
	   const int IZ,
	   const int IS,
	   int &KS,
	   pen_particleStack<pen_particleState>& stackE,
	   pen_particleStack<pen_state_gPol>& stackG,
	   pen_rand& random)
{
  //  This subroutine simulates the relaxation of a singly ionised atom of
  //  the element IZ with a vacancy in the IS shell (the K shell or an L, M
  //  or N subshell). This initial vacancy is filled by electrons from
  //  outer shells through radiative and non-radiative transitions, which
  //  may produce additional vacancies.
  
  //  We use the following notation to designate the possible transitions:
  //  *  Radiative: IS0-IS1 (an electron from the IS1 shell fills the
  //     vacancy in the IS0 shell, leaving a hole in the IS1 shell).
  //  *  Non-radiative: IS0-IS1-IS2 (an electron from the IS1 shell fills
  //     the vacancy in the IS0 shell, and the released energy is taken
  //     away by an electron in the IS2 shell; this process leaves two
  //     holes, in the IS1 and IS2 shells).
  //  The de-excitation cascade (i.e. the set of transitions that occur for
  //  a given initial vacancy) is sampled from the transition probabilities
  //  contained in the Livermore Evaluated Atomic Data Library (EADL). The
  //  energy of the radiation emitted in each transition is read from the
  //  PENELOPE database.

  //  The simulation of the de-excitation cascade is discontinued either
  //  when the K to N shells have been filled up or when there is not
  //  enough energy to produce 'active' radiation (with energy larger than
  //  EABS). The excitation energy of the residual ion is assumed to be
  //  deposited locally. We disregard the emission and transport of soft
  //  x-rays and slow electrons, whose energies are less than the binding
  //  energy of the N1 shell of the heavier element in the medium. This
  //  sets a lower limit for the energy interval that can be covered by the
  //  simulation program in a consistent way.
  
  //  De-excitation data for the loaded elements are stored in the common
  //  block /CRELAX/, in a form designed to minimise the amount of memory
  //  and to facilitate the random sampling. The quantities in the common
  //  block are the following:
  //  IFIRST(99,16) ... de-excitation data for a vacancy in the shell IS of
  //     the element IZ start at the position K=IFIRST(IZ,IS) in the
  //     storage arrays. The allowed values for IS are 1 to 16 (K shell
  //     and L, M and N subshells).
  //  ILAST(99,16) ... the de-excitation data for a vacancy in the shell
  //     IS of the element IZ end at the position K=ILAST(IZ,IS) in the
  //     storage arrays.
  //  IS1(K), IS2(K) ... shells that are active in the transition (see the
  //     shell label code below). For radiative transitions, IS2(K)=0.
  //  P(K) ... relative probability for the transition IS-IS1(K)-IS2(K).
  //  ET(K) ... energy of the secondary particle emitted in the transition.
  //  F(K), IAL(K) ... cutoff and alias values (Walker's sampling method).
  
  //  ---------------------------------------------------------------------
  //  Label code IS for electron shells:
  //      1 = K  (1s1/2),     11 = N2 (4p1/2),     21 = O5 (5d5/2),
  //      2 = L1 (2s1/2),     12 = N3 (4p3/2),     22 = O6 (5f5/2),
  //      3 = L2 (2p1/2),     13 = N4 (4d3/2),     23 = O7 (5f7/2),
  //      4 = L3 (2p3/2),     14 = N5 (4d5/2),     24 = P1 (6s1/2),
  //      5 = M1 (3s1/2),     15 = N6 (4f5/2),     25 = P2 (6p1/2),
  //      6 = M2 (3p1/2),     16 = N7 (4f7/2),     26 = P3 (6p3/2),
  //      7 = M3 (3p3/2),     17 = O1 (5s1/2),     27 = P4 (6d3/2),
  //      8 = M4 (3d3/2),     18 = O2 (5p1/2),     28 = P5 (6d5/2),
  //      9 = M5 (3d5/2),     19 = O3 (5p3/2),     29 = Q1 (7s1/2),
  //     10 = N1 (4s1/2),     20 = O4 (5d3/2),     30 = outer shells.
  //  ---------------------------------------------------------------------
  
  //const double TWOPI = PI+PI;
  double OCUP[30];
  int ISV[30];
  double PTIME[16][256];
  OCUP[0] = 2.0E0;
  OCUP[1] = 2.0E0;
  OCUP[2] = 2.0E0;
  OCUP[3] = 4.0E0;
  OCUP[4] = 2.0E0;
  OCUP[5] = 2.0E0;
  OCUP[6] = 4.0E0;
  OCUP[7] = 4.0E0;
  OCUP[8] = 6.0E0;
  OCUP[9] = 2.0E0;
  OCUP[10] = 2.0E0;
  OCUP[11] = 4.0E0;
  OCUP[12] = 4.0E0;
  OCUP[13] = 6.0E0;
  OCUP[14] = 6.0E0;
  OCUP[15] = 8.0E0;
  OCUP[16] = 2.0E0;
  OCUP[17] = 2.0E0;
  OCUP[18] = 4.0E0;
  OCUP[19] = 4.0E0;
  OCUP[20] = 6.0E0;
  OCUP[21] = 6.0E0;
  OCUP[22] = 8.0E0;
  OCUP[23] = 2.0E0;
  OCUP[24] = 2.0E0;
  OCUP[25] = 4.0E0;
  OCUP[26] = 4.0E0;
  OCUP[27] = 6.0E0;
  OCUP[28] = 2.0E0;
  OCUP[29] = 1.0E9;

  int KPARS;
  
  //  ****  Initialisation.

  if(IZ < 3 || IS > 16){ return;}
  //  ****  If the shell ionisation energy is less than ECUTR, the cascade
  //        is not followed.
  if(elements.EB[IZ-1][IS-1] < mat.ECUTR){ return;}
  
  double PAGE0 = state.PAGE;
  
  int NV = 0;
  for (int I = 0; I < 30; I++)
    {
      ISV[I] = 0;
    }
  int ISA = IS;
  state.PAGE = PAGE0;

  //  ****  Next transition.

  bool Eixir = false;
  while(!Eixir)
  {
    Eixir = true;
    int KF = elements.IFIRST[IZ-1][ISA-1];
    int KL = elements.ILAST[IZ-1][ISA-1];
    
    if (KL <= KF)
      {
        Eixir = true;
        break;
      }
      
    double RN, TST;
    int K1;
    int JS1, JS2;
    bool Eixir2 = false;
    while (!Eixir2)
      {
        //  ****  Walker's sampling algorithm.
        Eixir2 = true;
        RN = random.rand()*(double)(KL-KF+1);
	
        K1 = (int)RN;
        TST = RN-(double)K1;
        if(TST > elements.F[KF+K1-1])
          {
            KS = elements.IAL[KF+K1-1];
          }
        else
          {
            KS = KF+K1;
          }

        JS1 = elements.IS1[KS - 1];
        JS2 = elements.IS2[KS - 1];
//  ****  Vacancies in the intervening subshells.
        double FREJ;
        if (ISV[JS1 - 1] > 0)
          {
            FREJ = (OCUP[JS1 - 1] - ISV[JS1 - 1]) / OCUP[JS1 - 1];
          }
        else
          {
            FREJ = 1.0E0;
          }
        if (JS2 > 0)
          {
            FREJ = FREJ * ((OCUP[JS2 - 1] - ISV[JS2 - 1]) / OCUP[JS2 - 1]);
          }
	
        if (FREJ < 1.0E0)
          {
            if (random.rand() > FREJ)
              {
                Eixir2 = false;
                continue;
              }
          }
      }
    
    //  ****  If MODER=0, control is returned to the calling program after
    //  determining the first transition, KS. Useful for testing the random
    //  sampling. For normal operation, we can comment out the following
    //  statement.
    if(MODER == 0){ return;}
    //  ****  If LAGE=.TRUE., the particle age is recorded.
    if(state.LAGE)
    {
      if (ISV[ISA - 1] > 1)
        {
          state.PAGE = state.PAGE-(elements.ALW[IZ-1][ISA-1]/double(ISV[ISA - 1]))*log(random.rand());
        }
      else
        {
          state.PAGE = state.PAGE-elements.ALW[IZ-1][ISA-1]*log(random.rand());
        }
    }
     
    //  ****  Fluorescence radiation.
      
    if(JS2 == 0)
    {
      KPARS = PEN_PHOTON;
      if(JS1 < 17)
      {
        if(elements.EB[IZ-1][JS1-1] > mat.ECUTR)
        {
          NV = NV+1;
          ISV[JS1 - 1] = ISV[JS1 - 1] + 1;
          PTIME[JS1 - 1][ISV[JS1 - 1] - 1] = state.PAGE;
        }
      }
    }
    else
    {
      KPARS = PEN_ELECTRON;
      if(JS1 < 17)
      {
        if(elements.EB[IZ-1][JS1-1] > mat.ECUTR)
        {
          NV = NV+1;
          ISV[JS1 - 1] = ISV[JS1 - 1] + 1;
          PTIME[JS1 - 1][ISV[JS1 - 1] - 1] = state.PAGE;
        }
      }
      if(JS2 < 17)
      {
        if(elements.EB[IZ-1][JS2-1] > mat.ECUTR)
        {
          NV = NV+1;
          ISV[JS2 - 1] = ISV[JS2 - 1] + 1;
          PTIME[JS2 - 1][ISV[JS2 - 1] - 1] = state.PAGE;
        }
      }
    }
      
      //  ****  The emitted particle is stored in the secondary stack when
      //        its energy ET(K) is greater than EABS.
    
    double WS, SDTS, DF;
    if(elements.ET[KS-1] > mat.EABS[KPARS])
    {
      //  ****  Initial direction (isotropic).
      WS = -1.0+2.0*random.rand();
      
      SDTS = sqrt(1.0-WS*WS);
      DF = constants::TWOPI*random.rand();
      
      if(KPARS == PEN_ELECTRON)
	    {
	      // *** Electron emited
	      pen_particleState storeState;

	      storeState.E = elements.ET[KS-1];

	      storeState.X = state.X;
	      storeState.Y = state.Y;
	      storeState.Z = state.Z;
	  
	      storeState.U = cos(DF)*SDTS;
	      storeState.V = sin(DF)*SDTS;
	      storeState.W = WS;

	      storeState.WGHT = state.WGHT;
	      storeState.IBODY = state.IBODY;
	      storeState.MAT = state.MAT;	      

	      storeState.ILB[0] = state.ILB[0]+1;
	      storeState.ILB[1] = KPARS;
	      storeState.ILB[2] = ICOL;
        storeState.ILB[3] = IZ*1000000+ISA*10000+JS1*100+JS2;
	      storeState.ILB[4] = state.ILB[4];

	      storeState.LAGE = state.LAGE;
	      storeState.PAGE = state.PAGE;
	  
	      stackE.store(storeState);
	    }
      else
	    {
	      // *** Photon emited
	      pen_state_gPol storeState;

	      storeState.E = elements.ET[KS-1];

	      storeState.X = state.X;
	      storeState.Y = state.Y;
	      storeState.Z = state.Z;
	  
	      storeState.U = cos(DF)*SDTS;
	      storeState.V = sin(DF)*SDTS;
	      storeState.W = WS;

	      storeState.WGHT = state.WGHT;
	      storeState.IBODY = state.IBODY;
	      storeState.MAT = state.MAT;	      

	      storeState.ILB[0] = state.ILB[0]+1;
	      storeState.ILB[1] = KPARS;
	      storeState.ILB[2] = ICOL;
	      storeState.ILB[3] = IZ*1000000+ISA*10000+JS1*100+JS2;
	      storeState.ILB[4] = state.ILB[4];

	      storeState.LAGE = state.LAGE;
	      storeState.PAGE = state.PAGE;
	  
	      stackG.store(storeState);
	    }
      
      if(penGetError() != PEN_SUCCESS){ return;}
    }
      
//
//  ****  Are there any K-, L-, M- or N-vacancies unfilled?
//
    if (NV > 0)
      {
        Eixir2 = false;
        bool GOTO1, GOTO3, GOTO4;
        while (!Eixir2)
          {
            Eixir2 = true;
            GOTO1 = false;
            GOTO4 = false;
            GOTO3 = false;
            if (elements.EB[IZ - 1][ISA - 1] < mat.ECUTR)
              {
                GOTO4 = true;
                break;
              }
            if (ISV[ISA - 1] > 0) // Vacancies in the current subshell.
              {
                state.PAGE = PTIME[ISA - 1][ISV[ISA - 1] - 1];
                ISV[ISA - 1] = ISV[ISA - 1] - 1;
                NV = NV - 1;
                GOTO1 = true;
                break;
              }
            else
              {
                for (int IST = ISA + 1 - 1; IST < 16; IST++)  // Outer subshells.
                   {
                     if (ISV[IST] > 0)
                       {
                         ISA = IST + 1;
                         GOTO3 = true;
                         break;
                       }
                   }
              }
            if (GOTO3)
              {
                Eixir2 = false;
                continue;
              }
          }
        if (GOTO4)
          {
            break;
          }
        if (GOTO1)
          {
            Eixir = false;
            continue;
          }
      }
  }       // 1 CONTINUE. Eixir
  state.PAGE = PAGE0;

}

//  *********************************************************************
//                       SUBROUTINE GPHaT
//  *********************************************************************
void GPHaT(const double E, double &XS, const pen_material& mat, const pen_elementDataBase& elemDB)
{
  //  Delivers the photoelectric cross section XS (in cm**2) for photons of
  //  energy E in material M.


  double XEL = log(E);
  XS = 0.0;
  double DEE, PCSL;
  int IZZ, I, IU, IT;
  for(int IEL = 0; IEL < mat.NELEM; IEL++)
  {
    IZZ = mat.IZ[IEL];
    I = elemDB.IPHF[IZZ-1];
    IU = elemDB.IPHL[IZZ-1];

    while(IU-I > 1)
    {
      IT = (I+IU)/2;
      if(XEL > elemDB.EPH[IT])
      {
        I = IT;
      }
      else
      {
        IU = IT;
      }
    }
    DEE = elemDB.EPH[I+1]-elemDB.EPH[I];
    if(DEE > 1.0E-15)
    {
      PCSL = elemDB.XPH[I][0]+(elemDB.XPH[I+1][0]-elemDB.XPH[I][0])*(XEL-elemDB.EPH[I])/DEE;
    }
    else
    {
      PCSL = elemDB.XPH[I][0];
    }
    XS = XS+mat.STF[IEL]*exp(PCSL);
  }
}
