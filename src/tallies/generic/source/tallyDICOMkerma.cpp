
//
//
//    Copyright (C) 2021-2023 Universitat de València - UV
//    Copyright (C) 2021-2023 Universitat Politècnica de València - UPV
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
//        sanolgi@upvnet.upv.es (Sandra Oliver Gil)
//        vicent.gimenez.alventosa@gmail.com  (Vicent Giménez Alventosa)
//        vicente.gimenez@uv.es (Vicente Giménez Gómez)
//    
//


#ifdef _PEN_USE_DICOM_
#include "tallyDICOMkerma.hh"

pen_tallyDICOMkerma::pen_tallyDICOMkerma() : pen_genericTally(USE_JUMP | USE_STEP),
					     contVol(nullptr),
					     contourVox(nullptr),
					     ncontours(0)
{
  //Kerma
  setResultsGenerator<0>
    ([this](const unsigned long long nhists) -> penred::measurements::results<double, 3>{
      return this->tallyKerma.generateResults<0>(nhists);
    });

  //DVH
  setResultsGenerator<1>
    ([this](const unsigned long long nhists) -> penred::measurements::results<double, 2>{

      double invn = 1.0/static_cast<double>(nhists);

      //Create results
      const unsigned long nCont = static_cast<unsigned long>(ncontours);
      penred::measurements::results<double, 2> results;
      results.initFromLists
	({DVHnbins, nCont},
	 {penred::measurements::limitsType(0.0, DVHmaxDose),
	  penred::measurements::limitsType(0.0, static_cast<double>(nCont))
	 });
	  
      results.description =
	"PenRed: DVH report\n\n"
	"  Contours:\n";
      for(unsigned long j = 0; j < nCont; ++j){
	results.description += std::to_string(j) + " " + contNames[j] + "\n";
      }
      results.description += "\n";
  
      results.setDimHeader(0, "Dose (Gy)");
      results.setDimHeader(1, "Contour");
      results.setValueHeader("Volume (%)");

      //Read kerma value from tallyKermaTrackLength
      const double* kermaValue;
      kermaValue = nullptr;
  
      if(tallyKerma.enabledCart())
	{
	  kermaValue = tallyKerma.readCartesians();
	}
      else
	{
	  results.description +=
	    "Error at kerma tally reading. \n"
	    "Kerma value can not be read from tallyKermaTrackLength.\n";
	      
	    printf("pen_DICOMkerma: Error at kerma tally reading. "
		   "Kerma value can not be read from tallyKermaTrackLength.\n");
	  return results;
	}

      //Cartesian element volume
      const double volumeCart = dx*dy*dz;

      //Calculate number of voxels in each energy interval for each contour
      for(long int i = 0; i < nbin; ++i){
	//Get voxel contour index
	int icont = contourVox[i];
	if(icont < 0)
	  continue; //This voxel is not in any contour
          
	for(long int j = DVHnbins-1; j >= 0; --j){
	  if(kermaValue[i]*DVHfactor*invn/volumeCart >= static_cast<double>(j)*DVHbinWidth-1.0e-20){
	    results.data[icont*DVHnbins + j] += 1.0;
	    break;
	  }
	}
      }

      //Get the cummulative DVH
      for(unsigned icont = 0; icont < nCont; ++icont){
	for(long int j = DVHnbins-2; j >= 0; --j){
	  results.data[icont*DVHnbins + j] += results.data[icont*DVHnbins + j+1];
	}
      }
	
      for(long unsigned icont = 0; icont < nCont; ++icont){
	for(long unsigned j = 0; j < DVHnbins; ++j){
	  results.data[icont*DVHnbins + j] *= 100.0*voxVol/contVol[icont];
	}
      }
	
      return results;	
    });    
}

void pen_tallyDICOMkerma::flush(){
    tallyKerma.flush();
}

void pen_tallyDICOMkerma::tally_jump(const unsigned long long nhist,
					   const pen_KPAR kpar,
					   const pen_particleState& state,
					   const double ds){
    
    tallyKerma.tally_jump(nhist,kpar,state,ds);
}


void pen_tallyDICOMkerma::tally_step(const unsigned long long nhist,
					   const pen_KPAR kpar,
					   const pen_particleState& state,
					   const tally_StepData& stepData){
    
    tallyKerma.tally_step(nhist,kpar,state,stepData);
}


int pen_tallyDICOMkerma::sumTally(const pen_tallyDICOMkerma& tally){
    
    tallyKerma.sumTally(tally.tallyKerma);
    return 0;
}


int pen_tallyDICOMkerma::configure(const wrapper_geometry& geometry,
				    const abc_material* const materials[constants::MAXMAT],
				    const pen_parserSection& config,
				    const unsigned verbose){
    
    //Check if this geometry is a DICOM based geometry
  const pen_dicomGeo* pDICOMgeo = dynamic_cast<const pen_dicomGeo*>(&geometry);
  if(pDICOMgeo == nullptr){
    if(verbose > 0)
      printf("pen_DICOMkerma::configure: Error: DICOM kerma distribution "
	     "tally requires a DICOM based geometry.\n");
    return 1;
  }  
    
    
  ///////////////////////////////
  //  Read dicom information   //
  ///////////////////////////////  
    
  //Get DICOM
  const pen_dicom& dicom = pDICOMgeo->readDicom();
  
  nx = dicom.getNX();
  ny = dicom.getNY();
  nz = dicom.getNZ();

  nxy = nx*ny;
  nbin = dicom.getNVox();

  dx = dicom.getDX();
  dy = dicom.getDY();
  dz = dicom.getDZ();

  xmin = ymin = zmin = 0.0;

  voxVol = dicom.getVoxVol();
  

  //Get number of contours
  ncontours = dicom.nContours();
    
  //Allocate memory
  contVol       = static_cast<double*>(malloc(sizeof(double)*ncontours));
  
  //Save contour names
  for(long int i = 0; i < ncontours; ++i){
     contNames.push_back(dicom.contour(i).name);
  }

  //Get voxel contour information
  contourVox = dicom.readContour();
    
  for(long int i = 0; i < ncontours; ++i)
      contVol[i]   = 0.0;
      
  
  if(ncontours > 0){
      for(long int i = 0; i < nbin; ++i){
          int icont = contourVox[i];
          if(icont >= 0){
              contVol[icont] += voxVol;
          }
      }
  }
  ///////////////////////////////
  //   Read DVH information    //
  /////////////////////////////// 

  const double Gy2eV = 1.0/1.60217662E-16;
  //const double eV2Gy = 1.60217662E-16;

  int err = 0;  
    
  err = config.read("prescribedDose", prescribedDose);
  if(err != INTDATA_SUCCESS){
        prescribedDose = 1.0;
  }
  prescribedDose = fabs(prescribedDose);
  prescribedDose *= Gy2eV;
    
  err = config.read("MCkerma2Dose", DVHfactor);
  if(err != INTDATA_SUCCESS){
      DVHfactor = 1.0;
  }

  err = config.read("DVH-maxDose", DVHmaxDose);
  if(err != INTDATA_SUCCESS){
      DVHmaxDose = 3.0*prescribedDose;
  }else{
      DVHmaxDose = fabs(DVHmaxDose);
      //DVHmaxDose *= Gy2eV; 
  }


  int auxBins;
  err = config.read("DVH-bins", auxBins);
  if(err != INTDATA_SUCCESS){
      auxBins = 200;
  }
    
  if(auxBins <= 0){
      DVHnbins = 200;
  }else{
      DVHnbins = static_cast<unsigned long>(auxBins); 
  }
    
  DVHbinWidth = DVHmaxDose/static_cast<double>(DVHnbins);
  
  if(verbose > 1){
      printf("Number of contours: %ld\n",ncontours);
      if(ncontours > 0){
      printf("Contour names:\n");
      for(const std::string& contName : contNames)
      printf(" %s\n",contName.c_str());
      
      printf("Contour volumes (cm^3):\n");
      for(long int icont = 0; icont < ncontours; icont++){

        printf(" %15s:  %12.5E\n",contNames[icont].c_str(),contVol[icont]);
      }
      
        
      printf("\n");
      printf("DVH bins          : %lu\n",DVHnbins);
      printf("DVH bin width (Gy): %12.5E\n",DVHbinWidth);
      }
        
  }

    
      
    ///////////////////////////////
    //    Get kerma section      //
    /////////////////////////////// 
    
    if(verbose > 1)
    {
        printf("Configuring tally kerma track length from DICOM kerma tally\n");
    }
    
    pen_parserSection kermaSection;

    err = config.readSubsection("kerma",kermaSection);

    if(err != INTDATA_SUCCESS){
        if(verbose > 0){
            printf("pen_DICOMkerma: configure: Error: Configuration 'kerma' section does not exist.\n");
        }
        return -1;
    }

    //Only cartesian form is valid for DVH calculations of this tally.
    
    err = kermaSection.remove("cylindrical");
    
    err = kermaSection.remove("spherical");
    
    
    //Set cartesian variables to DICOM mesh
    err = kermaSection.set("cartesian/nx", static_cast<int>(nx));
    err = kermaSection.set("cartesian/ny", static_cast<int>(ny));
    err = kermaSection.set("cartesian/nz", static_cast<int>(nz));
    
    
    err = kermaSection.set("cartesian/xmin", xmin);
    err = kermaSection.set("cartesian/ymin", ymin);
    err = kermaSection.set("cartesian/zmin", zmin);
    
    err = kermaSection.set("cartesian/xmax", dx*nx);
    err = kermaSection.set("cartesian/ymax", dy*ny);
    err = kermaSection.set("cartesian/zmax", dz*nz);
    
    //Set sub tally name and thread
    tallyKerma.setName(readName());
    tallyKerma.setThread(getThread());
    
    if(tallyKerma.configure(geometry, materials, kermaSection, verbose) != 0){
        if(verbose > 0)
        {
            printf("pen_DICOMkerma: configure: Error at kerma tally configuration.\n");
        }
        return -2;
    }

    //Register kerma tally dump to be dumped
    addSubDump(tallyKerma);
    
 return 0;
}



void pen_tallyDICOMkerma::saveData(const unsigned long long nhist) const{

  // ** Save track length information
  tallyKerma.saveData(nhist);

  //*****************************
  //*     DVH information       *
  //*****************************

  FILE* out;
  double invn = 1.0/static_cast<double>(nhist);  
    
  out = fopen("DVH-info.dat","w");

  //Read kerma value from tallyKermaTrackLength
  const double* kermaValue;
  kermaValue = nullptr;
  
  if(tallyKerma.enabledCart())
  {
      kermaValue = tallyKerma.readCartesians();
  }
  else
  {
    printf("pen_DICOMkerma: saveData: Error at kerma tally reading. Kerma value can not be read from tallyKermaTrackLength.\n");
    return;
  }
  
  
  std::vector<std::vector<unsigned long>> DVHs;
  //Create one DVH for each contour
  DVHs.resize(ncontours);

  for(auto& DVH : DVHs){
      DVH.resize(DVHnbins+1); //+1 to store voxels with energy out of range 
      for(unsigned long& i : DVH)
          i = 0;
  }
  
  //Cartesian element volume
  const double volumeCart = dx*dy*dz;
  
  //Calculate number of voxels in each energy interval for each contour
  for(long int i = 0; i < nbin; ++i){
      //Get voxel contour index
      int icont = contourVox[i];
      if(icont < 0)
          continue; //This voxel is not in any contour
          
        for(long int j = DVHnbins; j >= 0; --j){
            if(kermaValue[i]*DVHfactor*invn/volumeCart >= static_cast<double>(j)*DVHbinWidth-1.0e-20){
                ++DVHs[icont][j];
                break;
            }
        }
  }
  
  
   //Get the cummulative DVH
  for(auto& DVH : DVHs){
      for(long int j = DVHnbins-1; j >= 0; --j){
          DVH[j] += DVH[j+1];
      }
  }
  
  
  //Print DVH
  fprintf(out,"#Contours:\n");
  for(int j = 0; j < ncontours; ++j)
      fprintf(out,"# %3d: %s\n", j, contNames[j].c_str());
  fprintf(out,"#\n");

  fprintf(out,"#    Dose low     |    Dose mid     | Volumes by contour (%%):\n");
  fprintf(out,"#      (Gy)       |      (Gy)       |");
  for(int j = 0; j < ncontours; ++j)
      fprintf(out,"       %3d       |",j);
  fprintf(out,"\n");
  
  for(long unsigned i = 0; i < DVHnbins; ++i){
      fprintf(out,"%15.5E   %15.5E  ",  
               static_cast<double>(i)*DVHbinWidth,
              (static_cast<double>(i)+0.5)*DVHbinWidth);
      for(int j = 0; j < ncontours; ++j)
        fprintf(out," %15.5E  ",100.0*static_cast<double>(DVHs[j][i])*voxVol/contVol[j]);
      
      fprintf(out,"\n");
  }
  
  fclose(out);

}

REGISTER_COMMON_TALLY(pen_tallyDICOMkerma)

#endif
