
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

#include "brachySource.hh"



void brachy_specificSampler::skip(const unsigned long long dhists){
  psf.skip(dhists);
}


int brachy_specificSampler::configure(double& Emax,
				      const pen_parserSection& config,
				      const unsigned nthreads,
				      const unsigned verbose){
    
  //First, initialize the phase space file sampler
  psf.setThread(getThread());
  
  //Get configuration subsections
  pen_parserSection psfSection;
  
  int err = config.readSubsection("psf", psfSection);
  if(err != INTDATA_SUCCESS){
    printf("brachySource:configure: Error: Unable to extract psf configuration. PSF configuration is needed for BRACHY source.\n");
    printf("error code: %d\n",err);
    return err;
  }
  
  err = psf.configure(Emax,psfSection,nthreads,verbose);
  if(err != 0){
    if(verbose > 0)
      printf("brachySource:configure: Error: Unable to configure psf.\n"
	     "               error code: %d\n",err);
    return err;
  }
  
  if(verbose > 1)
    printf("Phase space file source configured.\n");
  
  err = config.read("seedID", seedID);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0 && getThread() == 0){
      printf("brachySource:configure: Warning: Unable to read 'seedID' field. Integrer expected.\n"
    "                                 All kind of seeds will be read.\n");
    }
  }

  //Try to read seed rotation option
  if(config.read("seed-rotation", seedRotation) != INTDATA_SUCCESS)
    seedRotation = true;

  if(verbose > 1)
    printf("Seed rotation is %s.\n", seedRotation ? "enabled" : "disabled");
  
  return 0;
}

void brachy_specificSampler::sample(pen_particleState& state,
				 pen_KPAR& genKpar,
				 unsigned long long& dhist,
				 pen_rand& random){
  
        pen_particleState localstate;
        
        //PSF must be located at 0, 0, 0
        //then, the corresponding trnaslation must be applied in the PSF configure
        psf.sample(localstate,genKpar,dhist,random);
            
        double pos[3] = {localstate.X,localstate.Y,localstate.Z};
        double dir[3] = {localstate.U,localstate.V,localstate.W};
        
        //Get the seed information
        
        
        // Get random number
        double rand = random.rand();
        unsigned seedIndex = seeki(weights.data(),rand, weights.size());

	if(seedRotation){
	  matmul3D(rotations[seedIndex].data(),pos);
	  matmul3D(rotations[seedIndex].data(),dir);
	}
        
        
        localstate.X=pos[0];
        localstate.Y=pos[1];
        localstate.Z=pos[2];
        localstate.U=dir[0];
        localstate.V=dir[1];
        localstate.W=dir[2];

        
        localstate.X += positions[seedIndex].x;
        localstate.Y += positions[seedIndex].y;
        localstate.Z += positions[seedIndex].z;
        
        state = localstate;
}


void brachy_specificSampler::updateGeometry(const wrapper_geometry* geometryIn){
    
    //Check if the geometry is a DICOM
    const pen_dicomGeo* pDicomGeo = geometryIn->convertTo<pen_dicomGeo>();

    if(pDicomGeo == nullptr){
      //This geometry type is not a DICOM, check if it contains a DICOM
      size_t geoPos;
      pDicomGeo = geometryIn->getInternalGeoType<pen_dicomGeo>(geoPos);
    }
    
    if(pDicomGeo != nullptr){
        const pen_dicom& dicom = pDicomGeo->readDicom();
        unsigned long nSeeds = dicom.nSeeds();
        size_t nPositions = 0;

        for(size_t i=0;i<nSeeds;++i){
           const pen_seed& seed = dicom.seed(i);
           if(seedID < 0 || seedID == seed.ID){
               unsigned long seednPositions = seed.nPoints();
               for(unsigned long j = 0; j < seednPositions; ++j){
                   if(seed.getWeight(j) > 1.0E-8){
                        nPositions++;
                   }
                }
            }
        }
        if(nPositions == 0)
        {
            printf("brachySource:updateGeometry: Error: The DICOM geometry has not any seed.\n");
            return;
        }
        
        
        printf("Number of positions read: %li.\n", nPositions);
        
        positions.resize(nPositions);
        rotations.resize(nPositions);
        weights.resize(nPositions+1);
        
        
        size_t icp = 0; //Index of next control point
        double wghtSum = 0.0; //Index to sum weight values
        weights[0] = 0.0;
        for(size_t i=0;i<nSeeds;++i){
           const pen_seed& seed = dicom.seed(i);
           if(seedID < 0 || seedID == seed.ID){
               
               unsigned long seednPositions = seed.nPoints();
               
               double pos[3];
               double dir[3];
               
               for(unsigned long j = 0; j < seednPositions; ++j){
                   if(seed.getWeight(j) > 1.0E-8)
                   {
                        seed.getPos(pos,j);
                        seed.getDirection(dir,j);
                        positions[icp] = vector3D<double>(pos[0],pos[1],pos[2]);
                        weights[icp+1] = seed.getWeight(j);
                        wghtSum += weights[icp+1];
                        rollAlign(dir[0],dir[1],dir[2],0.0,rotations[icp].data());
                        ++icp;
                   }
               }
           }
        }
        
        if(icp != nPositions){
                printf("brachySource:updateGeometry: Error: Inconsistent DICOM readings"
                ". Please contact the developer.\n");
        }
        
        //Obtain cummulative function of weights
        //Loop for other elements
        for(size_t i=0; i<icp; ++i)
        {
            weights[i+1] = weights[i] + weights[i+1]/wghtSum;
        }
        
    }else{
        
        printf("brachySource:updateGeometry: Warning: DICOM seeds can not be read.\n"
               "Source position has been set to (0,0,0).\n");    
        positions.resize(1);
        rotations.resize(1);
        weights.resize(2);
        
        positions[0] = vector3D<double>(0.0,0.0,0.0);
        
        //Set rotation to identity matrix
        rotations[0][0] = 1.0; rotations[0][1] = 0.0; rotations[0][2] = 0.0;
        rotations[0][3] = 0.0; rotations[0][4] = 1.0; rotations[0][5] = 0.0;
        rotations[0][6] = 0.0; rotations[0][7] = 0.0; rotations[0][8] = 1.0;
        
        weights[0] = 0.0;
        weights[1] = 1.0;
        
        return;
    }
    
    
}

int brachy_specificSampler::sharedConfig(const brachy_specificSampler& o){
  //Share psf configuration
  return psf.sharedConfig(o.psf);
}

REGISTER_SPECIFIC_SAMPLER(brachy_specificSampler,pen_particleState,BRACHY)

#endif
