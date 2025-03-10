
//
//
//    Copyright (C) 2024-2025 Universitat de València - UV
//    Copyright (C) 2024-2025 Universitat Politècnica de València - UPV
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
//        sanolgi@upvnet.upv.es (Sandra Oliver Gil)
//    
//


#include "tallySingles.hh"

//Reader functions
int tallyReader_Singles::storeElement(const std::string& pathInSection,
				      const pen_parserData& element,
				      const unsigned){

  if(pathInSection.compare("energy/detection/min") == 0){
    emin = element;
  }
  else if(pathInSection.compare("energy/detection/max") == 0){
    emax = element;
  }
  else if(pathInSection.compare("energy/single/min") == 0){
    singleEmin = element;
  }
  else if(pathInSection.compare("energy/single/max") == 0){
    singleEmax = element;
  }
  else if(pathInSection.compare("time/min") == 0){
    tmin = element;
  }
  else if(pathInSection.compare("time/max") == 0){
    tmax = element;
  }
  else if(pathInSection.compare("time/join") == 0){
    tjoin = element;
  }
  else if(pathInSection.compare("clear") == 0){
    removeOnEnd = element;
  }
  else if(pathInSection.compare("binary") == 0){
    binary = element;
  }
  else{
    return errors::UNHANDLED;
  }
  return errors::SUCCESS;
  
}

int tallyReader_Singles::beginArray(const std::string& pathInSection,
				    const size_t,
				    const unsigned) {
  if(pathInSection.compare("detectors") == 0){
    return errors::SUCCESS;
  }else{
    return errors::UNHANDLED;
  }
}

int tallyReader_Singles::storeArrayElement(const std::string& pathInSection,
					   const pen_parserData& element,
					   const size_t,
					   const unsigned) {

  if(pathInSection.compare("detectors") == 0){
    const unsigned aux = element;
    kdets.emplace_back(aux);
    return errors::SUCCESS;
  }
  
  return errors::UNHANDLED;
}

void pen_Singles::flush(const unsigned actualKDet){

  const unsigned det = detInternalIndex[actualKDet];
  
  //Get the buffer
  singlesBuffer& buffer = buffers[det];
  
  //Flush the buffer
  std::vector<unsigned char> data = buffer.flush();

  //Avoid empty buffers
  if(data.empty())
    return;

  //Store data
  FILE* fdata = ::fopen(fDataPaths[det].c_str(), "ab");
  fwrite(data.data(), sizeof(unsigned char), data.size(), fdata);
  fclose(fdata);
    
  //Store data entry information
  const unsigned long nBytes = static_cast<unsigned long>(sizeof(unsigned char)*data.size());
  FILE* finfo = ::fopen(fInfoPaths[det].c_str(), "a");
  fprintf(finfo,"%lu %lu\n", offsets[det], nBytes);
  offsets[det] += nBytes;
  fclose(finfo);
}

void pen_Singles::flush(){
  
  for(unsigned i = 0; i < kdets.size(); ++i){
    if(kdets[i]){
      flush(i);
    }
  }
}

 
void pen_Singles::tally_localEdep(const unsigned long long nhist,
				  const pen_KPAR /*kpar*/,
				  const pen_particleState& state,
				  const double dE){
  
  //Energy deposited. Register a single only if is caused by an interaction
  if(knocked){

    knocked = false;
    
    //Check detector
    if(!toDetect)
      return;

    if(dE == 0.0)
      return; //Nothing to register
    
    //Check time
    if(state.PAGE < tmin || state.PAGE > tmax)
      return;
    

    //It is due an interaction in a sensible detector.
    //Remove energy of secondary created particles
    
    double secDE = 0.0;
    for(unsigned i = 0; i < constants::nParTypes; ++i){
      const abc_particleStack* stack = readStack(static_cast<pen_KPAR>(i));
      for(unsigned int j = nInStack[i]; j < stack->getNSec(); ++j){
	const pen_particleState secState = stack->readBaseState(j);
	secDE += secState.E*secState.W/state.W;
      }
    }
    
    count(dE-secDE, state.X, state.Y, state.Z, state.PAGE, state.WGHT, nhist);
  }
}

void pen_Singles::tally_beginPart(const unsigned long long /*nhist*/,
				  const unsigned kdet,
				  const pen_KPAR /*kpar*/,
				  const pen_particleState& /*state*/){

  //A particle starts. Save kdet and reset knock flag
  actualKdet = kdet;
  toDetect = activeDet();
  knocked = false;
}

void pen_Singles::tally_jump(const unsigned long long /*nhist*/,
			     const pen_KPAR /*kpar*/,
			     const pen_particleState& /*state*/,
			     const double /*ds*/){

  //Check detector
  if(!toDetect)
    return;
  
  //The particle is in an active detector.
  //Update the number of particles in each stack
  for(unsigned i = 0; i < constants::nParTypes; ++i){
    const abc_particleStack* stack = readStack(static_cast<pen_KPAR>(i));
    nInStack[i] = stack->getNSec();
  }
}

void pen_Singles::tally_step(const unsigned long long nhist,
			     const pen_KPAR /*kpar*/,
			     const pen_particleState& state,
			     const tally_StepData& stepData){

  //Register singles due continuous energy deposition

  //Check detector
  if(!toDetect)
    return;
  
  if(stepData.softDE == 0.0){return;}  //Nothing to deposit.

  //Check time
  if(state.PAGE < tmin || state.PAGE > tmax)
    return;
  
  //Energy deposited, register single
  count(stepData.softDE, stepData.softX, stepData.softY, stepData.softZ,
	state.PAGE, state.WGHT, nhist);
}


void pen_Singles::tally_knock(const unsigned long long /*nhist*/,
			      const pen_KPAR /*kpar*/,
			      const pen_particleState& /*state*/,
			      const int /*icol*/){
  //The particle suffered an interaction, flag it
  knocked = true;
}

void pen_Singles::tally_endPart(const unsigned long long nhist,
				const pen_KPAR /*kpar*/,
				const pen_particleState& state){

  //Register singles due absorption

  //Check detector
  if(!toDetect)
    return;

  if(state.MAT == 0){return;} //The particle escpaed from geometry system
  
  if(state.E == 0.0){return;}  //Nothing to deposit.

  //Check time
  if(state.PAGE < tmin || state.PAGE > tmax)
    return;

  //Register the remaining kinetic energy from the particle
  count(state.E, state.X, state.Y, state.Z, state.PAGE, state.WGHT, nhist);
  
}

void pen_Singles::tally_interfCross(const unsigned long long /*nhist*/,
				    const unsigned kdet,
				    const pen_KPAR /*kpar*/,
				    const pen_particleState& /*state*/){
  //An interface is crossed. Update kdet
  actualKdet = kdet;
  toDetect = activeDet();
  knocked = false;  
}

int pen_Singles::configure(const wrapper_geometry& geometry,
			   const abc_material* const /*materials*/[pen_geoconst::NB],
			   const pen_parserSection& config,
			   const unsigned verbose){

  geo = &geometry;

  //Flag simulation as unfinished
  simFinished = false;
  
  //Read material information from config section
  tallyReader_Singles reader;
  int err = reader.read(config,verbose);
  if(err != tallyReader_Singles::SUCCESS){
    return err;
  }

  binary = reader.binary;

  joinTime = reader.tjoin;

  removeOnEnd = reader.removeOnEnd;

  //Save energy limits
  emin = reader.emin;
  emax = reader.emax;

  singleEmin = reader.singleEmin;
  singleEmax = reader.singleEmax;

  //Save time limits
  tmin = reader.tmin;
  tmax = reader.tmax;

  //Set enabled detectors
  const unsigned maxkdet = *std::max_element(reader.kdets.cbegin(), reader.kdets.cend());
  kdets.clear();
  kdets.resize(maxkdet+1, false);

  detInternalIndex.clear();
  detInternalIndex.resize(maxkdet+1, maxkdet+1);

  //Flag used detectors and assign an internal index to them
  unsigned idet = 0;
  for(unsigned e : reader.kdets){
    kdets[e] = true;
    detInternalIndex[e] = idet++;
  }

  //Create a buffer for each detector
  buffers.resize(idet);

  //Init offsets
  offsets.resize(idet, 0);

  //Create information paths
  for(size_t i = 0; i < idet; ++i){
    //Create base filename
    std::string filename =
      "tmp_det_" + std::to_string(i) + "_info.dat";
    //Create final filename, appending thread and MPI rank information if needed
    filename = createFilename(filename.c_str());
      
    //Create the file (avoiding adding prefixes, as we already have the final name)
    FILE* aux = ::fopen(filename.c_str(), "w");
    fclose(aux);

    //Store information filename
    fInfoPaths.push_back(filename);
  }

  //Create data paths
  for(size_t i = 0; i < idet; ++i){
    //Create base filename
    std::string filename =
      "tmp_det_" + std::to_string(i) + "_data.dat";
    //Create final filename, appending thread and MPI rank information if needed
    filename = createFilename(filename.c_str());
      
    //Create the file (avoiding adding prefixes, as we already have the final name)
    FILE* aux = ::fopen(filename.c_str(), "w");
    fclose(aux);

    //Store information filename
    fDataPaths.push_back(filename);
  }
  
  //Register file paths to dump
  for(std::string& p : fInfoPaths)
    dump.toDumpFile(p);

  //Register file paths to dump
  for(std::string& p : fDataPaths)
    dump.toDumpFile(p);
  
  if(verbose > 1){
    printf("Detection time interval (s): %15.5E - %15.5E\n"
	   "Singles adition window  (s): %15.5E\n"
	   "Sensible detectors:\n"
	   ,tmin, tmax, joinTime);
    for(unsigned long i = 0; i < kdets.size(); ++i){
      if(kdets[i])
	printf("    %lu\n",i);
    }
  }
  
  return 0;
}

void pen_Singles::saveData(const unsigned long long /*nhist*/) const{

  //Process and save final data only on simulation finish
  if(!simFinished)
    return;
  
  //Create final data files for each detector
  std::vector<FILE*> singlesDetFiles(kdets.size());
  for(size_t idet = 0; idet < kdets.size(); ++idet){

    if(!kdets[idet]){
      //Disabled detector, skip
      continue;
    }

    //Create a results data file for this detector
    std::string filename = "singles_det_" + std::to_string(idet) + ".dat";
    singlesDetFiles[idet] = fopen(filename.c_str(), "wb");
  }
  
  // Iterate over all detectors
#ifdef _PEN_USE_THREADS_

  unsigned int nProcessThreads = 
    std::max(static_cast<unsigned int>(2),
	     std::thread::hardware_concurrency());
  std::vector<std::thread> processingThreads;

  std::atomic<unsigned int> atomicCount{0};

  for(unsigned ith = 0; ith < nProcessThreads; ++ith){
    processingThreads.push_back(std::thread([&, ith](){

      unsigned idet = atomicCount++;
      while(idet < this->kdets.size()){
	orderDetectorData(idet, singlesDetFiles[idet]);
	idet = atomicCount++;
      }
      
    }));
  }

  //Wait until all threads have been finished
  for(std::thread& t : processingThreads){
    t.join();
  }

#else
  
  for(unsigned idet = 0; idet < kdets.size(); ++idet){
    orderDetectorData(idet, singlesDetFiles[idet]);
  }
  
#endif

  //Close final results files
  for(size_t idet = 0; idet < kdets.size(); ++idet){

    if(!kdets[idet]){
      //Disabled detector, skip
      continue;
    }

    fclose(singlesDetFiles[idet]);
  }

  //Check if information files must be removed
  if(removeOnEnd){
    // Iterate over all partitions
    for(unsigned det = 0; det < fInfoPaths.size(); ++det){
      //Remove information file
      std::remove(fInfoPaths[det].c_str());
    }
  }
  
}

int pen_Singles::sumTally(const pen_Singles& tally){

  //Check detector number
  if(tally.fInfoPaths.size() != fInfoPaths.size())
    return -1;
  
  //Join information files
  for(unsigned det = 0; det < fInfoPaths.size(); ++det){

    //Open information files from both tallies
    FILE* finfo = ::fopen(fInfoPaths[det].c_str(), "a");    
    if(finfo == nullptr){
      printf("Error in 'sumTally' information file missing: '%s'\n",
	     fInfoPaths[det].c_str());
      return -2;
    }
    FILE* finfo2 = ::fopen(tally.fInfoPaths[det].c_str(), "r");
    if(finfo2 == nullptr){
      printf("Error in 'sumTally' information file missing: '%s'\n",
	     tally.fInfoPaths[det].c_str());
      return -2;
    }
    
    char line[1000];
    while(fgets(line, 1000, finfo2) != nullptr){
      unsigned long offset;
      unsigned long size;
      if(sscanf(line, " %lu %lu", &offset, &size) != 2){
	printf("Error: Corrupted Singles information file line:"
	       "       File: %s\n"
	       "       Line: %s",
	       tally.fInfoPaths[det].c_str(),
	       line);
	break;
      }

      //Append this chunk information
      fprintf(finfo,"%lu %lu\n", offsets[det], size);
      offsets[det] += size;
    }

    //Close both information files
    fclose(finfo);
    fclose(finfo2);

    //Open data files from both tallies
    FILE* fdata = ::fopen(fDataPaths[det].c_str(), "ab");    
    if(fdata == nullptr){
      printf("Error in 'sumTally' data file missing: '%s'\n",
	     fDataPaths[det].c_str());
      return -3;
    }
    FILE* fdata2 = ::fopen(tally.fDataPaths[det].c_str(), "rb");
    if(fdata2 == nullptr){
      printf("Error in 'sumTally' data file missing: '%s'\n",
	     tally.fDataPaths[det].c_str());
      return -2;
    }

    constexpr const size_t buffSize = 10*1024*1024; // 10MB
    std::vector<unsigned char> buffer(buffSize);
    size_t bytesRead;
    while((bytesRead = fread(buffer.data(), 1, buffSize, fdata2)) > 0){
      fwrite(buffer.data(), 1, bytesRead, fdata);
    }
    
    fclose(fdata);      
    fclose(fdata2);      

    //Remove old files
    if(removeOnEnd){
      std::remove(tally.fInfoPaths[det].c_str());
      std::remove(tally.fDataPaths[det].c_str());
    }
  }
  
  return 0;
}


void pen_Singles::orderDetectorData(const unsigned idet, FILE* fout) const {

  if(!kdets[idet]){
    //Disabled detector, skip
    return;
  }
    
  //Get local buffer detector index
  const unsigned det = detInternalIndex[idet];
    
  //Open the information file for this detector.
  //Avoid adding prefixes, as we have already the final name
  FILE* fInfo = ::fopen(fInfoPaths[det].c_str(), "r");
  if(fInfo == nullptr){
    printf("Error: Unable to read singles information"
	   " file for detector %u\n",det);
    return;
  }

  //Read the information of each data block
  std::vector<std::pair<unsigned long, unsigned long>> chunksInfo;

  char line[1000];
  while(fgets(line, 1000, fInfo) != nullptr){
    unsigned long offset;
    unsigned long size;
    if(sscanf(line, " %lu %lu", &offset, &size) != 2){
      printf("Error: Corrupted Singles information file line:"
	     "       File: %s\n"
	     "       Line: %s",
	     fInfoPaths[det].c_str(),
	     line);
      break;
    }

    if(size == 0) //Avoid empty chunks
      continue;

    if(chunksInfo.size() > 0){
      if(offset <= chunksInfo.back().first ||
	 chunksInfo.back().second != offset){
	printf("Error: Inconsistent offset and size values in file '%s'.\n"
	       "       Line: %s\n"
	       "       Offset: %lu\n"
	       "       Previous offset: %lu\n"
	       "       Previous size: %lu\n",
	       fInfoPaths[det].c_str(),
	       line, offset,
	       chunksInfo.back().first,
	       chunksInfo.back().second-chunksInfo.back().first);
	break;
      }
    }

    chunksInfo.emplace_back(offset, size + offset);
  }

  //Close information file
  fclose(fInfo);

  if(chunksInfo.size() == 0) //Empty data file
    return;
    
  //Open detector data file
  FILE* dataFile = ::fopen(fDataPaths[det].c_str(), "rb");
  if(dataFile == nullptr){
    printf("Error: Unable to open data file '%s'\n",
	   fDataPaths[det]);
    return;
  }

  //Read the first single in each data chunk
  std::vector<single> toProcess(chunksInfo.size());
  int iChunk = 0;
  while(iChunk < (int)toProcess.size()){
    if(toProcess[iChunk].read(dataFile, chunksInfo[iChunk].first)){

      //Check if the chunk limit has been reached
      if(chunksInfo[iChunk].first >= chunksInfo[iChunk].second){
	//Remove this chunk data
	chunksInfo.erase(chunksInfo.begin() + iChunk);
	toProcess.erase(toProcess.begin() + iChunk);
      }else{
	//Increase chunk index
	++iChunk;
      }
    }else{
      printf("Error: Unexpected end of file '%s' at process start\n",
	     fDataPaths[det].c_str());
      fflush(stdout);
      //Remove this chunk data
      chunksInfo.erase(chunksInfo.begin() + iChunk);
      toProcess.erase(toProcess.begin() + iChunk);
    }
  }

  //Process all singles
  while(chunksInfo.size() > 0){

    //Get the first single among all remaining chunks
    auto firstS = std::min_element(toProcess.cbegin(), toProcess.cend());
    single s = *firstS;
    size_t firstSChunk = std::distance(toProcess.cbegin(), firstS);

    //Load another single from first single file
    if(!toProcess[firstSChunk].read(dataFile, chunksInfo[firstSChunk].first)){
      printf("Error: Unexpected end of file '%s' replacing opening single\n",
	     fDataPaths[det]);
      fflush(stdout);
      //Remove this chunk data
      chunksInfo.erase(chunksInfo.begin() + firstSChunk);
      toProcess.erase(toProcess.begin() + firstSChunk);
    }else if(chunksInfo[firstSChunk].first >= chunksInfo[firstSChunk].second){
      //The chunk limit has been reached
      //Remove this chunk data
      chunksInfo.erase(chunksInfo.begin() + firstSChunk);
      toProcess.erase(toProcess.begin() + firstSChunk);
    }

    //Go chunk by chunk adding singles close enough
    iChunk = 0;
    while(iChunk < (int)chunksInfo.size()){
      single& nextSingle = toProcess[iChunk];
      //Add the singles in this file which are in join time range
      while(nextSingle.t - s.t <= joinTime){

	//Add single
	s.add(nextSingle);

	//Read the next single in the file
	if(!nextSingle.read(dataFile, chunksInfo[iChunk].first)){
	  printf("Error: Unexpected end of file '%s' searching singles in window\n",
		 fDataPaths[det]);
	  fflush(stdout);
	  //Remove this chunk data
	  chunksInfo.erase(chunksInfo.begin() + iChunk);
	  toProcess.erase(toProcess.begin() + iChunk);	  
	  --iChunk; //Compensate index increment
	  break;
	} else if(chunksInfo[iChunk].first >= chunksInfo[iChunk].second){
	  //The chunk limit has been reached
	  //Remove this chunk data
	  chunksInfo.erase(chunksInfo.begin() + iChunk);
	  toProcess.erase(toProcess.begin() + iChunk);

	  --iChunk; //Compensate index increment
	  break;
	}
      }
      //Increase data file index
      ++iChunk;
    }

    //All singles in join time added.
    //Save the final single in the corresponding detector data, if it is in energy range
    if(s.E >= singleEmin && s.E <= singleEmax){
      if(binary){
	unsigned char auxBuff[single::dataSize];
	size_t pos = 0;
	s.toBufferFinalB(auxBuff, pos);
	fwrite(auxBuff, sizeof(unsigned char), single::dataSize, fout);
      }
      else{
	char auxBuff[single::maxBuffSize];
	int nwrite = s.toBufferFinal(auxBuff, single::maxBuffSize);
	fwrite(auxBuff, sizeof(char), nwrite, fout);
      }
    }
  }

  //Close detector data file
  fclose(dataFile);      

  //Check if the data file must be removed
  if(removeOnEnd){
    std::remove(fDataPaths[det].c_str());
  }    
    
}

REGISTER_COMMON_TALLY(pen_Singles, SINGLES)
