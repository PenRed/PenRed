
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

  if(pathInSection.compare("energy/single/min") == 0){
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
  else if(pathInSection.compare("pileup") == 0){
    pileup = element;
  }
  else if(pathInSection.compare("scatter") == 0){
    scatter = element;
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

void pen_Singles::flush(const unsigned det){
  
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
  const unsigned long nBytes = static_cast<unsigned long>(data.size());
  FILE* finfo = ::fopen(fInfoPaths[det].c_str(), "a");
  fprintf(finfo,"%lu %lu\n", offsets[det], nBytes);
  offsets[det] += nBytes;
  fclose(finfo);
}

void pen_Singles::flush(){
  
  for(unsigned i = 0; i < buffers.size(); ++i){

    //Flush data
    flush(i);

    //Reset last history start position
    lastHistStart[i] = 0;
    
  }
}

 
void pen_Singles::tally_localEdep(const unsigned long long nhist,
				  const pen_KPAR kpar,
				  const pen_particleState& state,
				  const double dE){
  //Energy deposited. Register it
  count(dE, state.X, state.Y, state.Z, state.PAGE, state.WGHT, kpar, nhist);
}

void pen_Singles::tally_beginPart(const unsigned long long nhist,
				  const unsigned kdet,
				  const pen_KPAR kpar,
				  const pen_particleState& state){

  //A particle starts. Save kdet and register energy "extraction" from material
  actualKdet = kdet;
  toDetect = activeDet();

  //Set last icol to "absorbed after produced"
  lastICol = 255;

  //Reset boolean mask
  const int gen = state.ILB[0];
  actualMask = (gen == 1 ? single::FIRST_GENERATION : (gen == 2 ? single::SECOND_GENERATION : 0));
  
  if(state.ILB[1] == PEN_POSITRON && state.ILB[2] == BETAp_ANNIHILATION)
    actualMask |= single::FROM_ANNIHILATION;
  
  if(skipBeginPart){
    skipBeginPart = false;
  }
  else{
    if(toDetect)
      actualMask |= single::PRODUCED_IN_DETECTOR;
    count(-state.E, state.X, state.Y, state.Z, state.PAGE, state.WGHT, kpar, nhist);
  }
}

void pen_Singles::tally_sampledPart(const unsigned long long nhist,
				    const unsigned long long /*dhist*/,
				    const unsigned /*kdet*/,
				    const pen_KPAR /*kpar*/,
				    const pen_particleState& state){
  //Update last hist
  lastHist = nhist;
  
  //Add a pulse to compensate the substracted energy on beginPart call
  if(state.MAT > 0){
    //If the particle has been sampled in a non void region, prevent energy counting by beginPart call
    skipBeginPart = true;
  }
}

void pen_Singles::tally_step(const unsigned long long nhist,
			     const pen_KPAR kpar,
			     const pen_particleState& state,
			     const tally_StepData& stepData){

  //Set last icol to soft deposition
  if(kpar == PEN_ELECTRON)
    lastICol = BETAe_SOFT_INTERACTION;
  else if(kpar == PEN_POSITRON)
    lastICol = BETAp_SOFT_INTERACTION;    

  //Register singles due continuous energy deposition
  count(stepData.softDE, stepData.softX, stepData.softY, stepData.softZ,
	state.PAGE, state.WGHT, kpar, nhist);

  //Set last icol to "absorbed after moved" to handle
  //this kind of absorptions via tally_localEdep call
  lastICol = 254;
}

void pen_Singles::tally_move2geo(const unsigned long long /*nhist*/,
				 const unsigned /*kdet*/,
				 const pen_KPAR /*kpar*/,
				 const pen_particleState& state,
				 const double /*dsef*/,
				 const double /*dstot*/){
  
  //Particle has been created at void volume. Check if the geomtry system is reached  
  if(state.MAT > 0){
    //Non void volume reached, prevent energy counting by beginPart call
    skipBeginPart = true;
  }
}

void pen_Singles::tally_interfCross(const unsigned long long /*nhist*/,
				    const unsigned kdet,
				    const pen_KPAR /*kpar*/,
				    const pen_particleState& /*state*/){
  //An interface is crossed. Update kdet
  actualKdet = kdet;
  toDetect = activeDet();
}

void pen_Singles::tally_endHist(const unsigned long long /*nhist*/){

  for(unsigned i = 0; i < buffers.size(); ++i){

    //Reduce last history singles
    buffers[i].reduce(lastHistStart[i], scatter, joinTime, tmin, tmax);
    //Update the last history start position after the reduce step
    lastHistStart[i] = buffers[i].size();

    //Flush buffer, if necessary
    if(buffers[i].size() >= singlesBuffer::baseSize){
      flush(i);
      //Restart last history start position
      lastHistStart[i] = 0;
    }
  }
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

  pileup = reader.pileup;
  scatter = reader.scatter;
  
  removeOnEnd = reader.removeOnEnd;

  singleEmin = reader.singleEmin;
  singleEmax = reader.singleEmax;

  //Save time limits
  tmin = reader.tmin;
  tmax = reader.tmax;

  if(verbose > 1){
    printf("Output format: %s\n"
	   "Joining time : %.5E\n"
	   "Pileup       : %s\n"
	   "Scatter      : %s\n"
	   "Energy range : [%.5E, %.5E] eV\n"
	   "Time range   : [%.5E, %.5E] s\n",
	   binary ? "Binary" : "ASCII",
	   joinTime,
	   pileup ? "Enabled" : "Disabled",
	   scatter ? "Enabled" : "Disabled",
	   singleEmin, singleEmax,
	   tmin, tmax);
  }

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
  lastHistStart.resize(idet,0);

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

  skipBeginPart = false;
  
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


    //Append second tally data file to the first one
    constexpr const size_t buffSize = single::dataSize*singlesBuffer::baseSize;
    constexpr const size_t histPos = single::dataSize - sizeof(unsigned long long);
    std::vector<unsigned char> buffer(buffSize);
    size_t bytesRead;
    unsigned long long globalLastHist = lastHist;

    while((bytesRead = fread(buffer.data(), 1, buffSize, fdata2)) > 0){

      //Check read chunk
      if(bytesRead % single::dataSize != 0){
	printf("Error in 'sumTally'. Data size in both files missmatch "
	       "or the file is corrupted: '%s'\n",
	       tally.fDataPaths[det].c_str());	
      }

      //Iterate over read singles to correct history number
      size_t readSingles = bytesRead / single::dataSize;
      for(size_t iSing = 0; iSing < readSingles; ++iSing){
	//Calculate buffer initial position
	const size_t singHistPos = iSing*single::dataSize + histPos;
	//Extract history number
	unsigned long long localHist;
	std::memcpy(&localHist, &buffer[singHistPos], sizeof(unsigned long long));
      
	//Add the last history in the local tally
	localHist += lastHist;

	//Update buffer
	std::memcpy(&buffer[singHistPos], &localHist, sizeof(unsigned long long));

	//Update global last history
	if(globalLastHist < localHist) globalLastHist = localHist;
      }
      
      //Append data to the local data file
      fwrite(buffer.data(), 1, bytesRead, fdata);      
    }
    if(bytesRead != 0){
      printf("Error in 'sumTally'. Data size in both files missmatch "
	     "or the file is corrupted: '%s'\n",
	     tally.fDataPaths[det].c_str());
    }

    //Update last history
    lastHist = globalLastHist;
    
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
	   fDataPaths[det].c_str());
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

    //Load another single from the file "firstS" belongs
    if(!toProcess[firstSChunk].read(dataFile, chunksInfo[firstSChunk].first)){
      printf("Error: Unexpected end of file '%s' replacing opening single\n",
	     fDataPaths[det].c_str());
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

    //Check if pileup is enabled
    if(pileup){
      //Go chunk by chunk adding singles close enough
      iChunk = 0;
      while(iChunk < (int)chunksInfo.size()){
	single& nextSingle = toProcess[iChunk];
	//Add the singles in this file which are in join time range
	while(nextSingle.t - s.t <= joinTime){

	  //Add single
	  s.add(nextSingle);

	  //Flag pileup
	  s.info[0] |= single::PILEUP;
	  
	  //Read the next single in the file
	  if(!nextSingle.read(dataFile, chunksInfo[iChunk].first)){
	    printf("Error: Unexpected end of file '%s' searching singles in window\n",
		   fDataPaths[det].c_str());
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
    }

    //All singles in join time added.
    //Save the final single in the corresponding detector data, if it falls within energy range
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

void pen_Singles::singlesBuffer::reduce(const size_t start, const bool saveScatter,
					const double dt, const double tmin, const double tmax){

  if(n <= start+1)
    return;      

  //Sort specified data range
  std::sort(buffer.begin() + start, buffer.begin() + n);

  // Add singles within the same time window
  
  //Create auxiliary single
  single auxSing;
  size_t iAux = start;
  
  //Get the first single with positive energy
  size_t ifirst = iAux;
  auxSing = buffer[ifirst];
  while(auxSing.E < 0.0 && ifirst < n-1){
    //Get the next one
    auxSing = buffer[++ifirst]; 
  }
      
  if(ifirst != iAux){
    //The first single with positive energy is not at the first position, correct it
    buffer[ifirst] = buffer[iAux];
  }

  if(auxSing.E < 0.0){
    penred::logs::logger::printf(penred::logs::SIMULATION,
				 "Error: Only singles with negative energy found!\n"
				 "   First Energy: %15.5E. Pleae, report this case\n",
				 buffer[ifirst].E);
  }

  for(size_t i = iAux+1; i < n; ++i){
    if(buffer[i].E < 0.0){ //Check if the next pulse has a negative energy 	  
      auxSing.addEnergy(buffer[i].E, buffer[i].weight/buffer[i].E);
    }
    else if(buffer[i].t - auxSing.t < dt){ //Check time window
      auxSing.add(buffer[i]);	  
    }else{
      //Pulse with positive energy outside the join window

      //Save the current single
      if(saveScatter || !auxSing.isScattered()){
	if(auxSing.E > 0.0){
	  if(auxSing.t >= tmin && auxSing.t <= tmax)
	    buffer[iAux++] = auxSing;
	}
      }
	  
      //Get the next one
      auxSing = buffer[i];
    }
	
  }
      
  //Save the last single
  if(auxSing.E > 0.0){
    buffer[iAux++] = auxSing;
  }
      
  //Update the number of elements in the buffer
  n = iAux;      
}

REGISTER_COMMON_TALLY(pen_Singles, SINGLES)
