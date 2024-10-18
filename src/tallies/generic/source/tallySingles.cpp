
//
//
//    Copyright (C) 2024 Universitat de València - UV
//    Copyright (C) 2024 Universitat Politècnica de València - UPV
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


#include "tallySingles.hh"

//Reader functions
int tallyReader_Singles::storeElement(const std::string& pathInSection,
				      const pen_parserData& element,
				      const unsigned){

  if(pathInSection.compare("energy/min") == 0){
    emin = element;
  }
  else if(pathInSection.compare("energy/max") == 0){
    emax = element;
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

void pen_Singles::flush(const unsigned ipart, const unsigned det){

  //Get the buffer
  singlesBuffer& buffer = getBuffer(ipart, det);
  
  //Flush the buffer
  std::string data = buffer.flush();

  //Avoid empty buffers
  if(data.empty())
    return;

  //Get the data hash
  size_t hash = std::hash<std::string>{}(data);

  //Check if the information file has already been created
  if(fInfoPaths[ipart].empty()){
    //Create base name
    fInfoPaths[ipart] =
      "tmp_part_" + std::to_string(ipart) + "_info_" + std::to_string(hash) + ".dat";
    //Save final filename
    fInfoPaths[ipart] = createFilename(fInfoPaths[ipart].c_str());
      
    //Create the file (avoiding adding prefixes, as we have already the final name)
    FILE* aux = ::fopen(fInfoPaths[ipart].c_str(), "w");
    fclose(aux);
  }

  //Create the new data filename
  std::string dataFileName =
    "tmp_det_" + std::to_string(det) +"_part_" + std::to_string(ipart) + "_flush_" + std::to_string(buffer.flushes()) +
    "_data_" + std::to_string(hash) + ".dat";
  dataFileName = createFilename(dataFileName.c_str());

  //Store data
  FILE* fdata = ::fopen(dataFileName.c_str(), "w");
  fprintf(fdata, "%s", data.c_str());
  fclose(fdata);
    
  //Store data filename
  FILE* finfo = ::fopen(fInfoPaths[ipart].c_str(), "a");
  fprintf(finfo,"%u %s\n", det, dataFileName.c_str());
  fclose(finfo);  
}

void pen_Singles::flush(){
  
  for(unsigned i = 0; i < tpart; ++i){
    for(unsigned j = 0; j < kdets.size(); ++j)
      if(kdets[j]){
	flush(i, j);
      }
  }
}

 
void pen_Singles::tally_localEdep(const unsigned long long /*nhist*/,
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
    int it = static_cast<size_t>(state.PAGE / dt);
    if(it < 0 || it >= static_cast<int>(tpart))
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
    
    count(dE-secDE, state.X, state.Y, state.Z, state.PAGE, state.WGHT,
	  static_cast<unsigned>(it));
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

void pen_Singles::tally_step(const unsigned long long /*nhist*/,
			     const pen_KPAR /*kpar*/,
			     const pen_particleState& state,
			     const tally_StepData& stepData){

  //Register singles due continuous energy deposition

  //Check detector
  if(!toDetect)
    return;
  
  if(stepData.softDE == 0.0){return;}  //Nothing to deposit.

  //Check time
  int it = static_cast<size_t>(state.PAGE / dt);
  if(it < 0 || it >= static_cast<int>(tpart))
    return;
  
  //Energy deposited, register single
  count(stepData.softDE, stepData.softX, stepData.softY, stepData.softZ,
	state.PAGE, state.WGHT, static_cast<unsigned>(it));
}


void pen_Singles::tally_knock(const unsigned long long /*nhist*/,
			      const pen_KPAR /*kpar*/,
			      const pen_particleState& /*state*/,
			      const int /*icol*/){
  //The particle suffered an interaction, flag it
  knocked = true;
}

void pen_Singles::tally_endPart(const unsigned long long /*nhist*/,
				const pen_KPAR /*kpar*/,
				const pen_particleState& state){

  //Register singles due absorption

  //Check detector
  if(!toDetect)
    return;

  if(state.MAT == 0){return;} //The particle escpaed from geometry system
  
  if(state.E == 0.0){return;}  //Nothing to deposit.

  //Check time
  int it = static_cast<size_t>(state.PAGE / dt);
  if(it < 0 || it >= static_cast<int>(tpart))
    return;

  //Register the remaining kinetic energy from the particle
  count(state.E, state.X, state.Y, state.Z, state.PAGE, state.WGHT,
	static_cast<unsigned>(it));
  
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

  joinTime = reader.tjoin;

  removeOnEnd = reader.removeOnEnd;

  //Save energy limits
  emin = reader.emin;
  emax = reader.emax;

  //Save time limits
  tmin = reader.tmin;
  tmax = reader.tmax;
  dt = (tmax-tmin)/static_cast<double>(tpart);

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

  //Create also a buffer for each detector and partition
  for(size_t i = 0; i < tpart; ++i){
    buffers[i].resize(idet);
  }
  
  //Register file paths to dump
  for(std::string& p : fInfoPaths)
    dump.toDumpFile(p);

  if(verbose > 1){
    printf("Detection time interval (s): %15.5E - %15.5E\n"
	   "Singles adition window  (s): %15.5E\n"
	   "Sensible detectors:\n"
	   ,tmin, tmax, joinTime);
    for(unsigned long i = 0; i < kdets.size(); ++i){
      if(kdets[i])
	printf("    %lu\n", i);
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
    std::string filename = "singles_det_" + std::to_string(idet) + "_.dat";
    singlesDetFiles[idet] = fopen(filename.c_str(), "w");
  }
  
  // Iterate over all partitions
  for(unsigned ipart = 0; ipart < tpart; ++ipart){
    //Open information file
    if(fInfoPaths[ipart].empty()){
      continue;
    }

    //Open the information file for this partition.
    //Avoid adding prefixes, as we have already the final name
    FILE* fInfo = ::fopen(fInfoPaths[ipart].c_str(), "r");

    //Read each data file
    std::vector<std::pair<unsigned, std::string>> dataFilenames;

    char line[1000];
    while(fgets(line, 1000, fInfo) != nullptr){
      int auxDet;
      char auxFilename[1000];
      if(sscanf(line, " %d %s", &auxDet, auxFilename) != 2){
	printf("Error: Corrupted Singles information file line:"
	       "       File: %s\n"
	       "       Line: %s",
	       fInfoPaths[ipart].c_str(),
	       line);
	continue;
      }

      if(auxDet <= 0){
	printf("Error: Invalid Singles detector number in information file."
	       " Must be grater than 0:"
	       "       File: %s\n"
	       "       Line: %s"
	       "   Detector: %d\n",
	       fInfoPaths[ipart].c_str(),
	       line, auxDet);
	continue;      
      }

      dataFilenames.emplace_back(static_cast<unsigned>(auxDet), auxFilename);
    }

    //Close information file
    fclose(fInfo);

    //Sort data files by detector
    std::sort(dataFilenames.begin(), dataFilenames.end());

  
    //Iterate over detectors
    unsigned nextDataFile = 0;
    for(size_t idet = 0; idet < kdets.size(); ++idet){

      if(!kdets[idet]){
	//Disabled detector, skip
	continue;
      }

      //Open all detector data files
      std::vector<FILE*> dataFiles;
      for(size_t ifile = nextDataFile; ifile < dataFilenames.size(); ++ifile){
	if(dataFilenames[nextDataFile].first == idet){
	  dataFiles.push_back(::fopen(dataFilenames[ifile].second.c_str(), "r"));
	}else{
	  nextDataFile = ifile;
	  break;
	}
      }

      //Read the first single in each file
      std::vector<single> toProcess(dataFiles.size());
      size_t ifile = 0;
      while(ifile < dataFiles.size()){
	if(toProcess[ifile].read(dataFiles[ifile])){
	  ++ifile;
	}else{
	  fclose(dataFiles[ifile]);
	  dataFiles.erase(dataFiles.begin() + ifile);
	  toProcess.erase(toProcess.begin() + ifile);
	}
      }

      //Process singles of all data files
      while(dataFiles.size() > 0){

	//Get the first single among all remaining files
	auto firstS = std::min_element(toProcess.cbegin(), toProcess.cend());
	single s = *firstS;
	size_t firstSFile = std::distance(toProcess.cbegin(), firstS);

	//Load another single from first single file
	if(!toProcess[firstSFile].read(dataFiles[firstSFile])){
	  //Empty file, close it and remove this position
	  fclose(dataFiles[firstSFile]);
	  dataFiles.erase(dataFiles.begin() + firstSFile);
	  toProcess.erase(toProcess.begin() + firstSFile);
	}

	//Go file by file adding singles close enough
	ifile = 0;
	while(ifile < dataFiles.size()){
	  //Get next data file to process
	  FILE* dataFile = dataFiles[ifile];
	  //Get next non processed single in this file
	  single& nextSingle = toProcess[ifile];

	  //Add the singles in this file which are in join time range
	  while(nextSingle.t - s.t <= joinTime){

	    //Add single
	    s.add(nextSingle);

	    //Read the next single in the file
	    if(!nextSingle.read(dataFile)){
	      //End of file reached, remove data file
	      fclose(dataFile);
	      dataFiles.erase(dataFiles.begin() + ifile);
	      toProcess.erase(toProcess.begin() + ifile);
	      --ifile; //Compensate file index increment
	      break;
	    }
	  }
	  //Increase data file index
	  ++ifile;
	}

	//All singles in join time added.
	//Save the final single in the corresponding detector data
	char auxBuff[single::maxBuffSize];
	s.toBufferFinal(auxBuff, single::maxBuffSize);
	fprintf(singlesDetFiles[idet], "%s", auxBuff);
      }

      
    }
  }

  //Close detectors files
  for(size_t idet = 0; idet < kdets.size(); ++idet){

    if(!kdets[idet]){
      //Disabled detector, skip
      continue;
    }

    fclose(singlesDetFiles[idet]);
  }

  //Check if auxiliary files must be removed
  if(removeOnEnd){
    // Iterate over all partitions
    for(unsigned ipart = 0; ipart < tpart; ++ipart){
      //Open information file
      if(fInfoPaths[ipart].empty()){
	continue;
      }

      //Open the information file for this partition.
      //Avoid adding prefixes, as we have already the final name
      FILE* fInfo = ::fopen(fInfoPaths[ipart].c_str(), "r");

      //Read and remove each data file
      char line[1000];
      while(fgets(line, 1000, fInfo) != nullptr){
	int auxDet;
	char auxFilename[1000];
	if(sscanf(line, " %d %s", &auxDet, auxFilename) != 2){
	  printf("Error: Corrupted Singles information file line:"
		 "       File: %s\n"
		 "       Line: %s",
		 fInfoPaths[ipart].c_str(),
		 line);
	  continue;
	}
	std::remove(auxFilename);
      }

      //Close information file
      fclose(fInfo);
      //Remove information file
      std::remove(fInfoPaths[ipart].c_str());
    }
    
  }
  
}

int pen_Singles::sumTally(const pen_Singles& tally){

  //Join information files

  for(unsigned ipart = 0; ipart < tpart; ++ipart){

    if(tally.fInfoPaths[ipart].empty()){
      //Argument tally has not registered any data for this partition
      continue;
    }
    
    
    //Check if the information file has already been created
    if(fInfoPaths[ipart].empty()){

      //Generate a hash from info path
      size_t hash = std::hash<std::string>{}(tally.fInfoPaths[ipart]);
      
      //Create base name
      fInfoPaths[ipart] =
	"tmp_part_" + std::to_string(ipart) + "_info_" + std::to_string(hash) + ".dat";
      //Save final filename
      fInfoPaths[ipart] = createFilename(fInfoPaths[ipart].c_str());
      
      //Create the file (avoiding adding prefixes, as we have already the final name)
      FILE* aux = ::fopen(fInfoPaths[ipart].c_str(), "w");
      fclose(aux);
    }

    //Add the data files from the other tally
    FILE* finfo = ::fopen(fInfoPaths[ipart].c_str(), "a");    
    FILE* finfo2 = ::fopen(tally.fInfoPaths[ipart].c_str(), "r");

    char line[1000];
    while(fgets(line, 1000, finfo2) != nullptr){
      fprintf(finfo, "%s", line);
    }
    fclose(finfo);      
    fclose(finfo2);

    //Remove the old file
    std::remove(tally.fInfoPaths[ipart].c_str());
  }
  
  return 0;
}

REGISTER_COMMON_TALLY(pen_Singles, SINGLES)
