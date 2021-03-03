//
//
//    Copyright (C) 2020-2021 Universitat Politècnica de València - UPV
//
//    This file is part of RUPER-LB: Runtime Unpredictable Performance Load Balancer.
//
//    RUPER-LB is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    RUPER-LB is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with RUPER-LB.  If not, see <https://www.gnu.org/licenses/>. 
//
//
//    contact emails:
//
//        vicent.gimenez.alventosa@gmail.com (Vicent Giménez Alventosa)
//

#include "loadBalance.hh"

// *** MPI Request functions

// **************************** MPI ********************************* //
#ifdef _PEN_USE_MPI_
int LB::testAny(std::vector<LB::request>& reqs,
		int& index){
  index = -1;
  int cont = 0;
  for(auto& req : reqs){
    int flag;
    int err = req.test(flag);
    if(err != MPI_SUCCESS)
      return err;
    if(flag){
      index = cont;
      return MPI_SUCCESS;
    }
    ++cont;
  }
  return MPI_SUCCESS;
}

int LB::checkAny(std::vector<LB::request>& reqs,
		 int& index){
  index = -1;
  int cont = 0;
  for(auto& req : reqs){
    int flag;
    int err = req.check(flag);
    if(err != MPI_SUCCESS)
      return err;
    if(flag){
      index = cont;
      return MPI_SUCCESS;
    }
    ++cont;
  }
  return MPI_SUCCESS;
}
#endif
  // ************************** MPI END ******************************* //


//
// *** load balance worker class ***
//

float LB::worker::addMeasure(const std::chrono::steady_clock::time_point& actualTime,
			    const unsigned long long nIter){

  if(nIter < iterDone){
    //No sense
    return 1.0;
  }

  if(init > actualTime || actualTime <= lastTime){
    //No sense
    return 1.0;
  }

  //Calculate elapsed time since previous iteration update
  std::chrono::seconds::rep elapsedTime = elapsed(actualTime);

  //Calculate elapsed time since task beginning
  std::chrono::seconds::rep measureStamp = timeStamp(actualTime);
      
  if(elapsedTime <= 0){
    //Insufficient elapsed time
    return 1.0;
  }

  //Calculate the number of new iterations executed
  //since last report
  unsigned long long newIters = nIter-iterDone;
  
  //Update lasttime
  lastTime = actualTime;
  
  //Calculate new speed
  float newSpeed = static_cast<float>(newIters)/static_cast<float>(elapsedTime);

  //Get last speed
  float lastSpeed = speed();
  if(lastSpeed <= 0.0)
    lastSpeed = newSpeed;
  
  //Store measured timestamp and speed
  measures.push_back(LB::measure(measureStamp,newSpeed));
  
  //Update done iterations
  iterDone = nIter;  

  //Return the deviation respect previous measure
  return newSpeed/lastSpeed;
}

int LB::worker::start(const std::chrono::steady_clock::time_point& initIn,
		      const unsigned long long nIter){
  
  //Check if this worker is already doing the task
  if(started){
    return LB_ERROR_WORKER_ALREADY_STARTED;
  }
  //Check if this worker has finished this task
  if(finished){
    return LB_ERROR_WORKER_FINISHED;    
  }

  //Label this worker as started
  started = true;
  finished = false;
  
  //Set init timestamp
  init = initIn;
  //Set last report timestamp to allow adding the first measure
  lastTime = init-std::chrono::seconds(2);

  //Initialize iterations done
  iterDone = 0;
  //Set assigned iterations
  assigned = nIter;
  
  //Register the first dummy measure on worker start
  addMeasure(init,0); //iterations done at worker init = 0

  return LB_SUCCESS;
}

void LB::worker::end() {
  //Calculate elapsed time since previous iteration update
  const auto actualTime = std::chrono::steady_clock::now();
  //Add a final measure with the same speed of the last report
  std::chrono::seconds::rep measureStamp = timeStamp(actualTime);
  measures.push_back(LB::measure(measureStamp,speed()));
  finished = true;
}

unsigned long long LB::worker::predDone(std::chrono::steady_clock::time_point now) const {
  //Return the prediction of iterations done at timestamp "now"

  if(measures.size() < 1) return 0;
	    
  if(lastTime >= now)
    return iterDone;
	    
  //Calculate elapsed time since previous report
  std::chrono::seconds::rep timeInterval = elapsed(now);

  //Calculate non registered time contribution
  double doneSinceUpdate = measures.back().speed*timeInterval;

  return iterDone+static_cast<unsigned long long>(doneSinceUpdate);
}

int LB::worker::printReport(FILE* fout,
			    const std::chrono::steady_clock::time_point globStart) const{
  if(fout == nullptr)
    return LB_ERROR_NULL_POINTER;

  //Calculate elapsed time between global and worker starts
  std::chrono::seconds::rep startDiff = std::chrono::duration_cast<
    std::chrono::seconds>(init-globStart).count();
  
  fprintf(fout,"# timestamp (s) | speed (iter/s)\n");  
  for(const auto& measure: measures){
    long long int timestamp =
      static_cast<long long int>(measure.timestamp+startDiff);
    fprintf(fout," %15lld   %10.4f\n",timestamp,measure.speed);
  }
  return LB_SUCCESS;
}

int LB::worker::dump(FILE* fdump) const{
  if(fdump == nullptr)
    return LB_ERROR_NULL_POINTER;

  //Calculate seconds from init until lastTime
  std::chrono::seconds::rep timeInterval = std::chrono::duration_cast<
    std::chrono::seconds>(lastTime-init).count();

  //Print worker state
  fprintf(fdump,"%llu %d %d %llu %lld %lu\n",
	  assigned,started,finished,iterDone,
	  static_cast<long long int>(timeInterval),
	  static_cast<unsigned long>(measures.size()));  

  //Print worker measures
  for(const auto& measure: measures){
    long long int timestamp =
      static_cast<long long int>(measure.timestamp);
    fprintf(fdump," %lld  %10.8E\n",timestamp,measure.speed);
  }
  return LB_SUCCESS;
}

int LB::worker::load(FILE* fin,
		     const std::chrono::steady_clock::time_point taskBegin){
  if(fin == nullptr)
    return LB_ERROR_NULL_POINTER;

  reset();

  //Save initial time point
  init = taskBegin;

  //Read worker state
  unsigned long nmeasures = 0;
  long long int timeInterval;
  int istarted, ifinished;
  fscanf(fin," %llu %d %d %llu %lld %lu ",
	 &assigned,&istarted,&ifinished,&iterDone,
	 &timeInterval,&nmeasures);
  started = istarted;
  finished = ifinished;
  //Calculate last time
  lastTime = init + std::chrono::seconds(timeInterval);
  
  for(unsigned long i = 0; i < nmeasures; ++i){
    long long int timestamp;
    float speed;
    fscanf(fin," %lld %f ",&timestamp,&speed);
    measures.push_back(LB::measure(static_cast<std::chrono::seconds::rep>
				   (timestamp),speed));
  }
  
  return LB_SUCCESS;
}

//
// *** load balance guess worker class ***
//

float LB::guessWorker::addMeasure(const std::chrono::steady_clock::time_point& actualTime,
			    const unsigned long long nIter){

  //Check if we have some speed to correct
  if(measures.size() < 1 || measures.back().speed == 0.0){
    //We can't correct any speed, so get the measure as normal worker
    return worker::addMeasure(actualTime,nIter);
  }
  
  if(init > actualTime || actualTime <= lastTime){
    //No sense
    return 1.0;
  }

  //Calculate elapsed time since previous iteration update
  std::chrono::seconds::rep elapsedTime = elapsed(actualTime);

  //Calculate elapsed time since task beginning
  std::chrono::seconds::rep measureStamp = timeStamp(actualTime);

  if(elapsedTime <= 0){
    //Insufficient elapsed time
    return 1.0;
  }

  float speed;
  float Iratio;
  //Check if reported iterations are greater than previous report
  if(nIter <= iterDone){
    //Previous prediction has been too optimistic.
    //Correct previous measure with mean velocities
    float v1 = measures.back().speed;
    std::chrono::seconds::rep dt1 = measures.back().timestamp;
    if(dt1 <= std::chrono::seconds::rep(0))
      dt1 = 1;
    double v1m = iterDone/static_cast<double>(dt1);
    double v2m = nIter/static_cast<double>(measureStamp);
    Iratio = v2m/v1m;
    speed = v1*Iratio;
  }
  else{
    // dI1 = dt*v1
    // dI2 = dt*v2
    // dI1/dI2 = v1/v2 -> v2 = v1*(dI2/dI1)    
    float v1 = measures.back().speed;
    double dI1 = static_cast<double>(elapsedTime)*v1;
    double dI2 = static_cast<double>(nIter - iterDone);
    Iratio = dI2/dI1;
    speed = v1*Iratio;
  }
        
  //Update lasttime
  lastTime = actualTime;
  
  //Store measured timestamp and speed
  measures.push_back(LB::measure(measureStamp,speed));
  
  //Update done iterations
  iterDone = nIter;

  return Iratio;
}

//
// *** load balance task class *** 
//

int LB::task::extConnect(FILE* flog,
			 bool blockWrite,
			 const unsigned verbose){

  int err = PEN_TCP_SUCCESS;
  // Connect to the server
  //***********************
  err = client.connect();
  if(err != PEN_TCP_SUCCESS){
    if(verbose > 0){
      if(blockWrite) const std::lock_guard<std::mutex> lockTh(mtxfth);
      fprintf(flog,"%07ld s - LB::task::extConnect: Error: "
	      "Connection fail\n"
	      "      Error code: %d\n",
	      static_cast<long int>(timeStamp()),err);
      fflush(flog);
    }
    //Close connections
    client.close();  
  }
  return err;
}

int LB::task::extSendAndRecv(char message[pen_tcp::messageLen],
			     std::string& response,
			     FILE* flog,
			     bool blockWrite,
			     const unsigned verbose){

  int err = PEN_TCP_SUCCESS;

  // Send the request
  //***********************
  err = client.write(message);
  if(err != PEN_TCP_SUCCESS){
    if(verbose > 0){
      if(blockWrite) const std::lock_guard<std::mutex> lockTh(mtxfth);
      fprintf(flog,"%07ld s - LB::task::extSendAndRecv: Error: Send fail\n",
	      static_cast<long int>(timeStamp()));
      fflush(flog);
    }
    //Close connections
    client.close();  
    return err;    
  }

  // Wait for response
  //*******************
  err = client.receive(response);
  //Close connection
  client.close();
  if(err != PEN_TCP_SUCCESS){
    if(verbose > 0){
      if(blockWrite) const std::lock_guard<std::mutex> lockTh(mtxfth);
      fprintf(flog,"%07ld s - LB::task::extSendAndRecv: Error: Receive fail\n",
	      static_cast<long int>(timeStamp()));
      fflush(flog);
    }
    return err;
  }

  return PEN_TCP_SUCCESS;
}

int LB::task::parseAssigned(const std::string& response,
			    unsigned long long& newAssign,
			    int& errCode){
  //Extract error flag from the response
  int flag;
  errCode = LB_SUCCESS;
  int nread = sscanf(response.c_str(),"%d",&flag);
  if(nread != 1){
    return LB_ERROR_UNEXPECTED_FORMAT;    
  }

  if(flag != 0){
    //Error, get error code
    nread = sscanf(response.c_str(),"%*d %d",&errCode);
    if(nread != 1){
      return LB_ERROR_UNEXPECTED_FORMAT;    
    }
    return LB_REMOTE_ERROR;
  }

  //Success, get assigned iterations
  nread = sscanf(response.c_str(),"%*d %*s %llu",
		 &newAssign);
  if(nread != 1){
    return LB_ERROR_UNEXPECTED_FORMAT;    
  }

  return LB_SUCCESS;
}

int LB::task::extStart(int& TCPerr,
		       int& serverErr,
		       const unsigned verbose,
		       bool trusted){

  //Check if MPI is enabled
  FILE* flog = fth;
  bool blockWrite = true;
#ifdef _PEN_USE_MPI_
  (void)trusted; //Avoid unused warning
  int rank;
  MPI_Comm_rank(comm, &rank);
  if(rank != 0){
    return LB_ERROR_MPI_RANK_0_EXPECTED;
  }
  flog = fMPI;
  blockWrite = false;
#endif
  
  // Connect with the server
  //*************************
  TCPerr = extConnect(flog,blockWrite,verbose);
  if(TCPerr != PEN_TCP_SUCCESS)
    return LB_TCP_CONNECTION_FAIL;
  
  // Send a start petition
  //***********************

  //Get elapsed time since task init
  long int sinceStart = static_cast<long int>(timeStamp());
  
  char message[pen_tcp::messageLen];
  int w = snprintf(message,pen_tcp::messageLen,"2 %ld %ld\n",
		   ext_iworker,sinceStart);
  if(w <= 0 || static_cast<size_t>(w) >= pen_tcp::messageLen){
    if(verbose > 0){
      if(blockWrite) const std::lock_guard<std::mutex> lockTh(mtxfth);
      fprintf(flog,"%07ld s - LB::task::extStart: Error: Unable to "
	      "create request message\n",
	      static_cast<long int>(timeStamp()));
      fflush(flog);
    }
    //Close connection
    client.close();
    return LB_UNABLE_TO_CREATE_REQUEST;
  }
  
  if(verbose > 1){
    if(blockWrite) const std::lock_guard<std::mutex> lockTh(mtxfth);
    fprintf(flog,"%07ld s - LB::task::extStart: Send start request "
	    "to external balancer.\n",
	    static_cast<long int>(timeStamp()));
    fflush(flog);
  }

  std::string response;
  TCPerr = extSendAndRecv(message,response,flog,blockWrite,verbose);
  if(TCPerr != PEN_TCP_SUCCESS){
    return LB_TCP_COMMUNICATION_FAIL;
  }

  // Process the response
  //**********************
  response = pen_tcp::trim(response);  
  if(verbose > 1){
    if(blockWrite) const std::lock_guard<std::mutex> lockTh(mtxfth);
    fprintf(flog,"%07ld s - LB::task::extStart: Received response:\n%s\n",
	    static_cast<long int>(timeStamp()),response.c_str());
    fflush(flog);
  }

  unsigned long long ext_assigned;
  int errParse = parseAssigned(response,ext_assigned,serverErr);
  if(errParse != LB_SUCCESS){
    if(errParse == LB_ERROR_UNEXPECTED_FORMAT){
      if(verbose > 0){
	if(blockWrite) const std::lock_guard<std::mutex> lockTh(mtxfth);
	fprintf(flog,"%07ld s - LB::task::extStart: Error: Unexpected "
		"response format\n",
		static_cast<long int>(timeStamp()));
	fflush(flog);
      }
    }
    return errParse;
  }

  // Update assignation
  //********************
#ifdef _PEN_USE_MPI_
  setIterationsMPI(ext_assigned);
#else
  if(trusted)
    setIterationsTrusted(ext_assigned);
  else
    setIterations(ext_assigned);
#endif

  if(verbose > 1){
    if(blockWrite) const std::lock_guard<std::mutex> lockTh(mtxfth);
    fprintf(flog,"%07ld s - LB::task::extStart: Updated iterations to: %llu\n",
	    static_cast<long int>(timeStamp()),ext_assigned);
    fflush(flog);
  }
  
  return LB_SUCCESS;
}

int LB::task::extReport(int& TCPerr,
			int& serverErr,
			const unsigned long long nIter,
			const std::chrono::steady_clock::time_point& rep_time,
			const unsigned verbose,
			bool trusted){

  //Check if MPI is enabled
  FILE* flog = fth;
  bool blockWrite = true;
#ifdef _PEN_USE_MPI_
  (void)trusted; //Avoid unused warning
  int rank;
  MPI_Comm_rank(comm, &rank);
  if(rank != 0){
    return LB_ERROR_MPI_RANK_0_EXPECTED;
  }
  flog = fMPI;
  blockWrite = false;
#endif

  // Connect with the server
  //*************************
  TCPerr = extConnect(flog,blockWrite,verbose);
  if(TCPerr != PEN_TCP_SUCCESS)
    return LB_TCP_CONNECTION_FAIL;
    
  // Send a report
  //****************

  //Get elapsed time since task init
  long int sinceStart = static_cast<long int>(timeStamp(rep_time));
  
  char message[pen_tcp::messageLen];
  int w = snprintf(message,pen_tcp::messageLen,"1 %ld %llu %ld\n",
		   ext_iworker,nIter,sinceStart);
  if(w <= 0 || static_cast<size_t>(w) >= pen_tcp::messageLen){
    if(verbose > 0){
      if(blockWrite) const std::lock_guard<std::mutex> lockTh(mtxfth);
      fprintf(flog,"%07ld s - LB::task::extReport: Error: Unable to "
	      "create request message\n",
	      static_cast<long int>(timeStamp()));
      fflush(flog);
    }
    return LB_UNABLE_TO_CREATE_REQUEST;
  }
  
  if(verbose > 1){
    if(blockWrite) const std::lock_guard<std::mutex> lockTh(mtxfth);
    fprintf(flog,"%07ld s - LB::task::extReport: Send report "
	    "to external balancer.\n",
	    static_cast<long int>(timeStamp()));
    fflush(flog);
  }

  std::string response;
  TCPerr = extSendAndRecv(message,response,flog,blockWrite,verbose);
  if(TCPerr != PEN_TCP_SUCCESS){
    return LB_TCP_COMMUNICATION_FAIL;
  }

  // Process the response
  //**********************
  response = pen_tcp::trim(response);  
  if(verbose > 1){
    if(blockWrite) const std::lock_guard<std::mutex> lockTh(mtxfth);
    fprintf(flog,"%07ld s - LB::task::extReport: Received response:\n%s\n",
	    static_cast<long int>(timeStamp()),response.c_str());
    fflush(flog);
  }

  unsigned long long ext_assigned;
  int errParse = parseAssigned(response,ext_assigned,serverErr);
  if(errParse != LB_SUCCESS){
    if(errParse == LB_ERROR_UNEXPECTED_FORMAT){
      if(verbose > 0){
	if(blockWrite) const std::lock_guard<std::mutex> lockTh(mtxfth);
	fprintf(flog,"%07ld s - LB::task::extReport: Error: Unexpected "
		"response format\n",
		static_cast<long int>(timeStamp()));
	fflush(flog);
      }
    }
    return errParse;
  }

  // Update assignation
  //********************
#ifdef _PEN_USE_MPI_
  setIterationsMPI(ext_assigned);
#else
  if(trusted)
    setIterationsTrusted(ext_assigned);
  else
    setIterations(ext_assigned);
#endif

  if(verbose > 1){
    if(blockWrite) const std::lock_guard<std::mutex> lockTh(mtxfth);
    fprintf(flog,"%07ld s - LB::task::extReport: Updated iterations to: %llu\n",
	    static_cast<long int>(timeStamp()),ext_assigned);
    fflush(flog);
  }
  
  return LB_SUCCESS;
}

int LB::task::extFinish(int& TCPerr,
			int& serverErr,
			const unsigned long long nIter,
			const std::chrono::steady_clock::time_point& rep_time,
			const unsigned verbose){

  //Check if MPI is enabled
  FILE* flog = fth;
  bool blockWrite = true;
#ifdef _PEN_USE_MPI_
  int rank;
  MPI_Comm_rank(comm, &rank);
  if(rank != 0){
    return LB_ERROR_MPI_RANK_0_EXPECTED;
  }
  flog = fMPI;
  blockWrite = false;
#endif

  // Connect with the server
  //*************************
  TCPerr = extConnect(flog,blockWrite,verbose);
  if(TCPerr != PEN_TCP_SUCCESS)
    return LB_TCP_CONNECTION_FAIL;
  
  // Send a finish request
  //***********************

  //Get elapsed time since task init
  long int sinceStart = static_cast<long int>(timeStamp(rep_time));
  
  char message[pen_tcp::messageLen];
  int w = snprintf(message,pen_tcp::messageLen,"3 %ld %llu %ld\n",
		   ext_iworker,nIter,sinceStart);
  if(w <= 0 || static_cast<size_t>(w) >= pen_tcp::messageLen){
    if(verbose > 0){
      if(blockWrite) const std::lock_guard<std::mutex> lockTh(mtxfth);
      fprintf(flog,"%07ld s - LB::task::extFinish: Error: Unable to "
	      "create request message\n",
	      static_cast<long int>(timeStamp()));
      fflush(flog);
    }
    return LB_UNABLE_TO_CREATE_REQUEST;
  }
  
  if(verbose > 1){
    if(blockWrite) const std::lock_guard<std::mutex> lockTh(mtxfth);
    fprintf(flog,"%07ld s - LB::task::extFinish: Send finish request "
	    "to external balancer.\n",
	    static_cast<long int>(timeStamp()));
    fflush(flog);
  }

  std::string response;
  TCPerr = extSendAndRecv(message,response,flog,blockWrite,verbose);
  if(TCPerr != PEN_TCP_SUCCESS){
    return LB_TCP_COMMUNICATION_FAIL;
  }

  // Process the response
  //**********************
  response = pen_tcp::trim(response);  
  if(verbose > 1){
    if(blockWrite) const std::lock_guard<std::mutex> lockTh(mtxfth);
    fprintf(flog,"%07ld s - LB::task::extFinish: Received response:\n%s\n",
	    static_cast<long int>(timeStamp()),response.c_str());
    fflush(flog);
  }

  //Extract error flag from the response
  int flag;
  int nread = sscanf(response.c_str(),"%d",&flag);
  if(nread != 1){
    if(verbose > 0){
      if(blockWrite) const std::lock_guard<std::mutex> lockTh(mtxfth);
      fprintf(flog,"%07ld s - LB::task::extFinish: Error: Unexpected "
	      "response format\n",
	      static_cast<long int>(timeStamp()));
      fflush(flog);
    }
    return LB_ERROR_UNEXPECTED_FORMAT;    
  }

  if(flag != 0){
    nread = sscanf(response.c_str(),"%*d %d",&serverErr);
    if(nread != 1){
      if(verbose > 0){
	if(blockWrite) const std::lock_guard<std::mutex> lockTh(mtxfth);
	fprintf(flog,"%07ld s - LB::task::extFinish: Error: Unexpected "
		"response format\n",
		static_cast<long int>(timeStamp()));
	fflush(flog);
      }
      return LB_ERROR_UNEXPECTED_FORMAT;    
    }
    return LB_REMOTE_ERROR;
  }
  
  return LB_SUCCESS;
}

int LB::task::extAssign(int& TCPerr,
			int& serverErr,
			const unsigned verbose,
			bool trusted){

  //Check if MPI is enabled
  FILE* flog = fth;
  bool blockWrite = true;
#ifdef _PEN_USE_MPI_
  (void)trusted; //Avoid unused warning
  int rank;
  MPI_Comm_rank(comm, &rank);
  if(rank != 0){
    return LB_ERROR_MPI_RANK_0_EXPECTED;
  }
  flog = fMPI;
  blockWrite = false;
#endif
  
  // Send a assign petition
  //************************
  char message[pen_tcp::messageLen];
  snprintf(message,pen_tcp::messageLen,"5 %ld\n",ext_iworker);
  
  if(verbose > 1){
    if(blockWrite) const std::lock_guard<std::mutex> lockTh(mtxfth);
    fprintf(flog,"%07ld s - LB::task::extAssign: Send assign request "
	    "to external balancer.\n",
	    static_cast<long int>(timeStamp()));
    fflush(flog);
  }

  std::string response;
  TCPerr = extSendAndRecv(message,response,flog,blockWrite,verbose);
  if(TCPerr != PEN_TCP_SUCCESS){
    return LB_TCP_COMMUNICATION_FAIL;
  }

  // Process the response
  //**********************
  response = pen_tcp::trim(response);  
  if(verbose > 1){
    if(blockWrite) const std::lock_guard<std::mutex> lockTh(mtxfth);
    fprintf(flog,"%07ld s - LB::task::extAssign: Received response:\n%s\n",
	    static_cast<long int>(timeStamp()),response.c_str());
    fflush(flog);
  }

  unsigned long long ext_assigned;
  int errParse = parseAssigned(response,ext_assigned,serverErr);
  if(errParse != LB_SUCCESS){
    if(errParse == LB_ERROR_UNEXPECTED_FORMAT){
      if(verbose > 0){
	if(blockWrite) const std::lock_guard<std::mutex> lockTh(mtxfth);
	fprintf(flog,"%07ld s - LB::task::extAssign: Error: Unexpected "
		"response format\n",
		static_cast<long int>(timeStamp()));
	fflush(flog);
      }
    }
    return errParse;
  }

  // Update assignation
  //********************
#ifdef _PEN_USE_MPI_
  setIterationsMPI(ext_assigned);
#else
  if(trusted)
    setIterationsTrusted(ext_assigned);
  else
    setIterations(ext_assigned);
#endif

  if(verbose > 1){
    if(blockWrite) const std::lock_guard<std::mutex> lockTh(mtxfth);
    fprintf(flog,"%07ld s - LB::task::extAssign: Updated iterations to: %llu\n",
	    static_cast<long int>(timeStamp()),ext_assigned);
    fflush(flog);
  }
  
  return LB_SUCCESS;
}

// ********************** HTTP SUPPORT **********************
#ifdef _PEN_USE_LB_HTTP_

static size_t curlCallback(char* ptr,
			   size_t size,
			   size_t nmemb,
			   void* userdata){
  static_cast<std::string*>(userdata)->append(ptr,nmemb);
  return size*nmemb;
}

int LB::task::httpGet(const std::string& url,
		      std::string& bodyResponse){

  //Check if MPI is enabled
  FILE* flog = fth;
  bool blockWrite = true;
#ifdef _PEN_USE_MPI_
  int rank;
  MPI_Comm_rank(comm, &rank);
  if(rank != 0)
    return LB_ERROR_MPI_RANK_0_EXPECTED;
  flog = fMPI;
  blockWrite = false;
#endif
  
  CURL *curl;
  CURLcode res;
  //Init curl client
  curl = curl_easy_init();
  //Check curl status
  if(curl) {

    if(blockWrite) const std::lock_guard<std::mutex> lockTh(mtxfth);
    fprintf(flog,"%07ld s - LB::task::httpGet: "
	    "Send GET petition: %s\n",
	    static_cast<long int>(timeStamp()),url.c_str());
    fflush(flog);
	
    curl_easy_setopt(curl, CURLOPT_URL, url.c_str());

    //Create write callback function
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION,curlCallback);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, static_cast<void*>(&bodyResponse));
    
    //Perform the request
    res = curl_easy_perform(curl);
    if(res != CURLE_OK){
      if(blockWrite) const std::lock_guard<std::mutex> lockTh(mtxfth);
      fprintf(flog,"%07ld s - LB::task::httpGet: "
	      "Petition failed: %s\n",
	      static_cast<long int>(timeStamp()),curl_easy_strerror(res));
      fflush(flog);
      return LB_CURL_ERROR;
    }

    //Clean up
    curl_easy_cleanup(curl);

    //Print the response
    if(blockWrite) const std::lock_guard<std::mutex> lockTh(mtxfth);
    fprintf(flog,"%07ld s - LB::task::httpGet: "
	    "Response received: %s\n",
	    static_cast<long int>(timeStamp()),bodyResponse.c_str());
    fflush(flog);
	
    return LB_SUCCESS;
  }
  else{
    if(blockWrite) const std::lock_guard<std::mutex> lockTh(mtxfth);
    fprintf(flog,"%07ld s - LB::task::httpGet: "
	    "Error Unable to init curl\n",
	    static_cast<long int>(timeStamp()));
    fflush(flog);
    return LB_CURL_ERROR;
  }
}

int LB::task::extStartHTTP(int& serverErr,
			   const unsigned verbose,
			   bool trusted){
  //Check if MPI is enabled
  FILE* flog = fth;
  bool blockWrite = true;
#ifdef _PEN_USE_MPI_
  (void)trusted; //Avoid unused warning
  int rank;
  MPI_Comm_rank(comm, &rank);
  if(rank != 0)
    return LB_ERROR_MPI_RANK_0_EXPECTED;
  flog = fMPI;
  blockWrite = false;
#endif

  
  //Create parameters sufix
  char sufix[600];
  snprintf(sufix,600,"?worker=%ld&dt=%ld",
	   ext_iworker,static_cast<long int>(timeStamp()));
	
  //Set URL
  std::string url(startURL);
  url.append(sufix);

  //Send request
  std::string response;
  int ret = httpGet(url,response);
  if(ret != LB_SUCCESS){
    return ret;
  }
  
  // Process the response
  //**********************
  response = pen_tcp::trim(response);  
  unsigned long long ext_assigned;
  int errParse = parseAssigned(response,ext_assigned,serverErr);
  if(errParse != LB_SUCCESS){
    if(errParse == LB_ERROR_UNEXPECTED_FORMAT){
      if(verbose > 0){
	if(blockWrite) const std::lock_guard<std::mutex> lockTh(mtxfth);
	fprintf(flog,"%07ld s - LB::task::extStartHTTP: Error: Unexpected "
		"response format\n",
		static_cast<long int>(timeStamp()));
	fflush(flog);
      }
    }
    return errParse;
  }

  // Update assignation
  //********************
#ifdef _PEN_USE_MPI_
  setIterationsMPI(ext_assigned);
#else
  if(trusted)
    setIterationsTrusted(ext_assigned);
  else
    setIterations(ext_assigned);
#endif

  if(verbose > 1){
    if(blockWrite) const std::lock_guard<std::mutex> lockTh(mtxfth);
    fprintf(flog,"%07ld s - LB::task::extStartHTTP: Updated iterations "
	    "to: %llu\n",
	    static_cast<long int>(timeStamp()),ext_assigned);
    fflush(flog);
  }
  
  return LB_SUCCESS;  
}

int LB::task::extReportHTTP(int& serverErr,
			    const unsigned long long nIter,
			    const std::chrono::steady_clock::time_point& rep_time,
			    const unsigned verbose,
			    bool trusted){
  //Check if MPI is enabled
  FILE* flog = fth;
  bool blockWrite = true;
#ifdef _PEN_USE_MPI_
  (void)trusted; //Avoid unused warning
  int rank;
  MPI_Comm_rank(comm, &rank);
  if(rank != 0)
    return LB_ERROR_MPI_RANK_0_EXPECTED;
  flog = fMPI;
  blockWrite = false;
#endif
  
  //Get elapsed time since task init
  long int sinceStart = static_cast<long int>(timeStamp(rep_time));
  
  //Create parameters sufix
  char sufix[600];
  snprintf(sufix,600,"?worker=%ld&nIter=%llu&dt=%ld",
	   ext_iworker,nIter,sinceStart);
	
  //Set URL
  std::string url(reportURL);
  url.append(sufix);

  //Send request
  std::string response;
  int ret = httpGet(url,response);
  if(ret != LB_SUCCESS){
    return ret;
  }
  
  // Process the response
  //**********************
  response = pen_tcp::trim(response);  
  unsigned long long ext_assigned;
  int errParse = parseAssigned(response,ext_assigned,serverErr);
  if(errParse != LB_SUCCESS){
    if(errParse == LB_ERROR_UNEXPECTED_FORMAT){
      if(verbose > 0){
	if(blockWrite) const std::lock_guard<std::mutex> lockTh(mtxfth);
	fprintf(flog,"%07ld s - LB::task::extReportHTTP: Error: Unexpected "
		"response format\n",
		static_cast<long int>(timeStamp()));
	fflush(flog);
      }
    }
    return errParse;
  }

  // Update assignation
  //********************
#ifdef _PEN_USE_MPI_
  setIterationsMPI(ext_assigned);
#else
  if(trusted)
    setIterationsTrusted(ext_assigned);
  else
    setIterations(ext_assigned);
#endif

  if(verbose > 1){
    if(blockWrite) const std::lock_guard<std::mutex> lockTh(mtxfth);
    fprintf(flog,"%07ld s - LB::task::extReportHTTP: Updated iterations "
	    "to: %llu\n",
	    static_cast<long int>(timeStamp()),ext_assigned);
    fflush(flog);
  }
  
  return LB_SUCCESS;  
}

int LB::task::extFinishHTTP(int& serverErr,
			    const unsigned long long nIter,
			    const std::chrono::steady_clock::time_point& rep_time,
			    const unsigned verbose){
  //Check if MPI is enabled
  FILE* flog = fth;
  bool blockWrite = true;
#ifdef _PEN_USE_MPI_
  int rank;
  MPI_Comm_rank(comm, &rank);
  if(rank != 0)
    return LB_ERROR_MPI_RANK_0_EXPECTED;
  flog = fMPI;
  blockWrite = false;
#endif
  
  //Get elapsed time since task init
  long int sinceStart = static_cast<long int>(timeStamp(rep_time));
  
  //Create parameters sufix
  char sufix[600];
  snprintf(sufix,600,"?worker=%ld&nIter=%llu&dt=%ld",
	   ext_iworker,nIter,sinceStart);
	
  //Set URL
  std::string url(finishURL);
  url.append(sufix);

  //Send request
  std::string response;
  int ret = httpGet(url,response);
  if(ret != LB_SUCCESS){
    return ret;
  }
  
  // Process the response
  //**********************
  response = pen_tcp::trim(response);  

  //Extract error flag from the response
  int flag;
  int nread = sscanf(response.c_str(),"%d",&flag);
  if(nread != 1){
    if(verbose > 0){
      if(blockWrite) const std::lock_guard<std::mutex> lockTh(mtxfth);
      fprintf(flog,"%07ld s - LB::task::extFinishHTTP: Error: Unexpected "
	      "response format\n",
	      static_cast<long int>(timeStamp()));
      fflush(flog);
    }
    return LB_ERROR_UNEXPECTED_FORMAT;    
  }

  if(flag != 0){
    nread = sscanf(response.c_str(),"%*d %d",&serverErr);
    if(nread != 1){
      if(verbose > 0){
	if(blockWrite) const std::lock_guard<std::mutex> lockTh(mtxfth);
	fprintf(flog,"%07ld s - LB::task::extFinishHTTP: Error: Unexpected "
		"response format\n",
		static_cast<long int>(timeStamp()));
	fflush(flog);
      }
      return LB_ERROR_UNEXPECTED_FORMAT;    
    }
    return LB_REMOTE_ERROR;
  }
  
  return LB_SUCCESS;  
}

#endif
// **********************   HTTP END   **********************

void LB::task::extStartHandler(const unsigned retries,
			       std::chrono::seconds sleeptime,
			       const unsigned verbose,
			       const bool trusted){
  //Check if external balance has been configured
  if(ext_balance){

    // ********************** HTTP SUPPORT **********************
#ifdef _PEN_USE_LB_HTTP_
    // HTTP
    //*******
    if(ext_http){
      for(unsigned i = 0; i < retries; ++i){
	int serverErr;
	int ret = extStartHTTP(serverErr,verbose,trusted);
	if(ret == LB_SUCCESS || ret == LB_ERROR_MPI_RANK_0_EXPECTED)
	  break;
	
	std::this_thread::sleep_for(sleeptime);	
      }
      return;
    }
#endif
    // **********************   HTTP END   **********************
    
    // Raw TCP
    //**********
    //Send start request
    for(unsigned i = 0; i < retries; ++i){
      int TCPerr;
      int serverErr;
      int ret = extStart(TCPerr,serverErr,verbose,trusted);
      if(ret == LB_SUCCESS || ret == LB_ERROR_MPI_RANK_0_EXPECTED)
	break;
	
      std::this_thread::sleep_for(sleeptime);	
    }
    
  }
}

void LB::task::extReportHandler(const unsigned long long nIterDone,
				const std::chrono::steady_clock::time_point& rep_time,
				const std::chrono::seconds sleeptime,
				const unsigned verbose,
				const bool trusted){
  if(ext_balance){

    // ********************** HTTP SUPPORT **********************
#ifdef _PEN_USE_LB_HTTP_
    // HTTP
    //*******
    if(ext_http){
      int serverErr;
      int errRep = extReportHTTP(serverErr,nIterDone,rep_time,verbose,trusted);
      if(errRep != LB_SUCCESS){
	if(errRep != LB_REMOTE_ERROR){
	  //Connection error, try again
	  std::this_thread::sleep_for(sleeptime);	
	  errRep = extReportHTTP(serverErr,nIterDone,rep_time,verbose,trusted);
	  if(errRep == LB_SUCCESS)
	    return;
	}

	if(errRep == LB_REMOTE_ERROR &&
	   (serverErr == LB_ERROR_WORKER_OUT_OF_RANGE ||
	    serverErr == LB_ERROR_WORKER_NOT_STARTED)){
	  //Worker start has not been registered. Send a start petition
	  extStartHandler(2,sleeptime,verbose);
	  extReportHTTP(serverErr,nIterDone,rep_time,verbose,trusted);
	}
      }     
      return;
    }
#endif
    // **********************   HTTP END   **********************
    
    int TCPerr;
    int serverErr;
    int errRep = extReport(TCPerr,serverErr,nIterDone,rep_time,verbose,trusted);
    if(errRep != LB_SUCCESS){
      if(errRep == LB_TCP_CONNECTION_FAIL ||
	 errRep == LB_TCP_COMMUNICATION_FAIL){
	//Connection error, try again
	std::this_thread::sleep_for(sleeptime);	
	errRep = extReport(TCPerr,serverErr,nIterDone,rep_time,verbose,trusted);
	if(errRep == LB_SUCCESS)
	  return;
      }

      if(errRep == LB_REMOTE_ERROR &&
	 (serverErr == LB_ERROR_WORKER_OUT_OF_RANGE ||
	  serverErr == LB_ERROR_WORKER_NOT_STARTED)){
	//Worker start has not been registered. Send a start petition
	extStartHandler(2,sleeptime,verbose);
	extReport(TCPerr,serverErr,nIterDone,rep_time,verbose,trusted);
      }
    }
  }
}

void LB::task::extFinishHandler(const unsigned long long nIterDone,
				const std::chrono::steady_clock::time_point& rep_time,
				const std::chrono::seconds sleeptime,
				const unsigned verbose){
  if(ext_balance){

    // ********************** HTTP SUPPORT **********************
#ifdef _PEN_USE_LB_HTTP_
    // HTTP
    //*******
    if(ext_http){
      int serverErr;
      int errRep = extFinishHTTP(serverErr,nIterDone,rep_time,verbose);
      if(errRep != LB_SUCCESS){
	if(errRep != LB_REMOTE_ERROR){
	  //Connection error, try again
	  std::this_thread::sleep_for(sleeptime);	
	  errRep = extFinishHTTP(serverErr,nIterDone,rep_time,verbose);
	  if(errRep == LB_SUCCESS)
	    return;
	}

	if(errRep == LB_REMOTE_ERROR &&
	   (serverErr == LB_ERROR_WORKER_OUT_OF_RANGE ||
	    serverErr == LB_ERROR_WORKER_NOT_STARTED)){
	  //Worker start has not been registered. Send a start petition
	  extStartHandler(2,sleeptime,verbose);
	  extFinishHTTP(serverErr,nIterDone,rep_time,verbose);
	}
      }     
      return;
    }
#endif
    // **********************   HTTP END   **********************
    
    int TCPerr;
    int serverErr;
    int errRep = extFinish(TCPerr,serverErr,nIterDone,rep_time,verbose);
    if(errRep != LB_SUCCESS){
      if(errRep == LB_TCP_CONNECTION_FAIL ||
	 errRep == LB_TCP_COMMUNICATION_FAIL){
	//Connection error, try again
	std::this_thread::sleep_for(sleeptime);	
	errRep = extFinish(TCPerr,serverErr,nIterDone,rep_time,verbose);
	if(errRep == LB_SUCCESS)
	  return;
      }

      if(errRep == LB_REMOTE_ERROR  &&
	 (serverErr == LB_ERROR_WORKER_OUT_OF_RANGE ||
	  serverErr == LB_ERROR_WORKER_NOT_STARTED)){
	//Try to register the worker start and, then, finish it
	extStartHandler(2,sleeptime,verbose);
	extFinish(TCPerr,serverErr,nIterDone,rep_time,verbose);
      }
    }
  }
}

int LB::task::setWorkers(const size_t nw, const unsigned verbose){

  //Clear current configuration and prepare a new one
  //fort the specified number of workers.

  // NOTE: This function is NOT thread-safe. Ensure that
  // any other thread uses this class instance during
  // the setWorkers call.
  
  if(nw < 1){
    if(verbose > 0){
      const std::lock_guard<std::mutex> lockTh(mtxfth);      
      fprintf(fth,"%07ld s - LB:task:setworkers: Error: Require, at least, "
	      "one worker.\n",
	      static_cast<long int>(timeStamp()));
      fflush(fth);
    }
    return LB_ERROR_INVALID_WORKER_NUMBER;
  }

  lockWorkers();
  
  if(started && !finished){
    unlockWorkers();
    if(verbose > 0){
      const std::lock_guard<std::mutex> lockTh(mtxfth);      
      fprintf(fth,"%07ld s - LB:task:setworkers: Error: This task "
	      "is already running.\n",
	      static_cast<long int>(timeStamp()));
      fflush(fth);
    }
    return LB_ERROR_CURRENT_TASK_RUNNING;
  }

  unlockWorkers();
  
  //Clear previous configuration
  clear();
  //Create mutex locks
  std::vector<std::mutex> auxMutex(nw);
  locks.swap(auxMutex);
  //Create workers
  workers.resize(nw);

  nworkers = nw;

  return LB_SUCCESS;
}

int LB::task::workerStart(const size_t iw, const unsigned verbose){

  //Before tasks start, the number of workers
  //can change. So, lock all workers before check
  lockWorkers();

  if(iw >= nworkers){
    //Release all locks
    unlockWorkers();
    return LB_ERROR_WORKER_OUT_OF_RANGE;
  }

  const auto& actualTime = std::chrono::steady_clock::now();
  
  if(verbose > 1){
    const std::lock_guard<std::mutex> lockTh(mtxfth);      
    fprintf(fth,"%07ld s - LB:task:worker %lu: Enter the task.\n",
	    static_cast<long int>(timeStamp()),
	    static_cast<unsigned long>(iw));
    fflush(fth);
  }
  
  //Check if this task has started already
  if(!started){
    //Get the timestamp when the task begins
    if(verbose > 1){
      const std::lock_guard<std::mutex> lockTh(mtxfth);      
      fprintf(fth,"%07ld s - LB:task:worker %lu: Initialize the task.\n",
	      static_cast<long int>(timeStamp()),
	      static_cast<unsigned long>(iw));
      fflush(fth);
    }
    
    taskBegin = actualTime;
    lastCheckPoint = taskBegin;
    started = true; //Task begins now
    finished = false;

    // **** External LB **** //
    //If enabled, send a start request to external balancer
    extStartHandler(5,std::chrono::seconds(5),verbose);
    // **** External LB **** //
    
// **************************** MPI ********************************* //
#ifdef _PEN_USE_MPI_
    //Start MPI monitor in a thread
    std::thread th(&LB::task::monitor,this,3);
    th.detach();
// **************************** MPI ********************************* //
#endif
  }

  //Get worker
  worker& taskWorker = workers[iw];

  //Calculate iterations to assign
  unsigned long long toDo = remaining/nworkers;

  //Start the worker
  int ret = taskWorker.start(actualTime,toDo);

  //Release all worker locks
  unlockWorkers();
	    
  return ret;
}

bool LB::task::workerFinish(const size_t iw,
			    int& why,
			    const unsigned verbose){
  
  //Ask to task manager if the specified
  //worker can finish this task.
  
  if(iw >= nworkers){
    //Invalid worker, return false
    why = -1;
    return false; 
  }

  if(verbose > 1){
    const std::lock_guard<std::mutex> lockTh(mtxfth);      
    fprintf(fth,"%07ld s - LB:task:workerFinish: Worker %lu: Asking "
	    "to finish.\n",
	    static_cast<long int>(timeStamp()),iw);
    fflush(fth);
  }
  
  { //Create a local scope for lock_guard
    //Lock this worker
    const std::lock_guard<std::mutex> lock(locks[iw]);

    //Get worker reference
    worker& taskWorker = workers[iw];
  
    //Check if done iterations reached assigned
    if(taskWorker.done() < taskWorker.assigned){
      //This worker must do more iterations
      if(verbose > 1){
	const std::lock_guard<std::mutex> lockTh(mtxfth);      
	fprintf(fth,"%07ld s - LB:task:workerFinish: Worker %lu: Permission "
		"denied, insuficient histories done.\n"
		"       Done: %llu\n"
		"   Assigned: %llu\n",
		static_cast<long int>(timeStamp()),iw,
		taskWorker.done(),taskWorker.assigned);
	fflush(fth);
      }
      why = 0;
      return false;
    }

    //All requested iterations has been done,
    //check if is necessary to balance the workload
    if(remainingTime > threshold){
      //Is necessary to create a new checkpoint
      if(verbose > 1){
	const std::lock_guard<std::mutex> lockTh(mtxfth);
	fprintf(fth,"%07ld s - LB:task:workerFinish: Worker %lu: Permission "
		"denied, local balance required.\n",
		static_cast<long int>(timeStamp()),iw);
	fflush(fth);
      }
      why = 1;
      return false;
    }

// **************************** MPI ********************************* //
#ifdef _PEN_USE_MPI_
    //Finally, check if rank 0 allows us to finish the task
    if(!checkFinishMPI(verbose)){
      if(verbose > 1){
	const std::lock_guard<std::mutex> lockTh(mtxfth);
	fprintf(fth,"%07ld s - LB:task:workerFinish: Worker %lu: Permission "
		"denied, awaiting response from rank 0.\n",
		static_cast<long int>(timeStamp()),iw);
	fflush(fth);
      }
      why = 2;
      return false;
    }
#endif
// **************************** MPI ********************************* //
    
    //This worker can finish the task.
    //Set the finished flag to true.
    taskWorker.end();
  
  }
  //Check if all workers have finished
  lockWorkers();
  if(!finished){
    bool someDoing = false;
    for(const LB::worker& worker : workers){
      if(!worker.finishes()){
	someDoing = true;
	break;
      }
    }
    if(!someDoing){
      finished = true;
      if(verbose > 1){
	const std::lock_guard<std::mutex> lockTh(mtxfth);      
	fprintf(fth,"%07ld s - LB:task:workerFinish: Last worker has been "
		"finished the task.\n",
		static_cast<long int>(timeStamp()));
	fflush(fth);
      }

  // **** External LB **** //      
#ifndef _PEN_USE_MPI_
      //If MPI is not enabled, external balance finish requests
      //are handled here
      if(ext_balance){
	//Send finish request
	extFinishHandler(doneTrusted(),
			 std::chrono::steady_clock::now(),
			 std::chrono::seconds(5),verbose);
      }
#endif
  // **** External LB **** //      
    }
  }
  unlockWorkers();

  if(verbose > 1){
    const std::lock_guard<std::mutex> lockTh(mtxfth);      
    fprintf(fth,"%07ld s - LB:task:workerFinish: Worker %lu: Permission "
	    "granted.\n",
	    static_cast<long int>(timeStamp()),iw);
    fflush(fth);
  }
  
  return true;
  
}

std::chrono::seconds::rep LB::task::report(const size_t iw,
					   const unsigned long long nIter,
					   int* errRet,
					   const unsigned verbose){
  
  //Updates worker information according to the number of
  //iterations done and the last report

  //We suppose that this function is not called before the
  //task start (before the first call to "workerStart").
  //With this assumption, we only need to lock the worker
  //specified by "iw".

  if(errRet != nullptr) (*errRet) = LB_SUCCESS;
  
  if(iw >= nworkers){
    if(verbose > 0){
      const std::lock_guard<std::mutex> lockTh(mtxfth);
      fprintf(fth,"%07ld s - LB:task:report:Error: Worker %lu out of range.\n",
	      static_cast<long int>(timeStamp()),iw);
      fflush(fth);
    }
    if(errRet != nullptr) (*errRet) = LB_ERROR_WORKER_OUT_OF_RANGE;
    return std::chrono::seconds::rep(-1);
  }
  
  //Lock this worker
  const std::lock_guard<std::mutex> lock(locks[iw]);

  //Get worker reference
  worker& taskWorker = workers[iw];
  
  //Check if this worker is doing this task
  if(!taskWorker.begins()){
    if(verbose > 0){
      const std::lock_guard<std::mutex> lockTh(mtxfth);      
      fprintf(fth,"%07ld s - LB:task:report:Error: Worker %lu "
	      "doesn't started.\n",
	      static_cast<long int>(timeStamp()),iw);
      fflush(fth);
    }
    if(errRet != nullptr) (*errRet) = LB_ERROR_WORKER_NOT_STARTED;
    return std::chrono::seconds::rep(-1);
  }
  
  if(taskWorker.finishes()){
    if(verbose > 0){
      const std::lock_guard<std::mutex> lockTh(mtxfth);      
      fprintf(fth,"%07ld s - LB:task:report:Error: Worker %lu has finished.\n",
	      static_cast<long int>(timeStamp()),iw);
      fflush(fth);
    }
    if(errRet != nullptr) (*errRet) = LB_ERROR_WORKER_FINISHED;
    return std::chrono::seconds::rep(-1);
  }
  
  //Ensures that the number of iterations has increased
  if(taskWorker.done() > nIter){
    if(verbose > 0){
      const std::lock_guard<std::mutex> lockTh(mtxfth);      
      fprintf(fth,"%07ld s - LB:task:report:Error: Worker %lu unconsistent "
	      "iteration progress.\n",
	      static_cast<long int>(timeStamp()),iw);
      fflush(fth);
    }
    if(errRet != nullptr) (*errRet) = LB_ERROR_LESS_ITERATIONS_DONE;
    return std::chrono::seconds::rep(-1);
  }

  //Register the measure
  const auto actualTime = std::chrono::steady_clock::now();
  std::chrono::seconds::rep nextReport = taskWorker.elapsed(actualTime);
  float dev = taskWorker.addMeasure(actualTime,nIter);

  if(taskWorker.registers() < 2){
    nextReport = std::chrono::seconds::rep(minTimeCheck/1.5);
    //Notice that minTimeCheck can't be changed if this
    //function has the worker lock
  }
  
  //Recalculate the time until next report
  dev = fabs(dev-1.0);
  //Check if the deviation is too big
  if(dev > maxDeviation){
    //Decrease the time until next updates
    nextReport *= std::max((1.0-(dev-maxDeviation)),0.8);
  }else if(dev < 0.1*maxDeviation){ //Check if the deviation is too small
    //Increase the time between updates
    nextReport *= std::min((1.0+(0.5*maxDeviation-dev)),1.2);
  }

  //Check next report limits
  if(minTimeCheck <= nextReport){
    nextReport = std::max(0.8*minTimeCheck,10.0);
  }

  if(nextReport < static_cast<std::chrono::seconds::rep>(10))
    nextReport = static_cast<std::chrono::seconds::rep>(10);
  
  //Return the time until next report
  return nextReport;
}

int LB::task::checkPoint(const unsigned verbose){

  //Lock all workers
  lockWorkers();

  if(nworkers < 1){
    if(verbose > 0){
      const std::lock_guard<std::mutex> lockTh(mtxfth);      
      fprintf(fth,"%07ld s - LB:task:checkPoint:Error: No available workers.\n",
	      static_cast<long int>(timeStamp()));
      fflush(fth);
    }
    unlockWorkers();
    return LB_WARNING_NO_AVAILABLE_WORKERS;
  }
  
  if(!started || finished || iterations == 0){
    if(verbose > 1){
      const std::lock_guard<std::mutex> lockTh(mtxfth);      
      fprintf(fth,"%07ld s - LB:task:checkPoint:Warning: Nothing to do.\n",
	      static_cast<long int>(timeStamp()));
      fflush(fth);
    }
    unlockWorkers();
    return LB_WARNING_NOTHING_TO_DO;
  }
  
  //Get actual time
  std::chrono::steady_clock::time_point actualTime =
    std::chrono::steady_clock::now();
  
  //Calculate elapsed time since previous checkpoint
  std::chrono::seconds::rep timeInterval = std::chrono::duration_cast<
    std::chrono::seconds>(actualTime-lastCheckPoint).count();

  //Check if a checkpoint has been done recently
  if(timeInterval < 1){
    if(verbose > 2){
      const std::lock_guard<std::mutex> lockTh(mtxfth);      
      fprintf(fth,"%07ld s - LB:task:checkPoint: Warning: Insufficient time "
	      "since previous checkpoint.\n",
	      static_cast<long int>(timeStamp()));
      fflush(fth);
    }
    timeInterval = 1;
  }

  //Update last check point timestamp
  lastCheckPoint = actualTime;
  
  //Sum up the speeds of all started and not finished
  //workers to update iteration assignation 
  float taskSpeed = 0.0;
  unsigned long long iterDone = 0;
  unsigned long long iterPred = 0;
  for(const worker& taskWorker : workers){
    //Add iterations done of all workers
    iterDone += taskWorker.done();
    //Avoid finished or non started workers to
    //calculate the current task speed
    if(taskWorker.begins() && !taskWorker.finishes()){
      taskSpeed += taskWorker.speed();
      iterPred += taskWorker.predDone(actualTime);
    }
    else
      iterPred += taskWorker.done();
  }
  
  // **** External LB **** //
  bool forceUpdate = false;
#ifndef _PEN_USE_MPI_
  //If MPI is not enabled, external balance reports are handled at checkpoints
  unsigned long long prevIterations = iterations;
  extReportHandler(iterPred,actualTime,
		   std::chrono::seconds(5),verbose,true);
  if(prevIterations != iterations)
      forceUpdate = true;
#endif
  // **** External LB **** //

  if(iterations <= iterDone){
    
    for(worker& taskWorker : workers){
      //Assign the number of iterations
      //to do equal to the iterations done
      if(taskWorker.begins() && !taskWorker.finishes()){
	taskWorker.assigned = taskWorker.done();
      }
    }
    remaining = 0;
    remainingTime = 0;

    if(verbose > 1){
      const std::lock_guard<std::mutex> lockTh(mtxfth);      
      fprintf(fth,"%07ld s - LB:task:checkPoint: Objective "
	      "iterations reached.\n",
	      static_cast<long int>(timeStamp()));
      fflush(fth);
    }
    
    //Release all workers
    unlockWorkers();
    return LB_SUCCESS;
  }
  
  //Update remaining iterations
  if(iterPred < iterations)
    remaining = iterations - iterPred;
  else
    remaining = 0;

  //Estimate the required time to finish the task
  remainingTime = remaining/taskSpeed;

  if(verbose > 2){
    const std::lock_guard<std::mutex> lockTh(mtxfth);      
    fprintf(fth,"%07ld s - LB:task:checkPoint: Predicted remaining %llu "
	    "iterations, ETA: %llu s .\n",
	    static_cast<long int>(timeStamp()),
	    remaining,remainingTime);
    fflush(fth);
  }
  
  //Calculate remaining registered iterations
  remaining = iterations - iterDone;
  
  //Check if is necessary to balance the workload
// **************************** MPI ********************************* //
#ifdef _PEN_USE_MPI_
  mtxMPI.lock();
  if(remainingTime > threshold || !MPIfinished){
    if(verbose > 2){
      const std::lock_guard<std::mutex> lockTh(mtxfth);      
      if(!MPIfinished)
	fprintf(fth,"%07ld s - LB:task:checkPoint: MPI balance has not "
		"finished, perform a balancing process.\n",
		static_cast<long int>(timeStamp()));
      else
	fprintf(fth,"%07ld s - LB:task:checkPoint: Remaining time above "
		"the threshold, perform a balancing process.\n",
		static_cast<long int>(timeStamp()));
      fflush(fth);
    }
    mtxMPI.unlock();
#else
  // **************************** MPI ********************************* //
    if(remainingTime > threshold || forceUpdate){
      if(verbose > 2){
	const std::lock_guard<std::mutex> lockTh(mtxfth);      
	fprintf(fth,"%07ld s - LB:task:checkPoint: Remaining time above "
		"the threshold, perform a balancing process.\n",
		static_cast<long int>(timeStamp()));
	fflush(fth);
      }
#endif
      
      //Store the number of assigned iterations
      unsigned long long assigned = 0;  

      //Iterate over all workers to reassign the
      //number of iterations to do
      for(worker& taskWorker : workers){
	//Avoid finished or non started workers
	if(taskWorker.begins() && !taskWorker.finishes()){

	  //Ensure that we reassign only the workload of
	  //workers which have reported, as least, one
	  //speed measurement (First "measures" element is fake)
	  if(taskWorker.registers() > 1){
	    //Calculate new speed factor
	    const float localSpeed = taskWorker.speed();
	    const float speedFact = localSpeed/taskSpeed;

	    //Update assignment
	    unsigned long long toDo =
	      static_cast<unsigned long long>(std::max(remaining*speedFact,1.0f));
	    taskWorker.assigned = taskWorker.done() + toDo;
	    assigned += toDo;

	    if(verbose > 3){
	      const std::lock_guard<std::mutex> lockTh(mtxfth);      
	      fprintf(fth,"%07ld s - LB:task:checkPoint: Worker speed: %f i/s, "
		      "assigned iterations: %llu\n",
		      static_cast<long int>(timeStamp()),
		      localSpeed,taskWorker.assigned);
	      fflush(fth);
	    }
	  }
	}
      }

      //Add residual iterations to the first available worker
      if(remaining > assigned){
	for(worker& taskWorker : workers)
	  if(taskWorker.begins() && !taskWorker.finishes()){
	    taskWorker.assigned += remaining-assigned;
	  }
      }    
    }
  // **************************** MPI ********************************* //
#ifdef _PEN_USE_MPI_  
  else{
    mtxMPI.unlock();
    checkFinishMPI(verbose);
  }
#endif
  // **************************** MPI ********************************* //

  //Release all workers
  unlockWorkers();

  return LB_SUCCESS;
}

void LB::task::clear(){

  if(nworkers < 1) return;

  const std::lock_guard<std::mutex> lockTh(mtxfth);      
  if(fth != nullptr){
    fclose(fth);
    fth = nullptr;
  }
  
  iterations = 0;
  remaining = 0;
  locks.clear();
  workers.clear();

  nworkers = 0;
}

unsigned long long LB::task::predDone(unsigned long long& donePred){
  //Return the prediction of done iterations

  //Get actual time
  std::chrono::steady_clock::time_point actualTime =
    std::chrono::steady_clock::now();

  unsigned long long iterDone = 0;
  donePred = 0;
  unsigned iw = 0;
  for(const LB::worker& taskWorker : workers){
    //Lock this worker
    const std::lock_guard<std::mutex> lock(locks[iw]);
    //Add its prediction of iterations done
    iterDone += taskWorker.done();
    donePred += taskWorker.predDone(actualTime);
    ++iw;
  }
  return iterDone;
}

unsigned long long LB::task::done(){
  //Return the registered done iterations

  unsigned long long iterDone = 0;
  unsigned iw = 0;
  for(const LB::worker& taskWorker : workers){
    //Lock this worker
    const std::lock_guard<std::mutex> lock(locks[iw]);
    //Add its prediction of iterations done
    iterDone += taskWorker.done();
    ++iw;
  }
  return iterDone;
}

int LB::task::setIterationsTrusted(const unsigned long long nIter){

  double lastSpeed = -1.0;
  if(remainingTime > 0)
    lastSpeed=
      static_cast<double>(remaining)/static_cast<double>(remainingTime);
    
  if(nIter > iterations){
    remaining += nIter-iterations;
  }
  else{
    unsigned long long diff = iterations-nIter;
    if(diff > remaining)
      remaining = 0;
    else
      remaining -= diff;
  }
  if(lastSpeed > 0.0)
    remainingTime = remaining/lastSpeed;
  else if(remaining == 0)
    remainingTime = 0;
  else
    remainingTime = 10000000000;
  
  iterations = nIter;
  
  return LB_SUCCESS;
}

// **************************** MPI ********************************* //
#ifdef _PEN_USE_MPI_

int LB::task::setIterationsMPI(const unsigned long long nIter){
  if(commSet){
    //Get rank
    int rank;
    MPI_Comm_rank(comm, &rank);
    if(rank != 0){
      return LB_ERROR_HANDLED_BY_RANK0;
    }
    else{
      MPIiterations = nIter;
    }
    return LB_SUCCESS;
  }  
  else
    return LB_ERROR_MPI_COMM_NOT_SET;
}

int LB::task::setMPI(const MPI_Comm commIn,
		     const int tagReqests,
		     const int tagProcess,
		     const unsigned long long nIter,
		     const unsigned verbose){

  lockWorkers();
	    
  if(started && !finished){
    if(verbose > 0)
      printf("LB:task:setMPI: Error: This task is already running.\n");
    return LB_ERROR_TASK_ALREADY_RUNNING;
  }
  
  clearMPI();
  comm = commIn;
  tagReq = tagReqests;
  tagProc = tagProcess;
  //Get current process "rank" and the total number of processes
  int rank, mpiSize;
  MPI_Comm_rank(commIn, &rank);
  MPI_Comm_size(commIn, &mpiSize);

  //Only rank 0 requires a vector to store node speeds
  if(rank == 0){
    //Save global iterations to do
    MPIiterations = nIter;
    //Create MPI workers
    MPIworkers.resize(mpiSize);
  }

  //Set local iterations for each node. Initially,
  //the task is equitativelly splitted 
  iterations = nIter/mpiSize;
  remaining = iterations;
    
  commSet = true;

  unlockWorkers();

  return LB_SUCCESS;
}

void LB::task::sendReportMPI(std::chrono::steady_clock::time_point& lastReport,
			     const unsigned verbose){
  
  int rank;
  MPI_Comm_rank(comm, &rank);

  //Get actual time
  std::chrono::steady_clock::time_point actualTime =
    std::chrono::steady_clock::now();

  //Calculate elapsed time since last report
  std::chrono::seconds::rep elapsedTime = std::chrono::duration_cast<
    std::chrono::seconds>(actualTime-lastReport).count();

  //Save last report timepoint
  lastReport = actualTime;
  
  if(rank == 0){
    //Rank 0 can update directly its speed information

    if(MPIfinished) //No lock required
      return;
    
    //Perform an estimation of iterations done
    unsigned long long predIterDone;
    predDone(predIterDone);

    if(verbose > 1){
      fprintf(fMPI,"%07ld s - task:monitor:Rank 0: Perform a local report with "
	     "%llu iterations predicted and elapsed time %lu s\n",
	      static_cast<long int>(timeStamp()),
	      predIterDone,static_cast<long int>(elapsedTime));
      fflush(fMPI);
    }
    
    //Add a measure
    MPIworkers[0].addMeasure(elapsedTime,predIterDone);

    //Rebalance MPI workers
    unsigned long long MPIremainingTime = balanceMPI(verbose);

    if(verbose > 1){
      fprintf(fMPI,"%07ld s - task:monitor:Rank 0: New local iteration "
	      "assignment: %llu\n",
	      static_cast<long int>(timeStamp()),
	      MPIworkers[0].assigned);
      fflush(fMPI);
    }
    
    //Store new assignment
    setIterations(MPIworkers[0].assigned);

    if(verbose > 1){
      fprintf(fMPI,"%07ld s - task:monitor:Rank 0: ETA: %llu s\n",
	      static_cast<long int>(timeStamp()),MPIremainingTime);
      fflush(fMPI);
    }
    
    if(MPIremainingTime <= threshold){
      //Perform a local checkpoint
      checkPoint(verbose);
      //Lock MPI
      const std::lock_guard<std::mutex> lock(mtxMPI);
      MPIfinished = true;
    }
    
    return;
  }
  
  //Create a buffer to send report data
  //[0]->predicted iterations
  //[1]->elapsed time (s)
  unsigned long long iterDone[2];

  //Get actual speed and perform an estimation of iterations done
  predDone(iterDone[0]);
  iterDone[1] = static_cast<unsigned long long>(elapsedTime);

  if(verbose > 1){
    fprintf(fMPI,"%07ld s - task:monitor:Rank %d: Send a report with "
	    "%llu iterations predicted and elapsed time %llu s\n",
	    static_cast<long int>(timeStamp()),
	    rank,iterDone[0],iterDone[1]);
    fflush(fMPI);
  }
  
  //send to rank 0 the actual
  //number of iterations done and the elapsed time
  MPI_Send(static_cast<void*>(iterDone),2,
	   MPI_UNSIGNED_LONG_LONG,0,tagProc,comm);

  //Receive information from rank 0
  // buffer[0] -> assigned
  // buffer[1] -> MPI finished flag (0-> no, != 0 -> yes)
  unsigned long long buffer[2];
  MPI_Recv(static_cast<void*>(buffer),2,
	   MPI_UNSIGNED_LONG_LONG,0,tagProc,comm, MPI_STATUS_IGNORE);

  if(verbose > 1){
    fprintf(fMPI,"%07ld s - task:monitor:Rank %d: Received new "
	    "assignation: %llu\n",
	    static_cast<long int>(timeStamp()),
	    rank,buffer[0]);
    fflush(fMPI);
  }
  
  //Store new assignment
  setIterations(buffer[0]);

  //Check if MPI balance has finished
  if(buffer[1] != 0){
    //Perform a local checkpoint
    checkPoint(verbose);
    //Lock MPI
    const std::lock_guard<std::mutex> lock(mtxMPI);
    MPIfinished = true;
    if(verbose > 1){
      fprintf(fMPI,"%07ld s - task:monitor:Rank %d: Received MPI balance finish flag.\n",
	      static_cast<long int>(timeStamp()),
	      rank);
      fflush(fMPI);
    }
    return;
  }
}

float LB::task::receiveReportMPI(const unsigned inode,
				 const unsigned verbose){
  int rank;
  MPI_Comm_rank(comm, &rank);

  if(rank != 0){
    if(verbose > 0){
      fprintf(fMPI,"%07ld s - LB:task:receiveReportMPI: Error: Only rank 0 can "
	      "receive reports from workers.\n",
	      static_cast<long int>(timeStamp()));
      fflush(fMPI);
    }
    throw std::invalid_argument("Non rank 0 receiving reports");    
  }

  //Report incoming

  //Take mpi worker
  guessWorker& MPIworker = MPIworkers[inode];
	
  //Create a buffer to receive report data
  //[0]->prediction
  //[1]->elapsed time (s)
  unsigned long long iterDone[2];
  MPI_Recv(static_cast<void*>(iterDone),2,
	   MPI_UNSIGNED_LONG_LONG,inode,tagProc,comm, MPI_STATUS_IGNORE);

  if(verbose > 1){
    fprintf(fMPI,"%07ld s - task:monitor:Rank 0: Report received from rank %d. "
	    "%llu iterations predicted and elapsed time %llu s\n",
	    static_cast<long int>(timeStamp()),
	    inode,iterDone[0],iterDone[1]);
    fflush(fMPI);
  }
  
  //Check if MPI balancing has already finished
  float dev = 1.0;
  if(MPIfinished){  //No lock required
    
    unsigned long long buffer[2] = {MPIworker.assigned,1};
    //Send new assignment and the deviation between
    //calculated and received expected iterations
    MPI_Send(static_cast<void*>(buffer),2,
	     MPI_UNSIGNED_LONG_LONG,inode,tagProc,comm);

    if(verbose > 1){
      fprintf(fMPI,"%07ld s - task:monitor:Rank 0: Sent MPI finished flag with "
	      "Final assignment: %llu\n",
	      static_cast<long int>(timeStamp()),
	      buffer[0]);
      fflush(fMPI);
    }
  }
  else{
    //Register done iterations and the elapsed time
    dev = MPIworker.addMeasure(std::chrono::seconds::rep(iterDone[1]),
			       iterDone[0]);

    //Rebalance MPI workers
    unsigned long long MPIremainingTime = balanceMPI(verbose);
    unsigned long long buffer[2] = {MPIworker.assigned,0};

    //Update local assignation
    setIterations(MPIworkers[0].assigned);

    if(verbose > 1){
      fprintf(fMPI,"%07ld s - task:monitor:Rank 0: ETA: %llu s\n",
	      static_cast<long int>(timeStamp()),MPIremainingTime);
      fflush(fMPI);
    }
    
    if(MPIremainingTime <= threshold){
      buffer[1] = 1;
      if(verbose > 1){
	fprintf(fMPI,"%07ld s - task:monitor:Rank 0: Remaining time "
		"below threshold. Send MPI finished flag\n",
		static_cast<long int>(timeStamp()));
	fflush(fMPI);
      }
      //Perform a local checkpoint
      checkPoint(verbose);
      const std::lock_guard<std::mutex> lock(mtxMPI);
      MPIfinished = true;
      MPIworker.end();
    }
  
    //Send new assignment and the deviation between
    //calculated and received expected iterations
    MPI_Send(static_cast<void*>(buffer),2,
	     MPI_UNSIGNED_LONG_LONG,inode,tagProc,comm);
    
    if(verbose > 1){
      fprintf(fMPI,"%07ld s - task:monitor:Rank 0: Sent new "
	      "iteration assignment: %llu\n",
	      static_cast<long int>(timeStamp()),buffer[0]);
      fflush(fMPI);
    }
  }
  
  return dev;
}

unsigned long long LB::task::balanceMPI(const unsigned verbose){
  
  //Get speed and iterations done
  float totalSpeed;
  unsigned long long totalIterDone = speedMPI(totalSpeed);
  if(totalSpeed <= 0.0)
    return 1000000000000000;

  // **** External LB **** //
  if(ext_balance){
    //Send a report to the external load balancer
    //Get actual time
    std::chrono::steady_clock::time_point actualTime =
      std::chrono::steady_clock::now();
    extReportHandler(predDoneMPI(actualTime),actualTime,
		     std::chrono::seconds(5),verbose,true);
  }
  // **** External LB **** //

  unsigned long long MPIremainingTime = 0ull;
  //Check if there are something to do
  if(totalIterDone >= MPIiterations){
    for(guessWorker& worker : MPIworkers){
      if(worker.begins() && !worker.finishes()){
	worker.assigned = worker.done();
      }
    }
  }
  else{
    //Calculate remaining iterations
    unsigned long long MPIremaining = MPIiterations - totalIterDone;

    //Reassing iterations to each worker
    for(guessWorker& worker : MPIworkers){
      if(worker.begins() && !worker.finishes()){
	unsigned long long toDo = MPIremaining*(worker.speed()/totalSpeed);
	toDo = std::max(toDo,1llu);
	worker.assigned = worker.done() + toDo;
      }
    }
    MPIremainingTime = static_cast<double>(MPIremaining)/totalSpeed;
  }

  // **** External LB **** //
  if(ext_balance && MPIremainingTime <= threshold){
    //Send finish signal with all iterations completed
    std::chrono::steady_clock::time_point actualTime =
      std::chrono::steady_clock::now();
    extFinishHandler(MPIiterations,actualTime + std::chrono::seconds(2),
		     std::chrono::seconds(5),verbose);
  }
  // **** External LB **** //
  
  return MPIremainingTime;
}

void LB::task::monitor(const unsigned verbose){

  const auto oneSecond = std::chrono::seconds(1);
  
  if(!commSet){
    if(verbose > 0){
      fprintf(fMPI,"%07ld s - LB:task:monitor:Error: MPI communication "
	      "not initialized\n",
	      static_cast<long int>(timeStamp()));
      fflush(fMPI);
    }
    throw std::invalid_argument("setMPI not called");      
  }
  
  lockWorkers();
  
  if(!started || finished){
    if(verbose > 0){
      fprintf(fMPI,"%07ld s - LB:task:monitor:Error: task has not started "
	      "or has finished\n",static_cast<long int>(timeStamp()));
      fflush(fMPI);
    }
    throw std::invalid_argument("Non started or finished task");      
  }

  unlockWorkers();
  
  //Get rank and size
  int rank, mpiSize;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &mpiSize);

  if(mpiSize < 2){ //Nothing to balance
    const std::lock_guard<std::mutex> lock(mtxMPI);
    MPIfinished = true;
    return;
  }
  
  //Begin communication loop
  if(rank == 0){

    //
    //  Monitor initialization
    //**************************    
    
    //Create an array of requests, two for each rank. First one is for
    //requests from non 0 rank nodes. The second one to receive messages from
    //the rank 0 at the other ranks.

    //The first "mpiSize" requests correspond to messages sent by non
    //rank 0 nodes. On the other hand, the second half of the request
    //corresponds to messages sent by rank 0. The vector position
    //indicates the non 0 rank involved in the communication:

    // rank = index%mpiSize

    //Vector with schedule report information
    std::vector<LB::reportSchedule>
      reports(mpiSize,LB::reportSchedule(std::chrono::seconds(0),
					 MPIminTimeCheck)
	      );

    //Init count down for MPI worker 0 (local)
    reports[0].reset();
    //We set remaining time to 0 until the corresponding worker starts
    std::chrono::seconds timeUntilNextReport = MPIminTimeCheck;

    //Vector with end of communications sent flags
    std::vector<bool> finishSent;
    finishSent.resize(mpiSize,false);
    
    //Vector of requests
    std::vector<LB::request> requests(2*mpiSize);
      
    //Create an array of instruction identifiers, one for each request
    std::vector<int> instructions;
    instructions.resize(2*mpiSize,-1);

    //Create persistent requests for receive and send operations
    for(int i = 1; i < mpiSize; ++i){ //Avoid rank 0
      int err = requests[i].receiveInit(static_cast<void*>(&instructions[i]),
				     1,MPI_INT,i,tagReq,comm);
      if(err != MPI_SUCCESS){
	if(verbose > 0){
	  fprintf(fMPI,"%07ld s - Unable to create receive requests "
		  "on rank 0\n"
		  "       Error code: %d\n",
		  static_cast<long int>(timeStamp()),err);
	  fflush(fMPI);
	}
	throw std::invalid_argument("Unable to receive request");
      }
    }

    for(size_t i = mpiSize+1; i < requests.size(); ++i){ //Avoid rank 0
      int err = requests[i].sendInit(static_cast<void*>(&instructions[i]),
				     1,MPI_INT,i%mpiSize,tagReq,comm);
      if(err != MPI_SUCCESS){
	if(verbose > 0){
	  fprintf(fMPI,"%07ld s - Unable to create send requests on rank 0\n"
		  "       Error code: %d\n",
		  static_cast<long int>(timeStamp()),err);
	  fflush(fMPI);
	}	  
	throw std::invalid_argument("Unable to send request");
      }
    }
    
    //Beggin the receive requests, avoiding rank 0
    for(int i = 1; i < mpiSize; ++i){
      int err = requests[i].start();
      if(err != MPI_SUCCESS){
	if(verbose > 0){
	  fprintf(fMPI,"%07ld s - Unable to start receive "
		  "requests of rank %d on rank 0\n"
		  "       Error code: %d\n",
		  static_cast<long int>(timeStamp()),
		  i,err);
	  fflush(fMPI);
	}
	throw std::invalid_argument("Unable to start receive requests");
      }
    }

    //Get actual time as last report timepoint
    std::chrono::steady_clock::time_point lastReport =
      std::chrono::steady_clock::now();    

    //Start local worker
    {
      float speed;
      unsigned long long iterDone = speedMPI(speed);
      unsigned long long toDo =
	std::max((MPIiterations-iterDone)/MPIworkers.size(),1llu);
    
      MPIworkers[0].start(std::chrono::steady_clock::now(),toDo);
    }
    
    //
    //  Monitor loop rank 0
    //************************
    if(verbose > 1){
      fprintf(fMPI,"%07ld s - task:monitor:Rank 0: Monitor loop beggins\n",
	      static_cast<long int>(timeStamp()));
      fflush(fMPI);
    }
    {
      const std::lock_guard<std::mutex> lock(mtxMPI);
      MPIfinishReqFlag = false;
    }
    for(;;){

      // * Handle requests
      //---------------------
      
      int index; //Node to communicate with
      //Wait until any request has been completed
      std::chrono::seconds elapsed(10);
      int err = LB::noBusyWaitAny(requests,index,
				  std::chrono::seconds(sleepTime),
				  timeUntilNextReport,
				  elapsed,
				  mtxMPI,MPIfinishReqFlag);
      if(err != MPI_SUCCESS){
	if(verbose > 0){
	  fprintf(fMPI,"%07ld s - LB:task:monitor:Error: Rank 0 is unable to "
		  "test workers requests\n"
		  "           Error code: %d\n",
		  static_cast<long int>(timeStamp()),err);
	  fflush(fMPI);
	}
	throw std::invalid_argument("Unable to test requests on rank 0");
      }

      //Send report petitions if required
      timeUntilNextReport = std::chrono::seconds(1000000000000);
      //If MPI has finished, force a report for all remaining nodes
      if(MPIfinished)  //No lock required
	elapsed = timeUntilNextReport;
      
      for(int i = 0; i < mpiSize; ++i){

	//Skip reports already done and unstarted workers
	if(reports[i].remaining == std::chrono::seconds(0))
	  continue;

	//Check if this node require a report
	if(reports[i].remaining <= elapsed){
	  //Requires a report
	  if(verbose > 1){
	    fprintf(fMPI,"%07ld s - task:monitor:Rank 0: Request a report "
		    "to rank %d\n",
		    static_cast<long int>(timeStamp()),i);
	    fflush(fMPI);
	  }
	  if(i == 0){
	    sendReportMPI(lastReport,verbose);
	    reports[i].reset();
	  }
	  else{
	    reports[i].remaining = std::chrono::seconds(0);
	    //Send the report
	    instructions[mpiSize+i] = 1;
	    int errStart = requests[mpiSize+i].start();
	    if(errStart != MPI_SUCCESS){
	      if(verbose > 0){
		fprintf(fMPI,"%07ld s - LB:task:monitor:Error: Rank 0 is "
			"unable to start report request petition for rank %d.\n"
			"           Error code: %d\n",
			static_cast<long int>(timeStamp()),i,err);
		fflush(fMPI);
	      }
	      throw std::invalid_argument("Unable to test requests on rank 0");
	    }
	  }
	}
	else{
	  //Substract elapsed time
	  reports[i].remaining -= elapsed;
	  if(reports[i].remaining < timeUntilNextReport){
	    timeUntilNextReport = reports[i].remaining;
	  }
	}
      }

      //Check if timeout has been reached
      if(index < 0){
	if(index == -1){
	  //Local finish request received
	  if(verbose > 1){
	    fprintf(fMPI,"%07ld s - task:monitor:Rank 0: Local finish "
		    "request received. Balance MPI workers.\n",
		    static_cast<long int>(timeStamp()));
	    fflush(fMPI);
	  }
	  sendReportMPI(lastReport,verbose);

	  if(MPIfinished){ //No lock required
	    finishSent[0] = true;
	    reports[0].remaining = std::chrono::seconds(0);
	    timeUntilNextReport = std::chrono::seconds(1);
	    //Check if all MPI processes have finished the comunication
	    bool endLoop = true;
	    for(const bool sent : finishSent)
	      if(!sent){endLoop = false; break;}
	    if(endLoop)
	      break; //Exit from communication loop	    
	  }
	  else{
	    reports[0].reset();
	  }
	  //Disable finish request flat
	  const std::lock_guard<std::mutex> lock(mtxMPI);
	  MPIfinishReqFlag = false;	  
	}
	continue;
      }
      
      //Some request completed. Get the destination node
      int inode = index%mpiSize;

      //Send the request source (0 or inode) and the
      //instruction to handle to ensure that both nodes
      //process the very same request.
      int instruction = instructions[index];
      int sourceInstruction[2] = {index < mpiSize ? inode : 0,
				  instruction};
      
      err = noBusySsend(static_cast<void*>(sourceInstruction),
			2,MPI_INT,inode,tagProc,comm,
			oneSecond);
      if(err != MPI_SUCCESS){
	if(verbose > 0){
	  fprintf(fMPI,"%07ld s - LB:task:monitor:Error: Rank 0 is unable to "
		  "send the instruction to process to worker %d.\n"
		  "           Error code: %d\n",
		  static_cast<long int>(timeStamp()),inode,err);
	  fflush(fMPI);
	}
	throw std::invalid_argument("Rank 0 is unable to send the instruction");
      }
      
      //As worker will throw an exception if instructions don't coincide,
      //is not necessary to perform any check at this node.

      //At this point, both, rank 0 and the worker are synchronized
      //and they agree with the request and instruction to process. 
      //So, process the specified instruction. 	

      // * Proces the instruction
      //--------------------------
      if(instruction == 0){
	//Rank "inode" starts this task. 
	//Take mpi worker

	if(verbose > 1){
	  fprintf(fMPI,"%07ld s - task:monitor:Rank 0: Received start "
		  "petition from rank %d\n",
		  static_cast<long int>(timeStamp()),inode);
	  fflush(fMPI);
	}	  
	
	guessWorker& MPIworker = MPIworkers[inode];

	if(MPIworker.begins()){
	  //This node has already started
	  unsigned long long toDo = MPIworker.assigned;
	  MPI_Send(static_cast<void*>(&toDo),1,
		   MPI_UNSIGNED_LONG_LONG,inode,tagProc,comm);	  
	}
	else{
	  //Calculate iterations to assign
	  float speed;
	  unsigned long long iterDone = speedMPI(speed);
	  unsigned long long toDo =
	    std::max((MPIiterations-iterDone)/MPIworkers.size(),1llu);

	  //Try to start the worker
	  int ret = MPIworker.start(std::chrono::steady_clock::now(),toDo);
	  if(ret != LB_SUCCESS){
	    //This worker is doing the task or has finished
	    toDo = MPIworker.assigned;
	  }

	  //Return assigned iterations for this worker
	  MPI_Send(static_cast<void*>(&toDo),1,
		   MPI_UNSIGNED_LONG_LONG,inode,tagProc,comm);

	  //Reset report waiting time
	  reports[inode].reset();

	  if(verbose > 1){
	    fprintf(fMPI,"%07ld s - task:monitor:Rank 0: The initially "
		    "assignment is %llu iterations. Time until first "
		    "report: %ld\n",
		    static_cast<long int>(timeStamp()),toDo,
		    static_cast<long int>(reports[inode].dt()));
	    fflush(fMPI);
	  }
	}
	//Restart receive petition
	requests[index].start();
      }
      else if(instruction == 1){
	//Report incoming	
	float dev = receiveReportMPI(inode,verbose);

	if(MPIfinished){  //No lock required
	  finishSent[0] = true;
	  finishSent[inode] = true;
	  timeUntilNextReport = std::chrono::seconds(1);	  
	  //Check if all MPI processes have finished the comunication
	  bool endLoop = true;
	  for(const bool sent : finishSent)
	    if(!sent){endLoop = false; break;}
	  if(endLoop)
	    break; //Exit from communication loop
	}
	else{
	  dev = fabs(dev-1.0);
	  //Check if the deviation is too big
	  if(dev > maxDeviation){ //Deviation too big
	    //Decrease the time between updates
	    reports[inode].interval *= std::max((1.0-(dev-maxDeviation)),0.8);
	  }
	  else if(dev < 0.1*maxDeviation){ //Deviation too small
	    //Increase the time between updates
	    reports[inode].interval *= std::min((1.0+(0.5*maxDeviation-dev)),1.2);
	  }

	  //Check next report limits
	  if(MPIminTimeCheck > reports[inode].interval){
	    reports[inode].interval = MPIminTimeCheck;
	  }
	  
	  if(reports[inode].remaining == std::chrono::seconds(0))
	    reports[inode].reset();

	  if(verbose > 1){
	    fprintf(fMPI,"%07ld s - task:monitor:Rank 0: Time until "
		    "next rank %d report: %ld s\n",
		    static_cast<long int>(timeStamp()),inode,
		    static_cast<long int>(reports[inode].dt()));
	    fflush(fMPI);
	  }
	}
      }
      else if(instruction == 2){ //Worker want to finish
	//Receive a report from worker
	if(verbose > 1){
	  fprintf(fMPI,"%07ld s - task:monitor:Rank 0: Received finish "
		  "petition from rank %d\n",
		  static_cast<long int>(timeStamp()),inode);
	  fflush(fMPI);
	}
	
	receiveReportMPI(inode,verbose);
	if(MPIfinished){ //No lock required
	  finishSent[0] = true;
	  finishSent[inode] = true;
	  timeUntilNextReport = std::chrono::seconds(1);	  
	  //Check if we expect a report
	  if(reports[inode].remaining != std::chrono::seconds(0)){
	    //No report expected, disable reports
	    reports[inode].remaining = std::chrono::seconds(0);
	  }
	  else{
	    //Report expected, cancel the petition
	    requests[mpiSize+inode].cancel();
	  }
	  //Check if all MPI processes have finished the comunication
	  bool endLoop = true;
	  for(const bool sent : finishSent)
	    if(!sent){endLoop = false; break;}
	  if(endLoop)
	    break; //Exit from communication loop	  
	}
	else{
	  //Restart receive petition
	  requests[index].start();
	}
      }
    }

    if(verbose > 1){
      fprintf(fMPI,"%07ld s - task:monitor:Rank 0: Finishing monitor loop.\n",
	      static_cast<long int>(timeStamp()));
      fflush(fMPI);
    }	
  }
  else{

    //
    //  Monitor initialization
    //**************************    
    
    //First, create a request pair to send and receive
    //communications with rank 0
    std::vector<LB::request> requests(2); //[0]->receive [1]->send

    //Create an instruction identifier for each request.
    int instructions[2]; //[0]->receive [1]->send

    //Create persistent request for receive
    int errInit;
    errInit = requests[0].receiveInit(static_cast<void*>(&instructions[0]),
				      1,MPI_INT,0,tagReq,comm);
    if(errInit != MPI_SUCCESS){
      if(verbose > 0){
	fprintf(fMPI,"%07ld s - LB:task:monitor:Error: Unable to create "
		"the receive requests on rank %d\n"
		"       Error code: %d\n",
		static_cast<long int>(timeStamp()),rank,errInit);
	fflush(fMPI);
      }	
      throw std::invalid_argument("Unable to create receive request");
    }

    //Create persistent request to send with syncronization
    errInit = requests[1].sendInit(static_cast<void*>(&instructions[1]),
				   1,MPI_INT,0,tagReq,comm);
    if(errInit != MPI_SUCCESS){
      if(verbose > 0){
	fprintf(fMPI,"%07ld s - LB:task:monitor:Error: Unable to create "
		"the send request on rank %d\n"
		"       Error code: %d\n",
		static_cast<long int>(timeStamp()),rank,errInit);
	fflush(fMPI);
      }	
      throw std::invalid_argument("Unable to create send request");
    }
    
    //Once requests has been initializated, init the
    //receive process from rank 0
    errInit = requests[0].start();
    if(errInit != MPI_SUCCESS){
      if(verbose > 0){
	fprintf(fMPI,"%07ld s - LB:task:monitor:Error: Unable to start "
		"receive requests on rank %d\n"
		"       Error code: %d\n",
		static_cast<long int>(timeStamp()),rank,errInit);
	fflush(fMPI);
      }	
      throw std::invalid_argument("Unable to start receive request");
    }

    //Then, inform to the rank 0 that this rank
    //has started the task
    if(verbose > 1){
      fprintf(fMPI,"%07ld s - task:monitor:Rank %d: Send to rank 0 "
	      "a start petition\n",
	      static_cast<long int>(timeStamp()),rank);
      fflush(fMPI);
    }
    instructions[1] = 0;
    requests[1].start();
    errInit = noBusyWait(requests[1],oneSecond);
    if(errInit != MPI_SUCCESS){
      if(verbose > 0){
	fprintf(fMPI,"%07ld s - LB:task:monitor:Error: Worker %d is unable "
		"to send the 'start' request to rank 0\n"
		"       Error code: %d\n",
		static_cast<long int>(timeStamp()),rank,errInit);
	fflush(fMPI);
      }	
      throw std::invalid_argument("Unable to send 'start' request");
    }

    //Get, from rank 0, the request source and instruction to handle
    int sourceInstructionInit[2];
    errInit = LB::noBusyRecv(static_cast<void*>(&sourceInstructionInit),
			      2,MPI_INT,0,tagProc,comm,
			      oneSecond);
    if(errInit != MPI_SUCCESS){
      if(verbose > 0){
	fprintf(fMPI,"%07ld s - LB:task:monitor:Error: Worker is unable to "
		"receive request source from rank 0.\n"
		"           Error code: %d\n",
		static_cast<long int>(timeStamp()),errInit);
	fflush(fMPI);
      }	
      throw std::invalid_argument("Worker is unable to receive "
				  "request source");      
    }

    if(sourceInstructionInit[0] != rank || sourceInstructionInit[1] != 0){
      if(verbose > 0){
	fprintf(fMPI,"%07ld s - LB:task:monitor:Error: Rank 0 doesn't "
		"want to process the rank %d start petition.\n"
		"Received instruction to process: %d\n"
		"Instruction received from rank : %d\n",
		static_cast<long int>(timeStamp()),
		rank,sourceInstructionInit[1],sourceInstructionInit[0]);
	fflush(fMPI);
      }
      throw std::invalid_argument("Worker is unable to start communications");
    }
    
    //Get actual time as last report timepoint
    std::chrono::steady_clock::time_point lastReport =
      std::chrono::steady_clock::now();
    
    //Receive the number of iterations to do
    unsigned long long nIter;
    MPI_Recv(static_cast<void*>(&nIter),1,
	     MPI_UNSIGNED_LONG_LONG,0,tagProc,comm, MPI_STATUS_IGNORE);    

    if(verbose > 1){
      fprintf(fMPI,"%07ld s - task:monitor:Rank %d: Received "
	      "initial iteration assignation: %llu\n",
	      static_cast<long int>(timeStamp()),rank,nIter);
      fflush(fMPI);
    }
    
    //Update the number of iterations to do
    setIterations(nIter);

    //
    //  Monitor loop workers
    //************************
    if(verbose > 1){
      fprintf(fMPI,"%07ld s - task:monitor:Rank %d: Monitor loop beggins\n",
	      static_cast<long int>(timeStamp()),rank);
      fflush(fMPI);
    }
    for(;;){

      //Wait until any request has been completed
      int index;
      int err = LB::noBusyFreeWaitAny(requests,index,
				      std::chrono::seconds(sleepTime),
				      mtxMPI, MPIfinishReqFlag);
      if(err != MPI_SUCCESS){
	if(verbose > 0){
	  fprintf(fMPI,"%07ld s - LB:task:monitor:Error: Worker "
		  "is unable to test requests.\n"
		  "       Error code: %d\n",
		  static_cast<long int>(timeStamp()),err);
	  fflush(fMPI);
	}	  
	throw std::invalid_argument("Worker is unable to test "
				    "requests");      
      }

      //Check if some local worker set the finish request flag
      if(index < 0){
	const std::lock_guard<std::mutex> lock(mtxMPI);
	if(!MPIfinishReqFlag || MPIsent){
	  if(verbose > 0){
	    fprintf(fMPI,"%07ld s - LB:task:monitor:Error: This can't happens: "
		    "the monitor finish flag request is false but the "
		    "condition to finish wait has been accomplished!\n"
		    "      Request sent: %s"
		    "   Finish flag set: %s",
		    static_cast<long int>(timeStamp()),
		    MPIsent ? "true" : "false",
		    MPIfinishReqFlag ? "true" : "false");
	    fflush(fMPI);
	  }
	  throw std::invalid_argument("Woerker suffers an 'impossible' error "
				      "on finish request flag");
	}

	//Start the request to finish the task
	instructions[1] = 2; //Flag finish request
	requests[1].start();
	//Flag sent request
	MPIsent = true;
	//Unflag finish request to send
	MPIfinishReqFlag = false;

	if(verbose > 1){
	  fprintf(fMPI,"%07ld s - task:monitor:Rank %d: Sending finish "
		  "petition to rank 0.\n",
		  static_cast<long int>(timeStamp()),rank);
	  fflush(fMPI);
	}
	  
	continue;
      }

      //Get the expected instruction to handle
      int instruction = instructions[index];
      //Get expected request source
      int source = 0;
      if(index != 0) source = rank;
      //Get, from rank 0, the request source and instruction to handle
      int sourceInstruction[2];
      err = LB::noBusyRecv(static_cast<void*>(&sourceInstruction),
			   2,MPI_INT,0,tagProc,comm,
			   oneSecond);
      if(err != MPI_SUCCESS){
	if(verbose > 0){
	  fprintf(fMPI,"%07ld s - LB:task:monitor:Error: Worker is unable to "
		  "receive request source from rank 0.\n"
		  "           Error code: %d\n",
		  static_cast<long int>(timeStamp()),err);
	  fflush(fMPI);
	}
	throw std::invalid_argument("Worker is unable to receive "
				    "request source");      
      }
      
      //Check if the received request to process is the
      //same as the selected locally
      if(source != sourceInstruction[0]){
	//Rank 0 expects to process the other request, update the
	//index, wait it and free the request once received
	++index;
	index %= 2;
	err = LB::noBusyWait(requests[index],oneSecond);
	if(err != MPI_SUCCESS){
	  if(verbose > 0){
	    fprintf(fMPI,"%07ld s - LB:task:monitor:Error: Worker is unable to "
		    "receive the other request from from rank 0.\n"
		    "           Error code: %d\n",
		    static_cast<long int>(timeStamp()),err);
	    fflush(fMPI);
	  }	    
	  throw std::invalid_argument("Worker is unable to receive "
				      "the other request");      
	}
	//Change source and instruction
	source = sourceInstruction[0];
	instruction = instructions[index];
      }
      else{
	//Both expected requests coincide, complete this request
	requests[index].wait();
      }
      
      //Check if both instructions coincide
      if(instruction != sourceInstruction[1]){
	//Rank 0 expects another instruction. This should never happens
	if(verbose > 0){
	  fprintf(fMPI,"%07ld s - LB:task:monitor:Error: Worker and rank 0 "
		  "expects different instructions to process on the "
		  "same request.\n"
		  "           Worker expected: %d\n"
		  "           Rank 0 expected: %d\n"
		  "This should never happens, please, report "
		  "this bug to the developers.\n",
		  static_cast<long int>(timeStamp()),
		  instruction,sourceInstruction[1]);
	  fflush(fMPI);
	}
	throw std::invalid_argument("Unexpected error! Worker and rank0 "
				    "expectes different instructions "
				    "to process for the same request");      	
      }
      
      //At this point, both, rank 0 and this worker are synchronized
      //and they agree with the request and instruction to process. 
      //So, process the specified instruction. 

      if(index == 0){
	//Reset receive request
	requests[0].start();
      }
      
      // * Proces the instruction
      //--------------------------
      
      if(instruction == 1){ //Report
	
	sendReportMPI(lastReport,verbose);
	//Check if MPI balance has been finished
	if(MPIfinished){  //No lock required
	  //Cancel receive request
	  requests[0].cancel();
	  //Cancel send request
	  requests[1].cancel();
	  break; //Exit from communication loop
	}
      }
      else if(instruction == 2){ //Require to finish the task

	if(verbose > 1){
	  fprintf(fMPI,"%07ld s - task:monitor:Rank %d: Handle a "
		  "finish request.\n",
		  static_cast<long int>(timeStamp()),rank);
	  fflush(fMPI);
	}
	
	//First, send an updated report
	sendReportMPI(lastReport,verbose);
	
	if(MPIfinished){ //Check if MPI balance has been finished
	  //Cancel receive request
	  requests[0].cancel();
	  break; //Exit from communication loop
	}
	//Get MPI lock
	const std::lock_guard<std::mutex> lock(mtxMPI);
	//Change sent flag state if required
	if(MPIsent)
	  MPIsent = false;
      }
    }
  }

  if(verbose > 1){
    fprintf(fMPI,"%07ld s - Finishing MPI balancing.\n",
	    static_cast<long int>(timeStamp()));
    fflush(fMPI);
  }
}

int LB::task::init(const size_t nw,
		   const unsigned long long nIter,
		   const MPI_Comm commIn,
		   const int tagReqests,
		   const int tagProcess,
		   const char* logFileName,
		   const unsigned verbose){

  //Check log filename
  if(logFileName == nullptr){
    printf("LB:task:init: Error: Invalid log files.\n");
    return LB_ERROR_NULL_POINTER;
  }
  
  int err;
  err = setWorkers(nw,verbose);
  if(err != LB_SUCCESS){
    if(verbose > 0)
      printf("LB:task:init: Error creating workers. Error code: %d\n",err);
    return err;
  }
  
  err = setMPI(commIn,tagReqests,tagProcess,nIter,verbose);
  if(err != LB_SUCCESS){
    if(verbose > 0)
      printf("LB:task:init: Error at MPI parameters "
	     "initialization. Error code: %d\n",err);
    clear();
    clearMPI();
    return err;
  }

  //Create log files
  const std::lock_guard<std::mutex> lockTh(mtxfth);      
  int rank;
  MPI_Comm_rank(comm, &rank);
  char prefix[20];
  sprintf(prefix,"rank-%03d-",rank);
  std::string filenameTh(prefix);
  filenameTh += std::string(logFileName) + std::string("-local-LB.log");
  fth = fopen(filenameTh.c_str(),"w");
  std::string filenameMPI(prefix);
  filenameMPI += std::string(logFileName) + std::string("-MPI-LB.log");
  fMPI = fopen(filenameMPI.c_str(),"w");

  if(verbose > 0){
    if(fth == nullptr)
      printf("LB:task:init: Error: Unable to create local log file: %s\n",
	     filenameTh.c_str());
    if(fMPI == nullptr)
      printf("LB:task:init: Error: Unable to create MPI log file: %s\n",
	     filenameMPI.c_str());
  }

  if(fth == nullptr || fMPI == nullptr){
    clear();
    clearMPI();
    return LB_ERROR_INVALID_FILE;
  }
  
  return LB_SUCCESS;
}

int LB::task::printReportMPI(FILE* fout) const{

  if(fout == nullptr)
    return LB_ERROR_NULL_POINTER;

  int rank;
  MPI_Comm_rank(comm, &rank);

  if(rank != 0)
    return LB_ERROR_MPI_RANK_0_EXPECTED;
  
  fprintf(fout,"# Task load balance MPI report.\n");
  unsigned iw = 0;
  for(const auto& worker: MPIworkers){
    fprintf(fout,"# Worker %03u\n",iw);
    worker.printReport(fout);
    fprintf(fout,"\n\n");
    ++iw;
  }
  return LB_SUCCESS;
}

#endif
// **************************** MPI ********************************* //

int LB::task::init(const size_t nw,
		   const unsigned long long nIter,
		   const char* logFileName,
		   const unsigned verbose){

  //Check log filename
  if(logFileName == nullptr){
    printf("LB:task:init: Error: Invalid log files.\n");
    return LB_ERROR_NULL_POINTER;
  }
  
  int err;
  err = setWorkers(nw);
  if(err != LB_SUCCESS){
    if(verbose > 0)
      printf("LB:task:init: Error creating workers. Error code: %d\n",err);    
    return err;
  }
  
  err = setIterations(nIter);
  if(err != LB_SUCCESS){
    clear();
    if(verbose > 0)
      printf("LB:task:init: Error setting iterations. Error code: %d\n",err);
    return err;
  }

  //Create log files
  const std::lock_guard<std::mutex> lockTh(mtxfth);      
  std::string filenameTh(logFileName);
  filenameTh += std::string("-local-LB.log");
  fth = fopen(filenameTh.c_str(),"w");

  if(verbose > 0){
    if(fth == nullptr)
      printf("LB:task:init: Error: Unable to create local log file: %s\n",
	     filenameTh.c_str());
  }

  if(fth == nullptr){
    clear();
    return LB_ERROR_INVALID_FILE;
  }
  
  return LB_SUCCESS;
}

int LB::task::printReport(FILE* fout) const{

  if(fout == nullptr)
    return LB_ERROR_NULL_POINTER;

  fprintf(fout,"# Task load balance report.\n");
  unsigned iw = 0;
  for(const auto& worker: workers){
    fprintf(fout,"# Worker %03u\n",iw);
    worker.printReport(fout);
    fprintf(fout,"\n\n");
    ++iw;
  }
  return LB_SUCCESS;
}

 
void LB::task::reset(){
  clear();
// **************************** MPI ********************************* //
#ifdef _PEN_USE_MPI_
  clearMPI();
#endif
// **************************** MPI ********************************* //
}

// ********************** HTTP SUPPORT **********************
#ifdef _PEN_USE_LB_HTTP_
int LB::task::extHTTPserver(const unsigned extern_iw,
			    const char* url,
			    const unsigned verbose){

  FILE* flog = fth;
  int rank = 0;
#ifdef _PEN_USE_MPI_
  MPI_Comm_rank(comm, &rank);
  flog = fMPI;
#endif

  if(flog == nullptr){
    if(verbose > 0){
      printf("LB::task::initExternalBalance: Error: "
	     "Task has not initialised.\n");
    }
    return LB_ERROR_TASK_NOT_INITIALISED;
  }

  if(rank != 0){
    if(verbose > 0){
      fprintf(flog,"LB::task::initExternalBalance: Error: \n"
	      "Externam balance must be handled by rank 0.\n");
      fflush(flog);
    }
    return LB_ERROR_HANDLED_BY_RANK0;
  }  

  //Init CURL
  curl_global_init(CURL_GLOBAL_DEFAULT);  
  
  //Save extern worker number
  ext_iworker = extern_iw;
  //Enable external http balance
  ext_balance = true;
  ext_http = true;

  //Save base url
  baseURL.assign(url);

  //Create methods urls
  startURL.assign(baseURL);
  startURL.append("/start");

  reportURL.assign(baseURL);
  reportURL.append("/report");

  finishURL.assign(baseURL);
  finishURL.append("/finish");
  
  return LB_SUCCESS;  
}
#endif
// **********************   HTTP END   **********************

int LB::task::extLBserver(const unsigned extern_iw,
			  const char* host,
			  const char* port,
			  const char* CAfilename,
			  const char* certFilename,
			  const char* keyFilename,
			  const char* password,
			  const char* hostname,
			  const unsigned verbose){

  FILE* flog = fth;
  int rank = 0;
#ifdef _PEN_USE_MPI_
  MPI_Comm_rank(comm, &rank);
  flog = fMPI;
#endif

  if(flog == nullptr){
    if(verbose > 0){
      printf("LB::task::initExternalBalance: Error: "
	     "Task has not initialised.\n");
    }
    return LB_ERROR_TASK_NOT_INITIALISED;
  }
  
  if(rank != 0){
    if(verbose > 0){
      fprintf(flog,"LB::task::initExternalBalance: Error: \n"
	      "Externam balance must be handled by rank 0.\n");
      fflush(flog);
    }
    return LB_ERROR_HANDLED_BY_RANK0;
  }

  //Set log file
  client.setLogFile(flog);
  //Set verbose
  client.setVerbose(verbose);
  
  if(CAfilename != nullptr ||
     certFilename != nullptr ||
     keyFilename != nullptr ||
     password != nullptr ||
     hostname != nullptr){

    
#ifdef _PEN_USE_SSL_    
    // **** SSL enabled **** //
    int err = client.initSSL(CAfilename,certFilename,
			     keyFilename,password,hostname);
    if(err != PEN_TCP_SUCCESS){
      if(verbose > 0){
	fprintf(flog,"Error configuring SSL connection.\n"
		"       TCP Error code: %d\n",err);
	fflush(flog);
      }
      return LB_ERROR_CONFIGURING_CONNECTION;
    }
#else 
    // **** SSL disabled **** //
    if(verbose > 0){
      fprintf(flog,"Error configuring SSL connection.\n"
	            "SSL support has not enabled on compilation.\n"
	            "Recompile the code enabling SSL capabilities.\n");
    }
    return LB_SSL_DISABLED;
#endif
  }

  //Set host and port
  int err = client.setHost(host,port);
  if(err != PEN_TCP_SUCCESS){
    if(verbose > 0){
      fprintf(flog,"Error configuring host and port.\n"
	      "       TCP Error code: %d\n",err);
      fflush(flog);
    }
    return LB_ERROR_CONFIGURING_CONNECTION;
  }

  //Save extern worker number
  ext_iworker = extern_iw;
  //Enable external balance
  ext_balance = true;
  return LB_SUCCESS;
}
 
//********************
//**  Task server
//********************

int LB::taskServer::setLogFile(const char* logFileName,
			       const unsigned verbose){

  if(flog != nullptr)
    fclose(flog);
  flog = nullptr;

  if(logFileName != nullptr){
    flog = fopen(logFileName,"w");
    if(flog == nullptr){
      if(verbose > 0)
	printf("LB:taskServer:setLogFile: Error: Unable to create log "
	       "file: %s\n",logFileName);
      return LB_ERROR_UNABLE_TO_CREATE_FILE;
    }
  }

  return LB_SUCCESS;
}
  
int LB::taskServer::init(const size_t nw,
			 const unsigned long long nIter,
			 const char* logFileName,
			 const unsigned verbose){

  if(filterLog(verbose,1)){
    fprintf(flog,"%07ld s - Initialize task with %lu workers "
	    "and %llu iterations\n",
	    static_cast<long int>(timeStamp()),
	    static_cast<unsigned long>(nw),nIter);
    fflush(flog);
  }
  
  //Clear previous configuration
  clear();
      
  //Create workers
  workers.resize(nw);
  workerInits.resize(nw);
  workerResumes.resize(nw);
  
  //Save number of iterations
  iterations = nIter;
  remaining = nIter;
  
  if(logFileName != nullptr){
    return setLogFile(logFileName,verbose);
  }

  return LB_SUCCESS;
}

int LB::taskServer::workerStart(const size_t iw,
				const std::chrono::steady_clock::time_point& t,
				const unsigned verbose){

   if(filterLog(verbose,1)){
     fprintf(flog,"%07ld s - Worker %lu starting.\n",
	     static_cast<long int>(timeStamp()),
	     static_cast<unsigned long>(iw));
     fflush(flog);
   }
      
   if(iw >= workers.size()){
     //New worker added, resize vector of workers
     workers.resize(iw+1);
     workerInits.resize(iw+1);
     workerResumes.resize(iw+1);
     if(filterLog(verbose,1)){
       fprintf(flog,"%07ld s - Adding workers until worker %lu.\n",
	       static_cast<long int>(timeStamp()),
	       static_cast<unsigned long>(iw));
       fflush(flog);
     }
   }

   if(!started){
     started = true;
     finished = false;
   }

   //Get worker
   guessWorker& taskWorker = workers[iw];
      
   //Calculate iterations to assign
   unsigned long long toDo = remaining/workers.size();

   if(filterLog(verbose,1)){
     fprintf(flog,"   Initial iteration assign: %llu\n",toDo);
     fflush(flog);
   }
   
   //Start the worker
   int ret = taskWorker.start(t,toDo);
   if(ret == LB_SUCCESS){
     //Save worker init
     workerInits[iw] = t;
     //Save worker first "resume"
     workerResumes[iw] = t;
   } else if(ret == LB_ERROR_WORKER_ALREADY_STARTED){
     //Worker already started. We interpret this as a
     //crashed worker resuming its execution
     //So, add a dummy measure with the same iterations done
     //to avoid erroneous elapsed time measures
     taskWorker.addMeasure(t,taskWorker.done());
     //Ensure that the time is updated even if "addMeasure"
     //doesn't creates a new measure
     taskWorker.setLastTime(t);
     //Save this timepoint as the new worker resume point
     workerResumes[iw] = t;
     if(filterLog(verbose,1)){
       fprintf(flog,"   Resuming crashed worker\n");
       fflush(flog);
     }     
     return LB_SUCCESS;
   }
   return ret;
}

int LB::taskServer::workerFinish(const size_t iw,
				 const std::chrono::steady_clock::time_point& t,
				 const unsigned long long nIter,
				 const unsigned verbose){

   if(filterLog(verbose,1)){
     fprintf(flog,"%07ld s - Worker %lu finishing.\n",
	     static_cast<long int>(timeStamp()),
	     static_cast<unsigned long>(iw));
     fflush(flog);
   }

   //Perform a report without balance
   int err = report(iw,t,nIter,verbose,false);
   if(err == LB_SUCCESS){
     
     //Finish worker
     workers[iw].end();

     //Check if some worker is working
     for(const auto& worker : workers)
       if(!worker.finishes()){
	 //Some worker is workning, balance remaining workers
	 balance(verbose);
       }
   }
   else if(err == LB_ERROR_TASK_ALREADY_FINISHED){
     //Finish worker
     workers[iw].end();     
   }
   return err;
}
 
int LB::taskServer::report(const size_t iw,
			   const std::chrono::steady_clock::time_point& t,
			   const unsigned long long nIter,
			   const unsigned verbose,
			   const bool doBalance){
  //Updates worker information according to the number of
  //iterations done and the last report

  //We suppose that this function is not called before the
  //task start (before the first call to "workerStart").
  //With this assumption, we only need to lock the worker
  //specified by "iw".

  if(iw >= workers.size()){
    if(filterLog(verbose,0)){
      fprintf(flog,"%07ld s - Error: Worker index (%lu) out of range.\n",
	      static_cast<long int>(timeStamp()),
	      static_cast<unsigned long>(iw));
      fflush(flog);
    }
    return LB_ERROR_WORKER_OUT_OF_RANGE;
  }

  //Get worker reference
  guessWorker& taskWorker = workers[iw];
  
  if(filterLog(verbose,1)){
    fprintf(flog,"%07ld s - Reporting %llu iterations done by "
	    "worker %lu at local timestamp %07ld s.\n",
	    static_cast<long int>(timeStamp()),
	    nIter,static_cast<unsigned long>(iw),
	    static_cast<long int>(taskWorker.timeStamp(t)));
    fflush(flog);
  }

  //Check if this worker is doing this task
  if(!taskWorker.begins()){
    if(filterLog(verbose,0)){
      fprintf(flog,"%07ld s - Error: Worker %lu doesn't started.\n",
	      static_cast<long int>(timeStamp()),iw);
      fflush(flog);
    }
    return LB_ERROR_WORKER_NOT_STARTED;
  }

  if(taskWorker.finishes()){
    if(filterLog(verbose,0)){
      fprintf(flog,"%07ld s - Error: Worker %lu has finished.\n",
	      static_cast<long int>(timeStamp()),iw);
      fflush(flog);
    }
    return LB_ERROR_WORKER_FINISHED;
  }
      
  float dev = taskWorker.addMeasure(t,nIter);

  if(filterLog(verbose,0)){
    fprintf(flog,"%07ld s - Deviation since previous report: %.2f%%\n",
	    static_cast<long int>(timeStamp()),(dev-1.0)*100.0);	      
    fflush(flog);
  }
  
  if(finished){
    if(filterLog(verbose,0)){
      fprintf(flog,"%07ld s - Warning: Task already finished, "
	      "no balance will be performed.\n",
	      static_cast<long int>(timeStamp()));
      fflush(flog);
    }
    return LB_ERROR_TASK_ALREADY_FINISHED;
  }  
  if(doBalance && taskWorker.registers() > 2){
    //Balance workers
    balance(verbose);
  }
  
  return LB_SUCCESS;
}

void LB::taskServer::balance(const unsigned verbose){

  if(filterLog(verbose,1)){
    fprintf(flog,"%07ld s - Balancing iteration assignation\n",
	    static_cast<long int>(timeStamp()));	        
    fflush(flog);
  }
  
  if(!started){
    if(filterLog(verbose,1)){
      fprintf(flog,"          Task has not started\n");
      fflush(flog);
    }
    remainingTime = 1000000000000000ull;
    return;
  }

  if(finished){
    if(filterLog(verbose,1)){
      fprintf(flog,"          Task has finished\n");
      fflush(flog);
    }
    remainingTime = 0ull;
    return;
  }

  const auto now = std::chrono::steady_clock::now();
   
  float totalSpeed;
  unsigned long long totalIterDone = speed(totalSpeed,now);

  if(filterLog(verbose,1)){
    fprintf(flog,"          Iterations done: %llu\n"
	         "                    Speed: %12.5E\n",
	    totalIterDone,totalSpeed);
    fflush(flog);
  }
      
  //Check if there are something to do
  if(totalIterDone >= iterations){
    for(guessWorker& worker : workers){
      if(worker.begins() && !worker.finishes()){
	worker.assigned = worker.done();
      }
    }
    finished = true;
    remaining = 0ull;
    remainingTime = 0ull;
  }else{
    //Calculate remaining iterations
    remaining = iterations - totalIterDone;

    //Check global speed 
    if(totalSpeed <= 0.0){
      remainingTime = 1000000000000000ull;
      return;
    }

    unsigned long long pred = predDone();
    if(pred >= iterations)
      remainingTime = 0llu;
    else
      remainingTime = static_cast<double>(iterations-pred)/totalSpeed;
    if(filterLog(verbose,1)){
      fprintf(flog,"Predicted iterations done: %llu ETA: %llu s\n",
	      pred,remainingTime);
      fflush(flog);
    }

    //If remaining time is lesser than half threshold avoid rebalance
    if(remainingTime < threshold/2){
      return; //Avoid rebalance
    }
  
    //Assign iterations to each worker
    unsigned long long nworkers = workers.size();
    unsigned long cont = 0;
    for(guessWorker& worker : workers){
      if(worker.begins() && !worker.finishes()){
	//Check if the last report of this worker has been done recently.
	//If the server has not received reports from this worker
	//since "threshold", it will be considered as stopped
	//or crashed worker.
	unsigned long long telaps = worker.elapsed(now);
	if(telaps > threshold){
	  if(filterLog(verbose,1))
	    fprintf(flog,"    Worker %3lu is missing  -> done : %llu\n",
		    cont,worker.done());
	   
	  ++cont;
	  continue;
	}
	unsigned long long toDo = remaining*(worker.speed()/totalSpeed);
	toDo = std::max(toDo,nworkers);
	worker.assigned = worker.done() + toDo;
	 
	if(filterLog(verbose,1))
	  fprintf(flog,"    Worker %3lu is active   -> done : %llu/%llu\n",
		  cont,worker.done(),worker.assigned);
	 
      }else{
	if(filterLog(verbose,1))
	  fprintf(flog,"    Worker %3lu is inactive -> done : %llu\n",
		  cont,worker.done());
      }
      ++cont;
    }
    if(filterLog(verbose,1))
      fflush(flog);
  }
}

void LB::taskServer::clear(){
  iterations = remaining = 0llu;
  threshold = 300llu;
  started = finished = false;

  remainingTime = 100000000000;
  workers.clear();
  workerInits.clear();
  workerResumes.clear();
}
 
int LB::taskServer::load(FILE* fin,
			 const unsigned verbose){
  if(fin == nullptr)
    return LB_ERROR_NULL_POINTER;

  clear();

  unsigned long nworkers;
  int istarted, ifinished;
  fscanf(fin," %llu %llu %llu %llu %lu %d %d ",
	 &iterations, &remaining, &remainingTime,
	 &threshold,&nworkers,&istarted,&ifinished);

  started = istarted;
  finished = ifinished;
  
  //Allocate memory for each worker
  workers.resize(nworkers);
  workerInits.resize(nworkers);
  workerResumes.resize(nworkers);
  
  //Load worker inits
  for(auto& wInit : workerInits){
    long int aux;
    if(fscanf(fin,"%ld",&aux) != 1){
      return LB_ERROR_UNEXPECTED_FORMAT;
    }
    wInit = taskBegin + std::chrono::seconds(aux);
  }

  //Load worker resumes
  for(auto& wResume : workerResumes){
    long int aux;
    if(fscanf(fin,"%ld",&aux) != 1){
      return LB_ERROR_UNEXPECTED_FORMAT;
    }
    wResume = taskBegin + std::chrono::seconds(aux);
  }
  
  //Load worker information
  size_t cont = 0;
  for(auto& worker : workers){
    int err = worker.load(fin,workerInits[cont]);
    if(err != LB_SUCCESS){
      if(filterLog(verbose,0))
	fprintf(flog,"%07ld s - Unable to load worker %lu.\n",
		static_cast<long int>(timeStamp()),
		static_cast<long unsigned>(cont));
      worker.reset();
      fflush(flog);
    }
    ++cont;
  }
  return LB_SUCCESS;
}
 
int LB::taskServer::dump(const char* saveFilename,
			 const unsigned verbose){
  if(filterLog(verbose,1)){
    fprintf(flog,"%07ld s - Dumping recorded data to %s\n",
	    static_cast<long int>(timeStamp()),saveFilename);
    fflush(flog);
  }

  if(saveFilename == nullptr){
    if(filterLog(verbose,0)){
      fprintf(flog,"      Error: specified save filename is a null pointer\n");
      fflush(flog);
    }
    return LB_ERROR_NULL_POINTER;
  }

  FILE* fdump = nullptr;
  fdump = fopen(saveFilename,"w");
  if(fdump == nullptr){
    if(filterLog(verbose,0)){
      fprintf(flog,"     Error: Unable to create dump file '%s'\n",saveFilename);
      fflush(flog);
    }
    return LB_ERROR_UNABLE_TO_CREATE_FILE;
  }

  fprintf(fdump,"%llu %llu %llu %llu %lu %d %d\n",
	 iterations, remaining, remainingTime,
	 threshold,static_cast<unsigned long>(workers.size()),
	 started,finished);

  //Dump worker inits
  for(const auto& wInit : workerInits){
    fprintf(fdump,"%ld\n",static_cast<long int>(timeStamp(wInit)));
  }

  //Dump worker resumes
  for(const auto& wResume : workerResumes){
    fprintf(fdump,"%ld\n",static_cast<long int>(timeStamp(wResume)));
  }
  
  //Dump worker information
  for(const auto& worker : workers){
    worker.dump(fdump);
  }
      
  return LB_SUCCESS;
}


int LB::taskServer::receiveReport(const size_t iw,
				  const unsigned long long nIter,
				  const long int elapsed,
				  unsigned long long& newAssign,
				  const unsigned verbose){

  if(iw >= workers.size()){
    if(filterLog(verbose,0)){
      fprintf(flog,"%07ld s - Error: Worker index (%lu) out of range.\n",
	      static_cast<long int>(timeStamp()),iw);
      fflush(flog);
    }
    return LB_ERROR_WORKER_OUT_OF_RANGE;
  }
  else{
    //Calculate worker report steady point 
    auto reportTime = workerTP(iw,std::chrono::seconds(elapsed));

    //If the worker reports a elapsed time lesser than previous one
    //we require a new start petition to ajust the time. This could
    //happens on resumed crashed workers.
    if(workers[iw].elapsed(reportTime) <= 0){
      return LB_ERROR_WORKER_NOT_STARTED;
    }
    
    //Perform the report
    int ret = report(iw,reportTime,nIter,verbose);
    newAssign = workers[iw].assigned;
    return ret;
  }
}

int LB::taskServer::receiveStart(const size_t iw,
				 const long int elapsed,
				 unsigned long long& newAssign,
				 const unsigned verbose){
  //Calculate worker start steady point 
  auto workerInit =
    std::chrono::steady_clock::now() - std::chrono::seconds(elapsed);
	  
  //Start worker
  int ret = workerStart(iw,workerInit,verbose);
  if(ret == LB_SUCCESS)
    newAssign = workers[iw].assigned;
  return ret;
}

int LB::taskServer::receiveFinish(const size_t iw,
				  const unsigned long long nIter,
				  const long int elapsed,
				  const unsigned verbose){
  if(iw >= workers.size()){
    if(filterLog(verbose,0)){
      fprintf(flog,"%07ld s - Error: Worker index (%lu) out of range.\n",
	      static_cast<long int>(timeStamp()),iw);
      fflush(flog);
    }
    return LB_ERROR_WORKER_OUT_OF_RANGE;
  }
  else{
    //Calculate worker report steady point 
    auto finishTime = workerTP(iw,std::chrono::seconds(elapsed));    
    return workerFinish(iw,finishTime,nIter,verbose);
  }
}

 
void LB::taskServer::monitor(const bool local,
			     const unsigned short port,
			     const char* CAfilename,
			     const char* certFilename,
			     const char* keyFilename,
			     const char* dhFilename,
			     const char* password,
			     const unsigned verbose){

  
  //Create a sequential TCP server
  pen_tcp::server_seq server;  
  server.setLogFile(flog);
  server.setVerbose(verbose);
  server.setInit(taskBegin);

#ifdef _PEN_USE_SSL_    
  // **** SSL enabled **** //
  int errInit = server.init(local,port,
			    CAfilename,certFilename,keyFilename,
			    dhFilename,password);
  if(errInit != PEN_TCP_SUCCESS){
    if(filterLog(verbose,0)){
      fprintf(flog,"%07ld s - Error on server initialisation\n",
	      static_cast<long int>(timeStamp()));
      fflush(flog);
    }
    return;
  }
#else
  // **** SSL disabled **** //
  if(CAfilename != nullptr || certFilename != nullptr ||
     keyFilename != nullptr || dhFilename != nullptr ||
     password != nullptr){
    if(verbose > 0){
      fprintf(flog,"%07ld s - Error configuring SSL connection.\n"
	      "SSL support has not enabled on compilation.\n"
	      "Recompile the code enabling SSL capabilities.\n",
	      static_cast<long int>(timeStamp()));
    }
    return;
  }

  int errInit = server.init(local,port);
  if(errInit != PEN_TCP_SUCCESS){
    if(filterLog(verbose,0)){
      fprintf(flog,"%07ld s - Error on server initialisation\n",
	      static_cast<long int>(timeStamp()));
      fflush(flog);
    }
    return;
  }
#endif
  
  for(;;){

    //Wait until next incoming connection
    std::string message;
    int errList = server.listen();
    if(errList != PEN_TCP_SUCCESS){
      if(errList == PEN_TCP_ACCEPT_CONNECTION_FAIL)
	break; //Unable to accept connections

      //Close connection
      server.close();
      continue; //Possible error on handshake, skip this connection
    }
    int errRecv = server.receive(message);
    if(errRecv != PEN_TCP_SUCCESS){
      //Close connection
      server.close();      
      continue; //Error on read, skip this connection
    }
    //Message successfully read, extract the instruction to perform
    if(filterLog(verbose,1)){
      fprintf(flog,"%07ld s - Received message:\n  %s\n",
	      static_cast<long int>(timeStamp()),message.c_str());
      fflush(flog);
    }
	
    int instruction;
    if(sscanf(message.c_str(),"%d",&instruction) != 1){
      //Unable to read instruction
      if(filterLog(verbose,0)){
	fprintf(flog,"%07ld s - Unable to read instruction\n",
		static_cast<long int>(timeStamp()));
	fflush(flog);
      }
      //Send error message
      char errMessage[pen_tcp::messageLen];
      createError(1,LB_ERROR_UNEXPECTED_FORMAT,
		  "Unable to read instruction\n",errMessage,verbose);
      server.write(errMessage);
    }else{
      //Check the instruction
      long unsigned iw;
      unsigned long long nIter;
      long int elapsed;
      switch(instruction){
      case 1:{ //Report

	//Parse the instruction
	int nread = sscanf(message.c_str(),"%*d %lu %llu %ld",
			    &iw, &nIter, &elapsed);
	if(nread != 3){
	  //Unable to parse report instruction
	  if(filterLog(verbose,0)){
	    fprintf(flog,"%07ld s - Unable to parse report\n",
		    static_cast<long int>(timeStamp()));
	    fflush(flog);
	  }
	  //Send error message
	  char errMessage[pen_tcp::messageLen];
	  createError(1,LB_ERROR_UNEXPECTED_FORMAT,
		      "Unable to parse report\n",errMessage,verbose);
	  server.write(errMessage);
	} else{
	  //Perform the report
	  unsigned long long newAssign = 0;
	  int err = receiveReport(static_cast<size_t>(iw),nIter,
				  elapsed,newAssign,verbose);

	  //Send the response
	  if(err == LB_SUCCESS){
	    //Send new assignation
	    char response[pen_tcp::messageLen];
	    snprintf(response,pen_tcp::messageLen,
		     "0 \n Assigned: %llu\n ETA: %llu\n",
		     newAssign,ETA());
	    server.write(response);
	  }
	  else{
	    //Send error message
	    char errMessage[pen_tcp::messageLen];
	    createError(1,err,nullptr,errMessage,verbose);
	    server.write(errMessage);
	  }
	}
      }
	break;
      case 2:{ //Start petition

	//Parse the instruction
	int nread = sscanf(message.c_str(),"%*d %lu %ld",&iw,&elapsed);
	if(nread != 2){
	  //Unable to parse start petition
	  if(filterLog(verbose,0)){
	    fprintf(flog,"%07ld s - Unable to parse start request\n",
		    static_cast<long int>(timeStamp()));
	    fflush(flog);
	  }
	  //Send error message
	  char errMessage[pen_tcp::messageLen];
	  createError(2,LB_ERROR_UNEXPECTED_FORMAT,
		      "Unable to parse start request\n",errMessage,verbose);
	  server.write(errMessage);
	} else{
	  //Perform the worker start
	  unsigned long long newAssign = 0;
	  int err = receiveStart(static_cast<size_t>(iw),elapsed,
				 newAssign,verbose);

	  //Send the response
	  if(err == LB_SUCCESS){
	    //Send new assignation
	    char response[pen_tcp::messageLen];
	    snprintf(response,pen_tcp::messageLen,
		     "0 \n Assigned: %llu\n ETA: %llu\n",
		     newAssign,ETA());
	    server.write(response);
	  }
	  else{
	    char errMessage[pen_tcp::messageLen];
	    createError(2,err,nullptr,errMessage,verbose);
	    server.write(errMessage);
	  }
	}
      }
	break;
      case 3:{ //Finish petition

	//Parse the instruction
	int nread = sscanf(message.c_str(),"%*d %lu %llu %ld",
			   &iw,&nIter,&elapsed);
	if(nread != 3){
	  //Unable to parse finish instruction
	  if(filterLog(verbose,0)){
	    fprintf(flog,"%07ld s - Unable to parse finish request\n",
		    static_cast<long int>(timeStamp()));
	    fflush(flog);
	  }
	  //Send error message
	  char errMessage[pen_tcp::messageLen];
	  createError(3,LB_ERROR_UNEXPECTED_FORMAT,
		      "Unable to parse finish request\n",
		      errMessage,verbose);
	  server.write(errMessage);
	} else{
	  //Perform the finish request
	  int err = receiveFinish(static_cast<size_t>(iw),nIter,
				  elapsed,verbose);
	  //Send the response
	  if(err == LB_SUCCESS){
	    //Send new assignation
	    const char* res = "0\n Finish registered\n";
	    char response[pen_tcp::messageLen];
	    snprintf(response,pen_tcp::messageLen,"%s",res);
	    server.write(response);
	  }
	  else{
	    char errMessage[pen_tcp::messageLen];
	    createError(3,err,nullptr,errMessage,verbose);
	    server.write(errMessage);	    
	  }
	}
      }
	break;
      case 4:{ //Get worker status 
	if(sscanf(message.c_str(),"%*d %lu",&iw) != 1){
	  //Unable to parse petition
	  if(filterLog(verbose,0)){
	    fprintf(flog,"%07ld s - Unable to parse worker status request\n",
		    static_cast<long int>(timeStamp()));
	    fflush(flog);
	  }
	  //Send error message
	  char errMessage[pen_tcp::messageLen];
	  createError(4,LB_ERROR_UNEXPECTED_FORMAT,
		      "Unable to parse worker status request\n",
		      errMessage,verbose);
	  server.write(errMessage);	    
	} else{
	  if(iw < workers.size()){
	    //Send worker status
	    int status;
	    if(workers[iw].begins())
	      if(workers[iw].finishes())
		status = 2;
	      else
		status = 1;
	    else status = 0;
	    char response[pen_tcp::messageLen];
	    snprintf(response,pen_tcp::messageLen,"0 \n"
		     "status: %d\n"
		     "speed: %.2f\n"
		     "done: %llu\n"
		     "assigned: %llu\n",
		     status,workers[iw].speed(),
		     workers[iw].done(),workers[iw].assigned);
	    server.write(response);
	  }
	  else{
	    char errMessage[pen_tcp::messageLen];
	    createError(4,LB_ERROR_WORKER_OUT_OF_RANGE,nullptr,
			errMessage,verbose);
	    server.write(errMessage);	    
	  }
	}
      }
	break;
      case 5:{ //Get worker iteration assign
	if(sscanf(message.c_str(),"%*d %lu",&iw) != 1){
	  //Unable to parse petition
	  if(filterLog(verbose,0)){
	    fprintf(flog,"%07ld s - Unable to parse worker "
		    "assignation request\n",
		    static_cast<long int>(timeStamp()));
	    fflush(flog);
	  }
	  //Send error message
	  char errMessage[pen_tcp::messageLen];
	  createError(5,LB_ERROR_UNEXPECTED_FORMAT,
		      "Unable to parse worker assignation request\n",
		      errMessage,verbose);
	  server.write(errMessage);	    
	} else{
	  if(iw < workers.size()){
	    //Send worker assigned iterations
	    char response[pen_tcp::messageLen];
	    snprintf(response,pen_tcp::messageLen,"0 \n assigned: %llu\n",
		     workers[iw].assigned);
	    server.write(response);
	  }
	  else{
	    char errMessage[pen_tcp::messageLen];
	    createError(5,LB_ERROR_WORKER_OUT_OF_RANGE,nullptr,
			errMessage,verbose);
	    server.write(errMessage);	    
	  }
	}
      }
	break;
      case 6:{ //Dump information
	char filename[100];
	long int dt = static_cast<long int>(timeStamp());
	snprintf(filename,100,"dump-%lds.dat",dt);
	int err = dump(filename,verbose);
	if(err == LB_SUCCESS){
	  char response[pen_tcp::messageLen];
	  snprintf(response,pen_tcp::messageLen,"0 \n Dump done\n");
	  server.write(response);	  
	}
	else{
	  //Send error message
	  char errMessage[pen_tcp::messageLen];
	  createError(6,err,"Unable to create dump\n",
		      errMessage,verbose);
	  server.write(errMessage);	    	  	  
	}
      }
	break;
      case 7:{ //Print report information
	char filename[100];
	long int dt = static_cast<long int>(timeStamp());
	snprintf(filename,100,"report-%lds",dt);
	if(filterLog(verbose,1)){
	  fprintf(flog,"%07ld s - Print report to: %s\n",
		  static_cast<long int>(timeStamp()),filename);
	  fflush(flog);
	}
	FILE* frep = nullptr;
	frep = fopen(filename,"w");
	if(frep == nullptr){
	  if(filterLog(verbose,0)){
	    fprintf(flog,"%07ld s - Error: Unable to create file '%s' ",
		    static_cast<long int>(timeStamp()),filename);
	    fflush(flog);
	  }
	  //Send error message
	  char errMessage[pen_tcp::messageLen];
	  createError(7,LB_ERROR_UNABLE_TO_CREATE_FILE,
		      "Unable to create report file\n",
		      errMessage,verbose);
	  server.write(errMessage);	    	  
	}
	else{
	  char response[pen_tcp::messageLen];
	  snprintf(response,pen_tcp::messageLen,"0 \n Report done\n");
	  server.write(response);
	  printReport(frep);
	  fclose(frep);
	}
      }
	break;
      case 8:{ //Send server ID
	server.write(ID);
      }
	break;
      case 9:{ //Get global execution status
	const auto now = std::chrono::steady_clock::now();
	unsigned long wstarted, wfinished, wlost;
	wstarted = wfinished = wlost = 0;
	for(const LB::guessWorker& worker : workers){
	  if(workers[iw].begins()){
	    ++wstarted; //This worker has started
	    if(workers[iw].finishes()){
	      ++wfinished; //This worker has finished
	    }
	    else{
	      unsigned long long telaps = worker.elapsed(now);
	      if(telaps >= threshold)
		++wlost; //This worker is lost
	    }
	  }
	}
	//Send execution status
	char response[pen_tcp::messageLen];
	snprintf(response,pen_tcp::messageLen,"%lu %lu %lu %lu\n",
		 static_cast<unsigned long>(workers.size()),
		 wstarted,wfinished,wlost);
	server.write(response);
      }
	break;
      default:
	//Send error message
	char errMessage[pen_tcp::messageLen];
	createError(-1,LB_ERROR_UNEXPECTED_FORMAT,
		    "Invalid request\n",errMessage,verbose);
	server.write(errMessage);	    
      }
    }

    //Close connection
    server.close();
  }
}

void LB::taskServer::createError(const int prefix, const int errcode,
				 const char* errmessage,
				 char err[pen_tcp::messageLen],
				 const unsigned verbose){

  //Construct error message
  if(errmessage == nullptr){
    snprintf(err,pen_tcp::messageLen,"%d\n%d\n",prefix,errcode);
  }
  else{
    snprintf(err,pen_tcp::messageLen,"%d\n%d\n%s\n",
	     prefix,errcode,errmessage);
  }

  if(filterLog(verbose,2)){
    fprintf(flog,"%07ld s - Error message\n"
	    "                  Prefix: %d\n"
	    "              Error code: %d\n"
	    "           Error message: %s\n",
	    static_cast<long int>(timeStamp()),
	    prefix,errcode,
	    errmessage == nullptr ? "none" : errmessage);
    fflush(flog);
  }
}

int LB::taskServer::printReport(FILE* fout) const{

  if(fout == nullptr)
    return LB_ERROR_NULL_POINTER;

  fprintf(fout,"# Task load balance report.\n");
  unsigned iw = 0;
  for(const auto& worker: workers){
    fprintf(fout,"# Worker %03u\n",iw);
    worker.printReport(fout,taskBegin);
    fprintf(fout,"\n\n");
    ++iw;
  }
  return LB_SUCCESS;
}

LB::taskServer::~taskServer(){
  if(flog != nullptr)
    fclose(flog);
  flog = nullptr;
  clear();
}
