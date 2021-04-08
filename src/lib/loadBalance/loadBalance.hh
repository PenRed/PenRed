//
//
//    Copyright (C) 2020 Universitat Politècnica de València - UPV
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

#ifndef __PEN_LOAD_BALANCE__
#define __PEN_LOAD_BALANCE__

#include <chrono>
#include <ctime>
#include <vector>
#include <mutex>
#include <stdexcept>
#include <thread>
#include <cmath>
#include "tcp_cs.hh"

// ******************************* MPI ************************************ //
#ifdef _PEN_USE_MPI_
#include "mpi.h"
#endif
// ***************************** MPI END ********************************** //

// ********************** HTTP SUPPORT **********************
#ifdef _PEN_USE_LB_HTTP_
#include <curl/curl.h>
#endif
// **********************   HTTP END   **********************


enum LB_ERRCODE{
		LB_SUCCESS = 0,                 //0
		LB_ERROR_MPI_NOT_INITIALIZED,   //1
		LB_ERROR_THREAD_OUT_OF_RANGE,   //2
		LB_ERROR_WORKER_OUT_OF_RANGE,   //3
		LB_ERROR_LESS_ITERATIONS_DONE,  //4
		LB_ERROR_WORKER_NOT_STARTED,    //5
		LB_ERROR_WORKER_ALREADY_STARTED,//6
		LB_ERROR_WORKERS_NOT_SET,       //7
		LB_WARNING_ITERATION_OVERFLOW,  //8
		LB_WARNING_NO_AVAILABLE_WORKERS,//9
		LB_WARNING_NOTHING_TO_DO,       //10
		LB_WARNING_MORE_TIME_REQUIRED,  //11
		LB_ERROR_WORKER_FINISHED,       //12
		LB_ERROR_INVALID_WORKER_NUMBER, //13
		LB_ERROR_CURRENT_TASK_RUNNING,  //14
		LB_ERROR_TASK_ALREADY_RUNNING,  //15
		LB_ERROR_HANDLED_BY_RANK0,      //16
		LB_ERROR_MPI_COMM_NOT_SET,      //17
		LB_ERROR_MPI_UNEXPECTED_TAG,    //18
		LB_ERROR_MPI_ERROR,             //19
		LB_ERROR_MPI_UNEXPECTED_SOURCE, //20
		LB_ERROR_MPI_RANK_0_EXPECTED,   //21
		LB_ERROR_MPI_ON_UPDATE,         //22
		LB_ERROR_NULL_POINTER,          //23
		LB_ERROR_INVALID_FILE,          //24
		LB_ERROR_UNABLE_TO_CREATE_FILE, //25
		LB_ERROR_TASK_ALREADY_FINISHED, //26
		LB_ERROR_CONFIGURING_CONNECTION,//27
		LB_ERROR_TASK_NOT_INITIALISED,  //28
		LB_ERROR_UNEXPECTED_FORMAT,     //29
		LB_UNABLE_TO_CREATE_REQUEST,    //30
		LB_TCP_CONNECTION_FAIL,         //31
		LB_TCP_COMMUNICATION_FAIL,      //32
		LB_REMOTE_ERROR,                //33
		LB_SSL_DISABLED,                //34
		LB_CURL_ERROR,                  //35
};

namespace LB{

    // **************************** MPI ********************************* //
#ifdef _PEN_USE_MPI_

  struct request{
  private:
    bool active;
    bool init;
    MPI_Request req;
    void* buffer;
  public:
    request() : active(false), init(false), buffer(nullptr){}
    inline int sendInit(void* buff, int count, MPI_Datatype type,
			int dest, int tag, MPI_Comm comm){
      clear();
      int err = MPI_Ssend_init(buff,count,type,dest,tag,comm,&req);
      if(err == MPI_SUCCESS){
	buffer = buff;
	init = true;
      }
      return err;
    }
    inline int receiveInit(void* buff, int count, MPI_Datatype type,
			   int dest, int tag, MPI_Comm comm){
      clear();
      int err = MPI_Recv_init(buff,count,type,dest,tag,comm,&req);
      if(err == MPI_SUCCESS){
	buffer = buff;
	init = true;
      }
      return err;      
    }
    inline int start(){
      int err = MPI_Start(&req);
      if(err == MPI_SUCCESS)
	active = true;
      return err;
    }
    inline int test(int& flag){
      if(!active){
	flag = false;
	return MPI_SUCCESS;
      }
      else{
	int err = MPI_Test(&req,&flag,MPI_STATUS_IGNORE);
	if(flag)
	  active = false;
	return err;
      }
    }
    inline int check(int& flag){
      if(!active){
	flag = false;
	return MPI_SUCCESS;
      }
      else
	return MPI_Request_get_status(req,&flag,MPI_STATUS_IGNORE);
    }
    inline int wait(){
      if(!active){
	return MPI_SUCCESS;
      }
      else{
	int err = MPI_Wait(&req,MPI_STATUS_IGNORE);
	active = false;
	return err;
      }
    }
    inline void* read(){return buffer;}
    inline bool enabled(){return active;}
    inline int cancel(){
      if(active){
	active = false;
	return MPI_Cancel(&req);
      }
      return MPI_SUCCESS;
    }
    inline void clear(){
      if(active)
	MPI_Cancel(&req);
      if(init)
	MPI_Request_free(&req);
      init = false;
      active = false;
      buffer = nullptr;
    }
    ~request(){clear();}
  };
  
  template<class timeType>
  int iSendRecv(MPI_Request reqSendRecv[2],
		timeType sleepTime,
		int& MPIerr){
    //Start both requests
    MPIerr = MPI_Startall(2,reqSendRecv);
    if(MPIerr != MPI_SUCCESS)
      return LB_ERROR_MPI_ERROR;
  
    //Iterate until both has been completed
    for(;;){
      int flag;
      MPIerr = MPI_Testall(2,reqSendRecv,flag,MPI_STATUSES_IGNORE);
      if(MPIerr != MPI_SUCCESS)
	return LB_ERROR_MPI_ERROR;
      if(flag){
	return LB_SUCCESS;
      }
      std::this_thread::sleep_for(sleepTime);
    }
  
  }

  template<class timeType>
  int noBusyWait(LB::request& req,
		 const timeType sleepTime){
    int err;
    for(;;){
      int flag;
      err = req.test(flag);
      if(err != MPI_SUCCESS)
	return err;
      if(flag)
	return MPI_SUCCESS;	    
      std::this_thread::sleep_for(sleepTime);
    }
  }

  int testAny(std::vector<LB::request>& reqs,
	      int& index);

  int checkAny(std::vector<LB::request>& reqs,
	       int& index);
  
  template<class timeType>
  int noBusyWaitAny(std::vector<LB::request>& reqs,
		    int& index,
		    const timeType sleepTime){
    for(;;){
      int id;
      int err = testAny(reqs,id);
      if(id >= 0 || err != MPI_SUCCESS){
	index = id;
	return err;
      }
      std::this_thread::sleep_for(sleepTime);
    }
  }

  template<class timeType>
  int noBusyWaitAny(std::vector<LB::request>& reqs,
		    int& index,
		    const timeType sleepTime,
		    const timeType timeout,
		    timeType& elapsed){
    elapsed = timeType(0); //Set elapsed time to 0
    do{
      int id;
      int err = testAny(reqs,id);
      if(id >= 0 || err != MPI_SUCCESS){
	index = id;
	return err;
      }
      std::this_thread::sleep_for(sleepTime); //Sleep
      elapsed += sleepTime; //Accumulate elapsed time
    }while(elapsed < timeout); //Check timeout

    //No request completed
    index = -1;
    return MPI_SUCCESS;
  }

  template<class timeType>
  int noBusyFreeWaitAny(std::vector<LB::request>& reqs,
		    int& index,
		    const timeType sleepTime){
    for(;;){
      int id;
      int err = checkAny(reqs,id);
      if(id >= 0 || err != MPI_SUCCESS){
	index = id;
	return err;
      }
      std::this_thread::sleep_for(sleepTime);
    }
  }

  template<class timeType>
  int noBusyFreeWaitAny(std::vector<LB::request>& reqs,
			int& index,
			const timeType sleepTime,
			const timeType timeout,
			timeType& elapsed){
    elapsed = timeType(0); //Set elapsed time to 0
    do{
      int id;
      int err = checkAny(reqs,id);
      if(id >= 0 || err != MPI_SUCCESS){
	index = id;
	return err;
      }
      std::this_thread::sleep_for(sleepTime); //Sleep
      elapsed += sleepTime; //Accumulate elapsed time
    }while(elapsed < timeout); //Check timeout

    //No request completed
    index = -1;
    return MPI_SUCCESS;
  }
	
  template<class timeType>
  int noBusyFreeWaitAny(std::vector<LB::request>& reqs,
			int& index,
			const timeType sleepTime,
			std::mutex& mtx,
			bool& condition){
    for(;;){
      int id;
      int err = checkAny(reqs,id);
      if(id >= 0 || err != MPI_SUCCESS){
	index = id;
	return err;
      }
      //Check the condition
      {
	const std::lock_guard<std::mutex> lock(mtx);
	if(condition){
	  index = -1; //Flag condition accomplished
	  return MPI_SUCCESS;
	}
      }
      std::this_thread::sleep_for(sleepTime);
    }
  }

  template<class timeType>
  int noBusyWaitAny(std::vector<LB::request>& reqs,
		    int& index,
		    const timeType sleepTime,
		    const timeType timeout,
		    timeType& elapsed,
		    std::mutex& mtx,
		    bool& condition){
    elapsed = timeType(0); //Set elapsed time to 0
    do{
      int id;
      int err = testAny(reqs,id);
      if(id >= 0 || err != MPI_SUCCESS){
	index = id;
	return err;
      }
      //Check the condition
      {
	const std::lock_guard<std::mutex> lock(mtx);
	if(condition){
	  index = -1; //Flag condition accomplished
	  return MPI_SUCCESS;
	}
      }
      std::this_thread::sleep_for(sleepTime); //Sleep
      elapsed += sleepTime; //Accumulate elapsed time
    }while(elapsed < timeout); //Check timeout

    //No request completed
    index = -2;
    return MPI_SUCCESS;
  }
  
  template<class timeType>
  int noBusySsend(const void* buff, int count,
		  MPI_Datatype datatype, int dest,
		  int tag, MPI_Comm comm,
		  const timeType sleepTime){
    MPI_Request req;
    int err = MPI_Issend(buff,count,datatype,dest,tag,comm,&req);
    if(err != MPI_SUCCESS)
      return err;
    for(;;){
      int flag;
      err = MPI_Test(&req,&flag,MPI_STATUS_IGNORE);
      if(err != MPI_SUCCESS)
	return err;
      if(flag)
	return MPI_SUCCESS;
      std::this_thread::sleep_for(sleepTime);      
    }
  }

  template<class timeType>
  int noBusyRecv(void* buff, int count,
		 MPI_Datatype datatype, int source,
		 int tag, MPI_Comm comm,
		 const timeType sleepTime){
    MPI_Request req;
    int err = MPI_Irecv(buff,count,datatype,source,tag,comm,&req);
    if(err != MPI_SUCCESS)
      return err;
    for(;;){
      int flag;
      err = MPI_Test(&req,&flag,MPI_STATUS_IGNORE);
      if(err != MPI_SUCCESS)
	return err;
      if(flag)
	return MPI_SUCCESS;
      std::this_thread::sleep_for(sleepTime);      
    }
  }
  
#endif
  // ************************** MPI END ******************************* //
  
  struct reportSchedule{
    std::chrono::seconds remaining;
    std::chrono::seconds interval;

    reportSchedule(const std::chrono::seconds& arg1,
		   const std::chrono::seconds& arg2) : remaining(arg1),
						       interval(arg2){}
    inline void reset(){
      remaining = interval;
    }
    inline void correct(float dev){
      interval *= dev;
    }

    inline std::chrono::seconds::rep dt() const {return interval.count();}
    
  };
	
  struct measure{
    //Timestamp of each measure (in seconds).
    std::chrono::seconds::rep timestamp;
    //Speed measure
    float speed;

    measure(const std::chrono::seconds::rep& t,
	    const float s) : timestamp(t),
			     speed(s)
    {}
  };

  struct worker{
  public:
    //Number of iterations assigned
    unsigned long long assigned;
  protected:
    //Controls if this worker has started the task
    bool started;
    //Controls if this worker has finished the task
    bool finished;
    //Iterations done
    unsigned long long iterDone;
    //Timestamps:
    //LastTime stores the timestamp of the last report
    std::chrono::steady_clock::time_point lastTime;
    //Init time is set at "start" call
    std::chrono::steady_clock::time_point init;
    //Measures history.
    //The time origin (0s) is set to the first "start" call,
    //i.e. when the task begins.
    std::vector<LB::measure> measures;
	  
  public:
	  
    worker() : assigned(0), started(false),
	       finished(false), iterDone(0){
      measures.reserve(200); //Reserve space for 200 measures
    }

    inline void setLastTime(const std::chrono::steady_clock::time_point t){
      lastTime = t;
    }

    inline size_t registers() const {return measures.size();}

    
    inline LB::measure readMeasure() const {
      return measures.back();
    }
    inline LB::measure readMeasure(const size_t index) const {
      return measures[index];
    }
    
    float
    addMeasure(const std::chrono::steady_clock::time_point& actualTime,
	       const unsigned long long nIter);
	  
    inline float
    addMeasure(const unsigned long long nIter){
      return addMeasure(std::chrono::steady_clock::now(),nIter);
    }
	  
    inline float
    addMeasure(const std::chrono::seconds::rep elapsedTime,
	       const unsigned long long nIter){
	    
      //Calculate actual time
      std::chrono::steady_clock::time_point actualTime = lastTime;
      actualTime += std::chrono::seconds(elapsedTime);
	    
      //Call add measure
      return addMeasure(actualTime,nIter);
    }
	    
    inline unsigned long long done() const {return iterDone;}

    inline bool begins() const {return started;};
    inline bool finishes() const {return finished;};
	  
    unsigned long long predDone(std::chrono::steady_clock::time_point now) const;

    inline float speed() const {
      if(measures.size() > 0)
	return measures.back().speed;
      return 0.0;
    }

    inline std::chrono::seconds::rep elapsed(const std::chrono::steady_clock::time_point time) const {
      if(time < lastTime)
	return 0;
      return std::chrono::duration_cast<std::chrono::seconds>(time-lastTime).count();
    }

    inline std::chrono::seconds::rep timeStamp(const std::chrono::steady_clock::time_point time){
      return std::chrono::duration_cast<std::chrono::seconds>(time-init).count();
    }
    
    inline std::chrono::seconds::rep timeStamp(){
      return timeStamp(std::chrono::steady_clock::now());
    }    
    
    int start(const std::chrono::steady_clock::time_point& init,
	      const unsigned long long nIter);
    void end();

    inline void reset(){
      assigned = 0;
      started = finished = false;
      iterDone = 0;
      measures.clear();
    }

    int printReport(FILE* fout,
		    const std::chrono::steady_clock::time_point globStart) const;

    inline int printReport(FILE* fout) const{
      return printReport(fout,init);
    }
    
    int dump(FILE* fdump) const;
    int load(FILE* fin,
	     const std::chrono::steady_clock::time_point taskBegin);
    
  };

  struct guessWorker : public worker{

    float addMeasure(const std::chrono::steady_clock::time_point& actualTime,
		     const unsigned long long nIter);

    inline float
    addMeasure(const std::chrono::seconds::rep elapsedTime,
	       const unsigned long long nIter){
	    
      //Calculate actual time
      std::chrono::steady_clock::time_point actualTime = lastTime;
      actualTime += std::chrono::seconds(elapsedTime);
	    
      //Call add measure
      return addMeasure(actualTime,nIter);
    }
  };
	
  struct task{
  private:
    //Total number ot iterations to done
    unsigned long long iterations;
    //Locks
    std::vector<std::mutex> locks;

    //workers
    std::vector<worker> workers;

    //Starting task point. Is set when the
    //first worker begins the task
    std::chrono::steady_clock::time_point taskBegin;

    //Time point of last checkpoint
    std::chrono::steady_clock::time_point lastCheckPoint; 

    //Minimum time between checkpoints (in s)
    std::chrono::seconds::rep minTimeCheck;

    //Remaining iterations since last checkPoint
    unsigned long long  remaining;

    //Estimated Remaining time (s)
    unsigned long long remainingTime;

    //Number of workers
    size_t nworkers;

    //Task started flag
    bool started;
    //Task finished flag
    bool finished;
	  
    //Remaining time (in s) threshold to stop balancing
    unsigned long long  threshold;

    //Objective maximum deviation between updates
    float maxDeviation;

    //Local balancing log output file
    FILE* fth;
    //Create a mutex to use previous file
    std::mutex mtxfth;
    
    // **************************** MPI ********************************* //
#ifdef _PEN_USE_MPI_

    //MPI balancing log output file
    FILE* fMPI;

    MPI_Comm comm;
    //The task requires two tags, one for requests
    //and another one to process instructions
    int tagReq;
    int tagProc;
    //Vector with MPI workers
    std::vector<LB::guessWorker> MPIworkers;

    //Flag to control if MPI communication
    //has been finished. This variable should
    //be protected by mtxMPI lock. However,
    //notice that MPIfinished will be changed only
    //by "monitor" thread. So, reads on "monitor" thread
    //don't require any lock protection.
    bool MPIfinished;
	  
    //Flag which control if the "setMPI" method
    //has been called
    bool commSet;

    //Following variables are used only by the rank 0
    //----------------------------------------------------
    //Global iterations to do by all nodes and threads
    unsigned long long MPIiterations;
    //----------------------------------------------------

    //Minimum time between MPI checks
    std::chrono::seconds MPIminTimeCheck;
    //Sleep time between request checks
    std::chrono::seconds sleepTime;

    //Create a mutex to use the conditional variable
    std::mutex mtxMPI;
    //Create a flag for finish requests
    bool MPIfinishReqFlag;
    //Create a flag for request send 
    bool MPIsent;
    
#endif
    // **************************** MPI ********************************* //

    //External balance variables
    bool ext_balance;
    long int ext_iworker;
    pen_tcp::client client;

    bool ext_http;
    std::string baseURL;
    std::string startURL;
    std::string reportURL;
    std::string finishURL;
    
    inline void lockWorkers(){
      //Lock all workers
      for(std::mutex& mtx : locks)
	mtx.lock();	    
    }
    inline void unlockWorkers(){
      //Release all workers
      for(std::mutex& mtx : locks)
	mtx.unlock();
    }
	  
    void clear();

    int setIterationsTrusted(const unsigned long long nIter); 

    int extConnect(FILE* flog,
		   bool blockWrite,
		   const unsigned verbose);
    
    int extSendAndRecv(char message[pen_tcp::messageLen],
		       std::string& response,
		       FILE* flog,
		       bool blockWrite,
		       const unsigned verbose);
    int parseAssigned(const std::string& response,
		      unsigned long long& newAssign,
		      int& errCode);
    
    int extStart(int& TCPerr,
		 int& serverErr,
		 const unsigned verbose = 1,
		 bool trusted = false);
    int extReport(int& TCPerr,
		  int& serverErr,
		  const unsigned long long nIter,
		  const std::chrono::steady_clock::time_point& rep_time,
		  const unsigned verbose = 1,
		  bool trusted = true);
    int extFinish(int& TCPerr,
		  int& serverErr,
		  const unsigned long long nIter,
		  const std::chrono::steady_clock::time_point& rep_time,
		  const unsigned verbose = 1);

    int extAssign(int& TCPerr,
		  int& serverErr,
		  const unsigned verbose = 1,
		  bool trusted = true);

// ********************** HTTP SUPPORT **********************
#ifdef _PEN_USE_LB_HTTP_

    int httpGet(const std::string& url,
		std::string& bodyResponse);
    
    int extStartHTTP(int& serverErr,
		     const unsigned verbose,
		     bool trusted);

    int extReportHTTP(int& serverErr,
		      const unsigned long long nIter,
		      const std::chrono::steady_clock::time_point& rep_time,
		      const unsigned verbose,
		      bool trusted);

    int extFinishHTTP(int& serverErr,
		      const unsigned long long nIter,
		      const std::chrono::steady_clock::time_point& rep_time,
		      const unsigned verbose);
#endif
// **********************   HTTP END   **********************
    
    void extStartHandler(const unsigned retries = 5,
			 std::chrono::seconds sleeptime = std::chrono::seconds(10),
			 const unsigned verbose = 1,
			 const bool trusted = true);
    
    void extReportHandler(const unsigned long long nIterDone,
			  const std::chrono::steady_clock::time_point& rep_time,
			  const std::chrono::seconds sleeptime = std::chrono::seconds(10),
			  const unsigned verbose = 1,
			  const bool trusted = true);

    void extFinishHandler(const unsigned long long nIterDone,
			  const std::chrono::steady_clock::time_point& rep_time,
			  const std::chrono::seconds sleeptime = std::chrono::seconds(10),
			  const unsigned verbose = 1);

    unsigned long doneTrusted(){
      //Return the registered done iterations

      unsigned long long iterDone = 0;
      unsigned iw = 0;
      for(const LB::worker& taskWorker : workers){
	//Add its prediction of iterations done
	iterDone += taskWorker.done();
	++iw;
      }
      return iterDone;
    }
    
  public:
    task() : iterations(0),
	     minTimeCheck(300),
	     remaining(0),
	     remainingTime(100000000000000),
	     nworkers(0), started(false),
	     finished(false), threshold(300),
	     maxDeviation(0.1),ext_balance(false),
	     ext_iworker(-1),ext_http(false){
      
      fth = nullptr;
      // **************************** MPI ********************************* //
#ifdef _PEN_USE_MPI_
      fMPI = nullptr;
      commSet = false;
      MPIfinished = false;
      MPIsent = false;
      MPIfinishReqFlag = false;	  
      //Global iterations to do by all nodes and threads
      MPIiterations = 0;
      MPIminTimeCheck = std::chrono::seconds(100);
      sleepTime = std::chrono::seconds(10);
#endif
      // **************************** MPI ********************************* //
    }

    // ********************** HTTP SUPPORT **********************
#ifdef _PEN_USE_LB_HTTP_
    int extHTTPserver(const unsigned extern_iw,
		      const char* url,
		      const unsigned verbose);
#endif
    // **********************   HTTP END   **********************
    
    int extLBserver(const unsigned extern_iw,
		    const char* host,
		    const char* port,
		    const char* CAfilename,
		    const char* certFilename,
		    const char* keyFilename,
		    const char* password,
		    const char* hostname,
		    const unsigned verbose);    
    
    inline void setCheckTime(const std::chrono::seconds::rep t){
      lockWorkers();
      minTimeCheck = t;
      unlockWorkers();
    }
    inline std::chrono::seconds::rep getCheckTime(){
      lockWorkers();
      const std::chrono::seconds::rep ret = minTimeCheck;
      unlockWorkers();
      return ret;
    }
    
    inline void minTime(unsigned long long t){
      lockWorkers();
      threshold = t;
      unlockWorkers();
    }
	  
    inline size_t size() const {return nworkers;}

    inline std::chrono::seconds::rep elapsed(const std::chrono::steady_clock::time_point time) const {
      if(time < lastCheckPoint)
	return 0;
      return std::chrono::duration_cast<std::chrono::seconds>(time-lastCheckPoint).count();
    }

    inline std::chrono::seconds::rep timeStamp(const std::chrono::steady_clock::time_point time){
      return std::chrono::duration_cast<std::chrono::seconds>(time-taskBegin).count();
    }
    
    inline std::chrono::seconds::rep timeStamp(){
      return timeStamp(std::chrono::steady_clock::now());
    }
	  
    int setWorkers(const size_t nw,
		   const unsigned verbose = 1);

    int workerStart(const size_t iw, const unsigned verbose = 1);
    bool workerFinish(const size_t iw,
		      int& why,
		      const unsigned verbose = 1);

    std::chrono::seconds::rep 
    report(const size_t iw,
	   const unsigned long long nIter,
	   int* errRet,
	   const unsigned verbose = 3);

    inline void reportAsync(const size_t iw,
			    const unsigned long long nIter,
			    const unsigned verbose){
      std::thread th(&LB::task::report,this,
		     iw,nIter,nullptr,verbose);
      th.detach();
    }

    inline unsigned long long assigned(const size_t iw){

      if(iw >= nworkers){
	//Invalid worker, return 0 assigned iterations
	return 0; 
      }
  
      //Lock this worker
      const std::lock_guard<std::mutex> lock(locks[iw]);

      return workers[iw].assigned;
    }
    inline unsigned long long assigned(){
      //Return the number of iterations to do for this node
      lockWorkers();
      unsigned long long ret = iterations;
      unlockWorkers();
      return ret;
    }
    inline float speed(){
      float totalSpeed = 0;
      unsigned iw = 0;
      for(const LB::worker& worker : workers){
	const std::lock_guard<std::mutex> lock(locks[iw]);
	if(worker.begins() && !worker.finishes())
	  totalSpeed += worker.speed();
	++iw;
      }
      return totalSpeed;
    }
	  
    unsigned long long predDone(unsigned long long& donePred);
    unsigned long long done();
	  
    inline int setIterations(const unsigned long long nIter){
      lockWorkers();
      int ret = setIterationsTrusted(nIter);
      unlockWorkers();
      return ret;
    }
    
    int checkPoint(const unsigned verbose);
    
    void reset();    

    // **************************** MPI ********************************* //
#ifdef _PEN_USE_MPI_

  private:
    inline void clearMPI(){
      commSet = false;
      if(fMPI != nullptr){
	fclose(fMPI);
	fMPI = nullptr;
      }
      MPIworkers.clear();
      MPIiterations = 0;
      MPIminTimeCheck = std::chrono::seconds(100);
      MPIfinishReqFlag = false;
      MPIsent = false;
      MPIfinished = false;
    }
    void sendReportMPI(std::chrono::steady_clock::time_point& lastReport,
		       const unsigned verbose);
    float receiveReportMPI(const unsigned inode,
			  const unsigned verbose);
    unsigned long long balanceMPI(const unsigned verbose);

    inline unsigned long long speedMPI(float& speedOut){
      unsigned long long totalIterDone = 0;
      speedOut = 0;
      for(const LB::guessWorker& worker : MPIworkers){
	if(worker.begins()){
	  totalIterDone += worker.done();
	  if(!worker.finishes())
	    speedOut += worker.speed();
	}
      }
      return totalIterDone;
    }

    inline unsigned long long predDoneMPI(std::chrono::steady_clock::time_point actualTime){
      
      unsigned long long donePred = 0;
      for(const LB::worker& worker : MPIworkers){
	//Add its prediction of iterations done
	donePred += worker.predDone(actualTime);
      }
      return donePred;
    }
    
    inline bool checkFinishMPI(const unsigned verbose = 1){
      //Check if rank 0 allows us to finish the task
      //Lock MPI because this function is not
      //called from "monitor" loop
      const std::lock_guard<std::mutex> lock(mtxMPI);
      if(!MPIfinished){
	//Can't finish yet, check if a finish request
	//has been already sent or flaged to send
	int rank;
	MPI_Comm_rank(comm, &rank);
	if(rank == 0){
	  //Master node
	  if(!MPIfinishReqFlag){
	    if(verbose > 1){
	      fprintf(fMPI,"%07ld s - task:checkFinishMPI:Rank 0: MPI "
		      "finish petition enqueued.\n",
		      static_cast<long int>(timeStamp()));
	      fflush(fMPI);
	    }
	    MPIfinishReqFlag = true;	    
	  }
	}
	else{
	  //Worker node
	  if(!MPIfinishReqFlag && !MPIsent){
	    if(verbose > 1){
	      fprintf(fMPI,"%07ld s - task:checkFinishMPI:Rank %d: MPI "
		      "finish petition enqueued for send.\n",
		      static_cast<long int>(timeStamp()),rank);
	      fflush(fMPI);
	    }
	    MPIfinishReqFlag = true;
	  }
	  else if(verbose > 1){
	    fprintf(fMPI,"%07ld s - task:checkFinishMPI:Rank %d: MPI "
		    "finish petition already enqueued.\n",
		    static_cast<long int>(timeStamp()),rank);
	    fflush(fMPI);
	  }
	}
	return false;
      }
      return true;
    }
	  
  public:
	  
    int setMPI(const MPI_Comm commIn,
	       const int tagReqests,
	       const int tagProcess,
	       const unsigned long long nIter,
	       const unsigned verbose = 1);
    int setIterationsMPI(const unsigned long long nIter);

    inline void setCheckTimeMPI(const std::chrono::seconds::rep t){
      const std::lock_guard<std::mutex> lock(mtxMPI);
      MPIminTimeCheck = std::chrono::seconds(t);
    }
    inline std::chrono::seconds::rep getCheckTimeMPI(){
      const std::lock_guard<std::mutex> lock(mtxMPI);
      const std::chrono::seconds::rep ret = MPIminTimeCheck.count();
      return ret;
    }
    
    void monitor(const unsigned verbose);

    int init(const size_t nw,
	     const unsigned long long nIter,
	     const MPI_Comm commIn,
	     const int tagReqests,
	     const int tagProcess,
	     const char* logFileName,
	     const unsigned verbose = 1);

    int printReportMPI(FILE* fout) const;    
#endif
    // ************************** MPI END ******************************* //
    
    int init(const size_t nw,
	     const unsigned long long nIter,
	     const char* logFileName,
	     const unsigned verbose = 1);

    int printReport(FILE* fout) const;
    
    void skip(const unsigned long long nIter){

      //Skip "nIter" iterations removing they from the iterations to do counter
      
      // **************************** MPI ********************************* //
#ifdef _PEN_USE_MPI_
      
      //Get rank
      int rank;
      MPI_Comm_rank(comm, &rank);

      //Add skipp information of all ranks
      unsigned long long toSkip;
      MPI_Reduce(&nIter,&toSkip,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,0,comm);
      if(rank == 0 && toSkip > 0){
	//Recalculate the number of iterations to done
	if(MPIiterations > toSkip)
	  setIterationsMPI(MPIiterations-toSkip);
	else
	  setIterationsMPI(0);	  
      }
#else
      // **************************** MPI ********************************* //

      //Skip iterations locally
      if(nIter > 0){
	if(iterations > nIter)
	  setIterations(iterations-nIter);
	else
	  setIterations(0);
      }
#endif
    }

    ~task(){
      clear();
      // **************************** MPI ********************************* //
#ifdef _PEN_USE_MPI_
      clearMPI();
#endif
      // ************************** MPI END ******************************* // 
    }
  };

  class taskServer{

  private:
    //Total number ot iterations to done
    unsigned long long iterations;
    //Total number ot iterations done
    unsigned long long remaining;
    //Estimated remaining time 
    unsigned long long remainingTime;    

    //Init time is set at constructor call
    std::chrono::steady_clock::time_point taskBegin;

    //Vector of guess workers
    std::vector<guessWorker> workers;
    //Vector with workers init time points
    std::vector<std::chrono::steady_clock::time_point> workerInits;
    //Vector with workers resume time points
    std::vector<std::chrono::steady_clock::time_point> workerResumes;
    FILE* flog;

    //Remaining time (in s) threshold to stop balancing
    unsigned long long  threshold;

    //Task started flag
    bool started;
    //Task finished flag
    bool finished;
    
    void createError(const int prefix, const int errcode,
		     const char* errmessage,
		     char err[pen_tcp::messageLen],
		     const unsigned verbose);
    
  public:

    char ID[pen_tcp::messageLen];
    
    taskServer() : iterations(0), remaining(0),
		   remainingTime(100000000000),
		   flog(nullptr), threshold(300),
		   started(true), finished(false){
      taskBegin = std::chrono::steady_clock::now();
      ID[0] = '\0'; //Default empty ID
    }

    inline unsigned long long readThreshold() const {return threshold;}

    inline void setInit(){
      taskBegin = std::chrono::steady_clock::now();
    }

    inline void setInit(const std::chrono::steady_clock::time_point t){
      taskBegin = t;
    }
    
    inline void moveInit(std::chrono::seconds elapsed){
      taskBegin += elapsed;
    }
    
    inline std::chrono::steady_clock::time_point workerTP(const size_t iw,
							  const std::chrono::seconds elapsed) const {
      return workerResumes[iw] + elapsed;
    }
    
    inline void setID(const char* newID){
      snprintf(ID,pen_tcp::messageLen,"%s",newID);
    }

    inline std::chrono::seconds::rep timeStamp(const std::chrono::steady_clock::time_point& time) const{
      return std::chrono::duration_cast
	<std::chrono::seconds>(time-taskBegin).count();
    }
    
    inline std::chrono::seconds::rep timeStamp() const{
      return timeStamp(std::chrono::steady_clock::now());
    }

    inline bool filterLog(const unsigned verboseLevel,
			  const unsigned required) const {
      if(flog != nullptr && verboseLevel > required )
	return true;
      return false;
    }

    inline unsigned long long ETA() const{return remainingTime;}

    inline LB::measure readMeasure(const size_t iw) const{
      if(iw >= workers.size())
	return LB::measure(-1,-1.0);
     return workers[iw].readMeasure();
    }
    
    inline LB::measure readMeasure(const size_t iw, const size_t index) const{
      if(iw >= workers.size())
	return LB::measure(-1,-1.0);
      if(index >= workers[iw].registers())
	return LB::measure(-1,-1.0);	
      return workers[iw].readMeasure(index);
    }
    
    inline unsigned long long speed(float& globspeed,
				    const std::chrono::steady_clock::time_point& now) const{
      float totalSpeed = 0;
      unsigned long long done = 0;
      for(const LB::guessWorker& worker : workers){
	if(worker.begins() && !worker.finishes()){
	  //Check if the last report of this worker has been done recently.
	  //If the server has not received reports from this worker
	  //since "threshold", it will be considered as stopped
	  //or crashed worker.
	  unsigned long long telaps = worker.elapsed(now);
	  if(telaps <= threshold)
	    totalSpeed += worker.speed();
	}
	done += worker.done();
      }
      globspeed = totalSpeed;
      return done;
    }

    inline unsigned long long predDone() const{
      //Return the prediction of done iterations

      //Get actual time
      std::chrono::steady_clock::time_point actualTime =
	std::chrono::steady_clock::now();

      unsigned long long donePred = 0;
      for(const LB::guessWorker& worker : workers){
	//Add its prediction of iterations done
	donePred += worker.predDone(actualTime);
      }
      return donePred;
    }

    inline void setThreshold(unsigned long long t){threshold = t;}

    inline unsigned long long assigned(size_t iw){
      if(iw < workers.size())
	return workers[iw].assigned;
      else
	return 0ull;
    }
    
    void clear();

    int setLogFile(const char* logFileName,
		   const unsigned verbose = 1);
    
    int init(const size_t nw,
	     const unsigned long long nIter,
	     const char* logFileName,
	     const unsigned verbose = 1);

    int workerStart(const size_t iw,
		    const std::chrono::steady_clock::time_point& t,
		    const unsigned verbose = 1);

    int workerFinish(const size_t iw,
		     const std::chrono::steady_clock::time_point& t,
		     const unsigned long long nIter,
		     const unsigned verbose = 1);
    
    int report(const size_t iw,
	       const std::chrono::steady_clock::time_point& t,
	       const unsigned long long nIter,
	       const unsigned verbose = 1,
	       const bool toBalance = true);
    
    void balance(const unsigned verbose = 1);

    int load(FILE* fin,
	     const unsigned verbose = 1);
    
    int dump(const char* saveFilename,
	     const unsigned verbose = 1);

    int receiveReport(const size_t iw,
		      const unsigned long long nIter,
		      const long int elapsed,
		      unsigned long long& newAssign,
		      const unsigned verbose = 1);

    int receiveStart(const size_t iw,
		     const long int elapsed,
		     unsigned long long& newAssign,
		     const unsigned verbose = 1);
    
    int receiveFinish(const size_t iw,
		      const unsigned long long nIter,
		      const long int elapsed,
		      const unsigned verbose = 1);
    
    void monitor(const bool local,
		 const unsigned short port,
		 const char* CAfilename,
		 const char* certFilename,
		 const char* keyFilename,
		 const char* dhFilename,
		 const char* password,
		 const unsigned verbose = 1);

    int printReport(FILE* fout) const;
    inline int printReport(const char* filename) const{

      if(filename == nullptr)
	return LB_ERROR_NULL_POINTER;

      FILE* fout = nullptr;
      fout = fopen(filename,"w");
      if(fout == nullptr)
	return LB_ERROR_UNABLE_TO_CREATE_FILE;

      int ret = printReport(fout);
      fclose(fout);
      return ret;
    }

    ~taskServer();
  };

  std::string trim(const std::string& str,
		   const std::string& whitespace = " \t\n");  
}
#endif
