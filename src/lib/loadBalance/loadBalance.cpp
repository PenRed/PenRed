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
//        vicent.gimenez.alventosa@gmail.com
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
  std::chrono::seconds::rep elapsedTime = std::chrono::duration_cast<
    std::chrono::seconds>(actualTime-lastTime).count();

  //Calculate elapsed time since task beginning
  std::chrono::seconds::rep measureStamp = std::chrono::duration_cast<
    std::chrono::seconds>(actualTime-init).count();
      
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
  lastTime = init;

  //Initialize iterations done
  iterDone = 0;
  //Set assigned iterations
  assigned = nIter;
  
  //Register the first dummy measure on worker start
  addMeasure(0); //Elapsed time = 0, iterations done = 0

  return LB_SUCCESS;
}

void LB::worker::end() {
  //Calculate elapsed time since previous iteration update
  const auto actualTime = std::chrono::steady_clock::now();
  //Add a final measure with the same speed of the last report
  std::chrono::seconds::rep measureStamp = std::chrono::duration_cast<
    std::chrono::seconds>(actualTime-init).count();
  measures.push_back(LB::measure(measureStamp,speed()));
  finished = true;
}

unsigned long long LB::worker::predDone(std::chrono::steady_clock::time_point now) const {
  //Return the prediction of iterations done at timestamp "now"

  if(measures.size() < 1) return 0;
	    
  if(lastTime >= now)
    return iterDone;
	    
  //Calculate elapsed time since previous report
  std::chrono::seconds::rep timeInterval = std::chrono::duration_cast<
    std::chrono::seconds>(now-lastTime).count();

  //Calculate non registered time contribution
  double doneSinceUpdate = measures.back().speed*timeInterval;

  return iterDone+static_cast<unsigned long long>(doneSinceUpdate);
}

int LB::worker::printReport(FILE* fout) const{
  if(fout == nullptr)
    return LB_ERROR_NULL_POINTER;

  fprintf(fout,"# timestamp (s) | speed (iter/s)\n");  
  for(const auto& measure: measures){
    long long int timestamp =
      static_cast<long long int>(measure.timestamp);
    fprintf(fout," %15lld   %10.3f\n",timestamp,measure.speed);
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
  std::chrono::seconds::rep elapsedTime = std::chrono::duration_cast<
    std::chrono::seconds>(actualTime-lastTime).count();

  //Calculate elapsed time since task beginning
  std::chrono::seconds::rep measureStamp = std::chrono::duration_cast<
    std::chrono::seconds>(actualTime-init).count();

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

  if(verbose > 1){
    const std::lock_guard<std::mutex> lockTh(mtxfth);      
    fprintf(fth,"%07ld s - LB:task:worker %lu: Enter the task.\n",
	    static_cast<long int>(timeStamp()),
	    static_cast<unsigned long>(iw));
  }
  
  //Check if this task has started already
  if(!started){
    //Get the timestamp when the task begins
    if(verbose > 1){
      const std::lock_guard<std::mutex> lockTh(mtxfth);      
      fprintf(fth,"%07ld s - LB:task:worker %lu: Initialize the task.\n",
	      static_cast<long int>(timeStamp()),
	      static_cast<unsigned long>(iw));
    }
    
    taskBegin = std::chrono::steady_clock::now();
    lastCheckPoint = taskBegin;
    started = true; //Task begins now
    finished = false;
    //Sleep a few seconds to avoid an interval of 0 seconds
    //on first worker measure
    std::this_thread::sleep_for (std::chrono::seconds(3));

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
  int ret = taskWorker.start(taskBegin,toDo);
    
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
      }
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
      fprintf(fth,"%07ld s - LB:task:report:Error: Worker %lu doesn't started.\n",
	      static_cast<long int>(timeStamp()),iw);
    }
    if(errRet != nullptr) (*errRet) = LB_ERROR_WORKER_NOT_STARTED;
    return std::chrono::seconds::rep(-1);
  }
  
  if(taskWorker.finishes()){
    if(verbose > 0){
      const std::lock_guard<std::mutex> lockTh(mtxfth);      
      fprintf(fth,"%07ld s - LB:task:report:Error: Worker %lu has finished.\n",
	      static_cast<long int>(timeStamp()),iw);
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
    }
    unlockWorkers();
    return LB_WARNING_NO_AVAILABLE_WORKERS;
  }
  
  if(!started || finished || iterations == 0){
    if(verbose > 1){
      const std::lock_guard<std::mutex> lockTh(mtxfth);      
      fprintf(fth,"%07ld s - LB:task:checkPoint:Warning: Nothing to do.\n",
	      static_cast<long int>(timeStamp()));
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
    }
    mtxMPI.unlock();
#else
  // **************************** MPI ********************************* //
    if(remainingTime > threshold){
      if(verbose > 2){
	const std::lock_guard<std::mutex> lockTh(mtxfth);      
	fprintf(fth,"%07ld s - LB:task:checkPoint: Remaining time above "
		"the threshold, perform a balancing process.\n",
		static_cast<long int>(timeStamp()));            
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

int LB::task::setIterations(const unsigned long long nIter){

  lockWorkers();

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
  
  unlockWorkers();
  
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
    unsigned long long MPIremainingTime = balanceMPI();

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
    unsigned long long MPIremainingTime = balanceMPI();
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

unsigned long long LB::task::balanceMPI(){
  
  //Get speed and iterations done
  float totalSpeed;
  unsigned long long totalIterDone = speedMPI(totalSpeed);
  if(totalSpeed <= 0.0)
    return 1000000000000000;

  //Check if there are something to do
  if(totalIterDone >= MPIiterations){
    for(guessWorker& worker : MPIworkers){
      if(worker.begins() && !worker.finishes()){
	worker.assigned = worker.done();
      }
    }
    return 0ull;
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
    return static_cast<double>(MPIremaining)/totalSpeed;
  }
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
