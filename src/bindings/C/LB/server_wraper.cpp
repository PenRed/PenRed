
#include "loadBalance.hh"

extern "C"{
  LB::taskServer* taskServer_new(){return new LB::taskServer();}

  int taskServer_init(LB::taskServer* server,
		      unsigned long nw,
		      unsigned long long nIter,
		      char* logFileName,
		      unsigned verbose){
    return server->init(static_cast<size_t>(nw),nIter,logFileName,verbose);
  }
  
  unsigned long long taskServer_ETA(LB::taskServer* server){
    return server->ETA();
  }

  void taskServer_setThreshold(LB::taskServer* server,
			       unsigned long long t){
    server->setThreshold(t);
  }

  int taskServer_setLogFile(LB::taskServer* server,
			    char* logFileName){
    return server->setLogFile(logFileName);
  }

  int taskServer_receiveReport(LB::taskServer* server,
			       unsigned long iw,
			       unsigned long long nIter,
			       long int elapsed,
			       unsigned long long* newAssign,
			       unsigned verbose){
    return server->receiveReport(size_t(iw),nIter,elapsed,*newAssign,verbose);
  }

  int taskServer_receiveStart(LB::taskServer* server,
			      unsigned long iw,
			      long int elapsed,
			      unsigned long long* assigned,
			      unsigned verbose){
    return server->receiveStart(size_t(iw),elapsed,*assigned,verbose);
  }

  int taskServer_receiveFinish(LB::taskServer* server,
			       unsigned long iw,
			       unsigned long long nIter,
			       long int elapsed,
			       unsigned verbose){
    return server->receiveFinish(static_cast<size_t>(iw),nIter,elapsed,verbose);
  }
  
  void taskServer_monitor(LB::taskServer* server,
			  bool local,
			  unsigned short port,
			  char* CAfilename,
			  char* certFilename,
			  char* keyFilename,
			  char* dhFilename,
			  char* password,
			  unsigned verbose){
    server->monitor(local,port,CAfilename,certFilename,
		    keyFilename,dhFilename,password,verbose);
  }

  unsigned long long taskServer_speed(LB::taskServer* server,
				      float* globspeed){
    return server->speed(*globspeed,
			 std::chrono::steady_clock::now());
  }

  int taskServer_printReport(LB::taskServer* server,
			     char* filename){
    return server->printReport(filename);
  }

  void taskServer_setInit(LB::taskServer* server){
    server->setInit();
  }
  
  void taskServer_moveInit(LB::taskServer* server,
			   long int elapsed){
    server->moveInit(std::chrono::seconds(elapsed));
  }

  void taskServer_setID(LB::taskServer* server,
			char* newID){
    server->setID(newID);
  }
  
}
