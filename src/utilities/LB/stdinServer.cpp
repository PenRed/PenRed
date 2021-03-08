 
#include "loadBalance.hh"

int main(int argc, char** argv){
    
  if(argc < 4){
    printf("usage: %s nIter nworkers threshold [elapsed logfile verbose]\n"
	   "elapsed, logfile and verbose are optional values\n",argv[0]);
    return 1;
  }
    
  //Parse parameters
  long long int nIterG = atoll(argv[1]);
  if(nIterG < 1){
    printf("Invalid number of iterations: %lld\n",nIterG);
    return -2;
  }    
  int nw = atoi(argv[2]);
  if(nw < 1){
    printf("Invalid number of workers: %d\n",nw);
    return -3;
  }
  long long int threshold = atoll(argv[3]);
  if(threshold < 1){
    printf("Invalid threshold: %lld\n",threshold);
    return -4;
  }
  
  const char* logfilename = nullptr;
  long int elapsedG = 0;
  unsigned verbose = 2;
  if(argc >= 5)
    elapsedG = atol(argv[4]);
  if(argc >= 6)
    logfilename = argv[5];
  if(argc >= 7)
    verbose = static_cast<unsigned>(std::max(atoi(argv[6]),0));
 
  //Initialize server
  LB::taskServer server;
  if(elapsedG != 0)
    server.moveInit(std::chrono::seconds(elapsedG));
  server.init(nw,static_cast<unsigned long long>(nIterG),logfilename,verbose);
  server.setThreshold(static_cast<unsigned long long>(threshold));
  
  
  const size_t maxLineSize = 512;
  char message[maxLineSize];
  for(;;){
        
    //Wait until next input
    fgets(message,maxLineSize,stdin);
 
    //Get the instruction identifier
    int instruction = -1;
    if(sscanf(message,"%d",&instruction) != 1){
      printf("-1\n");
      continue;
    }
        
    long unsigned iw;
    unsigned long long nIter;
    long int elapsed;        
    switch(instruction){
    case 0:{ //Stop server
      return 0;
    }
    case 1:{ //Report
      //Parse the instruction
      int nread = sscanf(message,"%*d %lu %llu %ld",
			 &iw, &nIter, &elapsed);
      if(nread != 3){
	//Unable to parse report instruction
	printf("1 %d\n",LB_ERROR_UNEXPECTED_FORMAT);
      }
      else{
	//Perform the report
	unsigned long long newAssign = 0;
	int err = server.receiveReport(static_cast<size_t>(iw),nIter,
				       elapsed,newAssign,verbose);
	if(err == LB_SUCCESS)
	  printf("0\n Assigned: %llu\n ETA: %llu\n",
		 newAssign,server.ETA());
	else
	  printf("1 %d\n",err);
      }
    }
      break;
    case 2:{ //Start
      //Parse the instruction
      int nread = sscanf(message,"%*d %lu %ld",&iw,&elapsed);
      if(nread != 2){
	//Unable to parse start petition
	printf("2 %d\n",LB_ERROR_UNEXPECTED_FORMAT);
      }
      else{
	//Perform the worker start
	unsigned long long newAssign = 0;
	int err = server.receiveStart(static_cast<size_t>(iw),elapsed,
				      newAssign,verbose);
	//Write the response
	if(err == LB_SUCCESS){
	  printf("0 \n Assigned: %llu\n ETA: %llu\n",
		 newAssign,server.ETA());
	  fflush(stdout);
	}
	else{
	  printf("2 %d\n",err);
	  fflush(stdout);
	}
      }
    }
      break;
    case 3:{ //Finish
      //Parse the instruction
      int nread = sscanf(message,"%*d %lu %llu %ld",
			 &iw,&nIter,&elapsed);
      if(nread != 3){
	//Unable to parse finish instruction
	printf("3 %d\n",LB_ERROR_UNEXPECTED_FORMAT);
      } else{
	//Perform the finish request
	int err = server.receiveFinish(static_cast<size_t>(iw),nIter,
				       elapsed,verbose);
	//Send the response
	if(err == LB_SUCCESS)
	  printf("0\n");
	else
	  printf("3 %d\n",err);
      }
    }
      break;
    case 4:{ //Load
      //Parse the instruction
      char filename[200]; 
      int nread = sscanf(message,"%*d %200s",filename);
      if(nread != 1){
	//Unable to parse finish instruction
	printf("4 %d\n",LB_ERROR_UNEXPECTED_FORMAT);
      }else{
	if(strcmp(filename,"stdin") == 0){
	  int err = server.load(stdin,verbose);
	  if(err != LB_SUCCESS)
	    printf("4 %d\n",err);
	  else
	    printf("0\n");	  
	}
	else{
	  FILE* fin = nullptr;
	  fin = fopen(filename,"r");
	  if(fin == nullptr)
	    printf("4 -1\n");
	  else{
	    int err = server.load(fin,verbose);
	    if(err != LB_SUCCESS)
	      printf("4 %d\n",err);
	    else
	      printf("0\n");
	  }
	}
      }
    }
      break;
    case 5:{ //Save
      //Parse the instruction
      char filename[200]; 
      int nread = sscanf(message,"%*d %200s",filename);
      if(nread != 1){
	//Unable to parse finish instruction
	printf("5 %d\n",LB_ERROR_UNEXPECTED_FORMAT);
      }else{
	int err = server.dump(filename,verbose);
	if(err != LB_SUCCESS)
	  printf("5 %d\n",err);
	else
	  printf("0\n");	
      }
    }
      break;
    case 6:{ //Read measure
      //Parse worker number
      int nread = sscanf(message,"%*d %lu",&iw);
      if(nread != 1){
	printf("6 %d\n",LB_ERROR_UNEXPECTED_FORMAT);
      } else{
	LB::measure lastMeasure = server.readMeasure(iw);
	if(std::signbit(lastMeasure.speed))
	  printf("6 %d\n", LB_ERROR_WORKER_OUT_OF_RANGE);
	else
	  printf("0 \n timestamp: %ld\n speed: %.6E\n",
		 static_cast<long int>(lastMeasure.timestamp),lastMeasure.speed);
      }
    }
      break;
    default:
      printf("-1 %d\n",LB_ERROR_UNEXPECTED_FORMAT);
    }
  }
    
}
