 
#include "loadBalance.hh"

int main(const int argc, const char** argv){

  if(argc < 6){
    printf("usage: %s ID procedence port nworkers nIter CA-cert-file "
           "server-cert-file server-key-file dh-params-file key-password(optional)\n\n"
	   "procedence must be set to 0 to accept only connections "
	   "from localhost. Otherwise set it to 1\n\n"
       "NOTE: CA-cert-file, server-cert-file, server-key-file and dh-params-file "
       "are optional parameters. If ALL of them are provided, the server will "
       "enable SSL comunications. Otherwise the comunications will be performed "
       "with no security.\n",argv[0]);
    return 1;
  }

  bool local = atoi(argv[2]);
  local = !local;
  int port  = atoi(argv[3]);
  if(port < 0){
    printf("Invalid port: %d\n",port);
    return -1;
  }
  
  int nw    = atoi(argv[4]);
  if(nw < 1){
    printf("Invalid number of workers: %d\n",nw);
    return -2;
  }
  
  long long int nIter = atoll(argv[5]);
  if(nIter < 1){
    printf("Invalid number of iterations: %lld\n",nIter);
    return -3;
  }

  const char* caCertfile = nullptr;
  const char* serverCertfile = nullptr;
  const char* serverKeyfile = nullptr;
  const char* DHparamfile = nullptr;
  const char* keyPassword = nullptr;
  if(argc >= 10){
      caCertfile = argv[6];
      serverCertfile = argv[7];
      serverKeyfile = argv[8];
      DHparamfile = argv[9];
      if(argc > 10)
          keyPassword = argv[10];
  }

  unsigned verbose = 3;
    
  LB::taskServer server;
  server.setID(argv[1]);
  server.init(nw,static_cast<unsigned long long>(nIter),"serverLog.txt",verbose);

  server.monitor(local,port,caCertfile,serverCertfile,serverKeyfile,DHparamfile,keyPassword,verbose);
  
  return 0;
}
