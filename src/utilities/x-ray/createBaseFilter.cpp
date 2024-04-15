
#include "x-ray.hh"


int main(int argc, char** argv){

  //Check arguments
  if(argc < 5){
    printf("usage: %s dx dy dz n-x-vergex-groups", argv[0]);
    return 1;
  }

  //Get arguments
  double dx = std::atof(argv[1]);
  double dy = std::atof(argv[2]);
  double dz = std::atof(argv[3]);
  unsigned nVG = std::atoi(argv[4]);

  std::ofstream out("filter.msh", std::ofstream::out);
    
  penred::xray::createBaseFilter(dx,dy,dz,nVG,out,1,"filter","void",true);

  out.close();
  
  return 0;
}
