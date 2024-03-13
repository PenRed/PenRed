
#include "x-ray.hh"


int main(int argc, char** argv){

  //Check arguments
  if(argc < 8){
    printf("usage: %s dxOut dxInTop dxInBot "
	   "dyOut dyInTop dyInBot dz\n", argv[0]);
    return 1;
  }

  //Get arguments
  double dxOut    = std::atof(argv[1]);
  double dxInTop  = std::atof(argv[2]);
  double dxInBot  = std::atof(argv[3]);
  double dyOut    = std::atof(argv[4]);
  double dyInTop  = std::atof(argv[5]);
  double dyInBot  = std::atof(argv[6]);
  double dz       = std::atof(argv[7]);

  std::ofstream out("collimator.msh", std::ofstream::out);
    
  penred::xray::createBaseCollimator(dxOut, dxInTop, dxInBot,
				     dyOut, dyInTop, dyInBot,
				     dz,
				     out,
				     1,"filter","void",true);

  out.close();
  
  return 0;
}
