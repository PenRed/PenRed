 
#include "pen_contexts.hh"

int main(){

  const pen_readerSection& sMat =
    pen_readerSection::readObjectSection<pen_contextReaderMat>();

  printf("Material configuration: \n\n %s\n", sMat.stringify().c_str());

  printf("Material configuration YAML: \n\n%s\n", sMat.stringifyYAML().c_str());  

  const pen_readerSection& sVR =
    pen_readerSection::readObjectSection<pen_contextReaderVR>();

  printf("VR configuration: \n\n %s\n", sVR.stringify().c_str());

  printf("VR configuration YAML: \n\n%s\n", sVR.stringifyYAML().c_str());  
  
  return 0;
}
