 
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
//        vicent.gimenez.alventosa@gmail.com
//    
//

 
#include "pen_data.hh"
#include <array>

constexpr const char* readFormat = R"===(

# Test format reader

## Sources

# History number
sources/${subsection}/use/Isotope/reader-value False 
sources/${subsection}/use/Isotope/reader-required/type "optional"
sources/${subsection}/isotope/to/use/reader-value "I-192" 
sources/${subsection}/isotope/to/use/reader-required/type "required_if"  
sources/${subsection}/isotope/to/use/reader-required/value "use/Isotope"  
sources/${subsection}/nhists/reader-value 1.0e5
sources/${subsection}/nhists/reader-description "Number of particles to be simulated"
sources/${subsection}/nhists/reader-required/type "optional_if_exist"
sources/${subsection}/nhists/reader-required/value "use/Isotope"
sources/${subsection}/nhists/reader-conditions/positive/type "positive"
sources/${subsection}/nhists/reader-conditions/not-too-much/type "lesser"
sources/${subsection}/nhists/reader-conditions/not-too-much/value 1.0e10

# Particle type
sources/${subsection}/particle/reader-value "gamma"
sources/${subsection}/particle/reader-required/type "optional"

# Particle type
sources/${subsection}/positions/${subsection}/reader-description "Origin and orientation section"
sources/${subsection}/positions/${subsection}/origin/reader-value [0.0,0.0,0.0]
sources/${subsection}/positions/${subsection}/direction/reader-value [1.0,0.0,0.0]


## Tally parameters
tally/mesh/dx/reader-description "voxel 'x' size in cm"
tally/mesh/dx/reader-value 1.0
tally/mesh/dy/reader-value 1.0
tally/mesh/dz/reader-value 1.0
tally/mesh/nx/reader-value 10
tally/mesh/nx/reader-conditions/positive/type "positive"
tally/mesh/ny/reader-value 10
tally/mesh/ny/reader-conditions/positive/type "positive"
tally/mesh/nz/reader-value 10
tally/mesh/nz/reader-conditions/positive/type "positive"
tally/mesh/nz/reader-conditions/lesserThan/type "lesser"
tally/mesh/nz/reader-conditions/lesserThan/value "tally/mesh/nx"

## Geometry parameters
geometry/body/${subsection}/reader-value 1
geometry/body/${subsection}/reader-conditions/positive/type "positive"
)===";

const char* testConfig = R"===(

tally/mesh/dx    2.0
tally/mesh/dy    0.5
tally/mesh/dz    1
tally/mesh/nx 100
tally/mesh/ny 120.2
tally/mesh/nz 99

sources/gamma/nhists    98

sources/gamma/positions/1g/direction [   1.00000E+00,   1.00000E+00,   0.00000E+00]
sources/gamma/positions/1g/origin [   0.50000E+00,   0.30000E+00,   0.20000E+00]

sources/gamma/positions/2g/direction [   -1.00000E+00,   -1.00000E+00,   0.00000E+00]
sources/gamma/positions/2g/origin [   -0.50000E+00,   -0.30000E+00,   -0.20000E+00]

sources/electron/use/Isotope false
sources/electron/particle "electron"

sources/electron/positions/1e/direction [   1.30000E+00,   0.60000E+00,   0.07000E+00]
sources/electron/positions/1e/origin [   0.40000E+00,   0.20000E+00,   0.10000E+00]

sources/electron/positions/2e/direction [   -1.30000E+00,   -0.60000E+00,   -0.07000E+00]
sources/electron/positions/2e/origin [   -0.40000E+00,   -0.20000E+00,   -0.10000E+00]

geometry/body/detector 1
geometry/body/anode -2
geometry/body/filter 3
)===";

struct test_reader : public pen_readerStorage{

  // 0 root, 1 sources, 2 sources/${subsection}/positions, 3 bodies
  unsigned actualSection; 

  struct positions{
    std::string name;
    std::array<double,3> origin;
    std::array<double,3> direction;
  };

  struct sources{
    std::string name;
    unsigned nhists;
    std::string particleType;
    std::vector<positions> vpos;
  };

  std::vector<sources> vsources;
  double dx,dy,dz;
  unsigned nx,ny,nz;

  struct body{
    std::string name;
    int index;
  };
  std::vector<body> bodies;
  
  inline int beginSectionFamily(const std::string& pathInSection,
				const size_t /*size*/,
				const unsigned) override {
    if(actualSection == 0){
      if(pathInSection.compare("sources") == 0){
	actualSection = 1;
	return 0;
      }else if(pathInSection.compare("geometry/body") == 0){
	actualSection = 3;
	return 0;
      }
      return -1;
    }else if(actualSection == 1){
      if(pathInSection.compare("positions") == 0){
	actualSection = 2;
	return 0;
      }
      return -2;
    }
    return -3;
  }

  inline int endSectionFamily(const unsigned) override {
    if(actualSection == 1){
      actualSection = 0;
      return 0;
    }else if(actualSection == 2){
      actualSection = 1;
      return 0;
    }else if(actualSection == 3){
      actualSection = 0;
      return 0;
    }
    return -1;
  }

  inline int beginSection(const std::string& name,
			  const unsigned) override {
    if(actualSection == 1){
      vsources.emplace_back();
      vsources.back().name = name;
      return 0;
    }else if(actualSection == 2){
      vsources.back().vpos.emplace_back();
      vsources.back().vpos.back().name = name;
      return 0;
    }else if(actualSection == 3){
      bodies.emplace_back();
      bodies.back().name = name;
      return 0;
    }

    return -1;
  }

  inline int endSection(const unsigned) override { return 0; }

  inline int beginArray(const std::string& ,
			const size_t size,
			const unsigned) override {
    
    if(size != 3)
      return -2;
    
    return 0;
  }

  inline int endArray(const unsigned ) override { return 0; }

  inline int storeElement(const std::string& pathInSection,
			  const pen_parserData& element,
			  const unsigned) override {

    if(actualSection == 0){
      if(pathInSection.compare("tally/mesh/dx") == 0)
	dx = element;
      else if(pathInSection.compare("tally/mesh/dy") == 0)
	dy = element;
      else if(pathInSection.compare("tally/mesh/dz") == 0)
	dz = element;
      else if(pathInSection.compare("tally/mesh/nx") == 0)
	nx = element;
      else if(pathInSection.compare("tally/mesh/ny") == 0)
	ny = element;
      else if(pathInSection.compare("tally/mesh/nz") == 0)
	nz = element;
      else
	return -1;
      return 0;
    } else if (actualSection == 1){
      if(pathInSection.compare("nhists") == 0){
	vsources.back().nhists = element;
      }else if(pathInSection.compare("use/Isotope") == 0){
	return 0;
      }
      else{
	return -1;
      }
      return 0;
    } else if(actualSection == 3){
      if(pathInSection.empty()){
	bodies.back().index = element;
	return 0;
      }
      return -1;
    } else{
      return -1;
    }
    return -1;
  }

  inline int storeArrayElement(const std::string& pathInSection,
				const pen_parserData& element,
				const size_t pos,
				const unsigned ) override {

    if(actualSection == 2){
      if(pathInSection.compare("direction") == 0){
	if(pos < 3){
	  vsources.back().vpos.back().direction[pos] = element;
	  return 0;
	}
      }else if(pathInSection.compare("origin") == 0){
	if(pos < 3){
	  vsources.back().vpos.back().origin[pos] = element;
	  return 0;
	}
      }
    }
    
    return -1;
  }

  inline int storeString(const std::string& pathInSection,
			 const std::string& element,
			 const unsigned) override {
    if(actualSection == 1){
      if(pathInSection.compare("particle") == 0){
	vsources.back().particleType = element;
	return 0;
      }else if(pathInSection.compare("isotope/to/use") == 0)
	return 0;
    }

    return -1;
  }

  inline void print(){

    printf("Tally parameters: \n"
	   "dx = %E\n"
	   "dy = %E\n"
	   "dz = %E\n\n"
	   "nx = %u\n"
	   "ny = %u\n"
	   "nz = %u\n",	   
	   dx, dy, dz,
	   nx, ny, nz);
    
    for(const sources& s : vsources){

      printf("\nSource: %s\n", s.name.c_str());
      printf(" histories: %u\n", s.nhists);
      printf(" particle: %s \n", s.particleType.c_str());
      
      for(const positions& p : s.vpos){

	printf("\n   Position: %s\n", p.name.c_str());
	printf("       origin: (%E,%E,%E)\n", p.origin[0], p.origin[1], p.origin[2]);
	printf("    direction: (%E,%E,%E)\n", p.direction[0], p.direction[1], p.direction[2]);
	
      }
      
    }

    printf("\n Bodies: \n\n");
    for(const body& b : bodies){
      printf(" - %s %d\n", b.name.c_str(), b.index);
    }
  }
  
};

int main(){

  pen_readerSection rs;

  std::string errorString;
  unsigned long errorLine;  
  int err = rs.setFormat(readFormat, errorString, errorLine);
  if(err != INTDATA_SUCCESS){
    printf("Error on format set: %s\n", pen_parserError(err));
    printf("Error data:\n%s\n", errorString.c_str());
    return -1;
  }

  printf("Format:\n\n");
  printf("%s\n",rs.stringify().c_str());

  // ** Try to read the provided example

  //Parse provided configuration
  pen_parserSection data;

  err = parseString(testConfig, data, errorString, errorLine);
  if(err != INTDATA_SUCCESS){
    printf("Error parsing data at line %lu. Error code: %d\n"
	   " Error string: %s\n", errorLine, err, errorString.c_str());
    return -1;
  }

  //Fill the reader structure
  test_reader reader;
  std::string errorElement;
  std::string errorElementSpecs;
  err = rs.read(data, reader,
		errorElement,
		errorElementSpecs, 2);
  if(err != INTDATA_SUCCESS){
    printf("Error reading data. %s\n"
	   "Error at element: %s\n"
	   "Error element specs:\n%s\n",
	   pen_parserError(err),
	   errorElement.c_str(), errorElementSpecs.c_str());
    return -1;
  }

  //Print read data
  reader.print();

  //Set the element 'sources/electron/use/Isotope' to true.
  //Then, the read should crash because the existence requirements
  err = data.set("sources/electron/use/Isotope", true, true);
  if(err != INTDATA_SUCCESS){
    printf("Unable to change 'sources/electron/use/Isotope' value."
	   "Error code: %d\n"
	   "Error: %s\n",
	   err, pen_parserError(err));
  }

  printf("\nTry to force a controled error...\n");
  
  err = rs.read(data, reader,
		errorElement,
		errorElementSpecs, 2);
  if(err == INTDATA_SUCCESS){
    printf("The second parsing does not failed!\n\n%s\n",
	   data.stringify().c_str());
    return -1;
  }

  //Ensure the error is the one we expect
  if(errorElement.find("sources/electron/isotope/to/use") != 0){
    printf("Unexpected error on another element!");
    printf("Error reading data. %s\n"
	   "Error at element: %s\n"
	   "Error element specs:\n%s\n",
	   pen_parserError(err),
	   errorElement.c_str(), errorElementSpecs.c_str());    
  }

  printf("\n\nExpected error. Test passed\n\n");
  
  
  return 0;
}

