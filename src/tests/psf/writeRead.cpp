
//
//
//    Copyright (C) 2019 Universitat de València - UV
//    Copyright (C) 2019 Universitat Politècnica de València - UPV
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
//        vicente.gimenez@uv.es
//    
//

 
#include "pen_samplers.hh"
#include "pen_phaseSpaceFile.hh"

bool compare(const double a, const double b){
  if(a == 0.0 || b == 0.0){
    if(a != b)
      return false;
  }
  else{
    if(fabs(1.0-a/b) > 1.0e-14){
      return false;
    }
  }
  return true;
}

bool compareStates(const pen_particleState& a, const pen_particleState& b){

  if(!compare(a.E,b.E)) return false;
  if(!compare(a.X,b.X)) return false;
  if(!compare(a.Y,b.Y)) return false;
  if(!compare(a.Z,b.Z)) return false;
  if(!compare(a.U,b.U)) return false;
  if(!compare(a.V,b.V)) return false;
  if(!compare(a.W,b.W)) return false;
  if(!compare(a.WGHT,b.WGHT)) return false;
  if(!compare(a.PAGE,b.PAGE)) return false;

  for(unsigned i = 0; i < 5; i++){
    if(a.ILB[i] != b.ILB[i]) return false;
  }
  
  return true;
}

int main(int argc, char** argv){

  //Create phase space file writer
  pen_psfwriter psfWriter;
  //Create phase space file reader
  pen_psfreader psfReader;

  //Create a buffer of particle states
  pen_particleState* states1 = nullptr;
  size_t dimBuffer1 = 1351368;
  states1 = new pen_particleState[1351368];
  double* nhists1 = (double*) malloc(sizeof(double)*dimBuffer1);
  pen_KPAR* genKpars1 = (pen_KPAR*) malloc(sizeof(pen_KPAR)*dimBuffer1);

  size_t minChunks = dimBuffer1/pen_psfwriter::BUFFERSIZE;
  if(dimBuffer1 % pen_psfwriter::BUFFERSIZE > 0) minChunks++;

  size_t* particlesInChunk = (size_t*) malloc(sizeof(size_t)*minChunks*2);
  
  //Create random generator
  pen_rand random;

  if(argc > 1){
    random.rand0(atoi(argv[1]));
  }

  //Create random sampler
  random_specificSampler sampler;

  //Fill the buffer
  double prob0 = 0.8;
  double prob1 = 0.9;
  double prob2 = 0.96;
  double prob3 = 0.99;
  
  double probe = 0.4;
  double probg = 0.9;

  unsigned long long nhist = 1; //Ensure non zero history number
  for(size_t i = 0; i < dimBuffer1; i++){
    pen_KPAR kpar;
    unsigned long long dhist;
    sampler.sample(states1[i],kpar,dhist,random);

    double rand = random.rand();

    if(rand < prob0)
      dhist = 0;
    else if(rand < prob1)
      dhist = 1;
    else if(rand < prob2)
      dhist = 2;
    else if(rand < prob3)
      dhist = 3;
    else
      dhist = 4;

    nhist += dhist;
    nhists1[i] = nhist;
    
    rand = random.rand();
    if(rand < probe)
      genKpars1[i] = PEN_ELECTRON;
    else if(rand < probg)
      genKpars1[i] = PEN_PHOTON;
    else
      genKpars1[i] = PEN_POSITRON;
  }

  //Create a PSF file with generated data
  FILE* fout = fopen("random.psf","wb");
  size_t chunk = 0;
  for(size_t i = 0; i < dimBuffer1; i++){
    if(psfWriter.store(nhists1[i],genKpars1[i],states1[i]) == -1){

      particlesInChunk[chunk] = psfWriter.partConcluded();
      printf("Chunk %lu full (%lu)\n",chunk,psfWriter.stored());
      printf("Concluded particles (%lu)\n",psfWriter.partConcluded());
      printf("Last history particles (%lu)\n",psfWriter.partLastHist());
      
      //Buffer is full, write it to output file
      psfWriter.write(fout,1);
      //Clear buffer
      psfWriter.clearBuffer();
      //Save state
      psfWriter.store(nhists1[i],genKpars1[i],states1[i]);
      
      chunk++;
    }
  }

  printf("Remaining particles: %lu\n",psfWriter.stored());
  particlesInChunk[chunk] = psfWriter.stored();
  chunk++;
  //Write remaining states
  psfWriter.write(fout,1,true);
  psfWriter.clear();
  //Close output file
  fclose(fout);

  //Read the generated psf and compare with sampled data
  FILE* fin = fopen("random.psf","rb");
  double psfHist = 0.0;
  size_t nread = 0;
  size_t diffs = 0;
  size_t chunk2 = 0;
  while(psfReader.read(fin,1) == PEN_PSF_SUCCESS){

    if(psfReader.stored() != particlesInChunk[chunk2]){
      printf("Read particles at chunk %lu mismatch written.\n",chunk2);
      printf("            Expected: %lu\n",particlesInChunk[chunk2]);
      printf("                Read: %lu\n",psfReader.stored());
      diffs++;
    }

    chunk2++;
      
    unsigned long dhist;
    unsigned kpar;
    pen_particleState state;
    while(psfReader.get(dhist,kpar,state) > 0){
      psfHist += dhist;
      //Compare read state with generated
      if(psfHist != nhists1[nread]){
	printf("History number mismatch at particle %lu, chunk %lu\n",nread,chunk2);
	printf("   Expected: %.0f\n",nhists1[nread]);
	printf("       Read: %.0f\n",psfHist);
	diffs++;
	printf("%lu %.0f %u %s\n",dhist,psfHist,kpar,state.stringify().c_str());
	printf("%.0f %.0f %u %s\n", nread > 0 ? nhists1[nread]-nhists1[nread-1] : nhists1[nread],
	       nhists1[nread],
	       genKpars1[nread],
	       states1[nread].stringify().c_str());
      }
      if(kpar != genKpars1[nread]){
	printf("Kpar number mismatch at particle %lu, chunk %lu\n",nread,chunk2);
	printf("   Expected: %u\n",genKpars1[nread]);
	printf("       Read: %u\n",kpar);
	diffs++;
	printf("%lu %.0f %u %s\n",dhist,psfHist,kpar,state.stringify().c_str());
	printf("%.0f %.0f %u %s\n", nread > 0 ? nhists1[nread]-nhists1[nread-1] : nhists1[nread],
	       nhists1[nread],
	       genKpars1[nread],
	       states1[nread].stringify().c_str());
      }
      if(!compareStates(state,states1[nread])){
	printf("States mismatch at particle %lu, chunk %lu\n",nread,chunk2);
	printf("   Expected: %s\n",state.stringify().c_str());
	printf("       Read: %s\n",states1[nread].stringify().c_str());
	diffs++;
	printf("%lu %.0f %u %s\n",dhist,psfHist,kpar,state.stringify().c_str());
	printf("%.0f %.0f %u %s\n", nread > 0 ? nhists1[nread]-nhists1[nread-1] : nhists1[nread],
	       nhists1[nread],
	       genKpars1[nread],
	       states1[nread].stringify().c_str());
      }
      
      nread++;
    }
  }

  if(chunk2 != chunk){
    printf("Read chunks mismatch written.\n");
    printf("        Expected: %lu\n",chunk);
    printf("            Read: %lu\n",chunk2);
    diffs++;    
  }
  
  if(diffs == 0){
    printf("\n\nTest completed!\n");
  }
  else{
    printf("\n\nTest failed!\n");
  }
  
  delete [] states1;
  free(nhists1);
  free(genKpars1);
  
  return 0;
}
