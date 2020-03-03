
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
#include "pen_geometries.hh"

void printHist(const char* filename, double* hist, double low, double step, unsigned dim, unsigned nEvents);

int main(int argc, char** argv){

  if(argc < 16){
    printf("usage: %s config-filename spatial-sampler direction-sampler energy-sampler time-sampler\n",argv[0]);
    return -1;
  }
    
  //Create random generator
  pen_rand random;
  //Set seeds
  random.setSeeds(1,1);
    
  //Read configuration file
  pen_parserSection config;
  int err = parseFile(argv[1],config);
  
  //Print read configuration
  std::string strconfig;
  config.stringify(strconfig);
  printf("Configuration:\n");
  printf("%s",strconfig.c_str());

  if(err != INTDATA_SUCCESS){
    printf("Error parsing configuration.\n");
    printf("Error code: %d\n", err);
    return -2;
  }
  
  //Create a general particle generator
  pen_genericStateGen generator;

  //Get configuration subsections
  pen_parserSection spatialConf;
  pen_parserSection directionConf;
  pen_parserSection energyConf;
  pen_parserSection timeConf;

  //Extract spatial
  err = config.readSubsection("spatialSampling",spatialConf);
  if(err != INTDATA_SUCCESS){
    printf("Unable to extract spatial sampling configuration\n");
    printf("Error code: %d\n",err);
    return err;
  }

  printf("\nSpatial config:\n");
  std::string spatConfstr;
  spatialConf.stringify(spatConfstr);
  printf("%s\n",spatConfstr.c_str());

  //Extract direction
  err = config.readSubsection("directionSampling",directionConf);
  if(err != INTDATA_SUCCESS){
    printf("Unable to extract direction sampling configuration\n");
    printf("Error code: %d\n",err);
    return err;
  }  

  printf("Direction config:\n");
  std::string dirConfstr;
  directionConf.stringify(dirConfstr);
  printf("%s\n",dirConfstr.c_str());
  
  
  //Extract energy
  err = config.readSubsection("energySampling",energyConf);
  if(err != INTDATA_SUCCESS){
    printf("Unable to extract energy sampling configuration\n");
    printf("Error code: %d\n",err);
    return err;
  }

  printf("Energy config:\n");
  std::string energyConfstr;
  energyConf.stringify(energyConfstr);
  printf("%s\n",energyConfstr.c_str());
  
  //Extract time
  err = config.readSubsection("timeSampling",timeConf);
  if(err != INTDATA_SUCCESS){
    printf("Unable to extract time sampling configuration\n");
    printf("Error code: %d\n",err);
    return err;
  }

  printf("Time config:\n");
  std::string timeConfstr;
  timeConf.stringify(timeConfstr);
  printf("%s\n",timeConfstr.c_str());  

  
  // Select samplers
  ///////////////////////

  printf("\n\nAvailable samplers:\n");
  printf("%s\n",pen_genericStateGen::samplersList().c_str());
  
  //Select spatial sampler
  printf("Spatial sampler:\n");
  err = generator.selectSpatialSampler(argv[2],spatialConf,3);
  if(err == -1){
    printf("Unable to create spatial sampler: %s\n",argv[2]);
    return -3;
  }
  else if(err == -2){
    printf("Error at configuration file.\n");
    return -4;
  }
  printf("\nSelected '%s'\n",argv[2]);

  //Select Direction sampler
  printf("Direction sampler:\n");
  err = generator.selectDirectionSampler(argv[3],directionConf,3);
  if(err == -1){
    printf("Unable to create direction sampler: %s\n",argv[3]);
    return -3;
  }
  else if(err == -2){
    printf("Error at configuration file.\n");
    return -4;
  }  
  printf("\nSelected '%s'\n",argv[3]);

  //Select Energy sampler
  printf("Energy sampler:\n");
  err = generator.selectEnergySampler(argv[4],energyConf,3);
  if(err == -1){
    printf("Unable to create energy sampler: %s\n",argv[4]);
    return -3;
  }
  else if(err == -2){
    printf("Error at configuration file.\n");
    return -4;
  }
  printf("\nSelected '%s'\n",argv[4]);

  //Select time sampler
  printf("Time sampler:\n");
  err = generator.selectTimeSampler(argv[5],timeConf,3);
  if(err == -1){
    printf("Unable to create time sampler: %s\n",argv[5]);
    return -3;
  }
  else if(err == -2){
    printf("Error at configuration file.\n");
    return -4;
  }    
  printf("\nSelected '%s'\n",argv[5]);

  //Create a dummy geometry
  pen_dummyGeo geometry;
  generator.setGeometry(&geometry);
  
  FILE* fout = fopen("sampled.data","w");

  unsigned const histDim = 1000;
  double Ebins[histDim];
  const double Elow = atof(argv[6]);
  const double Estep = atof(argv[7]);

  printf("Energy range: %12.4E - %12.4E\n",Elow,Elow+Estep*histDim);
  
  double Xbins[histDim];
  const double Xlow = atof(argv[8]);
  const double Xstep = atof(argv[9]);

  printf("X range: %12.4E - %12.4E\n",Xlow,Xlow+Xstep*histDim);
  
  double Ybins[histDim];
  const double Ylow = atof(argv[10]);
  const double Ystep = atof(argv[11]);

  printf("Y range: %12.4E - %12.4E\n",Ylow,Ylow+Ystep*histDim);

  double Zbins[histDim];
  const double Zlow = atof(argv[12]);
  const double Zstep = atof(argv[13]);

  printf("Z range: %12.4E - %12.4E\n",Zlow,Zlow+Zstep*histDim);
  
  double Ubins[histDim];
  const double Ulow = 0.0;
  const double Ustep = 1.0/double(histDim);

  printf("U range: %12.4E - %12.4E\n",Ulow,Ulow+Ustep*histDim);
  
  double Vbins[histDim];
  const double Vlow = 0.0;
  const double Vstep = 1.0/double(histDim);

  printf("V range: %12.4E - %12.4E\n",Vlow,Vlow+Vstep*histDim);
  
  double Wbins[histDim];
  const double Wlow = 0.0;
  const double Wstep = 1.0/double(histDim);

  printf("W range: %12.4E - %12.4E\n",Wlow,Wlow+Wstep*histDim);
  
  double Tbins[histDim];
  const double Tlow = atof(argv[14]);
  const double Tstep = atof(argv[15]);

  printf("Time range: %12.4E - %12.4E\n",Tlow,Tlow+Tstep*histDim);

  unsigned nEvents = 10000;
  for(unsigned i = 0; i < nEvents; i++){
    pen_particleState state;

    generator.sample(state,random);

    fprintf(fout,"%d %12.4E %12.4E %12.4E %12.4E %12.4E %12.4E %12.4E %12.4E\n",
	    i, state.E, state.X, state.Y, state.Z,
	    state.U, state.V, state.W, state.PAGE);

    int Ebin = (state.E - Elow)/Estep;
    int Xbin = (state.X - Xlow)/Xstep;
    int Ybin = (state.Y - Ylow)/Ystep;
    int Zbin = (state.Z - Zlow)/Zstep;
    int Ubin = (state.U - Ulow)/Ustep;
    int Vbin = (state.V - Vlow)/Vstep;
    int Wbin = (state.W - Wlow)/Wstep;
    int Tbin = (state.PAGE - Tlow)/Tstep;
    
    //Energy
    if(Ebin >= 0 && Ebin < (int)histDim){
      Ebins[Ebin]++;
    }

    //Position
    if(Xbin >= 0 && Xbin < (int)histDim){
      Xbins[Xbin]++;
    }
    if(Ybin >= 0 && Ybin < (int)histDim){
      Ybins[Ybin]++;
    }
    if(Zbin >= 0 && Zbin < (int)histDim){
      Zbins[Zbin]++;
    }
    
    //Direction
    if(Ubin >= 0 && Ubin < (int)histDim){
      Ubins[Ubin]++;
    }
    if(Vbin >= 0 && Vbin < (int)histDim){
      Vbins[Vbin]++;
    }
    if(Wbin >= 0 && Wbin < (int)histDim){
      Wbins[Wbin]++;
    }

    //Time
    if(Tbin >= 0 && Tbin < (int)histDim){
      Tbins[Tbin]++;
    }
    
  }

  fclose(fout);

  //Energy
  printHist("energies.dat", Ebins, Elow, Estep, histDim, nEvents);
  //Position
  printHist("posX.dat", Xbins, Xlow, Xstep, histDim, nEvents);
  printHist("posY.dat", Ybins, Ylow, Ystep, histDim, nEvents);
  printHist("posZ.dat", Zbins, Zlow, Zstep, histDim, nEvents);
  //Direction
  printHist("dirU.dat", Ubins, Ulow, Ustep, histDim, nEvents);
  printHist("dirV.dat", Vbins, Vlow, Vstep, histDim, nEvents);
  printHist("dirW.dat", Wbins, Wlow, Wstep, histDim, nEvents);
  //Time
  printHist("time.dat", Tbins, Tlow, Tstep, histDim, nEvents);
  
  return 0;
}

void printHist(const char* filename, double* hist, double low, double step, unsigned dim, unsigned nEvents){

  FILE* df = fopen(filename, "w");

  fprintf(df,"#    Bin        low        prob\n");
  for(unsigned i = 0; i < dim; i++)
    fprintf(df,"%10d %12.5E %12.5E\n",i,low+step*double(i),(double)hist[i]/(double)nEvents);

  fclose(df);
}
