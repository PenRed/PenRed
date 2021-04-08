//
//
//    Copyright (C) 2021 Universitat de València - UV
//    Copyright (C) 2021 Universitat Politècnica de València - UPV
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


#include "PenRed.hh"
#include "pen_geometries.hh"

bool dumpResults(const unsigned long long nhist,
		 const double EDPT, const double EDPT2,
		 const double ETRT, const double ETRT2,
		 const double RMU0, const double E0, const double TOL,
		 const long int TSEC, const bool forceFinish,
		 FILE* fout);

void simulate(const pen_context& context,
	      const double E0,
	      const double dumpTime,
	      const double simTime,
	      const double tolerance,
	      int& seed1, int& seed2,
	      FILE* fout);

int main(int argc, const char** argv){

  //Check the arguments
  if(argc < 5){
    printf("usage: %s material-file energy-spectrum-file "
	   "tolerance sim-time\n",argv[0]);
    return 1;
  }

  //Get configuration parameters

  double tolerance = atof(argv[3]);
  double simTime = atof(argv[4]);

  if(tolerance <= 1.0e-6){
    printf("Error: Tolerance must be positive and"
	   " greater than zero: %E\n",tolerance);
    return -1;
  }

  double dumpTime = 10.0;

  const char* filenameSpectrum = argv[2];

  //Get initial random seeds
  int seed1, seed2;
  rand0(5,seed1,seed2);
  
  //Create a dummy geometry
  pen_dummyGeo geometry;

  //Configure dummy geometry
  pen_parserSection dummySection;
  geometry.configure(dummySection,2);
  
  //Create elements data base
  pen_elementDataBase* elementsDB = new pen_elementDataBase;
  
  //Create a context
  pen_context context(*elementsDB);

  //Set context geometry
  context.setGeometry(&geometry);
  
  //Set the number of materials to context (1)
  int errmat = context.setMats<pen_material>(1);
  if(errmat != 0){
    printf("Error at context material creation: %d.\n",errmat);
    return -1;
  }

  //Get the material
  pen_material& mat = context.getBaseMaterial(0);
  
  //Configure the material
  mat.C1=0.2;
  mat.C2=0.2;
  mat.WCC=1.0e3;
  mat.WCR=1.0e3;

  mat.EABS[PEN_ELECTRON] = 50.0E0;
  mat.EABS[PEN_PHOTON]   = 50.0E0;
  mat.EABS[PEN_POSITRON] = 50.0E0;
  double EABS1 = mat.EABS[PEN_ELECTRON];

  //Configure context
  printf("  \n");
  printf("Configuring context with material '%s'. Please, wait...\n",argv[1]);
  
  FILE* fcontext = nullptr;
  fcontext = fopen("context.rep","w");
  if(fcontext == nullptr){
    printf("Error: unable to create file 'context.rep'\n");
    return -2;
  }

  double EMAX=1.0E9;
  int INFO = 2;
  std::string PMFILEstr[constants::MAXMAT];
  PMFILEstr[0].assign(argv[1]);
  int err = context.init(EMAX,fcontext,INFO,PMFILEstr);
  if(err != 0){
    printf("Error: Unable to configure context. Check 'context.rep'.\n");
    return -3;
  }

  fclose(fcontext);
  
  //Open output file
  FILE* fout =  nullptr;
  fout =  fopen("mutren.dat","w");
  if(fout == nullptr){
    printf("Error: Unable to create file 'mutren.dat'");
    return -4;
  }

  //Print header information
  fprintf(fout,"#  Photon mass energy-absorption coefficients\n");
  fprintf(fout,"#  Material file: %s\n",argv[1]);
  fprintf(fout,"#  Mass density: %E g/cm^3\n",mat.RHO);
  fprintf(fout,"#  \n");
  fprintf(fout,"#  Tolerance: %E %%\n",tolerance);
  fprintf(fout,"#  Allotted simulation time: %E s\n",simTime);
  fprintf(fout,"#  Dump period: %E s\n",dumpTime);
  fprintf(fout,"#  \n");
  fprintf(fout,"#    Energy       mu/rho      mu_en/rho     "
	  "(1-g)*f        f            1-g      uncert.      T\n");
  fprintf(fout,"#     (eV)       (cm^2/g)     (cm^2/g)%45.45s(%%)       (s)\n"," ");
  fprintf(fout,"# \n");
  fflush(fout);

  //Try to open spectrum file
  FILE* fspec = nullptr;
  fspec = fopen(filenameSpectrum,"r");
  if(fspec == nullptr){
    printf("Error: Unable to open spectrum file: '%s'\n",filenameSpectrum);
    return -5;
  }
  
  char line[1000];
  unsigned long auxRead = 0;
  unsigned long nread = 0;
  while(pen_getLine(fspec,1000,line,auxRead) == 0){
    nread += auxRead;

    //Get and check next energy
    double E0;
    int nscan = sscanf(line,"%lE",&E0);
    if(nscan != 1){
      printf("\nWarning: Unexpected format at line %lu. Skipping.\n"
	     "         '%s'",nread,line);
      continue;
    }

    if(E0 < 50.0){
      printf("\nWarning: Energy too low at line %lu. Skipping.\n"
	     "         '%s'",nread,line);
      continue;
    }

    printf("\n# *******************\n");
    printf("#  E0= %E eV\n",E0);

    //Ensure to include electron significant radiative yield
    for(long int i = static_cast<long int>(constants::NEGP)-1; i >= 0; --i){
      if(exp(mat.EBRY[i])*context.grid.ET[i] < 1.0e-4*E0){
	mat.EABS[PEN_ELECTRON] = std::max(context.grid.ET[i],EABS1);
	break;
      }
    }

    mat.EABS[PEN_POSITRON] =
      std::min(std::max(1.0e-4*E0,5.0e3),mat.EABS[PEN_ELECTRON]);

    printf("#  EABS(el)= %E eV\n",mat.EABS[PEN_ELECTRON]);
    printf("#  EABS(ph)= %E eV\n",mat.EABS[PEN_PHOTON]);
    printf("#  EABS(po)= %E eV\n",mat.EABS[PEN_POSITRON]);

    simulate(context,E0,dumpTime,simTime,tolerance,seed1,seed2,fout);
  }

  fclose(fspec);
  
  printf("Done!\n");

  return 0;
}

bool dumpResults(const unsigned long long nhist,
		 const double EDPT, const double EDPT2,
		 const double ETRT, const double /*ETRT2*/,
		 const double RMU0, const double E0, const double TOL,
		 const long int TSEC, const bool forceFinish,
		 FILE* fout){
  double FNT  = 1.0/static_cast<double>(nhist);
  double F    = ETRT*FNT;
  double FG   = EDPT*FNT;
  double EFG  = FNT*sqrt(fabs(EDPT2-EDPT*EDPT*FNT)); // 1 sigma
  double PERR = 100.0*EFG/FG;

  //Standard output report
  printf("   %.5E  %.5E  %.5E  %.5E  %.1E  %ld\n",
	 static_cast<double>(nhist),FG*RMU0/E0,FG/E0,F/E0,PERR,TSEC);

  if(PERR < TOL || forceFinish){
    //Results file report 
    fprintf(fout,"   %.5E  %.5E  %.5E  %.5E  %.5E  %.5E  %.1E  %09ld\n",
	    E0,RMU0,FG*RMU0/E0,FG/E0,F/E0,FG/F,PERR,TSEC);
    fflush(fout);
    //finish this energy
    return true;
  }
  return false;
}

void simulate(const pen_context& context,
	      const double E0,
	      const double dumpTime,
	      const double simTime,
	      const double tolerance,
	      int& seed1, int& seed2,
	      FILE* fout){

  //Declare random generator
  pen_rand random;
  //Initialize random
  random.setSeeds(seed1,seed2);

  //Get material
  const pen_material& mat = context.readBaseMaterial(0);
  
  unsigned long long ntot = 1000000000;

  //Create particle stacks for electrons, positrons and gamma
  pen_particleStack<pen_particleState> stackE;
  pen_particleStack<pen_particleState> stackP;
  pen_particleStack<pen_state_gPol> stackG;
  
  //Create particle simulations
  pen_betaE betaE(context,stackE,stackG);
  pen_gamma gamma(context,stackE,stackP,stackG);
  pen_betaP betaP(context,stackE,stackG,stackP);
    
  double RMU0 = 1.0/context.range(E0,PEN_PHOTON,0);
  RMU0 = RMU0/mat.RHO;

  printf("# Initial seeds: %d %d\n",seed1,seed2);
  printf("#  mu/rho= %E cm^2/g\n\n",RMU0);
  printf("#       N        mu_en/rho     (1-g)*f         f       unc.(%%)     T(s)\n");

  //Create a common initial history state
  pen_state_gPol genState;
  genState.E = E0;
  genState.WGHT = 1.0;
  genState.MAT = 1;
  genState.IBODY = 0;
  genState.ILB[0] = 1;
  genState.ILB[1] = 0;
  genState.ILB[2] = 0;
  genState.ILB[3] = 0;
  genState.ILB[4] = 0;
  genState.U = 0.0;
  genState.V = 0.0;
  genState.W = 1.0;

  double DSMAX = 1.0e20;

  double EDPT = 0.0;
  double ETRT = 0.0;
  double EDPT2 = 0.0;
  double ETRT2 = 0.0;
    
  //Create timers
  long long int mili2dump = static_cast<long long int>(dumpTime)*1000;
  pen_stopWatch dumpWatch(mili2dump);
  dumpWatch.start();
  long long int mili2finish = static_cast<long long int>(simTime)*1000;
  pen_stopWatch finishWatch(mili2finish);
  finishWatch.start();

  //Get start time
  const auto tstart = std::chrono::steady_clock::now();

  bool breaked = false;
  for(unsigned long long nhist = 1; nhist <= ntot; ++nhist){

    //Clear stacks
    stackE.cleans();
    stackP.cleans();
    stackG.cleans();

    //Copy initial state
    stateCopy(gamma.getState(),genState);
    const pen_state_gPol& state = gamma.readState();

    //Update material and body information
    gamma.updateMat();
    gamma.updateBody();
    
    gamma.START();

    //Interact until get a non delta interaction
    int icol;
    do{
      double ds, de;
      gamma.JUMP(ds,random,DSMAX);
      gamma.KNOCK(de,icol,random);
    }while(icol == GAMMA_DELTA);

    double EDP = E0 - state.E;
    double ETR = EDP;

    pen_state_gPol auxGammaState;
    bool continueSim = true;
    while(continueSim){

      continueSim = false;
	
      //Count all photons
      while(stackG.getNSec() > 0){
	continueSim = true;
	stackG.get(auxGammaState);
	EDP -= auxGammaState.E;
	if(auxGammaState.ILB[0] == 2) ETR -= auxGammaState.E;
      }

      //Count electrons
      while(stackE.getNSec() > 0){
	continueSim = true;
	stackE.get(betaE.getState());
	if(betaE.readState().E >= mat.EABS[PEN_ELECTRON]){
	  //Update material and body information
	  betaE.updateMat();
	  betaE.updateBody();
	  
	  betaE.START();
	  do{
	    double ds, de, X, Y, Z;
	    int icolDummy;
	    betaE.JUMP(ds,random,DSMAX);  //Call Jump
	    betaE.move(ds,de,X,Y,Z,random); //Move tha particle due soft interactions
	    if(betaE.readState().E < mat.EABS[PEN_ELECTRON])
	      break;
	    betaE.KNOCK(de,icolDummy,random);  //Call nock
	  }while(betaE.readState().E >= mat.EABS[PEN_ELECTRON]);
	}
      }

      //Count positrons
      while(stackP.getNSec() > 0){
	continueSim = true;
	stackP.get(betaP.getState());
	if(betaP.readState().E >= mat.EABS[PEN_POSITRON]){
	  //Update material and body information
	  betaP.updateMat();
	  betaP.updateBody();

	  betaP.START();
	  do{
	    double ds, de, X, Y, Z;
	    int icolDummy;
	    betaP.JUMP(ds,random,DSMAX);  //Call Jump
	    betaP.move(ds,de,X,Y,Z,random); //Move tha particle due soft interactions
	    if(betaP.readState().E < mat.EABS[PEN_POSITRON])
	      break;
	    betaP.KNOCK(de,icolDummy,random);  //Call nock
	  }while(betaP.readState().E >= mat.EABS[PEN_POSITRON]);
	}
      }
    }

    //Add history results to global counters
    EDPT += EDP;  EDPT2 += EDP*EDP;
    ETRT += ETR;  ETRT2 += ETR*ETR;
    
    //Get actual time
    const auto tnow = std::chrono::steady_clock::now();
    //Check if is time to dump
    if(dumpWatch.check(tnow)){
      //Time to dump
      long int elaps =
	std::chrono::duration_cast<std::chrono::seconds>(tnow-tstart).count();
	
      if(dumpResults(nhist,EDPT,EDPT2,ETRT,ETRT2,
		     RMU0,E0,tolerance,elaps,false,fout)){
	//Finish this point
	breaked = true;
	break;
      }
      dumpWatch.start();
    }
    //Check if is time to finish
    if(finishWatch.check(tnow)){
      //Time to finish
      long int elaps =
	std::chrono::duration_cast<std::chrono::seconds>(tnow-tstart).count();
	
      if(dumpResults(nhist,EDPT,EDPT2,ETRT,ETRT2,
		     RMU0,E0,tolerance,elaps,true,fout)){
	//Finish this point
	breaked = true;
	break;
      }
    }
      
  }

  if(!breaked){
    //Number of histories reached, report results
    //Get actual time
    const auto tnow = std::chrono::steady_clock::now();
    //Time to finish
    long int elaps =
      std::chrono::duration_cast<std::chrono::seconds>(tnow-tstart).count();
    dumpResults(ntot,EDPT,EDPT2,ETRT,ETRT2,
		RMU0,E0,tolerance,elaps,true,fout);
  }

  random.getSeeds(seed1,seed2);
}
