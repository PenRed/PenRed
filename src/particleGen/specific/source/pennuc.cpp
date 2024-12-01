//
//
//    Copyright (C) 2020-2023 Universitat de València - UV
//    Copyright (C) 2020-2023 Universitat Politècnica de València - UPV
//    Copyright (C) 2024 Vicent Giménez Alventosa
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
//        vicent.gimenez.alventosa@gmail.com  (Vicent Giménez Alventosa)
//        vicente.gimenez@uv.es (Vicente Giménez Gómez)
//        sanolgi@upvnet.upv.es (Sandra Oliver Gil)
//    
// 

// The following routine has been translated and adapted from the original PENNUC
// developed by the CIEMAT, Laboratoire National Henri Becquerel, and the
// Universitat de Barcelona which includes the following license:

/*
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  PENNUC (version 2018)                                               C
C  Copyright 2016-2018                                                 C
C  CIEMAT, Laboratoire National Henri Becquerel, and                   C
C  Universitat de Barcelona                                            C
C                                                                      C
C  Permission to use, copy, modify, distribute and sell this software  C
C  and its documentation for any purpose is hereby granted without     C
C  fee, provided that the above copyright notice appears in all        C
C  copies and that both that copyright notice and this permission      C
C  notice appear in all supporting documentation. The Ciemat and the   C
C  Universitat de Barcelona make no representations about the suita-   C
C  bility of this software for any purpose. It is provided "as is"     C
C  without express or implied warranty.                                C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
*/

//*******
//** The description of each routine has been extracted from the original code
//*******

#include "pennuc.hh"

int pennuc_specificSampler::configure(double& Emax,
				      const pen_parserSection& config,
				      const unsigned /*nthreads*/,
				      const unsigned verbose){

  //  ****  Read nuclear decay data from the NUCLEIDE file and
  //        initialise the random sampling.

  int IER;

  //Read the filename of the pennuc log file.
  std::string PennuclogFilename ("pennuc.dat");
  int err = config.read("pennuclog_filename",PennuclogFilename);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("pennuc_specificSampler:configure: Error: Unable to read "
       "'pennuclog_filename' field from the configuration.\n "
       "The default filename %s will be used.\n",PennuclogFilename.c_str());
    }
  }
  else
  {
    if(verbose > 0){
      printf("pennuc_specificSampler:configure: 'pennuclog_filename' field read from the configuration.\n "
       "pennuclog_filename = %s.\n",PennuclogFilename.c_str());
    }
  }

  FILE* fout = nullptr;
  char straux[200];
  strcpy (straux, PennuclogFilename.c_str());
  //fout = fopen("pennuc.dat", "w");
  fout = fopen(straux, "w");
  if(fout == nullptr){
    if(verbose > 0){
      printf("pennuc_specificSampler:configure: Error: Unable "
	     "to open output file %s\n", PennuclogFilename.c_str());
    }
    return -1;
  }

  if(spatial() == nullptr){
    if(verbose > 0){
      printf("pennuc_specificSampler:configure: Error: Expected "
	     "spatial sampler not provided\n");
    }
    return -2;
  }
  
  //Read the filename of the NUCLEIDE data file of the considered radionuclide
  std::string nucleideFilename;
  err = config.read("nucleide_filename",nucleideFilename);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("pennuc_specificSampler:configure: Error: Unable to read "
	     "'nucleide_filename' field from the configuration. "
	     "String expected\n");
    }
    return -3;
  }

  //Read the filename of the atomic configuration and subshell binding energies data file.
  std::string atomicFilename ("pdatconf.p14");
  err = config.read("atomic_filename",atomicFilename);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("pennuc_specificSampler:configure: Error: Unable to read "
       "'atomic_filename' field from the configuration.\n "
       "The default filename %s will be used.\n",atomicFilename.c_str());
    }
  }
  else
  {
    if(verbose > 0){
      printf("pennuc_specificSampler:configure: 'atomic_filename' field read from the configuration.\n "
       "atomic_filename = %s.\n",atomicFilename.c_str());
    }
  }

  //Read the filename of the singly ionised atoms with the initial vacancy in one of the K, L,
  //  M and N shells data file.
  std::string relaxationFilename ("pdrelax.p11");
  err = config.read("relaxation_filename",relaxationFilename);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("pennuc_specificSampler:configure: Error: Unable to read "
       "'atomic_filename' field from the configuration.\n "
       "The default filename %s will be used.\n",relaxationFilename.c_str());
    }
  }
  else
  {
    if(verbose > 0){
      printf("pennuc_specificSampler:configure: 'relaxation_filename' field read from the configuration.\n "
       "relaxation_filename = %s.\n",relaxationFilename.c_str());
    }
  }

  //Read the filename of the atreli log file.
  std::string AtrelilogFilename ("atreli.dat");
  err = config.read("atrelilog_filename",AtrelilogFilename);
  if(err != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("pennuc_specificSampler:configure: Error: Unable to read "
       "'atrelilog_filename' field from the configuration.\n "
       "The default filename %s will be used.\n",AtrelilogFilename.c_str());
    }
  }
  else
  {
    if(verbose > 0){
      printf("pennuc_specificSampler:configure: 'atrelilog_filename' field read from the configuration.\n "
       "atrelilog_filename = %s.\n",AtrelilogFilename.c_str());
    }
  }

  sourceMaterial = 0;
  int auxInt;
  err = config.read("source-material",auxInt);
  if(err == INTDATA_SUCCESS){

    if(auxInt <= 0){
      if(verbose > 0){
        printf("pennuc_specificSampler:configure: Error: Material source index"
               " must be greater than 0.\n");
      }
      return -4;
    }

    sourceMaterial = auxInt;

    if(verbose > 1){
      printf("Material %u selected as source.\n",sourceMaterial);
    }
  }


  // ** DRTIME
  DRTIME = 5.0E-6;
  /*
  //Read the detector resolution time
  if(config.read("drtime",DRTIME) != INTDATA_SUCCESS){
    DRTIME = 5.0E-6;
    if(verbose > 1){
      printf(" Detector resolution time not specified. "
	     "A default value of %E s will be used. \n",DRTIME);
    }
  }
  else if(verbose > 1){
    printf(" Detector resolution time set to %E seconds\n",DRTIME);
  }
  */
  
  // ** ECNUC
  //Read the cutoff energy for the de-excitation of inner subshells
  if(config.read("ecut",ECNUC) != INTDATA_SUCCESS){
    ECNUC = 200.0E0;
    if(verbose > 1){
      printf(" Cutoff energy for the de-excitation of inner subshells "
	     "not specified. A default value of %E eV will be used. \n",ECNUC);
    }
  }
  else if(verbose > 1){
    printf(" Cutoff energy for the de-excitation of inner subshells "
	   "set to %E eV\n",ECNUC);
  }

  // ** AGE
  //Read the 'age' parameter to enable or disable particle age recording
  if(config.read("age",LAGE) != INTDATA_SUCCESS){
    if(verbose > 1){
      printf(" Particle age disabled\n");
    }
    LAGE = false;
  }
  else if(verbose > 1){
    printf(" Particle age %s\n",LAGE ? "enabled" : "disabled");
  }

  // ** ISOT
  ISOT_ = 1;
  /*
  //Read the isotope to simulate
  if(config.read("isot",ISOT_) != INTDATA_SUCCESS){
    ISOT_ = 1;
    if(verbose > 1){
      printf(" Active radioisotope not selected. "
	     "The first one will be used.\n");
    }
  }
  */
  if(ISOT_ < 1 || ISOT_ > NISM){
    //Invalid radioisotope number
    if(verbose > 0){
      printf(" pennuc_specificSampler:configure: Error: "
	     "Invalid radioisotope number selected (%d).\n",ISOT_);
    }
    return -4;
  }
  
  double auxEmax;
  NIR=0;
  NR = 0; //No particles in the storage
  lastMETAST = 0; //Nucleus is not in a metastable level 
  PENNUC0(nucleideFilename.c_str(),atomicFilename.c_str(),relaxationFilename.c_str(),AtrelilogFilename.c_str(),auxEmax,fout,IER);
  fclose(fout);
  Emax = 1.001E0*auxEmax;
  if(IER != 0)
    {
      if(verbose > 0){
	printf("pennuc_specificSampler:configure: SOURCE0: PENNUC "
	       "initialisation failed.\n");
      }
      return -4;
    }

  NR = 0; //No particles in the storage
  lastMETAST = 0; //Nucleus is not in a metastable level 
  
  return 0;
}

bool pennuc_specificSampler::getNext(pen_particleState& state,
				     pen_KPAR& genKpar){

  //Avoid alpha
  while(NR > 0 && KPNR[NR-1] == 4){
    --NR;
  }

  //Check if a non alpha particle remains in the buffer
  if(NR > 0){

    //Get the particle type
    switch(KPNR[NR - 1]){
    case 1:
      genKpar = PEN_ELECTRON;
      break;
    case 2:
      genKpar = PEN_PHOTON;
      break;
    case 3:
      genKpar = PEN_POSITRON;
      break;
    case 4:
      --NR; //Remove the particle
      return false;
      break;
    default:
      --NR; //Remove the particle
      return false;
    }
      
    //Check the energy
    if(ENR[NR - 1] > 50.0E0){
      //Set the particle position to decay position
      state.X = x0;
      state.Y = y0;
      state.Z = z0;

      state.E = ENR[NR - 1];
      state.U = UNR[NR - 1];
      state.V = VNR[NR - 1];
      state.W = WNR[NR - 1];

      state.ILB[0] = 1;
      state.ILB[1] = 0;
      state.ILB[2] = 0;
      state.ILB[3] = ITNR[NR - 1];
      state.ILB[4] = 0;

      if(LAGE){
	state.LAGE = true;
	state.PAGE = AGENR[NR - 1];
      }

      state.WGHT = 1.0;

      //Remove one particle from the stack
      --NR;
      return true;
    }

    //Invalid energy
    --NR;
    return false;
  }
  //No more particles in the stack
  return false;
}

void pennuc_specificSampler::sample(pen_particleState& state,
				    pen_KPAR& genKpar,
				    unsigned long long& dhist,
				    pen_rand& random){

  //Try to extract the next particle
  while(NR > 0){
    if(getNext(state,genKpar)){
      //Particle extracted successfully
      dhist = 0;
      return;
    }
  }

  //No more particles remaining in the stack.
  //Generate a new cascade or continue the previous
  //if it is in a metastable state 

  // ** Spatial sampling
  if(lastMETAST == 0){ //Nucleus is not in a metastable level
    //Sample a new position
    spatial()->sample(state,random);

    //Check for specified material source
    if(sourceMaterial > 0){
      geometry->locate(state);
      while(state.MAT != sourceMaterial){
         spatial()->sample(state,random);
         geometry->locate(state);
      }
    }

    //Save decay position
    x0 = state.X;
    y0 = state.Y;
    z0 = state.Z;
    //Flag a new history
    dhist = 1;
  }else{
    //The same history remains
    dhist = 0;
  }

  //Call pennuc until a valid cascade is generated
  for(size_t i = 0; i < 1000; ++i){
    // ** Pennuc call
    int ISOT = ISOT_;
    int IER;
    PENNUC(ISOT,lastMETAST,IER,random);

    if(IER != 0){
      printf ("The intended isotope, ISOT =%3d is not loaded.\n", ISOT);
      printf ("PENNUC sampling failed.\n");
      genKpar = ALWAYS_AT_END;
      return;
    }

    //Extract the first valid particle
    while(NR > 0){
      if(getNext(state,genKpar)){
	//Particle extracted successfully
	return;
      }
    }

    //Unable to get a valid particle in this cascade.
    //Increase the history number if the nucleus is not
    //in a metastable level and try again
    if(lastMETAST == 0)
      ++dhist;
  }

  //Unable to extract any valid particle after 1000 cascades
  genKpar = ALWAYS_AT_END;
  printf ("No valid particle sampled after 1000 PENNUC calls.\n");
  printf ("PENNUC sampling failed.\n");
  
}

//  *********************************************************************
//                       SUBROUTINE PENNUC
//  *********************************************************************
void pennuc_specificSampler::PENNUC(int& IS,
				    int& METAST,
				    int& IER,
				    pen_rand& random){
  /*
    C  This subroutine generates random decay pathways of the radioactive
    C  nuclide IS using information from the NUCLEIDE database.
    C
    C  Metastable levels (i.e., levels with life times longer than the
    C  detector resolution time) are treated as effective halts of the decay
    C  process. When a metastable level is reached, the output value of
    C  METAST is set to 1 and control is returned to the main program. At
    C  the next call to subroutine PENNUC, the simulation of the decay
    C  continues down to the following metastable level or to the ground
    C  level. The output value METAST=0 indicates that the nucleus has
    C  reached its ground level. This operation scheme allows simulating the
    C  response of real detectors, which may give several output counts as
    C  the result of the decay of a single nucleus through a metastable
    C  state.
    C
    C  Each call to PENNUC returns the initial states of the set of
    C  particles released in a cascade. The initial state of each particle
    C  is defined by its type (KPNR), energy (ENR), age (AGENR), and
    C  direction cosines (UNR,VNR,WNR). The quantity ITNR identifies the
    C  transition that originated the particle. The number NR of emitted
    C  particles and their initial-state variables are returned through the
    C  module PENNUC_mod.
    C
    C  The transition identifier, ITNR, is defined as follows:
    C  -- For atomic transitions ISA-IS1 or ISA-IS1-IS2 of an ion of atomic
    C     number IZ,
    C        ITNR=IZ*1000000+ISA*10000+IS1*100+IS2
    C     where ISi is the label code of the ISi subshell (see the next
    C     table). For radiative transitions, IS2=0.
    C
    C  ---------------------------------------------------------------------
    C  Label code IS for electron subshells:
    C      1 = K  (1s1/2),     11 = N2 (4p1/2),     21 = O5 (5d5/2),
    C      2 = L1 (2s1/2),     12 = N3 (4p3/2),     22 = O6 (5f5/2),
    C      3 = L2 (2p1/2),     13 = N4 (4d3/2),     23 = O7 (5f7/2),
    C      4 = L3 (2p3/2),     14 = N5 (4d5/2),     24 = P1 (6s1/2),
    C      5 = M1 (3s1/2),     15 = N6 (4f5/2),     25 = P2 (6p1/2),
    C      6 = M2 (3p1/2),     16 = N7 (4f7/2),     26 = P3 (6p3/2),
    C      7 = M3 (3p3/2),     17 = O1 (5s1/2),     27 = P4 (6d3/2),
    C      8 = M4 (3d3/2),     18 = O2 (5p1/2),     28 = P5 (6d5/2),
    C      9 = M5 (3d5/2),     19 = O3 (5p3/2),     29 = Q1 (7s1/2),
    C     10 = N1 (4s1/2),     20 = O4 (5d3/2),     30 = outer shells.
    C  ---------------------------------------------------------------------
    C
    C  -- For nuclear disintegrations and de-excitations,
    C        ITNR=-(IDAUGH*1000000+IBRANCH*1000+ITR)
    C     where IDAUGH and IBRANCH are the indices of the daughter and the
    C     branch, and ITR is the type of nuclear transition that released
    C     the particle.
    C
    C        ITR= 0, alpha disintegration,
    C             1, beta- disintegration,
    C             2, beta+ disintegration,
    C             3, gamma emission,
    C             4, internal conversion in the K shell,
    C             5, internal conversion in the L1 subshell,
    C             6, internal conversion in the L2 subshell,
    C             7, internal conversion in the L3 subshell,
    C             8, internal conversion in the M, N and outer subshells.
    C
    C
    C  Input argument:
    C    IS ........... active radioisotope.
    C    METAST ....... flag for metastable levels. Set METAST=0 to start
    C                   the simulation of the decay of a new isotope, and
    C                   keep the output values unchanged.
    C
    C  Output arguments:
    C    METAST ....... output flag for metastable levels,
    C                   = 0, the de-excitation cascade has been completed,
    C                   = 1, a metastable state has been reached.
    C    IER  ......... error code. When IER.NE.0 the isotope IS is not
    C                   defined.
    C  */

  // Uses PENNUC_mod, CPN1, CPN2, CPN3, CPN4 and CPN5 commons;

  //  ----  Minimum lifetime for metastable levels.
  const double TIMELIM = DRTIME * 1.0E-2;

  //  ****  Output of ATREL.
  const int NEMP = 500;
  double ES[NEMP], PAGES[NEMP];
  int KPARS[NEMP], ITR[NEMP];
  //
  int& IDAUGH = lastIDAUGH;
  int& IBRANCH = lastIBRANCH;
  int& ISA = lastISA;
  //
  //      EXTERNAL RAND
  //
  NR = 0;
  if (IS < 1 || IS > NIR)
    {
      printf ("NIR =%d\n", NIR);
      printf (" The intended isotope, IS =%3d is not loaded.\n", IS);
      IER = 1;
      LEVD = 9999;
      return;
    }
  else
    {
      IER = 0;
    }
  double DPAGE = 0.0E0;
  //
  //  ****  Sample the daughter nuclide.
  //
  if (METAST != 1)
    {
      ISA = IS;
      METAST = 0;		// Just for protection.
      //
      int NDHS = NDAUGH[IS - 1];
      if (NDHS == 1)
	{
	  IDAUGH = 1;
	}
      else			// Binary search.
	{
	  double RD = random.rand();
	  if (RD < PDAUGH[IS - 1][0])
	    {
	      IDAUGH = 1;
	    }
	  else
	    {
	      IDAUGH = 1;
	      int I1 = NDHS;
	      bool Exit = false;
	      while (!Exit)
		{
		  Exit = true;
		  int IT = (IDAUGH + I1) / 2;
		  if (RD > PDAUGH[IS - 1][IT - 1])
		    {
		      IDAUGH = IT;
		    }
		  else
		    {
		      I1 = IT;
		    }
		  if (I1 - IDAUGH > 1)
		    {
		      Exit = false;
		      continue;
		    }
		}
	      IDAUGH = IDAUGH + 1;
	    }
	}
      //
      //  ****  Select the disintegration branch to the daughter nuclide.
      //
      int NBRS = NBRANCH[IS - 1][IDAUGH - 1];
      if (NBRS == 1)
	{
	  IBRANCH = 1;
	}
      else
	{
	  double RD = random.rand();
	  if (RD < BR[IS - 1][IDAUGH - 1][0])
	    {
	      IBRANCH = 1;
	    }
	  else			// Binary search.
	    {
	      IBRANCH = 1;
	      int I1 = NBRS;
	      bool Exit = false;
	      while (!Exit)
		{
		  Exit = true;
		  int IT = (IBRANCH + I1) / 2;
		  if (RD > BR[IS - 1][IDAUGH - 1][IT - 1])
		    {
		      IBRANCH = IT;
		    }
		  else
		    {
		      I1 = IT;
		    }
		  if (I1 - IBRANCH > 1)
		    {
		      Exit = false;
		      continue;
		    }
		}
	      IBRANCH = IBRANCH + 1;
	    }
	}
      //
      //  ****  Simulate the disintegration branch to the daughter nuclide.
      //
      int MODEBR = KTYPE[IS - 1][IDAUGH - 1][IBRANCH - 1];
      IAD = IANUC[IS - 1];
      //
      //  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      switch (MODEBR)
	{
	  //  ----  Alpha emission.
	case 1:
	  IZD = IZNUC[IS - 1] - 2;
	  IAD = IAD - 4;
	  NR = NR + 1;		// Released particle.
	  ENR[NR - 1] = EBRANCH[IS - 1][IDAUGH - 1][IBRANCH - 1];
	  ISOTROP(UNR[NR - 1], VNR[NR - 1], WNR[NR - 1],random);
	  AGENR[NR - 1] = DPAGE;
	  KPNR[NR - 1] = 4;
	  ITNR[NR - 1] = -(IDAUGH*1000000+IBRANCH*1000+0);
	  break;
	  //  ----  Beta- emission, energy sampled from the beta spectrum.
	case 2:
	  {
	    IZD = IZNUC[IS - 1] + 1;
	    IDT = INDBETA[IS - 1][IDAUGH - 1][IBRANCH - 1];
	    double EBETAS = SBETAS(IDT,random);
	    NR = NR + 1;	// Released particle.
	    ENR[NR - 1] = EBETAS;
	    ISOTROP (UNR[NR - 1], VNR[NR - 1], WNR[NR - 1],random);
	    AGENR[NR - 1] = DPAGE;
	    KPNR[NR - 1] = 1;
	    ITNR[NR - 1] = -(IDAUGH*1000000+IBRANCH*1000+1);
	    break;
	  }
	  //  ----  Beta+ emission, energy sampled from the beta spectrum.
	case 3:
	  {
	    IZD = IZNUC[IS - 1] - 1;
	    IDT = INDBETA[IS - 1][IDAUGH - 1][IBRANCH - 1];
	    double EBETAS = SBETAS(IDT,random);
	    NR = NR + 1;	// Released particle.
	    ENR[NR - 1] = EBETAS;
	    ISOTROP (UNR[NR - 1], VNR[NR - 1], WNR[NR - 1],random);
	    AGENR[NR - 1] = DPAGE;
	    KPNR[NR - 1] = 3;
	    ITNR[NR - 1] = -(IDAUGH*1000000+IBRANCH*1000+2);
	    break;
	  }
	  //  ----  Electron capture K.
	case 4:
	  {
	    IZD = IZNUC[IS - 1] - 1;
	    int NS;
	    ATREL(IZD, 1, ES, PAGES, KPARS, ITR, NS, random);
	    if (NS > 0)
	      {
		double PAGES1 = 1.0E35;
		for (int I = 0; I < NS; I++)
		  {
		    PAGES1 = (PAGES1 < PAGES[I] ? PAGES1 : PAGES[I]);
		  }
		for (int I = 0; I < NS; I++)
		  {
		    NR = NR + 1;	// Released particle.
		    ENR[NR - 1] = ES[I];
		    ISOTROP (UNR[NR - 1], VNR[NR - 1], WNR[NR - 1],random);
		    AGENR[NR - 1] = DPAGE + PAGES[I] - PAGES1;
		    KPNR[NR - 1] = KPARS[I];
		    ITNR[NR - 1] = ITR[I];
		  }
	      }
	    break;
	  }
	  //  ----  Electron capture L1.
	case 5:
	  {
	    IZD = IZNUC[IS - 1] - 1;
	    int NS;
	    ATREL(IZD, 2, ES, PAGES, KPARS, ITR, NS, random);
	    if (NS > 0)
	      {
		double PAGES1 = 1.0E35;
		for (int I = 0; I < NS; I++)
		  {
		    PAGES1 = (PAGES1 < PAGES[I] ? PAGES1 : PAGES[I]);
		  }
		for (int I = 0; I < NS; I++)
		  {
		    NR = NR + 1;	// Released particle.
		    ENR[NR - 1] = ES[I];
		    ISOTROP (UNR[NR - 1], VNR[NR - 1], WNR[NR - 1],random);
		    AGENR[NR - 1] = DPAGE + PAGES[I] - PAGES1;
		    KPNR[NR - 1] = KPARS[I];
		    ITNR[NR - 1] = ITR[I];
		  }
	      }
	    break;
	  }
	  //  ----  Electron capture L2.
	case 6:
	  {
	    IZD = IZNUC[IS - 1] - 1;
	    int NS;
	    ATREL(IZD, 3, ES, PAGES, KPARS, ITR, NS, random);
	    if (NS > 0)
	      {
		double PAGES1 = 1.0E35;
		for (int I = 0; I < NS; I++)
		  {
		    PAGES1 = (PAGES1 < PAGES[I] ? PAGES1 : PAGES[I]);
		  }
		for (int I = 0; I < NS; I++)
		  {
		    NR = NR + 1;	// Released particle.
		    ENR[NR - 1] = ES[I];
		    ISOTROP (UNR[NR - 1], VNR[NR - 1], WNR[NR - 1],random);
		    AGENR[NR - 1] = DPAGE + PAGES[I] - PAGES1;
		    KPNR[NR - 1] = KPARS[I];
		    ITNR[NR - 1] = ITR[I];
		  }
	      }
	    break;
	  }
	  //  ----  Electron capture L3.
	case 7:
	  {
	    IZD = IZNUC[IS - 1] - 1;
	    int NS;
	    ATREL(IZD, 4, ES, PAGES, KPARS, ITR, NS, random);
	    if (NS > 0)
	      {
		double PAGES1 = 1.0E35;
		for (int I = 0; I < NS; I++)
		  {
		    PAGES1 = (PAGES1 < PAGES[I] ? PAGES1 : PAGES[I]);
		  }
		for (int I = 0; I < NS; I++)
		  {
		    NR = NR + 1;	// Released particle.
		    ENR[NR - 1] = ES[I];
		    ISOTROP (UNR[NR - 1], VNR[NR - 1], WNR[NR - 1],random);
		    AGENR[NR - 1] = DPAGE + PAGES[I] - PAGES1;
		    KPNR[NR - 1] = KPARS[I];
		    ITNR[NR - 1] = ITR[I];
		  }
	      }
	    break;
	  }
	  //  ----  Electron capture M.
	case 8:
	  {
	    IZD = IZNUC[IS - 1] - 1;
	    int NS;
	    ATREL(IZD, 5, ES, PAGES, KPARS, ITR, NS, random);
	    if (NS > 0)
	      {
		double PAGES1 = 1.0E35;
		for (int I = 0; I < NS; I++)
		  {
		    PAGES1 = (PAGES1 < PAGES[I] ? PAGES1 : PAGES[I]);
		  }
		for (int I = 0; I < NS; I++)
		  {
		    NR = NR + 1;	// Released particle.
		    ENR[NR - 1] = ES[I];
		    ISOTROP (UNR[NR - 1], VNR[NR - 1], WNR[NR - 1],random);
		    AGENR[NR - 1] = DPAGE + PAGES[I] - PAGES1;
		    KPNR[NR - 1] = KPARS[I];
		    ITNR[NR - 1] = ITR[I];
		  }
	      }
	    break;
	  }
	  //  ----  Metastable level of parent nuclide.
	case 9:
	  IZD = IZNUC[IS - 1];
	  IAD = IANUC[IS - 1];
	  METAST = 1;
	  break;
	}
      //  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      //
      //  ****  Level fed by the desintegration.
      LEVD = IFLB[IS - 1][IDAUGH - 1][IBRANCH - 1];
      if (LEVD == 0)
	{
	  return;
	}
    }
  else
    {
      IS = ISA;
    }
  //
  //  ****  Continue the de-excitation cascade of the daughter nucleus
  //        until the ground level, or until the next metastable level.
  //
  while (LEVD > 0)
    {
      //
      if (METAST == 1)
	{
	  METAST = 0;		// This is the start of a new cascade.
	}
      else
	{
	  //  ----  Sample the de-excitation time of the level, except for levels
	  //        with extremely short life times.
	  if (DTIME[IS - 1][IDAUGH - 1][LEVD - 1] > TIMELIM)
	    {
	      double TAU = -log (random.rand()) * 1.4426950408889634E0
		* DTIME[IS - 1][IDAUGH - 1][LEVD - 1];
	      if (TAU > DRTIME)	// The level is metastable.
		{
		  METAST = 1;
		  return;
		}
	      else
		{
		  DPAGE = DPAGE + TAU;
		  METAST = 0;
		}
	    }
	  else
	    {
	      METAST = 0;
	    }
	}
      //
      //  ****  Now, select the transition that leaves the level.
      //
      int ITRANS;
      int NTRS = NTRANS[IS - 1][IDAUGH - 1][LEVD - 1];
      if (NTRS == 1)
	{
	  ITRANS = 1;
	}
      else			// Binary search.
	{
	  double RD = random.rand();
	  if (RD < PTRANS[IS - 1][IDAUGH - 1][LEVD - 1][0])
	    {
	      ITRANS = 1;
	    }
	  else
	    {
	      ITRANS = 1;
	      int I1 = NTRS;
	      bool Exit = false;
	      while (!Exit)
		{
		  Exit = true;
		  int IT = (ITRANS + I1) / 2;
		  if (RD > PTRANS[IS - 1][IDAUGH - 1][LEVD - 1][IT - 1])
		    {
		      ITRANS = IT;
		    }
		  else
		    {
		      I1 = IT;
		    }
		  if (I1 - ITRANS > 1)
		    {
		      Exit = false;
		      continue;
		    }
		}
	      ITRANS = ITRANS + 1;
	    }
	}
      //  ----  A transition has been determined.
      double ENERTR = ETRANS[IS - 1][IDAUGH - 1][LEVD - 1][ITRANS - 1];
      int MODETR = MODE[IS - 1][IDAUGH - 1][LEVD - 1][ITRANS - 1];
      LEVD = IFLT[IS - 1][IDAUGH - 1][LEVD - 1][ITRANS - 1];
      //
      //  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
      switch (MODETR)
	{
	  //  ----  Gamma emission.
	case 1:
	  NR = NR + 1;		// Released particle.
	  ENR[NR - 1] = ENERTR;
	  ISOTROP (UNR[NR - 1], VNR[NR - 1], WNR[NR - 1],random);
	  AGENR[NR - 1] = DPAGE;
	  KPNR[NR - 1] = 2;
	  ITNR[NR - 1] = -(IDAUGH*1000000+IBRANCH*1000+3);
	  break;
	  //  ----  K conversion electron.
	case 2:
	  {
	    NR = NR + 1;	// Released particle.
	    ENR[NR - 1] = ENERTR;
	    ISOTROP (UNR[NR - 1], VNR[NR - 1], WNR[NR - 1],random);
	    AGENR[NR - 1] = DPAGE;
	    KPNR[NR - 1] = 1;
	    ITNR[NR - 1] = -(IDAUGH*1000000+IBRANCH*1000+4);
	    int NS;
	    ATREL(IZD, 1, ES, PAGES, KPARS, ITR, NS, random);
	    if (NS > 0)
	      {
		for (int I = 0; I < NS; I++)
		  {
		    NR = NR + 1;	// Released particle.
		    ENR[NR - 1] = ES[I];
		    ISOTROP (UNR[NR - 1], VNR[NR - 1], WNR[NR - 1],random);
		    AGENR[NR - 1] = DPAGE + PAGES[I];
		    KPNR[NR - 1] = KPARS[I];
		    ITNR[NR - 1] = ITR[I];
		  }
	      }
	    break;
	  }
	  //  ----  L1 conversion electron.
	case 3:
	  {
	    NR = NR + 1;	// Released particle.
	    ENR[NR - 1] = ENERTR;
	    ISOTROP (UNR[NR - 1], VNR[NR - 1], WNR[NR - 1],random);
	    AGENR[NR - 1] = DPAGE;
	    KPNR[NR - 1] = 1;
	    ITNR[NR - 1] = -(IDAUGH*1000000+IBRANCH*1000+5);
	    int NS;
	    ATREL(IZD, 2, ES, PAGES, KPARS, ITR, NS, random);
	    if (NS > 0)
	      {
		for (int I = 0; I < NS; I++)
		  {
		    NR = NR + 1;	// Released particle.
		    ENR[NR - 1] = ES[I];
		    ISOTROP (UNR[NR - 1], VNR[NR - 1], WNR[NR - 1],random);
		    AGENR[NR - 1] = DPAGE + PAGES[I];
		    KPNR[NR - 1] = KPARS[I];
		    ITNR[NR - 1] = ITR[I];
		  }
	      }
	    break;
	  }
	  //  ----  L2 conversion electron.
	case 4:
	  {
	    NR = NR + 1;	// Released particle.
	    ENR[NR - 1] = ENERTR;
	    ISOTROP (UNR[NR - 1], VNR[NR - 1], WNR[NR - 1],random);
	    AGENR[NR - 1] = DPAGE;
	    KPNR[NR - 1] = 1;
	    ITNR[NR - 1] = -(IDAUGH*1000000+IBRANCH*1000+6);
	    int NS;
	    ATREL(IZD, 3, ES, PAGES, KPARS, ITR, NS, random);
	    if (NS > 0)
	      {
		for (int I = 0; I < NS; I++)
		  {
		    NR = NR + 1;	// Released particle.
		    ENR[NR - 1] = ES[I];
		    ISOTROP (UNR[NR - 1], VNR[NR - 1], WNR[NR - 1],random);
		    AGENR[NR - 1] = DPAGE + PAGES[I];
		    KPNR[NR - 1] = KPARS[I];
		    ITNR[NR - 1] = ITR[I];
		  }
	      }
	    break;
	  }
	  //  ----  L3 conversion electron.
	case 5:
	  {
	    NR = NR + 1;	// Released particle.
	    ENR[NR - 1] = ENERTR;
	    ISOTROP (UNR[NR - 1], VNR[NR - 1], WNR[NR - 1],random);
	    AGENR[NR - 1] = DPAGE;
	    KPNR[NR - 1] = 1;
	    ITNR[NR - 1] = -(IDAUGH*1000000+IBRANCH*1000+7);
	    int NS;
	    ATREL(IZD, 4, ES, PAGES, KPARS, ITR, NS, random);
	    if (NS > 0)
	      {
		for (int I = 0; I < NS; I++)
		  {
		    NR = NR + 1;	// Released particle.
		    ENR[NR - 1] = ES[I];
		    ISOTROP (UNR[NR - 1], VNR[NR - 1], WNR[NR - 1],random);
		    AGENR[NR - 1] = DPAGE + PAGES[I];
		    KPNR[NR - 1] = KPARS[I];
		    ITNR[NR - 1] = ITR[I];
		  }
	      }
	    break;
	  }
	  //  ----  M conversion electron.
	case 6:
	  {
	    NR = NR + 1;	// Released particle.
	    ENR[NR - 1] = ENERTR;
	    ISOTROP (UNR[NR - 1], VNR[NR - 1], WNR[NR - 1],random);
	    AGENR[NR - 1] = DPAGE;
	    KPNR[NR - 1] = 1;
	    ITNR[NR - 1] = -(IDAUGH*1000000+IBRANCH*1000+8);
	    int NS;
	    ATREL(IZD, 5, ES, PAGES, KPARS, ITR, NS, random);
	    if (NS > 0)
	      {
		for (int I = 0; I < NS; I++)
		  {
		    NR = NR + 1;	// Released particle.
		    ENR[NR - 1] = ES[I];
		    ISOTROP (UNR[NR - 1], VNR[NR - 1], WNR[NR - 1],random);
		    AGENR[NR - 1] = DPAGE + PAGES[I];
		    KPNR[NR - 1] = KPARS[I];
		    ITNR[NR - 1] = ITR[I];
		  }
	      }
	    break;
	  }
	  //  ----  Undefined transition. Assumed gamma to ground state.
	case 7:
	  NR = NR + 1;		// Released particle.
	  ENR[NR - 1] = ENERTR;
	  ISOTROP (UNR[NR - 1], VNR[NR - 1], WNR[NR - 1],random);
	  AGENR[NR - 1] = DPAGE;
	  KPNR[NR - 1] = 2;
	  ITNR[NR - 1] = 0;
	  break;
	}
      //  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      //
    }
  //
  return;
}

//  *********************************************************************
//                       SUBROUTINE ISOTROP
//  *********************************************************************
void pennuc_specificSampler::ISOTROP(double &U,
				     double &V,
				     double &W,
				     pen_rand& random) const{
  //  Gives the direction cosines of a random isotropic direction.
  const double TWOPI = 2.0E0 * constants::PI;
  W = 2.0E0 * random.rand() - 1.0E0;
  double UV = sqrt (1.0E0 - W * W);
  double PHI = random.rand() * TWOPI;
  U = UV * cos (PHI);
  V = UV * sin (PHI);
}

//  *********************************************************************
//                       SUBROUTINE PENNUC0
//  *********************************************************************
void pennuc_specificSampler::PENNUC0(const char* NUCFNAME,
             const char* ATOMFNAME, const char* RELAXFNAME,
             const char* ATRELIFNAME,
				     double& EPMAX,
				     FILE* IWR,
				     int &IER){
  //
  //  Initialisation of decay simulation for the radionuclide defined by
  //  the NUCLEIDE data file NUCFNAME.
  //
  //  PENNUC can handle up to NISM different nuclides simultaneously. To
  //  simulate decays of several nuclides, call subroutine PENNUC0 for each
  //  of them. Each radionuclide is assigned a label IS defined by the
  //  order in which they are initialised. That is, IS=1, 2, ... for the
  //  first, second, ... loaded nuclides. The sampling of the active
  //  nuclide (according to the characteristics of the radioactive source)
  //  must be performed in the calling program.
  //
  //  Input arguments:
  //    NUCFNAME ... filename of the NUCLEIDE data file of the considered
  //                 radionuclide.
  //    IWR ........ input decay data are summarised in an output file,
  //                 which must be open in the calling program as UNIT=IWR
  //                 before the first call to PENNUC0.
  //  Output values:
  //    EPMAX ...... highest energy of photons and electrons released in
  //                 the decay of the radionuclide (and of any other that
  //                 has been previously loaded).
  //    IER  ....... error code. When IER.NE.0 the initialisation has not
  //                 succeeded; the cause of the problem is written in
  //                 the standard output unit (6).
  //
  //  ---------------------------------------------------------------------
  //  Quantities stored in memory:
  //     NUCLIDE(I) ... name of the nuclide.
  //     NDAUGH(I) ... number of daughters.
  //
  //  For each daughter ID,
  //     PDAUGH(I,ID) ... relative disintegration probability to ID.
  //     DNUCLIDE(I,ID) ... name of the daughter nuclide.
  //     Q(I,ID) ... total available energy of the disintegration.
  //     NLEVEL(I,ID) ... number of (excited) levels of the daughter
  //        nuclide.
  //     ELEVEL(I,ID,LEV) ... energies of the daughter nuclide levels,
  //        sorted in increasing order of energies. The index LEV of the
  //        ground level is equal to 0.
  //     DTIME(I,ID,LEV) ... half-life of the level LEV.
  //
  //     NBRANCH(I,ID) ... number of disintegration branches.
  //     KTYPE(I,ID,IBR) ... disintegration type (1-9)
  //        1=ALPHA, 2=BETA-, 3=BETA+, 4=ECK, 5=ECL1, 6=ECL2, 7=ECL3,
  //        8=ECLM, 9=META.
  //     IFLB(I,ID,IBR) ... level fed of the daughter nuclide.
  //     IPRBETA(I,ID,IBr) ... prohibition parameter defining the nature
  //        of BETA branches (=0 for ALPHA and EC branches).
  //           IPRBETA=0: Allowed and first forbidden transitions.
  //           IPRBETA=1: Unique-first forbidden and second forbidden.
  //           IPRBETA=2: Unique-second forbidden and third forbidden.
  //           IPRBETA=3: Unique-third forbidden and fourth forbidden.
  //     EBRANCH(I,ID,IBR) ... energy of the disintegration branch
  //        (kinetic energy for alpha branches, maximum kinetic energy of
  //        the electron/positron in beta branches, none for EC and META).
  //     BR(I,ID,IBR) ... branching ratios, sum normalised to unity.
  //
  //     NTRANS(I,ID,LEV) ... number of transitions depopulating the level
  //        LEV, including electron conversion,
  //     PTRANS(I,ID,LEV,IT) ... probability of the transition IT that
  //        depopulates the level LEV.
  //     MODE(I,ID,LEV,IT) ... kind of transition (1-6):
  //        1=GAMMA, 2=ECK, 3=ECL1, 4=ECL2, 5=ECL3, 6=ECM.
  //     ETRANS(I,ID,LEV,IT) ... energy of the transition ITR.
  //     IFLT(I,ID,LEV,IT) ... level fed by the transition IT.
  //
  //     The uncertainty of a quantity QQQQ is named UQQQQ. Uncertainties
  //     are read from file NUCFNAME, but not used in the simulation.
  //  ---------------------------------------------------------------------
  //
  //Uses PENNUC_mod, CPN1, CPN2, CPN3, CPN4 and CPN5 commons
  //
  char TEMPBUF[100], WCODE[4];
  //      CHARACTER NUCFNAME*45
  //
  char WTYPE[9][6];
  strcpy (WTYPE[0], "ALPHA");
  strcpy (WTYPE[1], "BETA-");
  strcpy (WTYPE[2], "BETA+");
  strcpy (WTYPE[3], "ECK  ");
  strcpy (WTYPE[4], "ECL1 ");
  strcpy (WTYPE[5], "ECL2 ");
  strcpy (WTYPE[6], "ECL3 ");
  strcpy (WTYPE[7], "ECM  ");
  strcpy (WTYPE[8], "METAS");	// EC = electron capture.

  char WMODE[7][6];		// IC = internal conversion.
  strcpy (WMODE[0], "GAMMA");
  strcpy (WMODE[1], "ICK  ");
  strcpy (WMODE[2], "ICL1 ");
  strcpy (WMODE[3], "ICL2 ");
  strcpy (WMODE[4], "ICL3 ");
  strcpy (WMODE[5], "ICMN ");
  strcpy (WMODE[6], "NONE ");

  //
  IER = 0;
  NIR = NIR + 1;
  if (NIR > NISM)
    {
      printf
	("\n PENNUC0: Too many radionuclides. NIR =%3d\n Please, increase the value of the parameter NISM.\n",
	 NIR);
      fprintf (IWR,
	       "\n PENNUC0: Too many radionuclides. NIR =%3d\n Please, increase the value of the parameter NISM.\n",
	       NIR);
      IER = 11;
      return;
    }
  for (int IDc = 0; IDc < NDM; IDc++)
    {
      NBRANCH[NIR - 1][IDc] = 0;
      NLEVEL[NIR - 1][IDc] = 0;
      PDAUGH[NIR - 1][IDc] = 0.0E0;
      for (int IT = 0; IT <= NTM; IT++)
	{
	  ELEVEL[NIR - 1][IDc][IT] = 0.0E0;
	}
      for (int IB = 0; IB < NBM; IB++)
	{
	  KTYPE[NIR - 1][IDc][IB] = 0;
	  BR[NIR - 1][IDc][IB] = 0.0E0;
	  EBRANCH[NIR - 1][IDc][IB] = 0.0E0;
	  IFLB[NIR - 1][IDc][IB] = 0;
	  IPRBETA[NIR - 1][IDc][IB] = 0;
	  INDBETA[NIR - 1][IDc][IB] = 0;
	}
      for (int IL = 0; IL < NLM; IL++)
	{
	  NTRANS[NIR - 1][IDc][IL] = 0;
	  DTIME[NIR - 1][IDc][IL] = 0.0E0;
	  for (int IT = 0; IT < NTM; IT++)
	    {
	      PTRANS[NIR - 1][IDc][IL][IT] = 0.0E0;
	      ETRANS[NIR - 1][IDc][IL][IT] = 0.0E0;
	      MODE[NIR - 1][IDc][IL][IT] = 0;
	      IFLT[NIR - 1][IDc][IL][IT] = 0;
	    }
	}
    }
  //
  if (NIR == 1)
    {
      ATREL0 ();		// Initialise the atomic relaxation routines.
      IDT = 0;
    }
  else
    {
      fputc ('\n', IWR);
      fputc ('\n', IWR);
      fputc (' ', IWR);
      for (int i = 0; i < 76; i++)
	fputc ('=', IWR);
      fputc ('\n', IWR);
    }
  EPMAX = 2.0E5;
  //
  FILE *NRD;
  char straux[200];
  strcpy (straux, NUCFNAME);
  NRD = fopen (straux, "r");
  if (NRD == NULL)
    {
      fprintf (IWR,
	       "\n PENNUC0: The program could not open the file %s\n Please, copy that file into the PENNUC data directory.\n",
	       straux);
      printf
	("\n PENNUC0: The program could not open the file %s\n Please, copy that file into the PENNUC data directory.\n",
	 straux);
      IER = 12;
      return;
    }
  //
  READSKP (NRD, WCODE, TEMPBUF);
  sscanf (TEMPBUF, "%s", NUCLIDE[NIR - 1]);
  READSKP (NRD, WCODE, TEMPBUF);
  sscanf (TEMPBUF, "%d , %d", &IANUC[NIR - 1], &IZNUC[NIR - 1]);
  READSKP (NRD, WCODE, TEMPBUF);
  sscanf (TEMPBUF, "%d", &NDAUGH[NIR - 1]);
  //
  if (NIR > 1)
    {
      int IZZ = IZNUC[NIR - 1];
      int IAA = IANUC[NIR - 1];
      for (int IS = 0; IS < NIR - 1; IS++)
	{
	  if (IZZ == IZNUC[IS] && IAA == IANUC[IS])
	    {
	      fprintf (IWR,
		       "\n PENNUC0: The radioisotope decay file %s\n          has already been loaded.\n",
		       straux);
	      printf
		("\n PENNUC0: The radioisotope decay file %s\n          has already been loaded.\n",
		 straux);
	      IER = 13;
	      return;
	    }
	}
    }
  //
  int IZACT = IZNUC[NIR - 1];
  ATRELI (ATRELIFNAME, ATOMFNAME, RELAXFNAME, IZACT, 1, IER);	// Parent nucleus.
  if (IER != 0)
    {
      return;
    }
  //
  if (IZNUC[NIR - 1] > 2)
    {
      IZACT = IZNUC[NIR - 1] - 2;
      ATRELI (ATRELIFNAME, ATOMFNAME, RELAXFNAME, IZACT, 1, IER);	// Daughter nucleus (alpha).
      if (IER != 0)
	{
	  return;
	}
    }
  //
  if (IZNUC[NIR - 1] > 1)
    {
      IZACT = IZNUC[NIR - 1] - 1;
      ATRELI (ATRELIFNAME, ATOMFNAME, RELAXFNAME, IZACT, 1, IER);	// Daughter (beta+ and el. capture).
      if (IER != 0)
	{
	  return;
	}
    }
  //
  IZACT = IZNUC[NIR - 1] + 1;
  ATRELI (ATRELIFNAME, ATOMFNAME, RELAXFNAME, IZACT, 1, IER);	// Daughter nucleus (beta-).
  if (IER != 0)
    {
      return;
    }
  //
  const char s[2] = ",";
  for (int IDAUGH = 1; IDAUGH <= NDAUGH[NIR - 1]; IDAUGH++)
    {
      READSKP (NRD, WCODE, TEMPBUF);
      sscanf (TEMPBUF, "%s", DNUCLIDE[NIR - 1][IDAUGH - 1]);
      READSKP (NRD, WCODE, TEMPBUF);
      sscanf (strtok (TEMPBUF, s), "%lf", &PDAUGH[NIR - 1][IDAUGH - 1]);
      sscanf (strtok (NULL, s), "%lf", &UPDAUGH[NIR - 1][IDAUGH - 1]);
      sscanf (strtok (NULL, s), "%d", &NLEVEL[NIR - 1][IDAUGH - 1]);
      sscanf (strtok (NULL, s), "%d", &NBRANCH[NIR - 1][IDAUGH - 1]);
      READSKP (NRD, WCODE, TEMPBUF);
      sscanf (TEMPBUF, "%lf , %lf", &Q[NIR - 1][IDAUGH - 1],
	      &UQ[NIR - 1][IDAUGH - 1]);
      Q[NIR - 1][IDAUGH - 1] = Q[NIR - 1][IDAUGH - 1] * 1.0E3;	// keV to eV.
      UQ[NIR - 1][IDAUGH - 1] = UQ[NIR - 1][IDAUGH - 1] * 1.0E3;
      //
      for (int IB = 1; IB <= NBRANCH[NIR - 1][IDAUGH - 1]; IB++)
	{
	  fflush (stdout);
	  int ITPRO;
	  READSKP (NRD, WCODE, TEMPBUF);
	  sscanf (strtok (TEMPBUF, s), "%lf",
		  &BR[NIR - 1][IDAUGH - 1][IB - 1]);
	  sscanf (strtok (NULL, s), "%lf", &UBR[NIR - 1][IDAUGH - 1][IB - 1]);
	  sscanf (strtok (NULL, s), "%d", &IFLB[NIR - 1][IDAUGH - 1][IB - 1]);
	  sscanf (strtok (NULL, s), "%lf",
		  &EBRANCH[NIR - 1][IDAUGH - 1][IB - 1]);
	  sscanf (strtok (NULL, s), "%lf",
		  &UEBRANCH[NIR - 1][IDAUGH - 1][IB - 1]);
	  sscanf (strtok (NULL, s), "%d", &ITPRO);
	  //  The sign of ITPRO distinguishes non-unique beta transitions.
	  IPRBETA[NIR - 1][IDAUGH - 1][IB - 1] = abs (ITPRO);
	  //
	  if (strcmp (WCODE, "ALP") == 0)
	    {
	      KTYPE[NIR - 1][IDAUGH - 1][IB - 1] = 1;
	    }
	  if (strcmp (WCODE, "BEM") == 0)
	    {
	      KTYPE[NIR - 1][IDAUGH - 1][IB - 1] = 2;
	    }
	  if (strcmp (WCODE, "BEP") == 0)
	    {
	      KTYPE[NIR - 1][IDAUGH - 1][IB - 1] = 3;
	    }
	  if (strcmp (WCODE, "CK ") == 0)
	    {
	      KTYPE[NIR - 1][IDAUGH - 1][IB - 1] = 4;
	    }
	  if (strcmp (WCODE, "CL1") == 0)
	    {
	      KTYPE[NIR - 1][IDAUGH - 1][IB - 1] = 5;
	    }
	  if (strcmp (WCODE, "CL2") == 0)
	    {
	      KTYPE[NIR - 1][IDAUGH - 1][IB - 1] = 6;
	    }
	  if (strcmp (WCODE, "CL3") == 0)
	    {
	      KTYPE[NIR - 1][IDAUGH - 1][IB - 1] = 7;
	    }
	  if (strcmp (WCODE, "CM ") == 0)
	    {
	      KTYPE[NIR - 1][IDAUGH - 1][IB - 1] = 8;
	    }
	  if (strcmp (WCODE, "MS ") == 0)
	    {
	      KTYPE[NIR - 1][IDAUGH - 1][IB - 1] = 9;
	    }			// Metastable level.
	  if (strcmp (WCODE, "CL ") == 0)
	    {
	      KTYPE[NIR - 1][IDAUGH - 1][IB - 1] = 5;
	    }			// EL ~ EL1.
	  if (strcmp (WCODE, "CN ") == 0)
	    {
	      KTYPE[NIR - 1][IDAUGH - 1][IB - 1] = 8;
	    }			// EN ~ EM.
	  //
	  EBRANCH[NIR - 1][IDAUGH - 1][IB - 1] =
	    EBRANCH[NIR - 1][IDAUGH - 1][IB - 1] * 1.0E3;
	  UEBRANCH[NIR - 1][IDAUGH - 1][IB - 1] =
	    UEBRANCH[NIR - 1][IDAUGH - 1][IB - 1] * 1.0E3;

	  IZD = IZNUC[NIR - 1];
	  if (KTYPE[NIR - 1][IDAUGH - 1][IB - 1] == 2)
	    {
	      IDT = IDT + 1;
	      INDBETA[NIR - 1][IDAUGH - 1][IB - 1] = IDT;
	      IZD = IZNUC[NIR - 1] + 1;
	      int IMASS = IANUC[NIR - 1];
	      SBETA0 (IMASS, IZD, EBRANCH[NIR - 1][IDAUGH - 1][IB - 1],
		      IPRBETA[NIR - 1][IDAUGH - 1][IB - 1], IDT, IER);
	      if (IER != 0)
		{
		  return;
		}
	    }
	  //
	  if (KTYPE[NIR - 1][IDAUGH - 1][IB - 1] == 3)
	    {
	      IDT = IDT + 1;
	      INDBETA[NIR - 1][IDAUGH - 1][IB - 1] = IDT;
	      IZD = -(IZNUC[NIR - 1] - 1);
	      int IMASS = IANUC[NIR - 1];
	      SBETA0 (IMASS, IZD, EBRANCH[NIR - 1][IDAUGH - 1][IB - 1],
		      IPRBETA[NIR - 1][IDAUGH - 1][IB - 1], IDT, IER);
	      if (IER != 0)
		{
		  return;
		}
	    }

	}
      //
      for (int IL = NLEVEL[NIR - 1][IDAUGH - 1]; IL >= 1; IL--)
	{
	  READSKP (NRD, WCODE, TEMPBUF);
	  sscanf (TEMPBUF, "%lf , %lf , %d , %lf , %lf",
		  &ELEVEL[NIR - 1][IDAUGH - 1][IL - 1],
		  &UELEVEL[NIR - 1][IDAUGH - 1][IL - 1],
		  &NTRANS[NIR - 1][IDAUGH - 1][IL - 1],
		  &DTIME[NIR - 1][IDAUGH - 1][IL - 1],
		  &UDTIME[NIR - 1][IDAUGH - 1][IL - 1]);
	  //
	  ELEVEL[NIR - 1][IDAUGH - 1][IL - 1] =
	    ELEVEL[NIR - 1][IDAUGH - 1][IL - 1] * 1.0E3;
	  UELEVEL[NIR - 1][IDAUGH - 1][IL - 1] =
	    UELEVEL[NIR - 1][IDAUGH - 1][IL - 1] * 1.0E3;
	  //  ****  Levels with no defined transitions.
	  bool Goto200 = true;
	  if (NTRANS[NIR - 1][IDAUGH - 1][IL - 1] < 1)
	    {
	      PTRANS[NIR - 1][IDAUGH - 1][IL - 1][0] = 1.0E0;
	      UPTRANS[NIR - 1][IDAUGH - 1][IL - 1][0] = 0.0E0;
	      ETRANS[NIR - 1][IDAUGH - 1][IL - 1][0] =
		ELEVEL[NIR - 1][IDAUGH - 1][IL - 1];
	      UETRANS[NIR - 1][IDAUGH - 1][IL - 1][0] =
		UELEVEL[NIR - 1][IDAUGH - 1][IL - 1];
	      NTRANS[NIR - 1][IDAUGH - 1][IL - 1] = 1;
	      MODE[NIR - 1][IDAUGH - 1][IL - 1][0] = 7;
	      IFLT[NIR - 1][IDAUGH - 1][IL - 1][0] = 0;
	      Goto200 = false;
	    }
	  //
	  if (Goto200)
	    {
	      for (int IT = NTRANS[NIR - 1][IDAUGH - 1][IL - 1]; IT >= 1;
		   IT--)
		{
		  READSKP (NRD, WCODE, TEMPBUF);
		  sscanf (TEMPBUF, "%lf , %lf , %lf , %lf , %d",
			  &PTRANS[NIR - 1][IDAUGH - 1][IL - 1][IT - 1],
			  &UPTRANS[NIR - 1][IDAUGH - 1][IL - 1][IT - 1],
			  &ETRANS[NIR - 1][IDAUGH - 1][IL - 1][IT - 1],
			  &UETRANS[NIR - 1][IDAUGH - 1][IL - 1][IT - 1],
			  &IFLT[NIR - 1][IDAUGH - 1][IL - 1][IT - 1]);
		  //
		  if (strcmp (WCODE, "GA ") == 0)
		    {
		      MODE[NIR - 1][IDAUGH - 1][IL - 1][IT - 1] = 1;
		    }
		  if (strcmp (WCODE, "EK ") == 0)
		    {
		      MODE[NIR - 1][IDAUGH - 1][IL - 1][IT - 1] = 2;
		    }
		  if (strcmp (WCODE, "EL1") == 0)
		    {
		      MODE[NIR - 1][IDAUGH - 1][IL - 1][IT - 1] = 3;
		    }
		  if (strcmp (WCODE, "EL2") == 0)
		    {
		      MODE[NIR - 1][IDAUGH - 1][IL - 1][IT - 1] = 4;
		    }
		  if (strcmp (WCODE, "EL3") == 0)
		    {
		      MODE[NIR - 1][IDAUGH - 1][IL - 1][IT - 1] = 5;
		    }
		  if (strcmp (WCODE, "EM ") == 0)
		    {
		      MODE[NIR - 1][IDAUGH - 1][IL - 1][IT - 1] = 6;
		    }
		  if (strcmp (WCODE, "EL ") == 0)
		    {
		      MODE[NIR - 1][IDAUGH - 1][IL - 1][IT - 1] = 3;
		    }		// EL ~ EL1.
		  if (strcmp (WCODE, "EN ") == 0)
		    {
		      MODE[NIR - 1][IDAUGH - 1][IL - 1][IT - 1] = 6;
		    }		// EN ~ EM.
		  //
		  ETRANS[NIR - 1][IDAUGH - 1][IL - 1][IT - 1] =
		    ETRANS[NIR - 1][IDAUGH - 1][IL - 1][IT - 1] * 1.E3;
		  UETRANS[NIR - 1][IDAUGH - 1][IL - 1][IT - 1] =
		    UETRANS[NIR - 1][IDAUGH - 1][IL - 1][IT - 1] * 1.E3;
		}
	    }
	}
    }
  fclose (NRD);
  //
  // ************  Write nuclear data.
  //
  fprintf (IWR,
	   " Radioisotope number %3d\n Parent nuclide:  %s (Z =%3d, A = %3d)\n",
	   NIR, NUCLIDE[NIR - 1], IZNUC[NIR - 1], IANUC[NIR - 1]);
  fprintf (IWR, " Number of daughters =%3d\n", NDAUGH[NIR - 1]);
  fprintf (IWR, "\n EC = electron capture,  IC = internal conversion.\n");
  fprintf (IWR, " BR = branching ratio,   PR = probability.\n");

  for (int IDc = 1; IDc <= NDAUGH[NIR - 1]; IDc++)
    {
      if (KTYPE[NIR - 1][IDc - 1][0] == 1)	// ALPHA
	{
	  IZD = IZNUC[NIR - 1] - 2;
	  IAD = IANUC[NIR - 1] - 2;
	}
      else if (KTYPE[NIR - 1][IDc - 1][0] == 2)	// BETA-
	{
	  IZD = IZNUC[NIR - 1] + 1;
	  IAD = IANUC[NIR - 1];
	}
      else if (KTYPE[NIR - 1][IDc - 1][0] == 9)	// Metastable
	{
	  IZD = IZNUC[NIR - 1];
	  IAD = IANUC[NIR - 1];
	}
      else			// BETA+ and EC
	{
	  IZD = IZNUC[NIR - 1] - 1;
	  IAD = IANUC[NIR - 1];
	}
      fputc (' ', IWR);
      for (int i = 0; i < 76; i++)
	fputc ('-', IWR);
      fputc ('\n', IWR);
      fprintf (IWR,
	       " Daughter nuclide number%4d:  %s (Z =%3d, A = %3d)\n   Q =%13.6E eV,  BR =%13.6E\n",
	       IDc, DNUCLIDE[NIR - 1][IDc - 1], IZD, IAD, Q[NIR - 1][IDc - 1],
	       PDAUGH[NIR - 1][IDc - 1]);
      //
      if (NLEVEL[NIR - 1][IDc - 1] > 0)
	{
	  fprintf (IWR, "\n Excited levels:\n");
	  for (int IL = NLEVEL[NIR - 1][IDc - 1]; IL >= 1; IL--)
	    {
	      fprintf (IWR, " %3d,  E =%13.6E eV,  half-life =%13.6E s\n", IL,
		       ELEVEL[NIR - 1][IDc - 1][IL - 1],
		       DTIME[NIR - 1][IDc - 1][IL - 1]);
	    }
	}
      //
      fprintf (IWR, "\n Disintegration branches:\n");
      for (int IB = 1; IB <= NBRANCH[NIR - 1][IDc - 1]; IB++)
	{
	  IFLB[NIR - 1][IDc - 1][IB - 1] = (IFLB[NIR - 1][IDc - 1][IB - 1] < NLEVEL[NIR - 1][IDc - 1] ? IFLB[NIR - 1][IDc - 1][IB - 1] : NLEVEL[NIR - 1][IDc - 1]);	// Caution
	  if (KTYPE[NIR - 1][IDc - 1][IB - 1] == 1)	// Alpha decay.
	    {
	      fprintf (IWR,
		       " %3d,  %s to level %3d,  BR =%13.6E,  E =%13.6E eV\n",
		       IB, WTYPE[KTYPE[NIR - 1][IDc - 1][IB - 1] - 1],
		       IFLB[NIR - 1][IDc - 1][IB - 1],
		       BR[NIR - 1][IDc - 1][IB - 1],
		       EBRANCH[NIR - 1][IDc - 1][IB - 1]);
	    }
	  else if (KTYPE[NIR - 1][IDc - 1][IB - 1] == 2)	// Beta- decay.
	    {
	      fprintf (IWR,
		       " %3d,  %s to level %3d,  BR =%13.6E,  E =%13.6E eV, Ipro =%2d\n",
		       IB, WTYPE[KTYPE[NIR - 1][IDc - 1][IB - 1] - 1],
		       IFLB[NIR - 1][IDc - 1][IB - 1],
		       BR[NIR - 1][IDc - 1][IB - 1],
		       EBRANCH[NIR - 1][IDc - 1][IB - 1],
		       IPRBETA[NIR - 1][IDc - 1][IB - 1]);
	      EPMAX =
		(EPMAX >
		 EBRANCH[NIR - 1][IDc - 1][IB - 1] ? EPMAX : EBRANCH[NIR -
								    1][IDc -
								       1][IB -
									  1]);
	    }
	  else if (KTYPE[NIR - 1][IDc - 1][IB - 1] == 3)	// Beta+ decay.
	    {
	      fprintf (IWR,
		       " %3d,  %s to level %3d,  BR =%13.6E,  E =%13.6E eV, Ipro =%2d\n",
		       IB, WTYPE[KTYPE[NIR - 1][IDc - 1][IB - 1] - 1],
		       IFLB[NIR - 1][IDc - 1][IB - 1],
		       BR[NIR - 1][IDc - 1][IB - 1],
		       EBRANCH[NIR - 1][IDc - 1][IB - 1],
		       IPRBETA[NIR - 1][IDc - 1][IB - 1]);
	      //  ----  Positrons eventually give annihilation gamma rays. The maximum
	      //        energy of annihilation photons is .lt. 1.21*(E0+me*c^2).
	      EPMAX =
		(EPMAX >
		 1.21E0 * (EBRANCH[NIR - 1][IDc - 1][IB - 1] +
			   5.12E5) ? EPMAX : 1.21E0 * (EBRANCH[NIR - 1][IDc -
									1][IB
									   -
									   1]
						       + 5.12E5));
	    }
	  else			// Electron capture and metastable nuclides.
	    {
	      fprintf (IWR, " %3d,  %s to level %3d,  BR =%13.6E\n",
		       IB, WTYPE[KTYPE[NIR - 1][IDc - 1][IB - 1] - 1],
		       IFLB[NIR - 1][IDc - 1][IB - 1],
		       BR[NIR - 1][IDc - 1][IB - 1]);
	    }
	}
      //
      if (NLEVEL[NIR - 1][IDc - 1] > 0)
	{
	  for (int IL = NLEVEL[NIR - 1][IDc - 1]; IL >= 1; IL--)
	    {
	      fprintf (IWR, "\n Transitions from level %3d\n", IL);
	      for (int IT = 1; IT <= NTRANS[NIR - 1][IDc - 1][IL - 1]; IT++)
		{
		  fprintf (IWR,
			   " %3d --> %3d %s,  E =%13.6E eV,  PR =%13.6E\n",
			   IL, IFLT[NIR - 1][IDc - 1][IL - 1][IT - 1],
			   WMODE[MODE[NIR - 1][IDc - 1][IL - 1][IT - 1] - 1],
			   ETRANS[NIR - 1][IDc - 1][IL - 1][IT - 1],
			   PTRANS[NIR - 1][IDc - 1][IL - 1][IT - 1]);
		  EPMAX =
		    (EPMAX >
		     ETRANS[NIR - 1][IDc - 1][IL - 1][IT -
						     1] ? EPMAX : ETRANS[NIR -
									 1][IDc
									    -
									    1]
		     [IL - 1][IT - 1]);
		}
	    }
	}
    }
  //
  //  ****  Accumulated daugther probabilities.
  //
  double PTOT = 0.0E0;
  for (int IDc = 1; IDc <= NDAUGH[NIR - 1]; IDc++)
    {
      PTOT = PTOT + PDAUGH[NIR - 1][IDc - 1];
    }
  for (int IDc = 1; IDc <= NDAUGH[NIR - 1]; IDc++)
    {
      PDAUGH[NIR - 1][IDc - 1] = PDAUGH[NIR - 1][IDc - 1] / PTOT;
      UPDAUGH[NIR - 1][IDc - 1] = UPDAUGH[NIR - 1][IDc - 1] / PTOT;
    }
  fputc (' ', IWR);
  for (int i = 0; i < 76; i++)
    fputc ('-', IWR);
  fputc ('\n', IWR);
  fprintf (IWR, " Branching ratios to daughters (normalised):\n");
  for (int IDc = 1; IDc <= NDAUGH[NIR - 1]; IDc++)
    {
      fprintf (IWR, "%13.6E", PDAUGH[NIR - 1][IDc - 1]);
    }
  fprintf (IWR, "\n");
  fprintf (IWR, "\n Maximum energy of released particles = %13.6E eV\n",
	   EPMAX);
  fprintf (IWR, "\n ---  END  ---\n");
  //
  if (NDAUGH[NIR - 1] > 1)
    {
      for (int IDc = 2; IDc <= NDAUGH[NIR - 1]; IDc++)
	{
	  PDAUGH[NIR - 1][IDc - 1] =
	    PDAUGH[NIR - 1][IDc - 1] + PDAUGH[NIR - 1][IDc - 1 - 1];
	}
      PDAUGH[NIR - 1][NDAUGH[NIR - 1] - 1] = 1.0E0;
    }
  //
  //  ****  Accumulated branch probabilities.
  //
  for (int IDc = 1; IDc <= NDAUGH[NIR - 1]; IDc++)
    {
      double BRTOT = 0.0E0;
      for (int IB = 1; IB <= NBRANCH[NIR - 1][IDc - 1]; IB++)
	{
	  BRTOT = BRTOT + BR[NIR - 1][IDc - 1][IB - 1];
	}
      for (int IB = 1; IB <= NBRANCH[NIR - 1][IDc - 1]; IB++)
	{
	  BR[NIR - 1][IDc - 1][IB - 1] = BR[NIR - 1][IDc - 1][IB - 1] / BRTOT;
	  UBR[NIR - 1][IDc - 1][IB - 1] = UBR[NIR - 1][IDc - 1][IB - 1] / BRTOT;
	}
      //
      if (NBRANCH[NIR - 1][IDc - 1] > 1)
	{
	  for (int IB = 2; IB <= NBRANCH[NIR - 1][IDc - 1]; IB++)
	    {
	      BR[NIR - 1][IDc - 1][IB - 1] =
		BR[NIR - 1][IDc - 1][IB - 1] + BR[NIR - 1][IDc - 1][IB - 1 - 1];
	    }
	}
      BR[NIR - 1][IDc - 1][NBRANCH[NIR - 1][IDc - 1] - 1] = 1.0E0;
    }
  //
  //  ****  Accumulated level depopulation probabilities.
  //
  for (int IDc = 1; IDc <= NDAUGH[NIR - 1]; IDc++)
    {
      for (int IL = 1; IL <= NLEVEL[NIR - 1][IDc - 1]; IL++)
	{
	  PTOT = 0.0E0;
	  for (int IT = 1; IT <= NTRANS[NIR - 1][IDc - 1][IL - 1]; IT++)
	    {
	      PTOT = PTOT + PTRANS[NIR - 1][IDc - 1][IL - 1][IT - 1];
	    }
	  for (int IT = 1; IT <= NTRANS[NIR - 1][IDc - 1][IL - 1]; IT++)
	    {
	      PTRANS[NIR - 1][IDc - 1][IL - 1][IT - 1] =
		PTRANS[NIR - 1][IDc - 1][IL - 1][IT - 1] / PTOT;
	      UPTRANS[NIR - 1][IDc - 1][IL - 1][IT - 1] =
		UPTRANS[NIR - 1][IDc - 1][IL - 1][IT - 1] / PTOT;
	    }
	  //
	  if (NTRANS[NIR - 1][IDc - 1][IL - 1] > 1)
	    {
	      for (int IT = 2; IT <= NTRANS[NIR - 1][IDc - 1][IL - 1]; IT++)
		{
		  PTRANS[NIR - 1][IDc - 1][IL - 1][IT - 1] =
		    PTRANS[NIR - 1][IDc - 1][IL - 1][IT - 1] + PTRANS[NIR -
								     1][IDc -
									1][IL
									   -
									   1]
		    [IT - 1 - 1];
		}
	    }
	  PTRANS[NIR - 1][IDc - 1][IL - 1][NTRANS[NIR - 1][IDc - 1][IL - 1] -
					  1] = 1.0E0;
	}
    }

  IZD = 0;
  IAD = 0;
  LEVD = 9999;
  return;
}

//  *********************************************************************
//                       SUBROUTINE READSKP
//  *********************************************************************
void pennuc_specificSampler::READSKP(FILE * NRD,
				     char WCODE[4],
				     char TEMPBUF[100]) const {
  //
  //  Reads and reformats an input record. Skips comment records, and in
  //  data records replaces any ';' by ','.
  //
  //
  char straux[200];
  strcpy (WCODE, "COM");
  while (strcmp (WCODE, "COM") == 0)
    {
      fgets (straux, sizeof (straux), NRD);
      straux[strlen (straux) - 1] = '\0';
      sscanf (straux, "%3c %[^\n]%*[^\n]", WCODE, TEMPBUF);
      WCODE[strlen (WCODE)] = '\0';
      TEMPBUF[strlen (TEMPBUF) - 1] = '\0';
    }
  //

  char *pch;
  pch = strchr (TEMPBUF, ';');
  while (pch != NULL)
    {
      int IPOS = pch - TEMPBUF;
      TEMPBUF[IPOS] = ',';
      pch = strchr (pch + 1, ';');
    }
  //
  return;
}

//  *********************************************************************
//                       SUBROUTINE ATREL
//  *********************************************************************
void pennuc_specificSampler::ATREL(int IZ,
				   int IS,
				   double* ES,
				   double* PAGES,
				   int* KPARS,
				   int* ITR,
				   int &NS,
				   pen_rand& random){
  //
  //  This subroutine simulates the relaxation of a singly ionised atom of
  //  the element IZ with a vacancy in the IS subshell (the K shell or an
  //  L, M or N subshell). This initial vacancy is filled by electrons from
  //  outer subshells through radiative and non-radiative transitions,
  //  which may produce additional vacancies.
  //
  //  Input arguments:
  //     IZ ......... atomic number of the ionised element.
  //     IS ......... electron subshell containing the initial vacancy.
  //     ECNUC ...... cutoff energy for atomic transitions.
  //
  //  Output arguments:
  //     ES(500) .... energies of emitted particles (in eV).
  //     KPARS(500) ... kind of emitted particle (1=electron; 2=photon).
  //     PAGES(500) ... emission times (seconds), time elapsed since the
  //                  generation of the initial vacancy.
  //     ITR(500) ... transition (coded as in PENELOPE's ILB(4) flag).
  //     NS ......... number of particles in the relaxation cascade.
  //
  //  We use the following notation to designate the possible transitions:
  //  *  Radiative: IS0-IS1 (an electron from the IS1 subshell fills the
  //     vacancy in the IS0 subshell, leaving a hole in the IS1 subshell).
  //  *  Non-radiative: IS0-IS1-IS2 (an electron from the IS1 subshell
  //     fills the vacancy in the IS0 subshell, and the released energy is
  //     taken away by an electron in the IS2 subshell; this process leaves
  //     holes, in the IS1 and IS2 subshells).
  //  The de-excitation cascade (i.e. the set of transitions that occur for
  //  a given initial vacancy) is sampled from the transition probabilities
  //  extracted from the Livermore Evaluated Atomic Data Library (EADL).
  //  The energies of the particles emitted in each transition are read
  //  from the PENELOPE database.
  //
  //  The simulation of the de-excitation cascade is discontinued either
  //  when the K to N subshells have been filled up or when there is not
  //  enough energy to produce 'active' radiation (with energy larger than
  //  ECNUC).
  //
  //  De-excitation data for the loaded elements are stored in the common
  //  block /CATREL/, in a form designed to minimise the amount of memory
  //  and to facilitate the random sampling. The quantities in the common
  //  block are the following:
  //  IFIRST(99,16) ... de-excitation data for a vacancy in the subshell IS
  //     of the element IZ start at the position K=IFIRST(IZ,IS) in the
  //     storage arrays. The allowed values for IS are 1 to 16 (K shell
  //     and L, M and N subshells).
  //  ILAST(99,16) ... the de-excitation data for a vacancy in the subshell
  //     IS of the element IZ end at the position K=ILAST(IZ,IS) in the
  //     storage arrays.
  //  IS1(K), IS2(K) ... subshells that are active in the transition (see
  //     the subshell label code below). For radiative transitions,
  //     IS2(K)=0.
  //  P(K) ... relative probability for the transition IS-IS1(K)-IS2(K).
  //  ET(K) ... energy of the secondary particle emitted in the transition.
  //  F(K), IAL(K) ... cutoff and alias values (Walker's sampling method).
  //
  //  ---------------------------------------------------------------------
  //  Label code IS for electron subshells:
  //      1 = K  (1s1/2),     11 = N2 (4p1/2),     21 = O5 (5d5/2),
  //      2 = L1 (2s1/2),     12 = N3 (4p3/2),     22 = O6 (5f5/2),
  //      3 = L2 (2p1/2),     13 = N4 (4d3/2),     23 = O7 (5f7/2),
  //      4 = L3 (2p3/2),     14 = N5 (4d5/2),     24 = P1 (6s1/2),
  //      5 = M1 (3s1/2),     15 = N6 (4f5/2),     25 = P2 (6p1/2),
  //      6 = M2 (3p1/2),     16 = N7 (4f7/2),     26 = P3 (6p3/2),
  //      7 = M3 (3p3/2),     17 = O1 (5s1/2),     27 = P4 (6d3/2),
  //      8 = M4 (3d3/2),     18 = O2 (5p1/2),     28 = P5 (6d5/2),
  //      9 = M5 (3d5/2),     19 = O3 (5p3/2),     29 = Q1 (7s1/2),
  //     10 = N1 (4s1/2),     20 = O4 (5d3/2),     30 = outer shells.
  //  ---------------------------------------------------------------------
  //
  //  ****  Atomic relaxation data.
  //Uses CATREL and CATRE2 commons
  //  ****  Subhell occupancies, and auxiliary arrays.
  double OCUP[30] =
    { 2.0E0, 2.0E0, 2.0E0, 4.0E0, 2.0E0, 2.0E0, 4.0E0, 4.0E0, 6.0E0,
      2.0E0, 2.0E0, 4.0E0, 4.0E0, 6.0E0, 6.0E0, 8.0E0, 2.0E0, 2.0E0, 4.0E0,
      4.0E0, 6.0E0, 6.0E0, 8.0E0, 2.0E0, 2.0E0, 4.0E0, 4.0E0, 6.0E0, 2.0E0,
      1.0E9
    };
  double PTIME[16][256];
  int ISV[30];
  //
  const int NEMP = 500;
  //
  //  ****  Initialisation.
  //
  NS = 0;
  if (IZ < 3 || IS > 16)
    {
      return;
    }
  //  ****  If the subshell ionisation energy is less than ECNUC, the
  //        cascade is not followed.
  if (EBR[IZ - 1][IS - 1] < ECNUC)
    {
      return;
    }
  double PAGE0 = 0.0E0;
  //
  int NV = 0;
  for (int I = 1; I <= 30; I++)
    {
      ISV[I - 1] = 0;
    }
  int ISA = IS;
  double PAGE = PAGE0;
  //
  //  ****  Next transition.
  //
  // 1    CONTINUE
  bool Exit = false;
  while (!Exit)
    {
      Exit = true;
      int JS1, JS2, KS;
      int KF = IFIRST[IZ - 1][ISA - 1];
      int KL = ILAST[IZ - 1][ISA - 1];
      if (KL <= KF)
	{
	  return;
	}
      //  ****  Walker's sampling algorithm.
      // 2    CONTINUE
      bool Exit2 = false;
      while (!Exit2)
	{
	  Exit2 = true;
	  double FREJ;
	  double RN = random.rand() * (KL - KF + 1);
	  int K1 = int (RN);
	  double TST = RN - K1;
	  if (TST > F[KF + K1 - 1])
	    {
	      KS = IAL[KF + K1 - 1];
	    }
	  else
	    {
	      KS = KF + K1;
	    }
	  //
	  JS1 = IS1[KS - 1];
	  JS2 = IS2[KS - 1];
	  //  ****  Vacancies in the intervening subshells.
	  if (ISV[JS1 - 1] > 0)
	    {
	      FREJ = (OCUP[JS1 - 1] - ISV[JS1 - 1]) / OCUP[JS1 - 1];
	    }
	  else
	    {
	      FREJ = 1.0E0;
	    }
	  if (JS2 > 0)
	    {
	      FREJ = FREJ * ((OCUP[JS2 - 1] - ISV[JS2 - 1]) / OCUP[JS2 - 1]);
	    }
	  if (FREJ < 1.0E0)
	    {
	      if (random.rand() > FREJ)
		{
		  Exit2 = false;
		  continue;
		}
	    }
	}
      //  ****  The particle age is evaluated.
      if (ISV[ISA - 1] > 1)
	{
	  PAGE =
	    PAGE -
	    (ALW[IZ - 1][ISA - 1] /
	     double (ISV[ISA - 1])) *log (random.rand());
	}
      else
	{
	  PAGE = PAGE - ALW[IZ - 1][ISA - 1] * log (random.rand());
	}
      //
      //  ****  Fluorescence radiation.
      //
      if (JS2 == 0)
	{
	  if (JS1 < 17)
	    {
	      if (EBR[IZ - 1][JS1 - 1] > ECNUC)
		{
		  NV = NV + 1;
		  ISV[JS1 - 1] = ISV[JS1 - 1] + 1;
		  PTIME[JS1 - 1][ISV[JS1 - 1] - 1] = PAGE;
		}
	    }
	  if (ET[KS - 1] > ECNUC)
	    {
	      NS = NS + 1;
	      if (NS > NEMP)
		{
		  printf ("NS =%d\n", NS);
		  printf ("PENNUC warning: Dimension NEMP exceeded (1).\n");
		  return;
		}
	      ES[NS - 1] = ET[KS - 1];
	      KPARS[NS - 1] = 2;
	      ITR[NS - 1] = IZ * 1000000 + ISA * 10000 + JS1 * 100;
	      PAGES[NS - 1] = PAGE;
	    }
	}
      else
	{
	  if (JS1 < 17)
	    {
	      if (EBR[IZ - 1][JS1 - 1] > ECNUC)
		{
		  NV = NV + 1;
		  ISV[JS1 - 1] = ISV[JS1 - 1] + 1;
		  PTIME[JS1 - 1][ISV[JS1 - 1] - 1] = PAGE;
		}
	    }
	  if (JS2 < 17)
	    {
	      if (EBR[IZ - 1][JS2 - 1] > ECNUC)
		{
		  NV = NV + 1;
		  ISV[JS2 - 1] = ISV[JS2 - 1] + 1;
		  PTIME[JS2 - 1][ISV[JS2 - 1] - 1] = PAGE;
		}
	    }
	  if (ET[KS - 1] > ECNUC)
	    {
	      NS = NS + 1;
	      if (NS > NEMP)
		{
		  printf ("NS =%d\n", NS);
		  printf ("PENNUC warning: Dimension NEMP exceeded (2).\n");
		  return;
		}
	      ES[NS - 1] = ET[KS - 1];
	      KPARS[NS - 1] = 1;
	      ITR[NS - 1] = IZ * 1000000 + ISA * 10000 + JS1 * 100 + JS2;
	      PAGES[NS - 1] = PAGE;
	    }
	}
      //
      //  ****  Are there any K-, L-, M- or N-vacancies unfilled?
      //
      bool Goto1 = false;
      if (NV > 0)
	{
	  bool Goto3 = false;
	  Goto1 = false;
	  bool Exit3 = false;
	  while (!Exit3)
	    {
	      Exit3 = true;
	      Goto1 = false;
	      Goto3 = false;
	      if (EBR[IZ - 1][ISA - 1] < ECNUC)
		{
		  return;
		}
	      if (ISV[ISA - 1] > 0)	// Vacancies in the current subshell.
		{
		  PAGE = PTIME[ISA - 1][ISV[ISA - 1] - 1];
		  ISV[ISA - 1] = ISV[ISA - 1] - 1;
		  NV = NV - 1;
		  Goto1 = true;
		  Exit3 = true;
		  break;
		}
	      else
		{
		  for (int IST = ISA + 1; IST <= 16; IST++)	// Outer subshells.
		    {
		      if (ISV[IST - 1] > 0)
			{
			  ISA = IST;
			  Goto3 = true;
			  break;
			}
		    }
		}
	      if (Goto3)
		{
		  Exit3 = false;
		  continue;
		}
	    }
	}
      if (Goto1)
	{
	  Exit = false;
	  continue;
	}
    }
  //
  return;
}

//  *********************************************************************
//                       SUBROUTINE ATREL0
//  *********************************************************************
void pennuc_specificSampler::ATREL0(){
  //
  //  This subroutine sets all static values in common /CATREL/ to zero.
  //

  //  ****  Atomic relaxation data.
  //Uses CATREL common
  //
  for (int I = 1; I <= 99; I++)
    {
      for (int J = 1; J <= 16; J++)
	{
	  IFIRST[I - 1][J - 1] = 0;
	  ILAST[I - 1][J - 1] = 0;
	}
      for (int J = 1; J <= NSHP1; J++)
	{
	  EBR[I - 1][J - 1] = 0.0E0;
	}
    }
  //
  for (int I = 1; I <= NRX; I++)
    {
      P[I - 1] = 0.0E0;
      ET[I - 1] = 0.0E0;
      F[I - 1] = 0.0E0;
      IAL[I - 1] = 0;
      IS0[I - 1] = 0;
      IS1[I - 1] = 0;
      IS2[I - 1] = 0;
    }
  NCUR = 0;
  return;
}

//  *********************************************************************
//                       SUBROUTINE ATRELI
//  *********************************************************************
void pennuc_specificSampler::ATRELI(const char* ATRELIFNAME,
            const char* ATOMFNAME, const char* RELAXFNAME,
            int IZ,
				    int IWR,
				    int &IER){
  //
  //  This subroutine produces a table of atomic relaxation data for the
  //  element IZ.
  //
  //  Input arguments:
  //    IZ ...... atomic number.
  //    IWR ..... relaxation data are written in a file (UNIT=IWR) named
  //              'atreli.dat' if IWR is greater than zero. It is assumed
  //              that unit IWR is not open in the calling program.
  //  Output value:
  //    IER  .... error code. When IER.NE.0 the initialisation has not
  //              succeeded; the cause of the problem is written in
  //              the standard output unit (6).
  //
  //  Data are read from file 'pdrelax.p11', which contains data pertaining
  //  to singly ionised atoms with the initial vacancy in one of the K, L,
  //  M and N shells. This file was prepared from the Livermore Evaluated
  //  Atomic Data Library (EADL). The energies of x-ray lines were replaced
  //  by more accurate experimental and theoretical values given by
  //  Deslattes et al. (2003) -K and L shells- and by Burr (1967) -M and N
  //  shells.
  //
  //Uses PENNUC_mod, CATREL and CATRE2 common;
  //
  //      const double HBAR=6.58211928E-16;  // Planck's constant (eV*s)
  char CH5[6], CH2[3]; //, FNAME[46];
  char CSH5[30][6];
  double EB[30];
  int IFI[30], IKS[30];
  //
  strcpy (CSH5[0], "1s1/2");
  strcpy (CSH5[1], "2s1/2");
  strcpy (CSH5[2], "2p1/2");
  strcpy (CSH5[3], "2p3/2");
  strcpy (CSH5[4], "3s1/2");
  strcpy (CSH5[5], "3p1/2");
  strcpy (CSH5[6], "3p3/2");
  strcpy (CSH5[7], "3d3/2");
  strcpy (CSH5[8], "3d5/2");
  strcpy (CSH5[9], "4s1/2");
  strcpy (CSH5[10], "4p1/2");
  strcpy (CSH5[11], "4p3/2");
  strcpy (CSH5[12], "4d3/2");
  strcpy (CSH5[13], "4d5/2");
  strcpy (CSH5[14], "4f5/2");
  strcpy (CSH5[15], "4f7/2");
  strcpy (CSH5[16], "5s1/2");
  strcpy (CSH5[17], "5p1/2");
  strcpy (CSH5[18], "5p3/2");
  strcpy (CSH5[19], "5d3/2");
  strcpy (CSH5[20], "5d5/2");
  strcpy (CSH5[21], "5f5/2");
  strcpy (CSH5[22], "5f7/2");
  strcpy (CSH5[23], "6s1/2");
  strcpy (CSH5[24], "6p1/2");
  strcpy (CSH5[25], "6p3/2");
  strcpy (CSH5[26], "6d3/2");
  strcpy (CSH5[27], "6d5/2");
  strcpy (CSH5[28], "7s1/2");
  strcpy (CSH5[29], " free");
  //  ****  Atomic relaxation data.
  //
  const int NTRAN = 2500;
  int JS0[NTRAN], JS1[NTRAN], JS2[NTRAN];
  double PR[NTRAN], ER[NTRAN], WW[NTRAN], FF[NTRAN];
  int KK[NTRAN], IORD[NTRAN], ISR[NTRAN];
  double EE[30];
  //
  char LSHELL[31][3];
  strcpy (LSHELL[0], "  ");
  strcpy (LSHELL[1], "K ");
  strcpy (LSHELL[2], "L1");
  strcpy (LSHELL[3], "L2");
  strcpy (LSHELL[4], "L3");
  strcpy (LSHELL[5], "M1");
  strcpy (LSHELL[6], "M2");
  strcpy (LSHELL[7], "M3");
  strcpy (LSHELL[8], "M4");
  strcpy (LSHELL[9], "M5");
  strcpy (LSHELL[10], "N1");
  strcpy (LSHELL[11], "N2");
  strcpy (LSHELL[12], "N3");
  strcpy (LSHELL[13], "N4");
  strcpy (LSHELL[14], "N5");
  strcpy (LSHELL[15], "N6");
  strcpy (LSHELL[16], "N7");
  strcpy (LSHELL[17], "O1");
  strcpy (LSHELL[18], "O2");
  strcpy (LSHELL[19], "O3");
  strcpy (LSHELL[20], "O4");
  strcpy (LSHELL[21], "O5");
  strcpy (LSHELL[22], "O6");
  strcpy (LSHELL[23], "O7");
  strcpy (LSHELL[24], "P1");
  strcpy (LSHELL[25], "P2");
  strcpy (LSHELL[26], "P3");
  strcpy (LSHELL[27], "P4");
  strcpy (LSHELL[28], "P5");
  strcpy (LSHELL[29], "Q1");
  strcpy (LSHELL[30], "X ");
  //
  //  ****  Check if this element's data have already been loaded.
  //
  IER = 0;
  if (IFIRST[IZ - 1][0] != 0)
    {
      printf ("ATRELI: The element has already been loaded.\n");
      return;
    }
  //
  FILE* fileIWR = nullptr;
  if (IWR > 0)
    {
      char straux[200];
      strcpy (straux, ATRELIFNAME);
      //fileIWR = fopen ("atreli.dat", "w");
      fileIWR = fopen (straux, "w");
    }
  //
  //  ****  Loads atomic configuration and subshell binding energies.
  //
  //strcpy (FNAME, "pdatconf.p14");
  FILE *MRD;
  char straux[200];
  //strcpy (straux, FNAME);
  strcpy (straux, ATOMFNAME);
  MRD = fopen (straux, "r");
  if (MRD == NULL)
    {
      if (fileIWR != nullptr)
	{
	  fprintf (fileIWR,
		   "\n ATRELI: The program could not open the file %s",
		   straux);
	  fprintf (fileIWR,
		   "\n Please, copy that file into the PENNUC data directory\n");
	}
      printf ("\n ATRELI: The program could not open the file %s", straux);
      printf ("\n Please, copy that file into the PENNUC data directory\n");
      IER = 21;
      return;
    }
  //
  for (int J = 1; J <= 22; J++)
    {
      fgets (straux, sizeof (straux), MRD);
      straux[strlen (straux) - 1] = '\0';
      sscanf (straux, "%s", CH5);
    }
  for (int IS = 1; IS <= 30; IS++)
    {
      IKS[IS - 1] = 0;
    }
  int NS = 0;
  int IZT = 0;
  for (int J = 1; J <= 150000; J++)
    {
      int IIZ, IS, IIF;
      double EIE, CCP, GA1, GA2;
      if (fgets (straux, sizeof (straux), MRD) == NULL)
	{
	  break;
	}
      straux[strlen (straux) - 1] = '\0';
      sscanf (straux, "%3d%4d %s %s%3d%lf%lf%lf%lf", &IIZ, &IS, CH2, CH5,
	      &IIF, &EIE, &CCP, &GA1, &GA2);
      if (IIZ == IZ)
	{
	  NS = NS + 1;
	  if (NS > 30)
	    {
	      printf ("\n NS =%4d\n", NS);
	      printf ("ATRELI: Too many shells.\n");
	      IER = 22;
	      return;
	    }
	  if (IS < 1 || IS > 30)
	    {
	      printf ("\n IS =%4d\n", IS);
	      printf ("ATRELI: Wrong shell number.\n");
	      IER = 23;
	      return;
	    }
	  IZT = IZT + abs (IIF);
	  EB[IS - 1] = EIE;
	  if (GA2 > 0.0E0)
	    {
	      ALW[IZ - 1][IS - 1] = constants::HBAR / GA2;
	    }
	  else if (GA1 > 0.0E0)
	    {
	      ALW[IZ - 1][IS - 1] = constants::HBAR / GA1;
	    }
	  else
	    {
	      ALW[IZ - 1][IS - 1] = 0.0E0;
	    }
	  IFI[IS - 1] = abs (IIF);
	  IKS[NS - 1] = IS;
	}
    }
  int NSHT = NS;
  if (IZ != IZT)
    {
      printf ("IZ, IZT=%d %d\n", IZ, IZT);
      printf ("ATRELI: Unbalanced charges (element).\n");
      IER = 24;
      return;
    }
  fclose (MRD);
  //
  //strcpy (FNAME, "pdrelax.p11");
  //strcpy (straux, FNAME);
  strcpy (straux, RELAXFNAME);
  MRD = fopen (straux, "r");
  if (MRD == NULL)
    {
      if (fileIWR != nullptr)
	{
	  fprintf (fileIWR,
		   "\n ATRELI: The program could not open the file %s\n",
		   straux);
	  fprintf (fileIWR,
		   "\n Please, copy that file into the PENNUC data directory.\n");
	}
      printf ("\n ATRELI: The program could not open the file %s\n", straux);
      printf ("\n Please, copy that file into the PENNUC data directory.\n");
      IER = 25;
      return;
    }
  //
  int IZR, IS0R;
  int NT = 0;
  if (fgets (straux, sizeof (straux), MRD) != NULL)
    {
      straux[strlen (straux) - 1] = '\0';
      sscanf (straux, "%4d%4d", &IZR, &IS0R);
      for (int J = 1; J <= 150000; J++)
	{
	  int IS1R, IS2R;
	  double PPR, EIN;
	  if (fgets (straux, sizeof (straux), MRD) == NULL)
	    {
	      break;
	    }
	  straux[strlen (straux) - 1] = '\0';
	  sscanf (straux, "%4d%4d%4d%4d%lf%lf", &IZR, &IS0R, &IS1R, &IS2R,
		  &PPR, &EIN);
	  if (IZR == IZ)
	    {
	      NT = NT + 1;
	      JS0[NT - 1] = IS0R;
	      JS1[NT - 1] = IS1R;
	      JS2[NT - 1] = IS2R;
	      PR[NT - 1] = PPR;
	      ER[NT - 1] = EIN;
	    }
	}
    }
  fclose (MRD);
  //
  int NSHR;
  if (fileIWR != nullptr)
    {
      fprintf (fileIWR,
	       " *** ATREL:  Z =%3d,  No. of shells =%3d,  no. of transitions =%4d\n",
	       IZ, NSHT, NT);
    }
  NSHR = NSHT;
  if (NSHR >= NSHP1 - 1)
    {
      if (fileIWR != nullptr)
	{
	  fprintf (fileIWR, "Insufficient memory storage in ATRELI.\n");
	}
      if (fileIWR != nullptr)
	{
	  fprintf (fileIWR,
		   "Increase the value of the parameter NSHP1 to %3d\n",
		   NSHR + 1);
	}
      printf ("Insufficient memory storage in ATRELI.\n");
      printf ("Increase the value of the parameter NSHP1 to %3d\n", NSHR + 1);
      IER = 26;
      printf ("ATRELI: Insufficient memory storage.\n");
      return;
    }
  //
  int IS = 0;
  for (int J = 1; J <= 30; J++)
    {
      int KS = IKS[J - 1];
      if (KS > 0)
	{
	  if (IFI[KS - 1] != 0)
	    {
	      if (fileIWR != nullptr)
		{
		  fprintf (fileIWR, " %3d %5s %1d%16.8E\n", KS, CSH5[KS - 1],
			   IFI[KS - 1], EB[KS - 1]);
		}
	      IS = IS + 1;
	      ISR[IS - 1] = KS;
	      EE[IS - 1] = EB[KS - 1];
	    }
	}
    }
  ///
  if (NT > 0)
    {
      double ETT;
      for (int I = 1; I <= NT; I++)
	{
	  if (JS2[I - 1] == 0)
	    {
	      if (ER[I - 1] < 1.0E0)
		{
		  ETT = EE[JS0[I - 1] - 1] - EE[JS1[I - 1] - 1];
		}
	      else
		{
		  ETT = ER[I - 1];
		}
	    }
	  else
	    {
	      if (ER[I - 1] < 1.0E0)
		{
		  ETT =
		    EE[JS0[I - 1] - 1] - EE[JS1[I - 1] - 1] - EE[JS2[I - 1] -
								 1];
		}
	      else
		{
		  ETT = ER[I - 1];
		}
	    }
	  ETT = (ETT > 1.0E0 ? ETT : 1.0E0);
	  if (fileIWR != nullptr)
	    {
	      fprintf (fileIWR, " %3d%3d%3d%16.8E%16.8E\n", JS0[I - 1],
		       JS1[I - 1], JS2[I - 1], PR[I - 1], ETT);
	    }
	  ER[I - 1] = ETT;
	}
    }
  //
  //  ****  Check if this element's data have already been loaded.
  //
  if (IFIRST[IZ - 1][1 - 1] != 0)
    {
      return;
    }
  //
  EBR[IZ - 1][31 - 1] = double (NSHR + 1.0E-12);
  for (IS = 1; IS <= NSHR; IS++)
    {
      EBR[IZ - 1][ISR[IS - 1] - 1] = EE[IS - 1];
    }
  if (NT == 0)
    {
      for (IS = 1; IS <= 16; IS++)
	{
	  NCUR = NCUR + 1;
	  IFIRST[IZ - 1][IS - 1] = NCUR;
	  ILAST[IZ - 1][IS - 1] = NCUR;
	  P[NCUR - 1] = 1.0E0;
	  ET[NCUR - 1] = 0.0E0;
	  F[NCUR - 1] = 1.0E0;
	  IAL[NCUR - 1] = NCUR + 1;
	  IS0[NCUR - 1] = 1;
	  IS1[NCUR - 1] = 1;
	  IS2[NCUR - 1] = 1;
	}
      return;
    }
  //
  //  ****  Walker's aliasing.
  //
  for (IS = 1; IS <= 16; IS++)
    {
      int N = 0;
      for (int J = 1; J <= NT; J++)
	{
	  if (JS0[J - 1] == IS)
	    {
	      N = N + 1;
	      IORD[N - 1] = J;
	      WW[N - 1] = PR[J - 1];
	    }
	}
      if (N > 1)
	{
	  IRND0 (WW, FF, KK, N);
	  IFIRST[IZ - 1][IS - 1] = NCUR + 1;
	  ILAST[IZ - 1][IS - 1] = NCUR + N;
	  for (int L = 1; L <= N; L++)
	    {
	      P[NCUR + L - 1] = WW[L - 1];
	      ET[NCUR + L - 1] = ER[IORD[L - 1] - 1];
	      F[NCUR + L - 1] = FF[L - 1];
	      IAL[NCUR + L - 1] = IFIRST[IZ - 1][IS - 1] + KK[L - 1] - 1;
	      IS0[NCUR + L - 1] = JS0[IORD[L - 1] - 1];
	      IS1[NCUR + L - 1] = JS1[IORD[L - 1] - 1];
	      IS2[NCUR + L - 1] = JS2[IORD[L - 1] - 1];
	    }
	  NCUR = NCUR + N;
	}
      else
	{
	  NCUR = NCUR + 1;
	  IFIRST[IZ - 1][IS - 1] = NCUR;
	  ILAST[IZ - 1][IS - 1] = NCUR;
	  P[NCUR - 1] = 1.0E0;
	  ET[NCUR - 1] = ER[1 - 1];
	  F[NCUR - 1] = 1.0E0;
	  IAL[NCUR - 1] = NCUR;
	  IS0[NCUR - 1] = JS0[1 - 1];
	  IS1[NCUR - 1] = JS1[1 - 1];
	  IS2[NCUR - 1] = JS2[1 - 1];
	}
    }
  //
  //  ****  Verify that transition probabilities are correctly reproduced.
  //
  double TST = 0.0E0;
  for (IS = 1; IS <= 16; IS++)
    {
      int I0 = IFIRST[IZ - 1][IS - 1];
      int iN = ILAST[IZ - 1][IS - 1];
      double PT = 0.0E0;
      for (int I = I0; I <= iN; I++)
	{
	  PT = PT + P[I - 1];
	}
      for (int I = I0; I <= iN; I++)
	{
	  double PPI = 0.0E0;
	  for (int J = I0; J <= iN; J++)
	    {
	      if (IAL[J - 1] == I)
		{
		  PPI = PPI + (1.0E0 - F[J - 1]);
		}
	    }
	  PPI = (PPI + F[I - 1]) / double (iN - I0 + 1);
	  if (TST < fabs (1.0E0 - PPI * PT / P[I - 1]))
	    {
	      TST = fabs (1.0E0 - PPI * PT / P[I - 1]);
	    }
	}
    }
  if (TST > 1.0E-12)
    {
      printf ("TST =%16.8E\n", TST);
      printf ("ATRELI. Rounding error is too large.\n");
      IER = 27;
      return;
    }
  //
  if (fileIWR != nullptr)
    {
      fclose (fileIWR);
    }
  return;
}

//  *********************************************************************
//                        FUNCTION SBETAS
//  *********************************************************************
double pennuc_specificSampler::SBETAS(int IDt,
				      pen_rand& random){
  //
  //  Random sampling of the initial beta ray energy for the IDT beta
  //  transition.
  //
  //Uses CBET1 common;
  //
  //  ****  Selection of the interval (Walker's aliasing).
  double RN = random.rand() * NPM1[IDt - 1] + 1.0E0;
  int K = int (RN);
  double TST = RN - K;
  double RR, D;
  int I;
  double SBETAS_RETURN;
  if (TST < F_[K - 1][IDt - 1])
    {
      I = K;
      RR = TST;
      D = F_[K - 1][IDt - 1];
    }
  else
    {
      I = IA[K - 1][IDt - 1];
      RR = TST - F_[K - 1][IDt - 1];
      D = 1.0E0 - F_[K - 1][IDt - 1];
    }
  //  ****  Sampling from the rational inverse cumulative distribution.
  if (RR > 1.0E-12)
    {
      SBETAS_RETURN =
	X[I - 1][IDt - 1] +
	((1.0E0 + A[I - 1][IDt - 1] + B[I - 1][IDt - 1]) * D * RR / (D * D +
								     (A[I - 1]
								      [IDt -
								       1] *
								      D +
								      B[I -
									1][IDt
									   -
									   1]
								      * RR) *
								     RR)) *
	(X[I + 1 - 1][IDt - 1] - X[I - 1][IDt - 1]);
    }
  else
    {
      SBETAS_RETURN =
	X[I - 1][IDt - 1] + random.rand() * (X[I + 1 - 1][IDt - 1] -
					     X[I - 1][IDt - 1]);
    }
  //
  return SBETAS_RETURN;
}

//  *********************************************************************
//                       SUBROUTINE SBETA0
//  *********************************************************************
void pennuc_specificSampler::SBETA0(int IMASS,
				    int IZ,
				    double EMAX,
				    int IPRO,
				    int IDt,
				    int &IER){
  //
  //    Initialisation of the algorithm for random sampling of the initial
  //  kinetic energy of a beta ray according to the theoretical Fermi spec-
  //  trum.
  //
  //  Input parameters:
  //    IMASS ... nuclear mass number of the daughter nucleus.
  //    IZ ...... nuclear charge number (atomic number) of the daughter
  //              nucleus, with negative sign for positron emitters.
  //    EMAX .... maximum kinetic energy (in eV) of beta rays.
  //    IPRO .... prohibition parameter
  //              IPRO=0: Allowed and first forbidden transitions.
  //              IPRO=1: Unique-first forbidden and second forbidden.
  //              IPRO=2: Unique-second forbidden and third forbidden.
  //              IPRO=3: Unique-third forbidden and fourth forbidden.
  //    IDT ..... transition identification number (.LE. 100)
  //
  //  Output value:
  //    IER  .... error code. When IER.NE.0 the initialisation has not
  //              succeeded; the cause of the problem is written in
  //              the standard output unit (6).
  //
  //  After invoking subroutine SBETA0 in the form
  //        CALL SBETA0(IMASS,ZI,EMAX,IPRO,IDT)
  //  the statement
  //        ERND=SBETAS(IDT)
  //  produces a random value ERND of the initial energy of the beta ray.
  //

  //Uses BETAS, CBET1 and CRITAA

  if (IDt > ID)
    {
      printf ("SBETA0. The input value of IDT is too large.\n");
      printf ("Increase the value of the paramete ID to %d\n", IDt);
      IER = 31;
      return;
    }
  //
  IER = 0;

  BETAS betas;
  
  betas.IMASSC = IMASS;
  betas.IZC = IZ;
  betas.EMAXC = EMAX;
  betas.IPROC = IPRO;
  //

  CRITAA ritaa;
  
  double XL = 0.0E0;
  double XU = EMAX;
  double ERRM;
  RITA0 (FERMI, &betas, XL, XU, NN, NN / 4, ERRM, 0, ritaa);
  if (penGetError() != PEN_SUCCESS){
    IER = -1;
    return;
  }
  NPM1[IDt - 1] = ritaa.NPM1A;
  for (int I = 1; I <= ritaa.NPM1A + 1; I++){
    X[I - 1][IDt - 1] = ritaa.XA[I - 1];
    A[I - 1][IDt - 1] = ritaa.AA[I - 1];
    B[I - 1][IDt - 1] = ritaa.BA[I - 1];
    F_[I - 1][IDt - 1] = ritaa.FA[I - 1];
    IA[I - 1][IDt - 1] = ritaa.IA[I - 1];
  }
}

//  *********************************************************************
//                       FUNCTION FERMI
//  *********************************************************************
double FERMI(double E0,
	     void* arg){
  //
  //  ****  Computation of the (unnormalised) Fermi beta spectrum.
  //        Written by E. Garcia-Torano. December 1990.
  //
  //  ****  Beta emitter data (input through common /BETAS/):
  //  IMASS ... mass number of the daughter nucleus.
  //  IZ ...... atomic number of the daughter nucleus, with negative sign
  //            for positron emitters.
  //  EMAX .... maximum kinetic energy (in eV) of beta rays.
  //  IPRO .... prohibition parameter.
  //            IPRO=0: allowed and first forbidden transitions.
  //            IPRO=1: unique-first forbidden and second forbidden.
  //            IPRO=2: unique-second forbidden and third forbidden.
  //            IPRO=3: unique-third forbidden and fourth forbidden.
  //

  //Uses BETAS common

  BETAS* betas = static_cast<BETAS*>(arg);
  
  double FERMI_RETURN;
  //
  double ALPHA = 1.0E0 / 137.036E0;
  double AZ = ALPHA * betas->IZ;
  double GAM = sqrt (1.0E0 - pow (AZ, 2));
  double A13 = pow (double (betas->IMASS), 0.333333333333333E0);
  double RAD = (1.123E0 * A13 - 0.941E0 / A13) / 386.144E0;
  double EXPN = 2.0E0 * (GAM - 1.0E0);
  double W0 = betas->EMAX / 511.0E3 + 1.0E0;
  double CTE = pow ((2.0E0 * RAD), EXPN);
  //
  double W = E0 / 511.0E3 + 1.0E0;
  double P = sqrt (pow (W, 2) - 1.0E0);
  double EXPN1 = 2.0E0 * GAM - 1.0E0;
  double RL0 =
    1.0E0 - 7.0E0 * pow (AZ,
			 2) / 20.0E0 - 28.0E0 * (W * RAD * AZ) / 15.0E0 -
    8.0E0 * (RAD * AZ) / (15.0E0 * W) - pow ((P * RAD), 2) / 3.0E0 - pow (AZ,
									  4) /
    5.0E0 + AZ * pow ((W * RAD),
		      3) * (1.0E0 + 4.5E0 * AZ) - W * RAD * pow (AZ,
								 3) / 2.0E0 +
    (pow (AZ, 6)) / 20.0E0 - 3.0E0 * pow (AZ,
					  8) / 8.0E0 -
    3.0E0 * W * RAD * pow (AZ, 4) / 8.0E0 + 10.0E0 * pow ((W * RAD),
							  2) * pow (AZ, 6);
  if (E0 <= 1.0E-6 * betas->EMAX)
    {
      //
      //  ************  Energy = 0.
      //
      //  ****  Positron emitters.
      if (betas->IZ <= 0)
	{
	  FERMI_RETURN = 0.0E0;
	  return FERMI_RETURN;
	}
      //  ****  Negatron emitters.
      double ARGIN = 2.0E0 * GAM + 1.0E0;
      double YY = EXPN1 * log (AZ) - 2.0E0 * GAMMAL(ARGIN);
      FERMI_RETURN = 8.0E0 * constants::PI * CTE * pow ((betas->EMAX / 511.0E3), 2) * RL0
	* PROHIB (betas->IPRO, W, W0, *betas) * exp (YY);
      return FERMI_RETURN;
    }
  //
  //  ************  Energy .GT. 0.
  //
  double YY = AZ * W / P;
  double XX = EXPN1 * log (P) + constants::PI * YY + RATIOG (GAM, YY);
  double FAFERM = 4.0E0 * RL0 * CTE * exp (XX);
  FERMI_RETURN = pow ((W0 - W), 2) * W * PROHIB (betas->IPRO, W, W0, *betas) * FAFERM;
  return FERMI_RETURN;
}

//  *********************************************************************
//                       FUNCTION PROHIB
//  *********************************************************************
double PROHIB(int IPRO,
	      double W,
	      double W0,
	      const BETAS& betas){
  //
  //  Computes the shape factor from the following expressions :
  //  - Allowed and first forbidden transitions (IPRO=0),
  //            F = 1
  //  - Second and unique-first forbidden (IPRO=1),
  //            F = p^2 + q^2
  //  - Third and unique-second forbidden (IPRO=2),
  //            F = p^4 + q^4 + 10(p^2*q^2)/3.
  //  - Fourth and unique-third forbidden (IPRO=3),
  //            F = p^6 + q^6 + 7 p^2 q^2 (p^2+q^2)
  //
  //  Improved shape factors for some nuclides that do not follow the
  //  standard expressions are defined explicitly.
  //
  //  Version 2.0                          E.García-Toraño, April 2016.
  //

  //Uses BETAS common
  
  double PROHIB_RETURN;
  //
  double P2 = pow (W, 2) - 1.0E0;
  double Q2 = pow ((W0 - W), 2);
  //
  //  ****  Nuclides with special shape factors.
  //
  // C-11
  if (betas.IZ == -5 && betas.IMASS == 11)
    {
      PROHIB_RETURN = 1.0E0 - 0.0074E0 * W;
      return PROHIB_RETURN;
    }
  // N-13
  if (betas.IZ == -6 && betas.IMASS == 13)
    {
      PROHIB_RETURN = 1.0E0 + 0.0014E0 * W;
      return PROHIB_RETURN;
    }
  // F-18
  if (betas.IZ == -8 && betas.IMASS == 18)
    {
      PROHIB_RETURN = 1.0E0 + 0.0034E0 * W;
      return PROHIB_RETURN;
    }
  // NA-22
  if (betas.IZ == -10 && betas.IMASS == 22)
    {
      PROHIB_RETURN = 1.0E0 - 0.005E0 * W;
      return PROHIB_RETURN;
    }
  // NA-24
  if (betas.IZ == 12 && betas.IMASS == 24)
    {
      PROHIB_RETURN = 1.0E0 - 0.011E0 * W;
      return PROHIB_RETURN;
    }
  // CL-36
  if (betas.IZ == 18 && betas.IMASS == 36)
    {
      PROHIB_RETURN =
	1.0E0 - 0.282715E0 * W - 0.045988E0 / W - 0.491739E0 * pow (W,
								    2) +
	0.438731E0 * pow (W, 3);
      return PROHIB_RETURN;
    }
  // P-32
  if (betas.IZ == 16 && betas.IMASS == 32)
    {
      PROHIB_RETURN = 1.0E0 - 0.019E0 * W;
      return PROHIB_RETURN;
    }
  // K-40
  if (betas.IZ == 20 && betas.IMASS == 40)
    {
      PROHIB_RETURN =
	1.05E0 * pow (Q2, 3) + 6.3E0 * P2 * pow (Q2,
						 2) + 6.25E0 * Q2 * pow (P2,
									 2) +
	0.95E0 * pow (P2, 3);
      return PROHIB_RETURN;
    }
  // SC-49
  if (betas.IZ == 22 && betas.IMASS == 49)
    {
      PROHIB_RETURN = 1.0E0 - 0.008E0 * W;
      return PROHIB_RETURN;
    }
  // NI-63
  if (betas.IZ == 22 && betas.IMASS == 63)
    {
      PROHIB_RETURN = 1.0E0 + 1213.49E0 * W - 232.009E0 / W
	- 1517.16E0 * pow (W, 2) + 532.69E0 * pow (W, 3);
      return PROHIB_RETURN;
    }
  // GA-68
  if (betas.IZ == -30 && betas.IMASS == 68)
    {
      PROHIB_RETURN = 1.0E0 - 0.01E0 * W;
      return PROHIB_RETURN;
    }
  // SR-89
  if (betas.IZ == 39 && betas.IMASS == 89)
    {
      PROHIB_RETURN = (P2 + Q2) * (1.0E0 - 0.0112E0 * W);
      return PROHIB_RETURN;
    }
  // SR-90
  if (betas.IZ == 39 && betas.IMASS == 90)
    {
      PROHIB_RETURN = (P2 + Q2) * (1.0E0 - 0.054E0 * W);
      return PROHIB_RETURN;
    }
  // Y-90
  if (betas.IZ == 40 && betas.IMASS == 90)
    {
      PROHIB_RETURN = (P2 + Q2) * (1.0E0 - 0.0114E0 * W);
      return PROHIB_RETURN;
    }
  // ZR-89
  if (betas.IZ == -39 && betas.IMASS == 89)
    {
      PROHIB_RETURN = 1.0E0 - 0.39E0 * W + 0.09E0 * pow (W, 2);
      return PROHIB_RETURN;
    }
  // TC-99
  if (betas.IZ == 44 && betas.IMASS == 99)
    {
      PROHIB_RETURN = 0.529E0 * P2 + Q2;
      return PROHIB_RETURN;
    }
  // I-129
  if (betas.IZ == 54 && betas.IMASS == 129)
    {
      PROHIB_RETURN = P2 + 3.16E0 * Q2;
      return PROHIB_RETURN;
    }
  // LA-138
  if (betas.IZ == 58 && betas.IMASS == 138)
    {
      PROHIB_RETURN =
	1.0E0 + 407.71E0 * W - 50.695E0 / W - 583.794E0 * pow (W,
							       2) +
	246.279E0 * pow (W, 3);
      return PROHIB_RETURN;
    }
  // CE-141
  if (betas.IZ == 59 && betas.IMASS == 141)
    {
      PROHIB_RETURN = 1.0E0 - 0.28E0 * W;
      return PROHIB_RETURN;
    }
  // RE-186
  if (betas.IZ == 56 && betas.IMASS == 186)
    {
      PROHIB_RETURN = 1.0E0 + 0.038E0 * (P2 + Q2);
      return PROHIB_RETURN;
    }
  // TL-204
  if (betas.IZ == 82 && betas.IMASS == 204)
    {
      PROHIB_RETURN = 1.097E0 * Q2 + P2;
      return PROHIB_RETURN;
    }
  // BI-210
  if (betas.IZ == 84 && betas.IMASS == 210)
    {
      PROHIB_RETURN = 1.0E0 - 0.47E0 * W + 0.065E0 * pow (W, 2);
      return PROHIB_RETURN;
    }
  // PU-241
  if (betas.IZ == 95 && betas.IMASS == 241)
    {
      PROHIB_RETURN = 1.0E0 - 1.9582E0 * W + 0.96078E0 * pow (W, 2);
      return PROHIB_RETURN;
    }
  //
  //  ****  Nuclides with standard shape factors.
  //
  if (IPRO == 0)
    {
      PROHIB_RETURN = 1.0E0;
    }
  else if (IPRO == 1)
    {
      PROHIB_RETURN = P2 + Q2;
    }
  else if (IPRO == 2)
    {
      PROHIB_RETURN = pow (P2, 2) + pow (Q2, 2) + (10.0E0 * P2 * Q2) / 3.;
    }
  else if (IPRO == 3)
    {
      PROHIB_RETURN = pow (P2, 3) + pow (Q2, 3) + 7.0E0 * P2 * Q2 * (P2 + Q2);
    }
  else
    {
      PROHIB_RETURN = 1.0E0;
    }
  //
  return PROHIB_RETURN;
}

//  *********************************************************************
//                       FUNCTION RATIOG
//  *********************************************************************
double RATIOG(double G,
	      double DNU){
  //
  //  This function calculates the square modulus of the ratio
  //              GAMMA(G+CI*DNU))/GAMMA(2*G+1),
  //  where GAMMA stands for the complex gamma function and CI is the
  //  imaginary unit.
  //
  //      const double PI=3.1415926535897932E0;
  double RATIOG_RETURN = 0.0E0;
  double Y1, S1, S2, DCOC;
  double A = 6.0E0 / (G + 5.0E0);
  for (int J = 1; J <= 5; J++)
    {
      Y1 = A * DNU;
      S1 =
	(pow (Y1, 2) + pow ((J - 1), 2)) / (pow ((G + J - 1), 2) +
					    pow (DNU, 2));
      RATIOG_RETURN = RATIOG_RETURN + log (S1);
    }
  Y1 = A * DNU;
  if ((constants::PI * Y1) > 30.0E0)
    {
      S1 = log (2.0E0 * constants::PI / Y1) - constants::PI * Y1;
    }
  else
    {
      S1 = log (constants::PI / (Y1 * sinh (constants::PI * Y1)));
    }
  S2 = 25.0E0 + pow (Y1, 2);
  RATIOG_RETURN = RATIOG_RETURN + S1 + log (S2);
  S1 = pow ((5.0E0 + G), 2) + pow (DNU, 2);
  DCOC = DNU / (5.0E0 + G);
  S1 = 2.0E0 - log (S1) + 2.0E0 * DNU * atan (DCOC) / (5.0E0 + G)
    + 1.0E0 / (S1 * 6.0E0 * A);
  RATIOG_RETURN =
    RATIOG_RETURN + (1.0E0 - G) * S1 - 11.0E0 * log (A) - 2.0E0 * GAMMAL (G +
									  G +
									  1.0E0);
  return RATIOG_RETURN;
}
  

//  *********************************************************************
//                       FUNCTION GAMMAL
//  *********************************************************************
double GAMMAL(double R){
  //
  //  This function gives LOG(GAMMA(R)) for real arguments R.
  //
  double GAMMAL_RETURN;
  GAMMAL_RETURN = 80.5E0;
  double RA = R;
  double FAC = 1.0E0;
  if (fabs (RA) <= 10.0E0)
    {
      bool Eixir = false;
      while (!Eixir)
	{
	  Eixir = true;
	  FAC = FAC / RA;
	  RA = RA + 1.0E0;
	  if (fabs (RA) < 1.0E-35)
	    {
	      return GAMMAL_RETURN;
	    }
	  if (fabs (RA) < 10.0E0)
	    {
	      Eixir = false;
	      continue;
	    }
	}
    }
  //  ****  Stirling expansion.
  double RI2 = 1.0E0 / (RA * RA);
  double RS = (43867.0E0 / 244188.0E0) * RI2;
  RS = (RS - 3617.0E0 / 122400.0E0) * RI2;
  RS = (RS + 1.0E0 / 156.0E0) * RI2;
  RS = (RS - 691.0E0 / 360360.0E0) * RI2;
  RS = (RS + 1.0E0 / 1188.0E0) * RI2;
  RS = (RS - 1.0E0 / 1680.0E0) * RI2;
  RS = (RS + 1.0E0 / 1260.0E0) * RI2;
  RS = (RS - 1.0E0 / 360.0E0) * RI2;
  RS = (RS + 1.0E0 / 12.0E0) / RA;
  GAMMAL_RETURN =
    ((RA - 0.5E0) * log (RA) - RA + 9.1893853320467274E-1 + RS) + log (FAC);
  return GAMMAL_RETURN;
}

REGISTER_SPECIFIC_SAMPLER(pennuc_specificSampler,pen_particleState, PENNUC)
