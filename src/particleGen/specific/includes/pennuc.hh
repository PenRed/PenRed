 
//
//
//    Copyright (C) 2020-2023 Universitat de València - UV
//    Copyright (C) 2020-2023 Universitat Politècnica de València - UPV
//    Copyright (C) 2025 Vicent Giménez Alventosa
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
// Universitat de Barcelona with the following disclaimer and license:

/*
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C            PPPPP   EEEEEE  N    N  N    N  U    U   CCCC             C
C            P    P  E       NN   N  NN   N  U    U  C    C            C
C            P    P  E       N N  N  N N  N  U    U  C                 C
C            PPPPP   EEEE    N  N N  N  N N  U    U  C                 C
C            P       E       N   NN  N   NN  U    U  C    C            C
C            P       EEEEEE  N    N  N    N   UUUU    CCCC             C
C                                                                      C
C                                                   (version 2018).    C
C                                                                      C
C  Subroutine package for Monte Carlo simulation of decay paths of     C
C  radioactive nuclei.                                                 C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
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
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  This subroutine package simulates the decay of radionuclides by
C  using information from the NUCLEIDE database. The decay parameters
C  of the source nuclides are read from a text file, in the dedicated
C  PenNuc format, which can be downloaded from the NUCLEIDE web site
C     http://www.nucleide.org/DDEP_WG/DDEPdata.htm
C
C  A radioactive nucleus may reach a metastable level (i.e., a level
C  with a life time longer than the detector resolution time) in the
C  course of its de-excitation. Such levels are treated as effective
C  halts of the decay process. The radiations emitted between successive
C  halts constitute a 'cascade' and are considered to define a single
C  emission event.
C
C  The decay simulation of each radioisotope is initialised by calling
C  subroutine PENNUC0, which reads the NUCLEIDE datafile of that isotope
C  and extracts atomic relaxation data from the PENELOPE database files.
C  PENNUC0 prepares tables of relevant quantities and defines the value
C  of the output parameter EPMAX (the highest energy of the emitted
C  radiations). The PENNUC subroutine package can handle up to NISM
C  radionuclides, which are identified by a label IS defined by the
C  order in which they are loaded. That is, IS=1, 2, ...  denotes the
C  first, second, ... loaded nuclides.
C
C  At initialisation, the present subroutines read information from
C  the files 'pdatconf.p14' and 'pdrelax.p11' of the PENELOPE database.
C  When running the program, these two files must be present in the
C  PENNUC data directory, whose relative path (RPATH) is defined in the
C  module PENNUC_mod.
C
*/

//*******
//** The description of each routine has been extracted from the original code
//*******

#ifndef __PEN_PENNUC_SPECIFIC_SAMPLER__
#define __PEN_PENNUC_SPECIFIC_SAMPLER__

#include <type_traits>
#ifdef _PEN_EMBEDDED_DATA_BASE_
#include "materialCreator.hh"
#endif

struct BETAS
{
  double EMAX;
  int IMASS,IZ,IPRO2;

  double& EMAXC=EMAX;
  int&    IMASSC=IMASS;
  int&    IZC=IZ;
  int&    IPROC=IPRO2;
  int&    IPRO=IPRO2;
};

class pennuc_specificSampler : public abc_specificSampler<pen_particleState>{
  DECLARE_SPECIFIC_SAMPLER(pennuc_specificSampler, pen_particleState)

  private:

  //
  //  ****  Simulation of nuclear decay. Global quantities and state
  //        variables of emitted particles. All energies in eV.
  //

  //  ----  Resolution time of the detector.
  double DRTIME;
  //  ----  Cutoff energy for the de-excitation of inner subshells.
  double ECNUC;

  int ISOT_;

  //  ----  Physical dimensions of arrays:
  static const int NISM=5;   // Max. no. of nuclides.
  static const int NDM=5;    // Max. no. of daughters.
  static const int NBM=200;  // Max. no. of branches.
  static const int NLM=200;  // Max. no. of levels.
  static const int NTM=200;  // Max. no. of transitions.
  int NIR;  // Number of loaded radionuclides.
  //
  //  ----  Daughter nucleus (IZD,IAD) and current level (LEVD).
  int IZD,IAD,LEVD;
  //  ----  Particles emitted in the last decay cascade.
  static const int NRM=500;  // Max. no. of particles.
  int NR;             // No. of particles in the cascade (<= NRM).
  //  ----  Particle types (1=electron, 2=photon, 3=positron, 4=alpha).
  int KPNR[NRM];
  //  ----  Transitions that originated the particles.
  int ITNR[NRM];
  //  ----  Particle energies, ages, and direction cosines.
  double ENR[NRM],AGENR[NRM],UNR[NRM],VNR[NRM],WNR[NRM];
  //  
  
  //CBET1 common
  static const int NN = 64;
  static const int ID = 200;
  
  double X[NN][ID],A[NN][ID],B[NN][ID],F_[NN][ID];
  int IA[NN][ID],NPM1[ID];

  // CATREL
  static const long int NRX = 60000;
  static const int NSHP1 = 31;

  double EBR[99][NSHP1],P[NRX],ET[NRX],F[NRX];
  int IAL[NRX],IS0[NRX],IS1[NRX],IS2[NRX],IFIRST[99][16],ILAST[99][16],NCUR;

  // CATRE2
  double ALW[99][30];

  //CPN1
  char NUCLIDE[NISM][8],DNUCLIDE[NISM][NDM][8];

  //CPN2
  int KTYPE[NISM][NDM][NBM];
  int NBRANCH[NISM][NDM];
  int NLEVEL[NISM][NDM];
  int IZNUC[NISM];
  int IANUC[NISM];
  int NDAUGH[NISM];

  //CPN3
  double Q[NISM][NDM],UQ[NISM][NDM],PDAUGH[NISM][NDM];
  double UPDAUGH[NISM][NDM],ELEVEL[NISM][NDM][NTM+1];

     // CARE WITH ELEVEL ELEVEL(NISM,NDM,0:NTM),UELEVEL(NISM,NDM,0:NTM)
  double UELEVEL[NISM][NDM][NTM+1],DTIME[NISM][NDM][NLM];
  double UDTIME[NISM][NDM][NLM];


  //CPN4
  double BR[NISM][NDM][NBM],UBR[NISM][NDM][NBM];
  double EBRANCH[NISM][NDM][NBM],UEBRANCH[NISM][NDM][NBM];
  int IFLB[NISM][NDM][NBM],IPRBETA[NISM][NDM][NBM],INDBETA[NISM][NDM][NBM],IDT;

  //CPN5
  double PTRANS[NISM][NDM][NLM][NTM],UPTRANS[NISM][NDM][NLM][NTM];
  double ETRANS[NISM][NDM][NLM][NTM],UETRANS[NISM][NDM][NLM][NTM];
  int NTRANS[NISM][NDM][NLM],MODE[NISM][NDM][NLM][NTM],IFLT[NISM][NDM][NLM][NTM];

  
  int lastMETAST;
  double x0,y0,z0; //Saves decay position
  double t0;       //Saves base sampled time
  bool LAGE;

  int lastIDAUGH, lastIBRANCH, lastISA;

  unsigned sourceMaterial;

  const wrapper_geometry* geometry;

  public:

  pennuc_specificSampler() : abc_specificSampler<pen_particleState>(USE_SPATIAL),
			     DRTIME(5.0E-6),
			     ECNUC(200.0E0),
			     NIR(0),
			     lastMETAST(0),LAGE(false),sourceMaterial(0),geometry(nullptr){}
  
  int configure(double& Emax,
		const pen_parserSection& config,
		const unsigned /*nthreads*/,
		const unsigned verbose);

  bool getNext(pen_particleState& state,
	       pen_KPAR& genKpar);
  
  void sample(pen_particleState& state,
	      pen_KPAR& genKpar,
	      unsigned long long& dhist,
	      pen_rand& random);

  inline void updateGeometry(const wrapper_geometry* geometryIn){ geometry=geometryIn; }
  
  void PENNUC(int& IS, int& METAST, int& IER, pen_rand& random);
  
  void ISOTROP(double& U, double& V, double& W, pen_rand& random) const;
  
  void PENNUC0(const char* NUCFNAME, const char* ATOMFNAME, const char* RELAXFNAME, const char* ATRELIFNAME, double& EPMAX, FILE* IWR, int& IER);
  
  void ATREL(int IZ, int IS, double* ES,
	      double* PAGES, int* KPARS, int* ITR, int &NS, pen_rand& random);  
  
  void ATREL0();
  
  void ATRELI(const char* ATRELIFNAME, const char* ATOMFNAME, const char* RELAXFNAME, int IZ, int IWR, int &IER);  

  double SBETAS(int IDT, pen_rand& random);
  
  void SBETA0(int IMASS, int IZ, double EMAX, int IPRO, int IDT, int &IER);

  //  *********************************************************************
  //                       SUBROUTINE READSKP
  //  *********************************************************************
  // Base case (empty)
  template<typename...>
  struct all_pointers : std::true_type {};

  // Recursive case
  template<typename T, typename... others>
  struct all_pointers<T, others...> :
    std::integral_constant<bool,
			   std::is_pointer<T>::value && all_pointers<others...>::value
			   > {};
  
  template<class... readTypes, typename = std::enable_if_t<all_pointers<readTypes...>::value>>
  bool READSKP(FILE * NRD,
	       char WCODE[4],
	       const char* format,
	       readTypes... vars) const {
    //
    //  Reads and reformats an input record. Skips comment records, and in
    //  data records replaces any ';' by ','.
    //
    //
    
    char straux[300];
    strcpy (WCODE, "COM");
    while (strcmp (WCODE, "COM") == 0)
      {
	if(fgets (straux, sizeof (straux), NRD) == nullptr){
	  return 0;
	}
	if(strlen(straux) < 3) return 0; //Detect invalid empty line
	sscanf(straux, "%3c", WCODE);
	WCODE[strlen (WCODE)] = '\0';      
      }
  
    //Get the remaining line without the code
    std::string line(&straux[3], strlen(straux)-3);

    //Create the resulting string to read
    std::string result;

    //Define white spaces
    const char* ws = " \t\n\v\f\r";

    //Init initial position for the first field
    size_t pi = 0;
    while(pi != std::string::npos){

      //Find next delimiter position
      const size_t pf = line.find(';', pi+1);

      //Check if this field contains non white characters
      const size_t pndel = line.find_first_not_of(ws, pi+1);
      if(pndel == pf){
	//Empty field, set a zero
	result.append(" 0 ");
      }else{
	//Field with data, save it
	result.append(" ");
	result.append(line, pi+1, pf-pi-1);
	result.append(" ");
      }
      //Update itial delimiter position
      pi = pf;
    }

    //Try to read data
    int nRead = sscanf(result.c_str(), format, vars...);
    if(nRead != sizeof...(readTypes)){
      printf("READSKP: PENNUC format read error at line:\n"
	     "         Original:  %s%s"
	     "         Processed: %s"
	     "         Expected %d values, %d read instead.\n",
	     WCODE, line.c_str(), result.c_str(),
	     static_cast<int>(sizeof...(readTypes)), nRead);
      return false;
    }
    return true;
  }
  
};

double FERMI(double E0, void* arg);

double PROHIB(int IPRO_, double W, double W0, const BETAS& betas);
  
double RATIOG(double G, double DNU);
  
double GAMMAL(double R);


#endif
