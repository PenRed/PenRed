
//
//
//    Copyright (C) 2019-2021 Universitat de València - UV
//    Copyright (C) 2019-2021 Universitat Politècnica de València - UPV
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


#include "quadric_geo.hh"

REGISTER_GEOMETRY(pen_quadricGeo, PEN_QUADRIC)

pen_quadricGeo::pen_quadricGeo() : NSURF(0),
  LVERB(false)
  
{
  for(unsigned int i = 0; i < NB; i++)
    {
      strcpy(bodies[i].BALIAS, "    ");
      bodies[i].MATER = 0;
      bodies[i].KDET = 0;
    }
  
}

unsigned pen_quadricGeo::getIBody(const char* elementName) const{

  //Construct corrected alias
  char auxAlias[5];
  sprintf(auxAlias,"%4.4s",elementName);
  auxAlias[4] = '\0';
  for(unsigned j = 0; j < getElements(); j++){
    //Check if body alias is the expected one
    if(strcmp(auxAlias,bodies[j].BALIAS) == 0){
      return j;
    }
  }
  return getElements();
}

std::string pen_quadricGeo::getBodyName(const unsigned ibody) const{

  if(ibody < getBodies()){
    return std::string(bodies[ibody].BALIAS);
  }else{
    return std::string("NONE");
  }
    
}


int pen_quadricGeo::configure(const pen_parserSection& config, const unsigned verbose){

  int err;
  //Read input and output file from configuration
  std::string IRDfilename;
  std::string IWRfilename("null");
  if(config.read("input-file",IRDfilename) != INTDATA_SUCCESS){
    if(verbose > 0){
      printf("quadricGeo:configure:Error: 'input-file' field missing at configuration section.\n");
    }
    configStatus = PEN_QUAD_GEO_INPUT_SECTION;
    return configStatus;
  }
  
  //Generate output processed file if verbose is greather than zero
  if(verbose > 0){
    if(config.read("processed-geo-file",IWRfilename) != INTDATA_SUCCESS){
      if(verbose > 0){
	printf("quadricGeo:configure:Error: 'processed-geo-file' field missing at configuration section\n");
      }
      configStatus = PEN_QUAD_GEO_OUTPUT_SECTION;
      return configStatus;
    }
  }

  //Open input file
  //*****************
  FILE* IRD = nullptr;
  IRD = fopen(IRDfilename.c_str(),"r");
  if(IRD == nullptr){
    if(verbose > 0){
      printf("pen_quadricGeo:configure:Error: unable to open input file '%s'.\n", IRDfilename.c_str());
    }
    configStatus = PEN_QUAD_GEO_INPUT_FILE;
    return configStatus;
  }
  
  //Open output file
  //*****************
  FILE* IWR = nullptr;
  if(verbose > 0){
    IWR = fopen(IWRfilename.c_str(),"w");
    if(IWR == nullptr){
      printf("pen_quadricGeo:configure:Error: unable to open output file '%s'.\n", IWRfilename.c_str());
      configStatus = PEN_QUAD_GEO_OUTPUT_FILE;
      fclose(IRD);
      return configStatus;
    }
  }

  if(verbose > 1){
    printf("pen_quadricGeo:configure: Input and output files:\n");
    printf("                 input  : %s\n",IRDfilename.c_str());
    printf("                 output : %s\n\n",IWRfilename.c_str());
  }
  
  //Load geometry
  //*****************
  err = GEOMIN(IRD,IWR,verbose);
  if(err != PEN_QUAD_GEO_SUCCESS){
    if(verbose > 0){
      printf("pen_quadricGeo:configure: Error loading geometry.\n");
      printf("                          Error code: %d\n",err);
    }
    fclose(IRD);
    if(IWR != nullptr){fclose(IWR);}
    return err;
  }
  fclose(IRD);
  if(IWR != nullptr){fclose(IWR);}

  //Load dsmax 
  //*****************
  std::vector<std::string> bodiesAlias;
  err = config.ls("dsmax",bodiesAlias);
  if(err != INTDATA_SUCCESS){
    if(verbose > 1){
      printf("No dsmax specified for any body\n");
    }
  }
  else{
    if(verbose > 1){
      printf("dsmax specified for %lu bodies:\n\n",bodiesAlias.size());
    }
    for(unsigned i = 0; i < bodiesAlias.size(); i++){
      bool found = false;
      std::string key("dsmax/");      
      key += bodiesAlias[i];
      
      //Construct expected alias
      char auxAlias[5];
      sprintf(auxAlias,"%4.4s",bodiesAlias[i].c_str());
      auxAlias[4] = '\0';
      for(unsigned j = 0; j < getElements(); j++){
	//Check if body alias is the expected one
	if(strcmp(auxAlias,bodies[j].BALIAS) == 0){
	  //Get DSMAX
	  double auxDSmax;
	  err = config.read(key,auxDSmax);
	  if(err != INTDATA_SUCCESS){
	    if(verbose > 0){
	      printf("pen_quadricGeo:configure: Error reading 'dsmax' of body %s\n",bodies[j].BALIAS);
	      printf("                     key: %s\n",key.c_str());
	    }
	    configStatus = PEN_QUAD_GEO_BAD_READ_DSMAX;
	    return PEN_QUAD_GEO_BAD_READ_DSMAX;
	  }

	  //Check dsmax value
	  if(auxDSmax <= 0.0){
	    if(verbose > 0){
	      printf("pen_quadricGeo:configure: Error: 'DSMAX' must be greater than zero.\n");
	      printf("            Specified for body %s: %12.4E\n",bodies[j].BALIAS,auxDSmax);
	    }
	    configStatus = PEN_QUAD_GEO_INVALID_DSMAX;
	    return PEN_QUAD_GEO_INVALID_DSMAX;
	  }
	  //Assign dsmax
	  bodies[j].DSMAX = auxDSmax;
	  if(verbose > 1){
	    printf("Set DSMAX %12.4E for body %s\n",bodies[j].DSMAX,bodies[j].BALIAS);
	  }
	  found = true;
	  break;
	}
      }
      if(!found && verbose > 1){
	printf("Warning: body '%s' not found (key %s).\n",
	       auxAlias,key.c_str());
      }
    }
  }

  if(verbose > 1)
    printf("\n\n");
  
  // Load Detectors
  //*****************
  bodiesAlias.clear();
  err = config.ls("kdet",bodiesAlias);
  if(err != INTDATA_SUCCESS){
    if(verbose > 1){
      printf("No detector specified for any body\n");
    }
  }
  else{
    if(verbose > 1){
      printf("Detector specified for %lu bodies:\n\n",bodiesAlias.size());
    }
    for(unsigned i = 0; i < bodiesAlias.size(); i++){
      std::string key("kdet/");      
      key += bodiesAlias[i];

      //Get bodi index
      unsigned bIndex = getIBody(bodiesAlias[i].c_str());
      if(bIndex < getElements()){
	
	  //Get DSMAX
	  int auxKDET;
	  err = config.read(key,auxKDET);
	  if(err != INTDATA_SUCCESS){
	    if(verbose > 0){
	      printf("pen_quadricGeo:configure: Error reading 'kdet' for body %s\n",bodies[bIndex].BALIAS);
	      printf("                     key: %s\n",key.c_str());
	    }
	    configStatus = PEN_QUAD_GEO_BAD_READ_KDET;
	    return PEN_QUAD_GEO_BAD_READ_KDET;
	  }

	  //Check kdet value
	  if(auxKDET <= 0){
	    if(verbose > 0){
	      printf("pen_quadricGeo:configure: Error: 'KDET' must be greater than zero.\n");
	      printf("            Specified for body %s: %d\n",bodies[bIndex].BALIAS,auxKDET);
	    }
	    configStatus = PEN_QUAD_GEO_INVALID_KDET;
	    return PEN_QUAD_GEO_INVALID_KDET;
	  }
	  //Assign dsmax
	  bodies[bIndex].KDET = (unsigned)auxKDET;
	  if(verbose > 1){
	    printf("Set KDET %u for body '%s'\n",bodies[bIndex].KDET,bodies[bIndex].BALIAS);
	  }
	}
      else if(verbose > 1){
	printf("Warning: body '%s' not found (key %s).\n",
	       bodiesAlias[i].c_str(),key.c_str());
      }
    }
  }

  // Load Absorption Energies
  //*************************
  bodiesAlias.clear();
  err = config.ls("eabs",bodiesAlias);
  if(err != INTDATA_SUCCESS){
    if(verbose > 1){
      printf("No absorption energies specified for any body\n");
    }
  }
  else{
    if(verbose > 1){    
      printf("Absorption energies specified for %lu bodies:\n\n",
	     bodiesAlias.size());
    }
    for(unsigned i = 0; i < bodiesAlias.size(); i++){
      std::string key("eabs/");      
      key += bodiesAlias[i];

      //Get body index
      unsigned bIndex = getIBody(bodiesAlias[i].c_str());
      if(bIndex < getElements()){
	
	//Get particle names
	std::vector<std::string> particleNames;
	config.ls(key.c_str(),particleNames);

	//Get particle eabs
	double bodyEABS[constants::nParTypes];
	for(unsigned j = 0; j < constants::nParTypes; j++)
	  bodyEABS[j] = -1.0;
	  
	for(unsigned j = 0; j < particleNames.size(); j++){

	  unsigned kpar = particleID(particleNames[j].c_str());
	  if(kpar >= ALWAYS_AT_END){
	    if(verbose > 0){
	      printf("pen_quadricGeo:configure: Error on 'eabs' field, unknown particle '%s' on body '%s'.\n",particleNames[j].c_str(),bodiesAlias[i].c_str());
	    }
	    return PEN_QUAD_GEO_UNKNOWN_PARTICLE;
	  }
	  
	  std::string key2 = key + std::string("/") + particleNames[j];
	  double eabs;
	  err = config.read(key2,eabs);
	  if(err != INTDATA_SUCCESS){
	    if(verbose > 0){
	      printf("pen_quadricGeo:configure: Error reading energy "
		     "absorption at field '%s'. Double expected.\n",
		     key2.c_str());
	    }
	    return PEN_QUAD_GEO_BAD_READ_EABS;
	  }

	  if(eabs <= 0.0){
	    if(verbose > 0){
	      printf("pen_quadricGeo:configure: Error: Invalid energy "
		     "absorption %12.4E for body '%s' particle '%s'. "
		     "Must be greater than zero.\n",
		     eabs,bodiesAlias[i].c_str(),particleNames[j].c_str());
	    }
	    return PEN_QUAD_GEO_INVALID_EABS;
	  }
	  
	  bodyEABS[kpar] = eabs;
	}

	//Set body eabs for each specified particle
	if(setBodyEabs(bIndex,bodyEABS) != 0){
	  if(verbose > 0){
	    printf("pen_quadricGeo:configure: Error on 'eabs' field, unknown body '%s'\n",bodiesAlias[i].c_str());
	  }
	  return PEN_QUAD_GEO_UNDEF_BODY_LABEL;	  
	}
	
	if(verbose > 1){
	  printf("Absotion energies (eV) for body '%s':\n",bodies[bIndex].BALIAS);
	  for(unsigned j = 0; j < constants::nParTypes; j++){
	    printf(" %20.20s: %14.5E\n",particleName(j),bodies[bIndex].localEABS[j]);
	  }
	}
      }
      else if(verbose > 1){
	printf("Warning: body '%s' not found (key %s).\n",
	       bodiesAlias[i].c_str(),key.c_str());
      }
    }
  }
  
  if(verbose > 1){
    printf("\nBodies information:\n\n");

    for(unsigned j = 0; j < getElements(); j++){
      printf("  Body '%s':\n",bodies[j].BALIAS);
      printf("    Material: %u\n",bodies[j].MATER);
      printf("    DSMAX   : %10.4E\n",bodies[j].DSMAX);
      printf("    KDET    : %u\n",bodies[j].KDET);
      printf("    EABS    :\n");
      for(unsigned k = 0; k < constants::nParTypes; k++)
	printf("  %20.20s: %14.5E\n",particleName(k),bodies[j].localEABS[k]);
    }
    printf("\n");
  }
  
  configStatus = PEN_QUAD_GEO_SUCCESS;
  return PEN_QUAD_GEO_SUCCESS;
}

//  *********************************************************************
//                       SUBROUTINE GEOMIN
//  *********************************************************************
int pen_quadricGeo::GEOMIN(FILE* IRD, FILE* IWR, const unsigned verbose)
{
  
  //     Reads the geometry-definition file and sets up the arrays used
  //  for tracking particles through the system.

  //  Input arguments:
  //     IRD ....... input unit of the geometry definition file (opened in
  //                 the main program).
  //     IWR ....... output unit (opened in the main program). At output
  //                 this file contains a duplicate of the definition file
  //                 with the effective parameter values and with elements
  //                 of the various kinds labeled in strictly increasing
  //                 order. This part of the output file describes the
  //                 actual geometry used in the simulation. After the END
  //                 line of the geometry definition block, subroutine
  //                 GEOMIN writes a detailed report with the structure of
  //                 the tree of modules and with information on redundant
  //                 (duplicated) surfaces.

  
  // ----  Improvements suggested by F. Tola and M. Moreau  ---------------
  char CHR[7];  // Used to check the format of surface indices.
  char CCR, CNL;  // Used to identify 'C\r' and 'C\n' comments.
  // ----------------------------------------------------------------------
  char ALIAS[NS][6], ALIAB[NB][6];
  char C5[6], C5C[6], C4[5], C1;
  char CA[36] = {'0','1','2','3','4','5','6','7','8','9','A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z'};

  char GFILE[13];
  char BLINE[100], DEFB[NB][100], DEFS[NS][100], DEF[100];
  const char* LNUL = "00000000";
  const char* LSUR = "SURFACE ";
  const char* LIND = "INDICES=";
  const char* LBOD = "BODY    ";
  const char* LMAT = "MATERIAL";
  const char* LMOD = "MODULE  ";
  const char* LOPEN = ")       ";
  const char* LEND = "END     ";
  const char* LONE = "11111111";
  const char* LXSC = "X-SCALE=";
  const char* LYSC = "Y-SCALE=";
  const char* LZSC = "Z-SCALE=";
  const char* LTHE = "  THETA=";
  const char* LPHI = "    PHI=";
  const char* LOME = "  OMEGA=";
  const char* LXSH = "X-SHIFT=";
  const char* LYSH = "Y-SHIFT=";
  const char* LZSH = "Z-SHIFT=";
  const char* LAXX = "    AXX=";
  const char* LAXY = "    AXY=";
  const char* LAXZ = "    AXZ=";
  const char* LAYY = "    AYY=";
  const char* LAYZ = "    AYZ=";
  const char* LAZZ = "    AZZ=";
  const char* LAX = "     AX=";
  const char* LAY = "     AY=";
  const char* LAZ = "     AZ=";
  const char* LA0 = "     A0=";
  const char* LDEG = ") DEG   ";
  const char* LRAD = ") RAD   ";
  const char* LCLON = "CLONE   ";
  const char* LSUA = "SURFACE*";
  const char* LINC = "INCLUDE ";
  const char* LINA = "INCLUDE*";
  const char* LFIL = "   FILE=";

  char LARRAY[8][9];
  int  KQ[5], KM[NS];

  int IDESC[NB];
  int IDONE[NB];
  int ISCL[NS];
  int IBCL[NB];
  int IBOR[NB];

  int KSFF[NS];
  
  char LKEYW[9];

  unsigned int KB;
  
  //Check if input file is open
  if(IRD == nullptr){
    if(verbose > 0){
      printf("pen_quadricGeo:configure:Error: The input file is "
	     "a null pointer. Does the file exist?.\n");
    }
    configStatus = PEN_QUAD_GEO_INPUT; return configStatus;    
  }

  //  ************  Initialise parameters.

  //  ****  Character constants used to identify comment lines of the forms
  //  'C\r' and 'C\n', where '\r'=CHAR(13) is the carriage return character
  //  and '\n'=CHAR(10) is the new line character.
  CCR=(char) 13;
  CNL=(char) 10;

  //  ****  Surface coefficients, alias and KPLANE.

  int NSFF = 0;  // Number of fixed (starred) surfaces.
  for(unsigned int KS = 0; KS < NS; KS++)
    {
      strcpy(ALIAS[KS], "    0"); // Surface aliases (user labels).
      KSFF[KS] = 0;  // Clonable-free surfaces.
      surfaces[KS].AXX = 0.0;
      surfaces[KS].AXY = 0.0;
      surfaces[KS].AXZ = 0.0;
      surfaces[KS].AYY = 0.0;
      surfaces[KS].AYZ = 0.0;
      surfaces[KS].AZZ = 0.0;
      surfaces[KS].AX = 0.0;
      surfaces[KS].AY = 0.0;
      surfaces[KS].AZ = 0.0;
      surfaces[KS].A0 = 0.0;
      surfaces[KS].KPLANE = 0;
      //  KPLANE=1 if the surface is a plane.
    }
  //  **** Alias, material, mother, daughters, surfaces and side pointers
  //  of bodies. The last second component of the double arrays is the
  //  number of used components. For example, KSURF(KB,NXG) is the number
  //  of surfaces that limit the body or module KB. The number of values
  //  stored in array KFLAG(KB,--), however, is given by KSURF(KB,NXG).
  for(KB = 0; KB < NB; KB++)
    {
      //  **** The array KDET(KB) is used to label bodies that are parts of
      //  impact detectors. Detectors must be defined in the main program,
      //  after the call to subroutine GEOMIN.
      bodies[KB].KDET = 0;  // KDET(KB).ne.0 if body KB is part of a detector.
      strcpy(ALIAB[KB], "    0");  // Body aliases (user labels).
      bodies[KB].MATER = 0;  // Material in body KB.
      bodies[KB].KBOMO = 0;  // 0 for bodies, 1 for modules.
      bodies[KB].KMOTH = 0;  // Mother of body KB (must be unique).

      for(unsigned int K2 = 0; K2 < NXG; K2++)
	{
	  bodies[KB].KBODY[K2] = 0;  // Limiting bodies of body KB (keep it small).
	  //  KBODY(KB,K2), is the K2-th limiting body of body KB.
	  bodies[KB].KDGHT[K2] = 0;  // Daughters of module KB (should be small).
	  //  KDGHT(KB,K2), is the K2-th daughter (body or module) of module KB.
	  //  KSURF(KB,K2), is the K2-th limiting surface of body KB.
	  //  KFLAG(KB,KS)=1, if KS is a limiting surface of KB and KB is inside KS
	  bodies[KB].setSurf(K2,0,5); // Limiting surfaces of body KB and side pointers of material bodies	  
	  //                  (i.e. side pointer =-1).
	  //  KFLAG(KB,KS)=2, if KS is a limiting surface of KB and KB is outside
	  //                  KS (i.e. side pointer =+1).
	  //  KFLAG(KB,KS)=3, if KB is a body and KS does not directly limit KB,
	  //                  but KS appears in the definition of a body that
	  //                  limits KB.
	  //  KFLAG(KB,KS)=4, if KB is a module and KS limits one of its daughters
	  //                  (bodies and submodules), but does not appear in the
	  //                  definition of KB.
	  //  KFLAG(KB,KS)=5, otherwise.
	}
    }
  NMATG = 0;
  NSURF = 0;
  NBODYS = 0;
  int ICLONE = 0;
  int KEEPL = 0;

  //  ************  Reading the geometry input file.

  FILE* IR = IRD;
  FILE* IW = IWR;
  FILE* IRI = 0;
  int NINCL = 0;
  C1 = CA[NINCL+1-1];
  
  if(IW == IR)
    {
      if(verbose > 0){
	fprintf(IW, "SUBROUTINE GEOMIN. Input arguments.\n");
	fprintf(IW, "IRD =%p,  IWR =%p\n", (void*)IRD, (void*)IWR);
	fprintf(IW, "*** The input and output units must be different.\n");
	if(verbose > 0){
	  printf("pen_quadricGeo:configure:Error: The input and output units must be different.\n");
	}
      }
      configStatus = PEN_QUAD_GEO_WR; return configStatus;
    }

  bool Eixir0 = false;   //GOTO 1
  while(!Eixir0)
    {
      Eixir0 = true;
      char AUXBL[100];
      memset(AUXBL,0,sizeof(AUXBL));
      fscanf(IR,"%[^\r\n]%*[^\n]",AUXBL);
      if(fgetc(IR)=='\r'){fgetc(IR);}
      BLINE[0]='\0';
      if(strlen(AUXBL)>72){sprintf(BLINE,"%-.72s",AUXBL);}else{sprintf(BLINE,"%-72s",AUXBL);}
      sscanf(BLINE, "%8c", LKEYW);

      //Append end of string
      LKEYW[8] = '\0';
      
      if(strcmp(LKEYW, LNUL) != 0)
	{
	  if(IR == IRD)
	    {
	      if(verbose > 0){	      
		fprintf(IW, "%s\n", BLINE);
	      }
	    }
	  else
	    {
	      if(verbose > 0){	      
		fprintf(IW, "C %s\n", BLINE);
	      }
	    }
	  Eixir0 = false;
	  continue;
	}
      else
	{
	  if(IR == IRD)
	    {
	      if(verbose > 0){	      
		fprintf(IW, "%0*d\n",64,0);
	      }
	    }
	}
    }  

  bool Eixir = false;      //GOTO 2
  while(!Eixir)
    {
      Eixir = true;
      BLINE[0]='\0';
      fscanf(IR,"%[^\r\n]%*[^\n]",BLINE);
      if(fgetc(IR)=='\r'){fgetc(IR);}
      if(((BLINE[0] == 'C' || BLINE[0] == 'c') && (BLINE[1] == '\0' || BLINE[1] == CCR || BLINE[1] == CNL || BLINE[1] == ' ')) || BLINE[0] == '#')
	{
	  char AUXBL[100];
	  strcpy(AUXBL,BLINE);
	  if(strlen(AUXBL)>72){sprintf(BLINE,"%-.72s",AUXBL);}else{sprintf(BLINE,"%-72s",AUXBL);}
	  if(verbose > 0){	      
	    fprintf(IW, "%s\n", BLINE);
	  }
	  Eixir = false;
	  continue;
	}     
      for(int ii=0; ii<8; ii++){memset(LARRAY[ii],0,sizeof(LARRAY[ii]));}
      int Num_Elements_Llegits = sscanf(BLINE, "%8c%*c%4c%8c%8c%8c%8c%8c%8c%8c", LKEYW, C4, LARRAY[0], LARRAY[1], LARRAY[2], LARRAY[3], LARRAY[4], LARRAY[5], LARRAY[6]);

      //Append end of string chars
      LKEYW[8] = '\0';
      C4[4] = '\0';
      LARRAY[0][8] = '\0';
      LARRAY[1][8] = '\0';
      LARRAY[2][8] = '\0';
      LARRAY[3][8] = '\0';
      LARRAY[4][8] = '\0';
      LARRAY[5][8] = '\0';
      LARRAY[6][8] = '\0';
	  
      if(Num_Elements_Llegits == 0)
	{
	  if(verbose > 0){	      
	    if(strlen(BLINE)>72){fprintf(IW, "%-.72s\n", BLINE);}else{fprintf(IW, "%-72s\n", BLINE);}
	    fprintf(IW, "*** Wrong input format.\n");
	    if(verbose > 1){
	      printf("pen_quadricGeo:configure:Error: Wrong input format.\n");
	    }
	  }
	  configStatus = PEN_QUAD_GEO_INPUT; return configStatus;
	}
      
      char AuxStr[9];strcpy(AuxStr,"        ");
      for(int ii=0; ii<8; ii++)
	{
	  sprintf(AuxStr, "%-8.8s", LARRAY[ii]);
	  strcpy(LARRAY[ii],AuxStr);

	}

      if(strcmp(LKEYW, LSUR) == 0 || strcmp(LKEYW, LSUA) == 0)   //GOTO 100
	{
	  //	  
	  //  ************  Surfaces.
	  //	 
	  if(BLINE[8] != '(' || BLINE[13] != ')')
	    {
	      if(verbose > 0){	      
		if(strlen(BLINE)>72){fprintf(IW, "%-.72s\n", BLINE);}else{fprintf(IW, "%-72s\n", BLINE);}
		fprintf(IW, "*** Incorrect label format.\n");
		if(verbose > 1){
		  printf("pen_quadricGeo:configure:Error: Incorrect label format.\n");
		}
	      }
	      configStatus = PEN_QUAD_GEO_LABEL_FORMAT; return configStatus;
	    }
	  Num_Elements_Llegits = sscanf(BLINE, "%*9c%4c%8c%8c%8c%8c%8c%8c%8c", C4, LARRAY[0], LARRAY[1], LARRAY[2], LARRAY[3], LARRAY[4], LARRAY[5], LARRAY[6]);

	  //Append end of string chars
	  C4[4] = '\0';
	  LARRAY[0][8] = '\0';
	  LARRAY[1][8] = '\0';
	  LARRAY[2][8] = '\0';
	  LARRAY[3][8] = '\0';
	  LARRAY[4][8] = '\0';
	  LARRAY[5][8] = '\0';
	  LARRAY[6][8] = '\0';
	    
	  if(Num_Elements_Llegits == 0)
	    {
	      if(verbose > 0){	      
		if(strlen(BLINE)>72){fprintf(IW, "%-.72s\n", BLINE);}else{fprintf(IW, "%-72s\n", BLINE);}
		fprintf(IW, "*** Wrong input format.\n");
		if(verbose > 1){
		  printf("pen_quadricGeo:configure:Error: Wrong input format.\n");
		}
	      }
	      configStatus = PEN_QUAD_GEO_INPUT; return configStatus;
	    }
	  if(IR == IRD)
	    {
	      strcpy(C5, C4);
	      strcat(C5, "0");
	    }
	  else if(KEEPL == 1)
	    {
	      strcpy(C5, C4);
	      strcat(C5, "0");
	    }
	  else
	    {
	      strcpy(C5, C4);
	      strcat(C5, &C1);
	    }
	  if(NSURF > 0)
	    {
	      for(unsigned int KS0 = 0; KS0 < NSURF; KS0++)
		{
		  if(strcmp(C5,ALIAS[KS0]) == 0)
		    {
		      if(verbose > 0){	      
			if(strlen(BLINE)>72){fprintf(IW, "%-.72s\n", BLINE);}else{fprintf(IW, "%-72s\n", BLINE);}
			fprintf(IW, "*** Same label for two surfaces.\n");
			if(verbose > 1){
			  printf("pen_quadricGeo:configure:Error: Same label for two surfaces.\n");
			}
		      }
		      configStatus = PEN_QUAD_GEO_SAME_LABEL_SURF; return configStatus;
		    }
		}
	    }
	  NSURF = NSURF+1;
	  unsigned int KS = NSURF;
	  if(verbose > 0){	      
	    fprintf(IW, "%s(%4d%s%s%s%s%s%s%s\n", LKEYW, KS, LARRAY[0], LARRAY[1], LARRAY[2], LARRAY[3], LARRAY[4], LARRAY[5], LARRAY[6]);
	  }
	  if(KS > NS)
	    {
	      if(verbose > 0){	      
		fprintf(IW, "*** The parameter NS must be increased.\n");
		if(verbose > 1){
		  printf("pen_quadricGeo:configure:Error: The parameter NS must be increased.\n");
		}
	      }
	      configStatus = PEN_QUAD_GEO_NS; return configStatus;	      
	    }
	  sprintf(DEF, "%s(%4d%s%s%s%s%s%s%s", LKEYW, KS, LARRAY[0], LARRAY[1], LARRAY[2], LARRAY[3], LARRAY[4], LARRAY[5], LARRAY[6]);
	  strcpy(ALIAS[KS-1], C5);
	  if(strcmp (LKEYW, LSUA) == 0)
	    {
	      KSFF[KS-1] = 1;
	      NSFF = NSFF+1;
	    }
	  //  ****  Indices.
	  bool Eixir2 = false;
	  while(!Eixir2)
	    {
	      Eixir2 = true;

        BLINE[0]='\0';
        fscanf(IR,"%[^\r\n]%*[^\n]",BLINE);
        if(fgetc(IR)=='\r'){fgetc(IR);}
	    if(((BLINE[0] == 'C' || BLINE[0] == 'c') && (BLINE[1] == '\0' || BLINE[1] == CCR || BLINE[1] == CNL || BLINE[1] == ' ')) || BLINE[0] == '#')
		{
		  if(verbose > 0){	      
		    if(strlen(BLINE)>72){fprintf(IW, "%-.72s\n", BLINE);}else{fprintf(IW, "%-72s\n", BLINE);}
		  }
		  Eixir2 = false;
		  continue;
		}
	    }
	  Num_Elements_Llegits = sscanf(BLINE, "%8c%c%2d%c%2d%c%2d%c%2d%c%2d%c", LKEYW, &CHR[0], &KQ[0], &CHR[1], &KQ[1], &CHR[2], &KQ[2], &CHR[3], &KQ[3], &CHR[4], &KQ[4], &CHR[5]);

	  //Append end of string chars
	  LKEYW[8] = '\0';
	  CHR[6] = '\0';
	  
	  if(Num_Elements_Llegits == 0)
	    {
	      if(verbose > 0){	      
		if(strlen(BLINE)>72){fprintf(IW, "%-.72s\n", BLINE);}else{fprintf(IW, "%-72s\n", BLINE);}
		fprintf(IW, "*** Wrong input format.\n");
		if(verbose > 1){
		  printf("pen_quadricGeo:configure:Error: Wrong input format.\n");
		}
	      }
	      configStatus = PEN_QUAD_GEO_INPUT; return configStatus;
	    }

	  if(verbose > 0){	      
	    fprintf(IW, "%s%c%2d%c%2d%c%2d%c%2d%c%2d%c\n", LKEYW, CHR[0], KQ[0], CHR[1], KQ[1], CHR[2], KQ[2], CHR[3], KQ[3], CHR[4], KQ[4], CHR[5]);
	  }
	  //  ****  Check parentheses and commas.
	  if(CHR[0] != '(' || CHR[1] != ',' || CHR[2] != ',' || CHR[3] != ',' || CHR[4] != ',' || CHR[5] != ')')
	    {
	      if(verbose > 0){	      
		fprintf(IW, "*** Incorrect format of surface indices.\n");
		if(verbose > 1){
		  printf("pen_quadricGeo:configure:Error: Incorrect format of surface indices.\n");
		}
	      }
	      configStatus = PEN_QUAD_GEO_SURF_FORMAT; return configStatus;
	    }
	  //  ****  Now check the values of the indices:

	  int IMODE = abs(KQ[0]);
	  if(IMODE < abs(KQ[1])){ IMODE = abs(KQ[1]);} 
	  if(IMODE < abs(KQ[2])){ IMODE = abs(KQ[2]);}
	  if(IMODE < abs(KQ[3])){ IMODE = abs(KQ[3]);} 	  
	  if(IMODE < abs(KQ[4])){ IMODE = abs(KQ[4]);}
	  
	  if(IMODE > 1 || strcmp(LKEYW, LIND) != 0)
	    {
	      if(verbose > 0){	      
		fprintf(IW, "*** Incorrect surface indices.\n");
		if(verbose > 1){
		  printf("pen_quadricGeo:configure:Error: Incorrect surface indices.\n");
		}
	      }
	      configStatus = PEN_QUAD_GEO_SURF_IND; return configStatus;
	    }
	  if(IMODE != 0)
	    {

	      //  ****  Reduced form.

	      double XSCALE = 1.0;
	      double YSCALE = 1.0;
	      double ZSCALE = 1.0;
	      double OMEGA = 0.0;
	      double THETA = 0.0;
	      double PHI = 0.0;
	      double XSHIFT = 0.0;
	      double YSHIFT = 0.0;
	      double ZSHIFT = 0.0;
	      //  ****  Parameters of the quadric.

	      Eixir2 = false;
	      while(!Eixir2)        //GOTO 101
		{
		  Eixir2 = true;
		  memset(BLINE,0,sizeof(BLINE));
          fscanf(IR,"%[^\r\n]%*[^\n]",BLINE);
          if(fgetc(IR)=='\r'){fgetc(IR);}
		  if(((BLINE[0] == 'C' || BLINE[0] == 'c') && (BLINE[1] == '\0' || BLINE[1] == CCR || BLINE[1] == CNL || BLINE[1] == ' ')) || BLINE[0] == '#')
		    {
		      if(verbose > 0){	      
			if(strlen(BLINE)>72){fprintf(IW, "%-.72s\n", BLINE);}else{fprintf(IW, "%-72s\n", BLINE);}
		      }
		      Eixir2 = false;
		      continue;
		    }
		  double VALUE;
		  int ICHPAR;
		  char LANGLE[9]; 
		  char AUXSTR[9];strcpy(AUXSTR,"        ");
		  Num_Elements_Llegits = sscanf(BLINE, "%8c%*c%lf%*c%4d%8c", LKEYW, &VALUE, &ICHPAR, AUXSTR);
		  
		  //Append end of string chars
		  LKEYW[8] = '\0';
		  AUXSTR[8] = '\0';
		    
		  sprintf(LANGLE,"%-8s",AUXSTR);		  
		  if(strcmp(LKEYW, LNUL) == 0){ break;}   //GOTO 102
		  if(Num_Elements_Llegits == 0 )
		    {
		      if(verbose > 0){	      
			fprintf(IW, "%72s\n", BLINE);
			fprintf(IW, "*** Wrong input format.\n");
			if(verbose > 1){
			  printf("pen_quadricGeo:configure:Error: Wrong input format.\n");
			}
		      }
		      configStatus = PEN_QUAD_GEO_INPUT; return configStatus;
		    }
		  
		  //  ****  Transformation parameters.
		  if(strcmp(LKEYW, LXSC) == 0)
		    {
		      if(verbose > 0){	      
			fprintf(IW, "%s(%22.15E,%4d%s\n", LKEYW, VALUE, ICHPAR, LOPEN);
		      }
		      if(VALUE < 1.0E-15)
			{
			  if(verbose > 0){	      
			    fprintf(IW, "*** Scale factor less than 1.0E-15.\n");
			    if(verbose > 1){
			      printf("pen_quadricGeo:configure:Error: X-scale factor less than 1.0E-15.\n");
			    }
			  }
			  configStatus = PEN_QUAD_GEO_XSCALE; return configStatus;
			}
		      XSCALE = VALUE;
		    }
		  else if(strcmp(LKEYW, LYSC) == 0)
	  	    {
		      if(verbose > 0){	      
			fprintf(IW, "%s(%22.15E,%4d%s\n", LKEYW, VALUE, ICHPAR, LOPEN);
		      }
		      if(VALUE < 1.0E-15)
			{
			  if(verbose > 0){	      
			    fprintf(IW, "*** Scale factor less than 1.0E-15.\n");
			    if(verbose > 1){
			      printf("pen_quadricGeo:configure:Error: Y-scale factor less than 1.0E-15.\n");
			    }
			  }
			  configStatus = PEN_QUAD_GEO_YSCALE; return configStatus;
		        }
		      YSCALE = VALUE;
		    }
		  else if(strcmp(LKEYW, LZSC) == 0)
		    {
		      if(verbose > 0){	      
			fprintf(IW, "%s(%22.15E,%4d%s\n", LKEYW, VALUE, ICHPAR, LOPEN);
		      }
			  
		      if(VALUE < 1.0E-15)
			{
			  if(verbose > 0){	      
			    fprintf(IW, "*** Scale factor less than 1.0E-15.\n");
			    if(verbose > 1){
			      printf("pen_quadricGeo:configure:Error: Z-scale factor less than 1.0E-15.\n");
			    }
			  }
			  configStatus = PEN_QUAD_GEO_ZSCALE; return configStatus;
			}
		      ZSCALE = VALUE;
		    }
		  else if(strcmp(LKEYW, LOME) == 0)
		    {
		      if(strcmp(LANGLE, LRAD) == 0)
			{
			  if(verbose > 0){	      
			    fprintf(IW, "%s(%22.15E,%4d%s\n", LKEYW, VALUE, ICHPAR, LRAD);
			  }
			  OMEGA = VALUE;
			}
		      else
		        {
			  if(verbose > 0){	      
			    fprintf(IW, "%s(%22.15E,%4d%s\n", LKEYW, VALUE, ICHPAR, LDEG);
			  }
			  OMEGA = VALUE*constants::PI/180.0;			
			}
		    }
		  else if(strcmp(LKEYW, LTHE) == 0)
		    {
		      if(strcmp(LANGLE, LRAD) == 0)
			{
			  if(verbose > 0){	      
			    fprintf(IW, "%s(%22.15E,%4d%s\n", LKEYW, VALUE, ICHPAR, LRAD);
			  }
			  THETA = VALUE;
			}
		      else
		        {
			  if(verbose > 0){	      
			    fprintf(IW, "%s(%22.15E,%4d%s\n", LKEYW, VALUE, ICHPAR, LDEG);
			  }
			  THETA = VALUE*constants::PI/180.0;
			}
		    }
		  else if(strcmp(LKEYW, LPHI) == 0)
		    {
		      if(strcmp(LANGLE, LRAD) == 0)
			{
			  if(verbose > 0){	      
			    fprintf(IW, "%s(%22.15E,%4d%s\n", LKEYW, VALUE, ICHPAR, LRAD);
			  }
			  PHI = VALUE;
			}
		      else
		        {
			  if(verbose > 0){	      
			    fprintf(IW, "%s(%22.15E,%4d%s\n", LKEYW, VALUE, ICHPAR, LDEG);
			  }
			  PHI = VALUE*constants::PI/180.0;
			}
		    }
		  else if(strcmp(LKEYW, LXSH) == 0)
		    {
		      if(verbose > 0){	      
			fprintf(IW, "%s(%22.15E,%4d%s\n", LKEYW, VALUE, ICHPAR, LOPEN);
		      }
		      XSHIFT = VALUE;
		    }
		  else if(strcmp(LKEYW, LYSH) == 0)
		    {
		      if(verbose > 0){	      
			fprintf(IW, "%s(%22.15E,%4d%s\n", LKEYW, VALUE, ICHPAR, LOPEN);
		      }
		      YSHIFT = VALUE;
		    }
		  else if(strcmp(LKEYW, LZSH) == 0)
		    {
		      if(verbose > 0){	      
			fprintf(IW, "%s(%22.15E,%4d%s\n", LKEYW, VALUE, ICHPAR, LOPEN);
		      }
		      ZSHIFT = VALUE;
		    }
		  else
		    {
		      if(verbose > 0){	      
			fprintf(IW, "%72s\n", BLINE);
			fprintf(IW, "*** What do you mean?\n");
			if(verbose > 1){

			  printf("pen_quadricGeo:configure:Error: %72s\n",BLINE);
			  printf("pen_quadricGeo:configure:Error: What do you mean?\n");
			}
		      }
		      configStatus = PEN_QUAD_GEO_MEAN; return configStatus;
		    }
		  Eixir2 = false;
		  continue;
		}

	      //  ****  Expanded quadric.
	      
	      double QXX = KQ[0]/pow(XSCALE,2);
	      double QXY = 0.0;
	      double QXZ = 0.0;
	      double QYY = KQ[1]/pow(YSCALE,2);
	      double QYZ = 0.0;
	      double QZZ = KQ[2]/pow(ZSCALE,2);
	      double QX = 0.0;
	      double QY = 0.0;
	      double QZ = KQ[3]/ZSCALE;
	      double Q0 = KQ[4];
	      //  ****  Rotated and shifted quadric.
	      ROTSHF(OMEGA,THETA,PHI,XSHIFT,YSHIFT,ZSHIFT,QXX,QXY,QXZ,QYY,QYZ,QZZ,QX,QY,QZ,Q0);

	      double Aux = (QXX < QXY ? QXX : QXY);
	      if(Aux > QXZ){ Aux = QXZ;}
	      if(Aux > QYY){ Aux = QYY;}
	      if(Aux > QYZ){ Aux = QYZ;}
	      if(Aux > QZZ){ Aux = QZZ;}
	      Aux = fabs(Aux);

	      double Aux2 = (QXX > QXY ? QXX : QXY);
	      if(Aux2 < QXZ){ Aux2 = QXZ;}
	      if(Aux2 < QYY){ Aux2 = QYY;}
	      if(Aux2 < QYZ){ Aux2 = QYZ;}
	      if(Aux2 < QZZ){ Aux2 = QZZ;}
	      Aux2 = fabs(Aux2);

	      if(Aux < Aux2){ Aux = Aux2;}

	      if(Aux < 1.0E-30){ surfaces[KS-1].KPLANE = 1;}

	      surfaces[KS-1].AXX = QXX;
	      surfaces[KS-1].AXY = QXY;
	      surfaces[KS-1].AXZ = QXZ;
	      surfaces[KS-1].AYY = QYY;
	      surfaces[KS-1].AYZ = QYZ;
	      surfaces[KS-1].AZZ = QZZ;
	      surfaces[KS-1].AX = QX;
	      surfaces[KS-1].AY = QY;
	      surfaces[KS-1].AZ = QZ;
	      surfaces[KS-1].A0 = Q0;
	      if(verbose > 0){	      
		fprintf(IW, "%0*d\n",64,0);
	      }
	      strcpy(DEFS[KS-1],DEF);
	      Eixir = false;
	      continue;       //GOTO 2
	    }
	  //  ****  Implicit form.

  	  double QXX = 0.0;
	  double QXY = 0.0;
	  double QXZ = 0.0;
  	  double QYY = 0.0;
	  double QYZ = 0.0;
	  double QZZ = 0.0;
  	  double QX = 0.0;
	  double QY = 0.0;
	  double QZ = 0.0;
  	  double Q0 = 0.0;

	  bool Goto104 = false;
	  bool Goto107 = false;
  	  Eixir2 = false;
	  while(!Eixir2)       //GOTO 193
	    {
	      Eixir2 = true;
	      memset(BLINE,0,sizeof(BLINE));
          fscanf(IR,"%[^\r\n]%*[^\n]",BLINE);
          if(fgetc(IR)=='\r'){fgetc(IR);}
	      if(((BLINE[0] == 'C' || BLINE[0] == 'c') && (BLINE[1] == '\0' || BLINE[1] == CCR || BLINE[1] == CNL || BLINE[1] == ' ')) || BLINE[0] == '#')
		{
		  if(verbose > 0){	      
		    if(strlen(BLINE)>72){fprintf(IW, "%-.72s\n", BLINE);}else{fprintf(IW, "%-72s\n", BLINE);}
		  }
		  Eixir2 = false;
		  continue;
		}
	      double VALUE;
	      int ICHPAR;
	      char LANGLE[9];
	      char AUXSTR[9];strcpy(AUXSTR,"        ");
	      int Nombre_Elements_Escrits = sscanf(BLINE, "%8c%*c%lf%*c%4d%8c", LKEYW, &VALUE, &ICHPAR, AUXSTR);

	      //Append end of string chars
	      LKEYW[8] = '\0';
	      AUXSTR[8] = '\0';
	      
	      sprintf(LANGLE,"%-8s",AUXSTR);
	      if(strcmp(LKEYW, LNUL) == 0){ Eixir2 = true; Goto104 = false; Goto107 = true; break;}
	      if(strcmp(LKEYW, LONE) == 0){ Eixir2 = true; Goto104 = true ; Goto107 = true; break;}
	      if(Nombre_Elements_Escrits == 0)
		{
		  if(verbose > 0){	      
		    fprintf(IW, "%72s\n", BLINE);
		    fprintf(IW, "*** Wrong input format.\n");
		    if(verbose > 1){		    
		      printf("pen_quadricGeo:configure:Error: %72s\n",BLINE);
		      printf("pen_quadricGeo:configure:Error: Wrong input format.\n");
		    }
		  }
		  configStatus = PEN_QUAD_GEO_INPUT; return configStatus;
		}

	      
	      if(strcmp(LKEYW, LAXX) == 0)
		{
		  QXX = VALUE;
		  if(verbose > 0){	      
		    fprintf(IW, "%s(%22.15E,%4d%s\n", LAXX, QXX, ICHPAR, LOPEN);
		  }
		}
	      else if(strcmp(LKEYW, LAXY) == 0)
		{
		  QXY = VALUE;
		  if(verbose > 0){	      
		    fprintf(IW, "%s(%22.15E,%4d%s\n", LAXY, QXY, ICHPAR, LOPEN);
		  }
		}
	      else if(strcmp(LKEYW, LAXZ) == 0)
		{
		  QXZ = VALUE;
		  if(verbose > 0){	      
		    fprintf(IW, "%s(%22.15E,%4d%s\n", LAXZ, QXZ, ICHPAR, LOPEN);
		  }
		}
	      else if(strcmp(LKEYW, LAYY) == 0)
		{
		  QYY = VALUE;
		  if(verbose > 0){	      
		    fprintf(IW, "%s(%22.15E,%4d%s\n", LAYY, QYY, ICHPAR, LOPEN);
		  }
		}
	      else if(strcmp(LKEYW, LAYZ) == 0)
		{
		  QYZ = VALUE;
		  if(verbose > 0){	      
		    fprintf(IW, "%s(%22.15E,%4d%s\n", LAYZ, QYZ, ICHPAR, LOPEN);
		  }
		}
	      else if(strcmp(LKEYW, LAZZ) == 0)
		{
		  QZZ = VALUE;
		  if(verbose > 0){	      
		    fprintf(IW, "%s(%22.15E,%4d%s\n", LAZZ, QZZ, ICHPAR, LOPEN);
		  }
		}
	      else if(strcmp(LKEYW, LAX) == 0)
		{
		  QX = VALUE;
		  if(verbose > 0){	      
		    fprintf(IW, "%s(%22.15E,%4d%s\n", LAX, QX, ICHPAR, LOPEN);
		  }
		}
	      else if(strcmp(LKEYW, LAY) == 0)
		{
		  QY = VALUE;
		  if(verbose > 0){	      
		    fprintf(IW, "%s(%22.15E,%4d%s\n", LAY, QY, ICHPAR, LOPEN);
		  }
		}
	      else if(strcmp(LKEYW, LAZ) == 0)
		{
		  QZ = VALUE;
		  if(verbose > 0){	      
		    fprintf(IW, "%s(%22.15E,%4d%s\n", LAZ, QZ, ICHPAR, LOPEN);
		  }
		}
	      else if(strcmp(LKEYW, LA0) == 0)
		{
		  Q0 = VALUE;
		  if(verbose > 0){	      
		    fprintf(IW, "%s(%22.15E,%4d%s\n", LA0, Q0, ICHPAR, LOPEN);
		  }
		}
	      else
		{
		  if(verbose > 0){	      
		    fprintf(IW, "%72s\n", BLINE);
		    fprintf(IW, "*** What do you mean?\n");
		    if(verbose > 1){		    
		      printf("pen_quadricGeo:configure:Error: %72s\n",BLINE);
		      printf("pen_quadricGeo:configure:Error: What do you mean?\n");
		    }
		  }
		  configStatus = PEN_QUAD_GEO_MEAN; return configStatus;
		}
	      Eixir2 = false;
	      continue;
	    }    //END GOTO 193			    
		  
	  //  ****  Transformation parameters.
	  if(Goto104)
	    {
	      if(verbose > 0){	      
		fprintf(IW, "1111111111111111111111111111111111111111111111111111111111111111\n");
	      }
	      double OMEGA = 0.0;
	      double THETA = 0.0;
	      double PHI = 0.0;
	      double XSHIFT = 0.0;
	      double YSHIFT = 0.0;
	      double ZSHIFT = 0.0;

	      bool Eixir3 = false;
	      while(!Eixir3)
		{
		  Eixir3 = true;
		  memset(BLINE,0,sizeof(BLINE));
          fscanf(IR,"%[^\r\n]%*[^\n]",BLINE);
          if(fgetc(IR)=='\r'){fgetc(IR);}
		  if(((BLINE[0] == 'C' || BLINE[0] == 'c') && (BLINE[1] == '\0' || BLINE[1] == CCR || BLINE[1] == CNL || BLINE[1] == ' ')) || BLINE[0]  == '#')
		    {
		      if(verbose > 0){	      
			if(strlen(BLINE)>72){fprintf(IW, "%-.72s\n", BLINE);}else{fprintf(IW, "%-72s\n", BLINE);}
		      }
		      Eixir3 = false;
		      continue;
		    }
		  double VALUE;
		  int ICHPAR;
		  char LANGLE[9];
  	          char AUXSTR[9];strcpy(AUXSTR,"        ");
		  int Nombre_Elements_Escrits = sscanf(BLINE, "%8c%*c%lf%*c%4d%8c", LKEYW, &VALUE, &ICHPAR, AUXSTR);

		  //Append end of string chars
		  LKEYW[8] = '\0';
		  AUXSTR[8] = '\0';
		  
	    	  sprintf(LANGLE,"%-8s",AUXSTR);
		  if(Nombre_Elements_Escrits == 0)
		    {
		      if(verbose > 0){	      
			fprintf(IW, "%72s\n", BLINE);
			fprintf(IW, "*** Wrong input format.\n");
			if(verbose > 1){		    
			  printf("pen_quadricGeo:configure:Error: %72s\n",BLINE);
			  printf("pen_quadricGeo:configure:Error: Wrong input format\n");
			}
		      }
		      configStatus = PEN_QUAD_GEO_INPUT; return configStatus;
		    }
		  if(strcmp(LKEYW, LNUL) == 0){ break;}		      

		  

		  if(strcmp(LKEYW, LOME) == 0)
		    {
		      if(strcmp(LANGLE, LRAD) == 0)
		        {
			  if(verbose > 0){	      
			    fprintf(IW, "%s(%22.15E,%4d%s\n", LKEYW, VALUE, ICHPAR, LRAD);
			  }
			  OMEGA = VALUE;
			}
		      else
			{
			  if(verbose > 0){	      
			    fprintf(IW, "%s(%22.15E,%4d%s\n", LKEYW, VALUE, ICHPAR, LDEG);
			  }
			  OMEGA = VALUE*constants::PI/180.0;
			}
  		    }
		  else if(strcmp(LKEYW, LTHE) == 0)
	  	    {
		      if(strcmp(LANGLE, LRAD) == 0)
			{
			  if(verbose > 0){	      
			    fprintf(IW, "%s(%22.15E,%4d%s\n", LKEYW, VALUE, ICHPAR, LRAD);
			  }
			  THETA = VALUE;
			}
		      else
			{
			  if(verbose > 0){	      
			    fprintf(IW, "%s(%22.15E,%4d%s\n", LKEYW, VALUE, ICHPAR, LDEG);
			  }
			  THETA = VALUE*constants::PI/180.0;
			}
  		    }
		  else if(strcmp(LKEYW, LPHI) == 0)
		    {
		      if(strcmp(LANGLE, LRAD) == 0)
			{
			  if(verbose > 0){	      
			    fprintf(IW, "%s(%22.15E,%4d%s\n", LKEYW, VALUE, ICHPAR, LRAD);
			  }
			  PHI = VALUE;
			}
		      else
			{
			  if(verbose > 0){	      
			    fprintf(IW, "%s(%22.15E,%4d%s\n", LKEYW, VALUE, ICHPAR, LDEG);
			  }
			  PHI = VALUE*constants::PI/180.0;
			}
  		    }
		  else if(strcmp(LKEYW, LXSH) == 0)
	  	    {
		      if(verbose > 0){	      
			fprintf(IW, "%s(%22.15E,%4d%s\n", LKEYW, VALUE, ICHPAR, LOPEN);
		      }
		      XSHIFT = VALUE;
		    }
		  else if(strcmp(LKEYW, LYSH) == 0)
	  	    {
		      if(verbose > 0){	      
			fprintf(IW, "%s(%22.15E,%4d%s\n", LKEYW, VALUE, ICHPAR, LOPEN);
		      }
		      YSHIFT = VALUE;
		    }
		  else if(strcmp(LKEYW, LZSH) == 0)
	  	    {
		      if(verbose > 0){	      
			fprintf(IW, "%s(%22.15E,%4d%s\n", LKEYW, VALUE, ICHPAR, LOPEN);
		      }
		      ZSHIFT = VALUE;
		    }
		  else
	  	    {
		      if(verbose > 0){	      
			fprintf(IW, "%72s\n", BLINE);
			fprintf(IW, "*** What do you mean?\n");
			if(verbose > 1){		    
			  printf("pen_quadricGeo:configure:Error: %72s\n",BLINE);
			  printf("pen_quadricGeo:configure:Error: What do you mean?\n");
			}
		      }
		      configStatus = PEN_QUAD_GEO_MEAN; return configStatus;
		    }
		  Eixir3 = false;
		  continue;
		}  //END GOTO 105
	      //  ****  Rotation and translation of the surface.
	      ROTSHF(OMEGA,THETA,PHI,XSHIFT,YSHIFT,ZSHIFT,QXX,QXY,QXZ,QYY,QYZ,QZZ,QX,QY,QZ,Q0);

	    }

	  if(Goto107)
	    {
	      double TSTL = (QXX < QXY ? QXX : QXY);
	      if(TSTL > QXZ){ TSTL = QXZ;}
	      if(TSTL > QYY){ TSTL = QYY;}
	      if(TSTL > QYZ){ TSTL = QYZ;}
	      if(TSTL > QZZ){ TSTL = QZZ;}

	      double TSTU = (QXX > QXY ? QXX : QXY);
	      if(TSTU < QXZ){ TSTU = QXZ;}
	      if(TSTU < QYY){ TSTU = QYY;}
	      if(TSTU < QYZ){ TSTU = QYZ;}
	      if(TSTU < QZZ){ TSTU = QZZ;}

	      double Aux = fabs(TSTL);
	      double Aux2 = fabs(TSTU);
	      if(Aux < Aux2){ Aux = Aux2;}
	  
	      if(Aux < 1.0E-30){ surfaces[KS-1].KPLANE = 1;}
	      surfaces[KS-1].AXX = QXX;
	      surfaces[KS-1].AXY = QXY;
	      surfaces[KS-1].AXZ = QXZ;
	      surfaces[KS-1].AYY = QYY;
	      surfaces[KS-1].AYZ = QYZ;
	      surfaces[KS-1].AZZ = QZZ;
	      surfaces[KS-1].AX = QX;
	      surfaces[KS-1].AY = QY;
	      surfaces[KS-1].AZ = QZ;
	      surfaces[KS-1].A0 = Q0;
	      if(verbose > 0){	      
		fprintf(IW, "%0*d\n",64,0);
	      }
	      strcpy(DEFS[KS-1],DEF);
	      Eixir = false;
	      continue;
	    }
	}
      else if(strcmp(LKEYW, LBOD) == 0) //GOTO 200
	{
	  //	  
	  //  ************  Bodies.
	  //	  
	  if(BLINE[8] != '(' || BLINE[13] != ')')
	    {
	      if(verbose > 0){	      
		if(strlen(BLINE)>72){fprintf(IW, "%-.72s\n", BLINE);}else{fprintf(IW, "%-72s\n", BLINE);}
		fprintf(IW,"*** Incorrect label format.\n");
		if(verbose > 1){
		  printf("pen_quadricGeo:configure:Error: Incorrect label format.\n");
		}
	      }
	      configStatus = PEN_QUAD_GEO_LABEL_FORMAT; return configStatus;
	    }
	  int Nombre_Elements_Escrits = sscanf(BLINE, "%*9c%4c%8c%8c%8c%8c%8c%8c%8c", C4, LARRAY[0], LARRAY[1], LARRAY[2], LARRAY[3], LARRAY[4], LARRAY[5], LARRAY[6]);

	  //Append end of string chars
	  C4[4] = '\0';
	  LARRAY[0][8] = '\0';
	  LARRAY[1][8] = '\0';
	  LARRAY[2][8] = '\0';
	  LARRAY[3][8] = '\0';
	  LARRAY[4][8] = '\0';
	  LARRAY[5][8] = '\0';
	  LARRAY[6][8] = '\0';
	  
	  if(Nombre_Elements_Escrits == 0)
	    {
	      if(verbose > 0){	      
		if(strlen(BLINE)>72){fprintf(IW, "%-.72s\n", BLINE);}else{fprintf(IW, "%-72s\n", BLINE);}
		fprintf(IW, "*** Wrong input format.\n");
		if(verbose > 1){
		  printf("pen_quadricGeo:configure:Error: Wrong input format.\n");
		}
	      }
	      configStatus = PEN_QUAD_GEO_INPUT; return configStatus;
	    } 
	  if(IR == IRD)
	    {
	      strcpy(C5, C4);
	      strcat(C5, "0");
	    }
	  else if(KEEPL == 1)
	    {
	      strcpy(C5, C4);
	      strcat(C5, "0");
	    }
	  else
	    {
	      strcpy(C5, C4);
	      strcat(C5, &C1);
	    }
	  if(NBODYS > 0)
	    {
	      for(unsigned int KB0 = 0; KB0 < NBODYS; KB0++)
		{
		  if(strcmp(C5, ALIAB[KB0]) == 0)
		    {
		      if(verbose > 0){	      
			if(strlen(BLINE)>72){fprintf(IW, "%-.72s\n", BLINE);}else{fprintf(IW, "%-72s\n", BLINE);}
			fprintf(IW, "*** Same label for two bodies (or modules).\n");
			if(verbose > 1){
			  printf("pen_quadricGeo:configure:Error: Same label for two bodies (or modules).\n");
			}
		      }
		      configStatus = PEN_QUAD_GEO_SAME_LABEL_BODY; return configStatus;		      
		    }
		}
	    }
	  NBODYS = NBODYS+1;
	  
	  if(verbose > 0){	      
	    fprintf(IW, "%s(%4d%s%s%s%s%s%s%s\n", LKEYW, NBODYS, LARRAY[0], LARRAY[1], LARRAY[2], LARRAY[3], LARRAY[4], LARRAY[5], LARRAY[6]);
	  }
	  if(NBODYS > NB)
	    {
	      if(verbose > 0){	      
		fprintf(IW, "*** The parameter NB must be increased.\n");
		if(verbose > 1){
		  printf("pen_quadricGeo:configure:Error: The parameter NB must be increased.\n");
		}
	      }
	      configStatus = PEN_QUAD_GEO_NB; return configStatus;
	    }
	  
	  sprintf(DEF, "%s(%4d%s%s%s%s%s%s%s", LKEYW, NBODYS, LARRAY[0], LARRAY[1], LARRAY[2], LARRAY[3], LARRAY[4], LARRAY[5], LARRAY[6]);
	  strcpy(ALIAB[NBODYS-1],C5);
	  
	  
	  bool Eixir2 = false;
	  while(!Eixir2)    //GOTO 295
	    {
	      Eixir2 = true;
	      memset(BLINE,0,sizeof(BLINE));
          fscanf(IR,"%[^\r\n]%*[^\n]",BLINE);
          if(fgetc(IR)=='\r'){fgetc(IR);}
	      if(((BLINE[0] == 'C' || BLINE[0] == 'c') && (BLINE[1] == '\0' || BLINE[1] == CCR || BLINE[1] == CNL || BLINE[1] == ' ')) || BLINE[0] == '#')
		{
		  if(verbose > 0){	      
		    if(strlen(BLINE)>72){fprintf(IW, "%-.72s\n", BLINE);}else{fprintf(IW, "%-72s\n", BLINE);}
		  }

		  Eixir2 = false; continue;
		}
	    }
	  int IMAT;
	  Nombre_Elements_Escrits = sscanf(BLINE, "%8c%*c%4d", LKEYW, &IMAT);

	  //Append end of string chars
	  LKEYW[8] = '\0';
	  
	  if(Nombre_Elements_Escrits == 0)
	    {
	      if(verbose > 0){	      
		if(strlen(BLINE)>72){fprintf(IW, "%-.72s\n", BLINE);}else{fprintf(IW, "%-72s\n", BLINE);}
		fprintf(IW, "*** Wrong input format.\n");
		if(verbose > 1){
		  printf("pen_quadricGeo:configure:Error: Wrong input format.\n");
		}
	      }
	      configStatus = PEN_QUAD_GEO_INPUT; return configStatus;
	    }
	  if(strcmp(LKEYW, LMAT) != 0)
	    {
	      if(verbose > 0){	      
		if(strlen(BLINE)>72){fprintf(IW, "%-.72s\n", BLINE);}else{fprintf(IW, "%-72s\n", BLINE);}
		fprintf(IW, "*** Incorrect material definition line.\n");
		if(verbose > 1){
		  printf("pen_quadricGeo:configure:Error: Incorrect material definition line.\n");
		}
	      }
	      configStatus = PEN_QUAD_GEO_MAT; return configStatus;
	    }
	  if(verbose > 0){	      
	    fprintf(IW, "%s(%4d)\n", LKEYW, IMAT);
	  }
	  if(IMAT < 0){ IMAT = 0;}
	  bodies[NBODYS-1].MATER = IMAT;
	  if((int)NMATG < IMAT){ NMATG = IMAT;}
	  
	  Eixir2 = false;
	  while(!Eixir2)   //GOTO 201
	    {
	      Eixir2 = true;
	      memset(BLINE,0,sizeof(BLINE));
          fscanf(IR,"%[^\r\n]%*[^\n]",BLINE);
          if(fgetc(IR)=='\r'){fgetc(IR);}
	      if(((BLINE[0] == 'C' || BLINE[0] == 'c') && (BLINE[1] == '\0' || BLINE[1] == CCR || BLINE[1] == CNL || BLINE[1] == ' ')) || BLINE[0] == '#')
		{
		  if(verbose > 0){	      
		    if(strlen(BLINE)>72){fprintf(IW, "%-.72s\n", BLINE);}else{fprintf(IW, "%-72s\n", BLINE);}
		  }
		  Eixir2 = false;
		  continue;
		}
	      sscanf(BLINE, "%8c", LKEYW);

	      //Append end of string chars
	      LKEYW[8] = '\0';
	      
	      if(strcmp(LKEYW, LNUL) == 0){ break;}
	      Nombre_Elements_Escrits = sscanf(BLINE, "%*9c%4c", C4);

	      //Append end of string chars
	      C4[4] = '\0';
	      
	      if(Nombre_Elements_Escrits == 0)
		{
		  if(verbose > 0){	      
		    if(strlen(BLINE)>72){fprintf(IW, "%-.72s\n", BLINE);}else{fprintf(IW, "%-72s\n", BLINE);}
		    fprintf(IW, "*** Wrong input format.\n");
		    if(verbose > 1){
		      printf("pen_quadricGeo:configure:Error: Wrong input format.\n");
		    }
		  }
		  configStatus = PEN_QUAD_GEO_INPUT; return configStatus;
		}
	      if(IR == IRD)
		{
		  strcpy(C5, C4);
		  strcat(C5, "0");
		}
	      else if(KEEPL == 1)
		{
		  strcpy(C5, C4);
		  strcat(C5, "0");
		}
	      else
		{
		  strcpy(C5, C4);
		  strcat(C5, &C1);		  
		}
	      if(strcmp(LKEYW, LSUR) == 0 || strcmp(LKEYW, LSUA) == 0)
		{
		  int INDS;
		  Nombre_Elements_Escrits = sscanf(BLINE, "%*30c%2d", &INDS);
		  if(Nombre_Elements_Escrits == 0)
		    {
		      if(verbose > 0){	      
			if(strlen(BLINE)>72){fprintf(IW, "%-.72s\n", BLINE);}else{fprintf(IW, "%-72s\n", BLINE);}
			fprintf(IW, "*** Wrong input format.\n");
			if(verbose > 1){
			  printf("pen_quadricGeo:configure:Error: Wrong input format.\n");
			}
		      }
		      configStatus = PEN_QUAD_GEO_INPUT; return configStatus;
		    }
		  //  ****  Surface.
		  unsigned int KS;
		  bool Eixir3 = false;
		  for(unsigned int KS0 = 0; KS0 < NSURF; KS0++)
		    {        
		      if(strcmp(C5, ALIAS[KS0]) == 0)
			{
			  KS = KS0+1;
			  Eixir3 = true;
			  break;
			}
		    }
		  if(!Eixir3)
		    {
		      if(verbose > 0){	      
			fprintf(IW, "%72s\n", BLINE);
			fprintf(IW, "*** Undefined surface label.\n");
			if(verbose > 1){
			  printf("pen_quadricGeo:configure:Error: %72s\n", BLINE);
			  printf("pen_quadricGeo:configure:Error: Undefined surface label.\n");
			}
		      }
		      configStatus = PEN_QUAD_GEO_UNDEF_SURF_LABEL; return configStatus;
		    }
		      
		  if(verbose > 0){	      
		    fprintf(IW, "%s(%4d), SIDE POINTER=(%2d)\n", LKEYW, KS, INDS);
		  }
		    
		  unsigned int KST = bodies[NBODYS-1].getKSURF(NXG-1);
		  Eixir3 = false;
		  if(KST > 0)
		    {
		      int KST0 = KST;
		      for(int K = 0; K < KST0; K++)			
		        {
			  if(KS == bodies[NBODYS-1].getKSURF(K))
			    {
			      if(bodies[NBODYS-1].getKFLAG(K) < 3)
				{
				  if(verbose > 0){	      
				    fprintf(IW, "*** The last limiting surface has been defined twice.\n");
				    if(verbose > 1){
				      printf("pen_quadricGeo:configure:Error: The last limiting surface has been defined twice.\n");
				    }
				  }
				  configStatus = PEN_QUAD_GEO_LIMIT_SURF_DEF_TWICE; return configStatus;
				}
			      else
				{
				  KST = K+1;  // KS limits a limiting body.
				  Eixir3 = true;
				  break;
				}
			    }
			}
		    }
		  if(!Eixir3)
		    {
		      KST = KST+1;
		      if(KST >= NXG)
			{
			  if(verbose > 0){	      
			    fprintf(IW, "*** The number of limiting surfaces is too large.\n");
			    if(verbose > 1){
			      printf("pen_quadricGeo:configure:Error: The number of limiting surfaces is too large.\n");
			    }
			  }
			  configStatus = PEN_QUAD_GEO_MANY_LIMIT_SURFACE; return configStatus;
			}
		      
		      bodies[NBODYS-1].setKSURF(NXG-1,KST);
		      bodies[NBODYS-1].setKSURF(KST-1,KS);
		    }
		      
		  if(INDS == -1)
		    {
		      bodies[NBODYS-1].setKFLAG(KST-1,1);
		    }
		  else if(INDS == 1)
		    {
		      bodies[NBODYS-1].setKFLAG(KST-1,2);
		    }
		  else
		    {
		      if(verbose > 0){	      
			if(strlen(BLINE)>72){fprintf(IW, "%-.72s\n", BLINE);}else{fprintf(IW, "%-72s\n", BLINE);}
			fprintf(IW, "*** Check side pointer value.\n");
			if(verbose > 1){
			  printf("pen_quadricGeo:configure:Error: Check side pointer value.\n");
			}
		      }
		      configStatus = PEN_QUAD_GEO_SIDE_POINTER; return configStatus;
		    }
		}
	      else if(strcmp(LKEYW, LBOD) == 0)
		{
		  //  ****  Body.
		  if(NBODYS > 1)
		    {
		      bool Eixir3 = false;
		      for(unsigned int KB0 = 0; KB0 < NBODYS-1; KB0++)
		        {
			  if(strcmp(C5, ALIAB[KB0]) == 0)
			    {
			      KB = KB0+1;
			      Eixir3 = true;
			      break;
			    }
			}
		      if(!Eixir3)
			{
			  if(verbose > 0){	      
			    if(strlen(BLINE)>72){fprintf(IW, "%-.72s\n", BLINE);}else{fprintf(IW, "%-72s\n", BLINE);}
			    fprintf(IW, "*** Undefined body label.\n");
			    if(verbose > 1){
			      printf("pen_quadricGeo:configure:Error: Undefined body label.\n");
			    }
			  }
			  configStatus = PEN_QUAD_GEO_UNDEF_BODY_LABEL; return configStatus;
			}
		    }
		      
		  if(verbose > 0){	      
		    fprintf(IW, "%s(%4d)\n", LKEYW, KB);
		  }
		  if(bodies[KB-1].KBOMO != 0)
		    {
		      if(verbose > 0){	      
			fprintf(IW, "*** This body is a module.\n");
			if(verbose > 1){
			  printf("pen_quadricGeo:configure:Error: This body is a module.\n");
			}
		      }
		      configStatus = PEN_QUAD_GEO_BODY_IS_MODULE; return configStatus;
		    }
		      
		  int KN1 = bodies[KB-1].getKSURF(NXG-1);
		  for(int KS1 = 0; KS1 < KN1; KS1++)
		    {
		      bool Eixir3 = false;

		      int KSURF1 = bodies[KB-1].getKSURF(KS1);
		      int KN2 = bodies[NBODYS-1].getKSURF(NXG-1);
		      
		      for(int KS2 = 0; KS2 < KN2; KS2++)
			{
			  int KSURF2 = bodies[NBODYS-1].getKSURF(KS2);
			  if(KSURF2 == KSURF1){ Eixir3 = true; break;}
			}
		      if(Eixir3){ continue;}
		      unsigned int KST = KN2+1;
		      if(KST >= NXG)
			{
			  if(verbose > 0){	      
			    if(strlen(BLINE)>72){fprintf(IW, "%-.72s\n", BLINE);}else{fprintf(IW, "%-72s\n", BLINE);}
			    fprintf(IW, "*** The number of limiting surfaces is too large.\n");
			    if(verbose > 1){
			      printf("pen_quadricGeo:configure:Error: The number of limiting surfaces is too large.\n");
			    }
			  }
			  configStatus = PEN_QUAD_GEO_MANY_LIMIT_SURFACE; return configStatus;
			}

		      bodies[NBODYS-1].setKSURF(NXG-1,KST);
		      bodies[NBODYS-1].setSurf(KST-1,KSURF1,3);
		      
		    }
		      
		  bodies[NBODYS-1].KBODY[NXG-1] = bodies[NBODYS-1].KBODY[NXG-1]+1;
		  bodies[NBODYS-1].KBODY[bodies[NBODYS-1].KBODY[NXG-1]-1] = KB;
		}
	      else if(strcmp(LKEYW, LMOD) == 0)
		{
		  //  ****  Module.
		  if(NBODYS > 1)
		    {
		      bool Eixir3 = false;
		      for(unsigned int KB0 = 0; KB0 < NBODYS-1; KB0++)
			{
			  if(strcmp(C5, ALIAB[KB0]) == 0)
			    {
			      KB = KB0+1;
			      Eixir3 = true;
			      break;
			    }
			}
		      if(!Eixir3)
			{
			  if(verbose > 0){	      
			    if(strlen(BLINE)>72){fprintf(IW, "%-.72s\n", BLINE);}else{fprintf(IW, "%-72s\n", BLINE);}
			    fprintf(IW, "*** Undefined body label.\n");
			    if(verbose > 1){
			      printf("pen_quadricGeo:configure:Error: Undefined body label.\n");
			    }
			  }
			  configStatus = PEN_QUAD_GEO_UNDEF_BODY_LABEL; return configStatus;
			}
		    }
		      
		  if(verbose > 0){	      
		    fprintf(IW, "%s(%4d)\n", LKEYW, KB);
		  }
		  if(bodies[KB-1].KBOMO != 1)
		    {
		      if(verbose > 0){	      
			fprintf(IW, "*** This module is a body.\n");
			if(verbose > 1){
			  printf("pen_quadricGeo:configure:Error: This module is a body.\n");
			}
		      }
		      configStatus = PEN_QUAD_GEO_MODULE_IS_BODY; return configStatus;
		    }
		      
		  int KN1 = bodies[KB-1].getKSURF(NXG-1);
		  for(int KS1 = 0; KS1 < KN1; KS1++)
		    {
		      int KSURF1 = bodies[KB-1].getKSURF(KS1);
		      if(bodies[KB-1].getKFLAG(KS1) > 2){ continue;}
		      int KN2 = bodies[NBODYS-1].getKSURF(NXG-1);
		      bool Eixir3 = false;
		      for(int KS2 = 0; KS2 < KN2; KS2++)
			{
			  int KSURF2 = bodies[NBODYS-1].getKSURF(KS2);
			  if(KSURF2 == KSURF1){ Eixir3 = true; break;}
			}
		      if(Eixir3){ continue;}
		      unsigned int KST = KN2+1;
		      if(KST >= NXG)
			{
			  if(verbose > 0){	      
			    if(strlen(BLINE)>72){fprintf(IW, "%-.72s\n", BLINE);}else{fprintf(IW, "%-72s\n", BLINE);}
			    fprintf(IW, "*** The number of limiting surfaces is too large.\n");
			    if(verbose > 1){
			      printf("pen_quadricGeo:configure:Error: The number of limiting surfaces is too large.\n");
			    }
			  }
			  configStatus = PEN_QUAD_GEO_MANY_LIMIT_SURFACE; return configStatus;
			}

		      bodies[NBODYS-1].setKSURF(NXG-1,KST);
		      bodies[NBODYS-1].setSurf(KST-1,KSURF1,3);
		      
		    }
		      
		  bodies[NBODYS-1].KBODY[NXG-1] = bodies[NBODYS-1].KBODY[NXG-1]+1;
		  bodies[NBODYS-1].KBODY[bodies[NBODYS-1].KBODY[NXG-1]-1] = KB;
		}
	      else
		{
		  if(verbose > 0){	      
		    if(strlen(BLINE)>72){fprintf(IW, "%-.72s\n", BLINE);}else{fprintf(IW, "%-72s\n", BLINE);}
		    fprintf(IW, "*** What do you mean?\n");
		    if(verbose > 1){
		      printf("pen_quadricGeo:configure:Error: What do you mean?\n");
		    }
		  }
		  configStatus = PEN_QUAD_GEO_MEAN; return configStatus;
		}
	      Eixir2 = false;
	      continue;
	    }
	  if(verbose > 0){	      
	    fprintf(IW, "%0*d\n",64,0);
	  }
	  strcpy(DEFB[NBODYS-1],DEF);
	  Eixir = false;      //GOTO 2
	  continue;
	}
      else if(strcmp(LKEYW, LMOD) == 0)    //GOTO 300
	{
	  //		
	  //  ************  Modules.
	  //	    
	  if(BLINE[8] != '(' || BLINE[13] != ')')
	    {
	      if(verbose > 0){	      
		if(strlen(BLINE)>72){fprintf(IW, "%-.72s\n", BLINE);}else{fprintf(IW, "%-72s\n", BLINE);}
		fprintf(IW, "*** Incorrect label format.\n");
		if(verbose > 1){
		  printf("pen_quadricGeo:configure:Error: Incorrect label format.\n");
		}
	      }
	      configStatus = PEN_QUAD_GEO_LABEL_FORMAT; return configStatus;
	    }
	  int Nombre_Elements_Escrits = sscanf(BLINE, "%*9c%4c%8c%8c%8c%8c%8c%8c%8c", C4, LARRAY[0], LARRAY[1], LARRAY[2], LARRAY[3], LARRAY[4], LARRAY[5], LARRAY[6]);

	  //Append end of string chars
	  C4[4] = '\0';
	  LARRAY[0][8] = '\0';
	  LARRAY[1][8] = '\0';
	  LARRAY[2][8] = '\0';
	  LARRAY[3][8] = '\0';
	  LARRAY[4][8] = '\0';
	  LARRAY[5][8] = '\0';
	  LARRAY[6][8] = '\0';
	  
	  if(Nombre_Elements_Escrits == 0)
	    {
	      if(verbose > 0){	      
		fprintf(IW, "%72s\n", BLINE);
		fprintf(IW, "*** Wrong input format.\n");
		if(verbose > 1){
		  printf("pen_quadricGeo:configure:Error: %72s\n", BLINE);
		  printf("pen_quadricGeo:configure:Error: Wrong input format.\n");
		}
	      }
	      configStatus = PEN_QUAD_GEO_INPUT; return configStatus;
	    }
	  if(IR == IRD)
	    {
	      strcpy(C5, C4);
	      strcat(C5, "0");
	    }
	  else if(KEEPL == 1)
	    {
	      strcpy(C5, C4);
	      strcat(C5, "0");
	    }
	  else
	    {
	      strcpy(C5, C4);
	      strcat(C5, &C1);	  
	    }
	  if(NBODYS > 0)
	    {
	      for(unsigned int KB0 = 0; KB0 < NBODYS; KB0++)
		{
		  if(strcmp(C5, ALIAB[KB0]) == 0)
		    {
		      if(verbose > 0){	      
			fprintf(IW, "%72s\n", BLINE);
			fprintf(IW, "*** Same label for two bodies (or modules).\n");
			if(verbose > 1){
			  printf("pen_quadricGeo:configure:Error: %72s\n", BLINE);
			  printf("pen_quadricGeo:configure:Error: Same label for two bodies (or modules).\n");
			}
		      }
		      configStatus = PEN_QUAD_GEO_SAME_LABEL_BODY; return configStatus;
		    }
		}
	    }
	  NBODYS = NBODYS+1;
	  if(verbose > 0){	      
	    fprintf(IW, "%s(%4d%s%s%s%s%s%s%s\n", LKEYW, NBODYS, LARRAY[0], LARRAY[1], LARRAY[2], LARRAY[3], LARRAY[4], LARRAY[5], LARRAY[6]);
	  }
	  if(NBODYS > NB)
	    {
	      if(verbose > 0){	      
		fprintf(IW, "*** The parameter NB must be increased.\n");
		if(verbose > 1){
		  printf("pen_quadricGeo:configure:Error: The parameter NB must be increased.\n");
		}
	      }
	      configStatus = PEN_QUAD_GEO_NB; return configStatus;
	    }
	  sprintf(DEF, "%s(%4d%s%s%s%s%s%s%s", LKEYW, NBODYS, LARRAY[0], LARRAY[1], LARRAY[2], LARRAY[3], LARRAY[4], LARRAY[5], LARRAY[6]);      

	  strcpy(ALIAB[NBODYS-1], C5);

	  bool Eixir2 = false;
	  while(!Eixir2)
	    {
	      Eixir2 = true;
	      memset(BLINE,0,sizeof(BLINE));
          fscanf(IR,"%[^\r\n]%*[^\n]",BLINE);
          if(fgetc(IR)=='\r'){fgetc(IR);}
	      if(((BLINE[0] == 'C' || BLINE[0] == 'c') && (BLINE[1] == '\0' || BLINE[1] == CCR || BLINE[1] == CNL || BLINE[1] == ' ')) || BLINE[0] == '#')
		{
		  if(verbose > 0){	      
		    if(strlen(BLINE)>72){fprintf(IW, "%-.72s\n", BLINE);}else{fprintf(IW, "%-72s\n", BLINE);}
		  }
		  Eixir2 = false;
		  continue;
		}
	    }
	  int IMAT;
	  Nombre_Elements_Escrits = sscanf(BLINE, "%8c%*c%4d", LKEYW, &IMAT);

	  //Append end of string chars
	  LKEYW[8] = '\0';
	  
	  if(Nombre_Elements_Escrits == 0)
	    {
	      if(verbose > 0){	      
		fprintf(IW, "%72s\n", BLINE);
		fprintf(IW, "*** Wrong input format.\n");
		if(verbose > 1){
		  printf("pen_quadricGeo:configure:Error: %72s\n", BLINE);
		  printf("pen_quadricGeo:configure:Error: Wrong input format.\n");
		}
	      }
	      configStatus = PEN_QUAD_GEO_INPUT; return configStatus;
	    }
	  if(verbose > 0){	      
	    fprintf(IW, "%s(%4d)\n", LKEYW, IMAT);
	  }
	  if(strcmp(LKEYW, LMAT) != 0)
	    {
	      if(verbose > 0){	      
		fprintf(IW, "*** Incorrect material definition line.\n");
		if(verbose > 1){
		  printf("pen_quadricGeo:configure:Error: Incorrect material definition line.\n");
		}
	      }
	      configStatus = PEN_QUAD_GEO_MAT; return configStatus;
	    }
	  if(IMAT < 0){ IMAT = 0;}
	  bodies[NBODYS-1].MATER = IMAT;
	  if((int)NMATG < IMAT){ NMATG = IMAT;}

	  bodies[NBODYS-1].KDGHT[NXG-1] = 1;
	  bodies[NBODYS-1].KDGHT[0] = NBODYS;

	  //  ****  Limiting surfaces and components.

	  bool Goto308=false;
	  bool Goto312=false;

	  Eixir2 = false;
	  while(!Eixir2)  //GOTO 301
	    {
	      Eixir2 = true;
	      memset(BLINE,0,sizeof(BLINE));
          fscanf(IR,"%[^\r\n]%*[^\n]",BLINE);
          if(fgetc(IR)=='\r'){fgetc(IR);}
	      if(((BLINE[0] == 'C' || BLINE[0] == 'c') && (BLINE[1] == '\0' || BLINE[1] == CCR || BLINE[1] == CNL || BLINE[1] == ' ')) || BLINE[0] == '#')
		{
		  if(verbose > 0){	      
		    if(strlen(BLINE)>72){fprintf(IW, "%-.72s\n", BLINE);}else{fprintf(IW, "%-72s\n", BLINE);}
		  }
		  Eixir2 = false;
		  continue;
		}
	      sscanf(BLINE, "%8c", LKEYW);

	      //Append end of string chars
	      LKEYW[8] = '\0';
	      
	      if(strcmp(LKEYW, LNUL) == 0 || strcmp(LKEYW, LONE) == 0)
		{		
		  int KDT = bodies[NBODYS-1].KDGHT[NXG-1];
		  //  ****  Sort daughters in increasing order.
		  if(KDT > 1)
		    {
		      for(int KI = 0; KI < KDT-1; KI++)
			{
			  unsigned int KBMIN = bodies[NBODYS-1].KDGHT[KI];
			  int KMIN = KI+1;
			  for(int KJ = KI+1; KJ < KDT; KJ++)
			    {
			      if(bodies[NBODYS-1].KDGHT[KJ] < KBMIN)
				{
				  KBMIN = bodies[NBODYS-1].KDGHT[KJ];
				  KMIN = KJ+1;
				}
			    }
			  if(KMIN != KI+1)
			    {
			      int KSAVE = bodies[NBODYS-1].KDGHT[KI];
			      bodies[NBODYS-1].KDGHT[KI] = bodies[NBODYS-1].KDGHT[KMIN-1];
			      bodies[NBODYS-1].KDGHT[KMIN-1] = KSAVE;
			    }
			}
		    }
		  if(strcmp(LKEYW, LONE) == 0){Goto308=true ; Goto312=true; Eixir2=true; break;}
		  if(strcmp(LKEYW, LNUL) == 0){Goto308=false; Goto312=true; Eixir2=true; break;}
		}

	      //  ****  Limiting surface.
	      if(strcmp(LKEYW, LSUR) == 0 || strcmp(LKEYW, LSUA) == 0)
		{
		  int INDS;
		  Nombre_Elements_Escrits = sscanf(BLINE, "%*9c%4c%*17c%2d", C4, &INDS);

		  //Append end of string chars
		  C4[4] = '\0';
		  
		  if(Nombre_Elements_Escrits == 0)
		    {
		      if(verbose > 0){	      
			fprintf(IW, "%72s\n", BLINE);
			fprintf(IW, "*** Wrong input format.\n");
			if(verbose > 1){
			  printf("pen_quadricGeo:configure:Error: %72s\n", BLINE);
			  printf("pen_quadricGeo:configure:Error: Wrong input format.\n");
			}
		      }
		      configStatus = PEN_QUAD_GEO_INPUT; return configStatus;
		    }
		  if(IR == IRD)
		    {
		      strcpy(C5, C4);
		      strcat(C5, "0");
		    }
		  else if(KEEPL == 1)
		    {
		      strcpy(C5, C4);
		      strcat(C5, "0");
		    }
		  else
		    {
		      strcpy(C5, C4);
		      strcat(C5, &C1);
		    }
		  bool Eixir3 = false; //GOTO 302
		  unsigned int KS = 0;
		  for(unsigned int KS0 = 0; KS0 < NSURF; KS0++)
		    {
		      if(strcmp(C5, ALIAS[KS0]) == 0)
		        {
			  KS = KS0+1;
			  Eixir3 = true;
			  break;
			}
		    }
		  if(!Eixir3)
		    {
		      if(verbose > 0){	      
			fprintf(IW, "%s(%4s), SIDE POINTER=(%2d)\n", LKEYW, C5, INDS);
			fprintf(IW, "*** Undefined surface label.\n");
		      }
		      //PEN_QUAD_GEO_UNDEF_SURF_LABEL
		    }
			  
		  if(verbose > 0){	      
		    fprintf(IW, "%s(%4d), SIDE POINTER=(%2d)\n", LKEYW, KS, INDS);  
		  }
	      
		  unsigned int KST = bodies[NBODYS-1].getKSURF(NXG-1);
		  Eixir3 = false;    //GOTO 303
		  if(KST > 0) // Check whether KS is in the list.
		    {
		      unsigned int KST0 = KST;
		      for(unsigned int K = 0; K < KST0; K++)
			{
			  if(KS == bodies[NBODYS-1].getKSURF(K))
			    {
			      if(bodies[NBODYS-1].getKFLAG(K) < 3)
				{
				  if(verbose > 0){	      
				    fprintf(IW, "*** The last limiting surface has been defined twice.\n");
				    if(verbose > 1){
				      printf("pen_quadricGeo:configure:Error: The last limiting surface has been defined twice.\n");
				    }
				  }
				  configStatus = PEN_QUAD_GEO_LIMIT_SURF_DEF_TWICE; return configStatus;
				}
			      else
				{
				  KST = K+1;
				  Eixir3 = true;
				  break;			
				}
			    }
			}
		    }
		  if(!Eixir3)
		    {
		      KST = KST+1;
		      if(KST >= NXG)
			{
			  if(verbose > 0){	      
			    fprintf(IW, "*** The number of limiting surfaces is too large.\n");
			    if(verbose > 1){
			      printf("pen_quadricGeo:configure:Error: The number of limiting surfaces is too large.\n");
			    }
			  }
			  configStatus = PEN_QUAD_GEO_MANY_LIMIT_SURFACE; return configStatus;			  
			}

		      bodies[NBODYS-1].setKSURF(NXG-1,KST);
		      bodies[NBODYS-1].setKSURF(KST-1,KS);
		      
		    }
			  
		  if(INDS == -1)
		    {
		      bodies[NBODYS-1].setKFLAG(KST-1,1);
		    }
		  else if(INDS == 1)
		    {
		      bodies[NBODYS-1].setKFLAG(KST-1,2);
		    }
		  else
		    {
		      if(verbose > 0){	      
			fprintf(IW, "*** Check side pointer value.\n");
			if(verbose > 1){
			  printf("pen_quadricGeo:configure:Error: Check side pointer value.\n");
			}
		      }
		      configStatus = PEN_QUAD_GEO_SIDE_POINTER; return configStatus;
		    }
			  
		  //  ****  Body.
		}
	      else if(strcmp(LKEYW, LBOD) == 0)
		{
		  Nombre_Elements_Escrits = sscanf(BLINE, "%*9c%4c", C4);

		  //Append end of string chars
		  C4[4] = '\0';
		  
		  if(Nombre_Elements_Escrits == 0)
		    {
		      if(verbose > 0){	      
			fprintf(IW, "%72s\n", BLINE );			
			fprintf(IW, "*** Wrong input format.\n");
			if(verbose > 1){
			  printf("pen_quadricGeo:configure:Error: %72s\n", BLINE);
			  printf("pen_quadricGeo:configure:Error: Wrong input format.\n");
			}
		      }
		      configStatus = PEN_QUAD_GEO_INPUT; return configStatus;
		    }
		  if(IR == IRD)
		    {
		      strcpy(C5, C4);
		      strcat(C5, "0");
		    }
		  else if(KEEPL == 1)
		    {
		      strcpy(C5, C4);
		      strcat(C5, "0");			      
		    }
		  else
		    {
		      strcpy(C5, C4);
		      strcat(C5, &C1);			      
		    }
		  if(NBODYS > 1)  // Check whether KB is in the list.
		    {
		      bool Eixir3 = false;  //GOTO 304
		      for(unsigned int KB0 =  0; KB0 < NBODYS-1; KB0++)
			{
			  if(strcmp(C5, ALIAB[KB0]) == 0)
			    {
			      KB = KB0+1;
			      Eixir3 = true;
			      break;
			    }
			}
		      if(!Eixir3)
			{
			  if(verbose > 0){	      
			    fprintf(IW, "%s(%s)\n", LKEYW, C4);
			    fprintf(IW, "*** Undefined body label.\n");
			    if(verbose > 1){
			      printf("pen_quadricGeo:configure:Error: %s(%s)\n", LKEYW, C4);
			      printf("pen_quadricGeo:configure:Error: Undefined body label.\n");
			    }
			  }
			  configStatus = PEN_QUAD_GEO_UNDEF_BODY_LABEL; return configStatus;
			}
		    }

		  if(verbose > 0){	      
		    fprintf(IW, "%s(%4d)\n", LKEYW, KB);
		  }
		  if(bodies[KB-1].KBOMO != 0)
		    {
		      if(verbose > 0){	      
			fprintf(IW, "*** This body is a module.\n");
			if(verbose > 1){
			  printf("pen_quadricGeo:configure:Error: This body is a module.\n");
			}
		      }
		      configStatus = PEN_QUAD_GEO_BODY_IS_MODULE; return configStatus;
		    }
		  if(bodies[KB-1].KMOTH > 0 && bodies[KB-1].KMOTH != NBODYS)
		    {
		      if(verbose > 0){	      
			fprintf(IW, "*** You are trying to assign two mothers to the last body.\n");
			if(verbose > 1){
			  printf("pen_quadricGeo:configure:Error: You are trying to assign two mothers to the last body.\n");
			}
		      }
		      configStatus = PEN_QUAD_GEO_TWO_MOTHERS_LAST; return configStatus;
		    }
		  bodies[KB-1].KMOTH = NBODYS;
		  unsigned int KDT = bodies[NBODYS-1].KDGHT[NXG-1]+1;
		  bodies[NBODYS-1].KDGHT[NXG-1] = KDT;
		  bodies[NBODYS-1].KDGHT[KDT-1] = KB;
		  for(unsigned int K2 = 0; K2 < bodies[KB-1].KBODY[NXG-1]; K2++)
		    {
		      if(bodies[bodies[KB-1].KBODY[K2]-1].KMOTH == 0){ bodies[bodies[KB-1].KBODY[K2]-1].KMOTH = NBODYS;}
		    }
		  //  ****  Assign genealogical attributes to the sisters of body KB.
		  int K2M = bodies[KB-1].KBODY[NXG-1];
		  for(int K2 = 0; K2 < K2M; K2++)
		    {
		      unsigned int KBD = bodies[KB-1].KBODY[K2];
		      if(bodies[KBD-1].KMOTH == 0){ bodies[KBD-1].KMOTH = NBODYS;}
		      int IDGHT = 0;
		      for(unsigned int K = 0; K < KDT; K++)
			{
			  if(bodies[NBODYS-1].KDGHT[K] == KBD){ IDGHT = K+1;}
			}
		      if(IDGHT == 0)
			{
			  bodies[NBODYS-1].KDGHT[NXG-1] = bodies[NBODYS-1].KDGHT[NXG-1]+1;
			  bodies[NBODYS-1].KDGHT[bodies[NBODYS-1].KDGHT[NXG-1]-1] = KB;
			}
		    }
		  //  ****  Surfaces of the sister bodies.
		  unsigned int KN1 = bodies[KB-1].getKSURF(NXG-1);
		  //unsigned int KN1 = bodies[KB-1].KSURF[NXG-1];
		  for(unsigned int KS1 = 0; KS1 < KN1; KS1++)  //GOTO 317
		    {
		      if(bodies[KB-1].getKFLAG(KS1) > 3){ continue;}
		      unsigned int KSURF1 = bodies[KB-1].getKSURF(KS1);
		      unsigned int KN2 = bodies[NBODYS-1].getKSURF(NXG-1);
		      bool Eixir3 = false;
		      for(unsigned int KS2 = 0; KS2 < KN2; KS2++)
			{
			  unsigned int KSURF2 = bodies[NBODYS-1].getKSURF(KS2);
			  if(KSURF2 == KSURF1){ Eixir3 = true; break;}
			}
		      if(Eixir3){ continue;}
		      KN2 = KN2+1;
		      if(KN2 >= NXG)
			{
			  if(verbose > 0){	      
			    fprintf(IW, "*** The number of limiting surfaces is too large.\n");
			    if(verbose > 1){
			      printf("pen_quadricGeo:configure:Error: The number of limiting surfaces is too large.\n");
			    }
			  }
			  configStatus = PEN_QUAD_GEO_MANY_LIMIT_SURFACE; return configStatus;
			}

		      bodies[NBODYS-1].setKSURF(NXG-1,KN2);
		      bodies[NBODYS-1].setSurf(KN2-1,KSURF1,4);
		      
		    }
		  //
		  //  ****  Module.
		  //
		}
	      else if(strcmp(LKEYW, LMOD) == 0)  // Check whether KB is in the list.
		{
		  Nombre_Elements_Escrits = sscanf(BLINE, "%*9c%4c", C4);

		  //Append end of string chars
		  C4[4] = '\0';
		  
		  if(Nombre_Elements_Escrits == 0)
		    {
		      if(verbose > 0){	      
			fprintf(IW, "%72s\n", BLINE);
			fprintf(IW, "*** Wrong input format.\n");
			if(verbose > 1){
			  printf("pen_quadricGeo:configure:Error: %72s\n", BLINE);
			  printf("pen_quadricGeo:configure:Error: Wrong input format.\n");
			}
		      }
		      configStatus = PEN_QUAD_GEO_INPUT; return configStatus;
		    }
			  
		  if(IR == IRD)
	  	    {
		      strcpy(C5, C4);
		      strcat(C5, "0");
		    }
		  else if(KEEPL == 1)
	  	    {
		      strcpy(C5, C4);
		      strcat(C5, "0");
		    }
		  else
	  	    {
		      strcpy(C5, C4);
		      strcat(C5, &C1);
		    }
		  if(NBODYS > 0)
		    {
		      bool Eixir3 = false;  //GOTO 305
		      for(unsigned int KB0 = 0; KB0 < NBODYS-1; KB0++)
			{
			  if(strcmp(C5, ALIAB[KB0]) == 0)
			    {
			      KB = KB0+1;
			      Eixir3 = true;
			      break;
			    }
			}
		      if(!Eixir3)
			{
			  if(verbose > 0){	      
			    fprintf(IW, "%s(%s)\n", LKEYW, C4);
			    fprintf(IW, "*** Undefined body label.\n");
			    if(verbose > 1){
			      printf("pen_quadricGeo:configure:Error: %s(%s)\n", LKEYW, C4);
			      printf("pen_quadricGeo:configure:Error: Undefined body label.\n");
			    }
			  }
			  configStatus = PEN_QUAD_GEO_UNDEF_BODY_LABEL; return configStatus;
			}
		    }

		  if(verbose > 0){	      
		    fprintf(IW, "%s(%4d)\n", LKEYW, KB);
		  }
		  if(bodies[KB-1].KBOMO != 1)
		    {
		      if(verbose > 0){	      
			fprintf(IW, "*** This module is a body.\n");
			if(verbose > 1){
			  printf("pen_quadricGeo:configure:Error: This module is a body.\n");
			}
		      }
		      configStatus = PEN_QUAD_GEO_MODULE_IS_BODY; return configStatus;
		    }
		  if(bodies[KB-1].KMOTH > 0 && bodies[KB-1].KMOTH != NBODYS)
  		    {
		      if(verbose > 0){	      
			fprintf(IW, "*** You are trying to assign two mothers to the last module.\n");
			if(verbose > 1){
			  printf("pen_quadricGeo:configure:Error: You are trying to assign two mothers to the last module.\n");
			}
		      }
		      configStatus = PEN_QUAD_GEO_TWO_MOTHERS_LAST; return configStatus;
		    }
		  bodies[KB-1].KMOTH = NBODYS;
		  int KDT = bodies[NBODYS-1].KDGHT[NXG-1]+1;
		  bodies[NBODYS-1].KDGHT[NXG-1] = KDT;
		  bodies[NBODYS-1].KDGHT[KDT-1] = KB;
		  //  ****  Limiting surfaces.
		  unsigned int KN1 = bodies[KB-1].getKSURF(NXG-1);
		  for(unsigned int KS1 = 0; KS1 < KN1; KS1++)
		    {
		      if(bodies[KB-1].getKFLAG(KS1) > 2){ continue;}
		      unsigned int KSURF1 = bodies[KB-1].getKSURF(KS1);
		      unsigned int KN2 = bodies[NBODYS-1].getKSURF(NXG-1);

		      bool Eixir3 = false; //GOTO 307
		      for(unsigned int KS2 = 0; KS2 < KN2; KS2++)
			{
			  unsigned int KSURF2 = bodies[NBODYS-1].getKSURF(KS2);
			  if(KSURF2 == KSURF1){ Eixir3 = true; break;}
			}
		      if(Eixir3){continue;}
		      KN2 = KN2+1;
		      if(KN2 >= NXG)
			{
			  if(verbose > 0){	      
			    fprintf(IW, "*** The number of limiting surfaces is too large.\n");
			    if(verbose > 1){
			      printf("pen_quadricGeo:configure:Error: The number of limiting surfaces is too large.\n");
			    }
			  }
			  configStatus = PEN_QUAD_GEO_MANY_LIMIT_SURFACE; return configStatus;
			}
		      bodies[NBODYS-1].setKSURF(NXG-1,KN2);
		      bodies[NBODYS-1].setSurf(KN2-1,KSURF1,4);

		    }
		}
	      else
		{
		  if(verbose > 0){	      
		    fprintf(IW, "%72s\n", BLINE);
		    fprintf(IW, "*** What do you mean?\n");
		    if(verbose > 1){
		      printf("pen_quadricGeo:configure:Error: %72s\n", BLINE);
		      printf("pen_quadricGeo:configure:Error: What do you mean?\n");
		    }
		  }
		  configStatus = PEN_QUAD_GEO_MEAN; return configStatus;
		}
	      Eixir2 = false;         //GOTO 301
	      continue;
	    }
	  //
	  //  ****  Transformation parameters.
	  //
	  if(Goto308)
	    {
	      if(verbose > 0){	      
		fprintf(IW, "1111111111111111111111111111111111111111111111111111111111111111\n");
	      }
	      double OMEGA = 0.0;
	      double THETA = 0.0;
	      double PHI = 0.0;
	      double XSHIFT = 0.0;
	      double YSHIFT = 0.0;
	      double ZSHIFT = 0.0;

	      bool Eixir3 = false;   //GOTO 309
	      while(!Eixir3)
		{
		  Eixir3 = true;
		  memset(BLINE,0,sizeof(BLINE));
          fscanf(IR,"%[^\r\n]%*[^\n]",BLINE);
          if(fgetc(IR)=='\r'){fgetc(IR);}
		  if(((BLINE[0] == 'C' || BLINE[0] == 'c') && (BLINE[1] == '\0' || BLINE[1] == CCR || BLINE[1] == CNL || BLINE[1] == ' ')) || BLINE[0] == '#')
		    {
		      if(verbose > 0){	      
			if(strlen(BLINE)>72){fprintf(IW, "%-.72s\n", BLINE);}else{fprintf(IW, "%-72s\n", BLINE);}
		      }
		      Eixir3 = false;
		      continue;
		    }
		  double VALUE;
		  int ICHPAR;
		  char LANGLE[9];
		  char AUXSTR[9];strcpy(AUXSTR,"        ");
		  Nombre_Elements_Escrits = sscanf(BLINE, "%8c%*c%lf%*c%4d%8c", LKEYW, &VALUE, &ICHPAR, AUXSTR);

		  //Append end of string chars
		  LKEYW[8] = '\0';
		  AUXSTR[8] = '\0';
		  
		  if(Nombre_Elements_Escrits == 0)
		    {
		      if(verbose > 0){	      
			fprintf(IW, "%72s\n", BLINE);
			fprintf(IW, "*** Wrong input format.\n");
			if(verbose > 1){
			  printf("pen_quadricGeo:configure:Error: %72s\n", BLINE);
			  printf("pen_quadricGeo:configure:Error: Wrong input format.\n");
			}
		      }
		      configStatus = PEN_QUAD_GEO_INPUT; return configStatus;
		    }
		  if(strcmp(LKEYW, LNUL) == 0){ break;} //GOTO 310
		  sprintf(LANGLE,"%-8s",AUXSTR);		  

		  
		  if(strcmp(LKEYW, LOME) == 0)
		    {
		      if(strcmp(LANGLE, LRAD) == 0)
			{
			  if(verbose > 0){	      
			    fprintf(IW, "%s(%22.15E,%4d%s\n", LKEYW, VALUE, ICHPAR, LRAD);
			  }
			  OMEGA = VALUE;
			}
		      else
			{
			  if(verbose > 0){	      
			    fprintf(IW, "%s(%22.15E,%4d%s\n", LKEYW, VALUE, ICHPAR, LDEG);
			  }
			  OMEGA = VALUE*constants::PI/180.0;
			}
		    }
		  else if(strcmp(LKEYW, LTHE) == 0)
		    {
		      if(strcmp(LANGLE, LRAD) == 0)
			{
			  if(verbose > 0){	      
			    fprintf(IW, "%s(%22.15E,%4d%s\n", LKEYW, VALUE, ICHPAR, LRAD);
			  }
			  THETA = VALUE;
			}
		      else
			{
			  if(verbose > 0){	      
			    fprintf(IW, "%s(%22.15E,%4d%s\n", LKEYW, VALUE, ICHPAR, LDEG);
			  }
			  THETA = VALUE*constants::PI/180.0;
			}
		    }
		  else if(strcmp(LKEYW, LPHI) == 0)
		    {
		      if(strcmp(LANGLE, LRAD) == 0)
			{
			  if(verbose > 0){	      
			    fprintf(IW, "%s(%22.15E,%4d%s\n", LKEYW, VALUE, ICHPAR, LRAD);
			  }
			  PHI = VALUE;
			}
		      else
			{
			  if(verbose > 0){	      
			    fprintf(IW, "%s(%22.15E,%4d%s\n", LKEYW, VALUE, ICHPAR, LDEG);
			  }
			  PHI = VALUE*constants::PI/180.0;
			}
		    }
		  else if(strcmp(LKEYW, LXSH) == 0)
		    {
		      if(verbose > 0){	      
			fprintf(IW, "%s(%22.15E,%4d%s\n", LKEYW, VALUE, ICHPAR, LOPEN);
		      }
		      XSHIFT = VALUE;
		    }
		  else if(strcmp(LKEYW, LYSH) == 0)
		    {
		      if(verbose > 0){	      
			fprintf(IW, "%s(%22.15E,%4d%s\n", LKEYW, VALUE, ICHPAR, LOPEN);
		      }
		      YSHIFT = VALUE;
		    }
		  else if(strcmp(LKEYW, LZSH) == 0)
		    {
		      if(verbose > 0){	      
			fprintf(IW, "%s(%22.15E,%4d%s\n", LKEYW, VALUE, ICHPAR, LOPEN);
		      }
		      ZSHIFT = VALUE;
		    }
		  else
		    {
		      if(verbose > 0){	      
			fprintf(IW, "%72s\n", BLINE);
			fprintf(IW, "*** What do you mean?\n");
			if(verbose > 1){
			  printf("pen_quadricGeo:configure:Error: %72s\n", BLINE);
			  printf("pen_quadricGeo:configure:Error: What do you mean?\n");
			}
		      }
		      configStatus = PEN_QUAD_GEO_MEAN; return configStatus;
		    }
		  Eixir3 = false;   //GOTO 309
		  continue;
		}
	      //
	      //  ****  Rotation and translation of the module contents (its surfaces).
	      //
	      for(unsigned int KS = 0; KS < NSURF; KS++)
		{
		  KM[KS] = 0;
		}
	      for(KB = 0; KB < NBODYS; KB++)
		{
		  unsigned int KBMOTH = KB+1;  // We need to transform all the descendants.
		  Eixir3 = false;     // GOTO 311
		  while(!Eixir3)
		    {
		      Eixir3 = true;
		      if(KBMOTH == NBODYS)
			{
			  int KNS = bodies[KB].getKSURF(NXG-1);
			  for(int KSS = 0; KSS < KNS; KSS++)
			    {
			      int KS = bodies[KB].getKSURF(KSS);
			      if(bodies[KB].getKFLAG(KSS) < 5 && KM[KS-1] == 0 && KSFF[KS-1] == 0)				
				{
				  double QXX = surfaces[KS-1].AXX;
				  double QXY = surfaces[KS-1].AXY;
				  double QXZ = surfaces[KS-1].AXZ;
				  double QYY = surfaces[KS-1].AYY;
				  double QYZ = surfaces[KS-1].AYZ;
				  double QZZ = surfaces[KS-1].AZZ;
				  double QX = surfaces[KS-1].AX;
				  double QY = surfaces[KS-1].AY;
				  double QZ = surfaces[KS-1].AZ;
				  double Q0 = surfaces[KS-1].A0;
				  ROTSHF(OMEGA,THETA,PHI,XSHIFT,YSHIFT,ZSHIFT,QXX,QXY,QXZ,QYY,QYZ,QZZ,QX,QY,QZ,Q0);
				  surfaces[KS-1].AXX = QXX;
				  surfaces[KS-1].AXY = QXY;
				  surfaces[KS-1].AXZ = QXZ;
				  surfaces[KS-1].AYY = QYY;
				  surfaces[KS-1].AYZ = QYZ;
				  surfaces[KS-1].AZZ = QZZ;
				  surfaces[KS-1].AX = QX;
				  surfaces[KS-1].AY = QY;
				  surfaces[KS-1].AZ = QZ;
				  surfaces[KS-1].A0 = Q0;
				  KM[KS-1] = 1;
				}
			    }
			}
		      else
			{
			  KBMOTH = bodies[KBMOTH-1].KMOTH;  // Moves a generation up (grandmother).
			  if(KBMOTH > 0){ Eixir3 = false; continue;}
			}
		    }
		}
		
	    }
	  if(Goto312)
	    {
	      if(verbose > 0){	      
		fprintf(IW, "%0*d\n",64,0);
	      }
	      bodies[NBODYS-1].KBOMO = 1;
	      strcpy(DEFB[NBODYS-1],DEF);
	      Eixir = false;
	      continue;
	    }
	}
      else if(strcmp(LKEYW, LCLON) == 0)      //GOTO 400
	{ 
	  //
	  //  ************  Clone a module.
	  //
	  if(BLINE[8] != '(' || BLINE[13] != ')')
	    {
	      if(verbose > 0){	      
		if(strlen(BLINE)>72){fprintf(IW, "%-.72s\n", BLINE);}else{fprintf(IW, "%-72s\n", BLINE);}
		fprintf(IW, "*** Incorrect label format.\n");
		if(verbose > 1){
		  printf("pen_quadricGeo:configure:Error: Incorrect label format.\n");
		}
	      }
	      configStatus = PEN_QUAD_GEO_LABEL_FORMAT; return configStatus;
	    }
	  if(verbose > 0){	      
	    fprintf(IW, "C \n");
	    fprintf(IW, "C >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n");
	    fprintf(IW, "C ************  Cloned module: \n");
	    fprintf(IW, "C \n");
	  }
	  ICLONE = ICLONE+1;
	  int Nombre_Elements_Escrits = sscanf(BLINE, "%*9c%4c%8c%8c%8c%8c%8c%8c%8c", C4, LARRAY[0], LARRAY[1], LARRAY[2], LARRAY[3], LARRAY[4], LARRAY[5], LARRAY[6]);

	  //Append end of string chars
	  C4[4] = '\0';
	  LARRAY[0][8] = '\0';
	  LARRAY[1][8] = '\0';
	  LARRAY[2][8] = '\0';
	  LARRAY[3][8] = '\0';
	  LARRAY[4][8] = '\0';
	  LARRAY[5][8] = '\0';
	  LARRAY[6][8] = '\0';
	  
	  if(Nombre_Elements_Escrits == 0)
	    {
	      if(verbose > 0){	      
		fprintf(IW, "%72s\n", BLINE);
		fprintf(IW, "*** Wrong input format.\n");
		if(verbose > 1){
		  printf("pen_quadricGeo:configure:Error: %72s\n", BLINE);
		  printf("pen_quadricGeo:configure:Error: Wrong input format.\n");
		}
	      }
	      configStatus = PEN_QUAD_GEO_INPUT; return configStatus;
	    }
	  if(IR == IRD)
	    {
	      strcpy(C5, C4);
	      strcat(C5, "0");
	    }
	  else if(KEEPL == 1)
	    {
	      strcpy(C5, C4);
	      strcat(C5, "0");
	    }
	  else
	    {
	      strcpy(C5, C4);
	      strcat(C5, &C1);
	    }
	  strcpy(C5C, C5);
	  if(((BLINE[0] == 'C' || BLINE[0] == 'c') && (BLINE[1] == '\0' || BLINE[1] == CCR || BLINE[1] == CNL || BLINE[1] == ' ')) || BLINE[0] == '#')
	    {
	      if(verbose > 0){	      
		if(strlen(BLINE)>72){fprintf(IW, "%-.72s\n", BLINE);}else{fprintf(IW, "%-72s\n", BLINE);}
	      }
	    }
	  else
	    {
	      if(verbose > 0){	      
		if(strlen(BLINE)>72){fprintf(IW, "# %-.72s\n", BLINE);}else{fprintf(IW, "# %-72s\n", BLINE);}
	      }
	    }
	  if(NBODYS > 0)
	    {
	      for(unsigned int KB0 = 0; KB0 < NBODYS; KB0++)
		{
		  if(strcmp(C5, ALIAB[KB0]) == 0)
		    {
		      if(verbose > 0){	      
			if(strlen(BLINE)>72){fprintf(IW, "%-.72s\n", BLINE);}else{fprintf(IW, "%-72s\n", BLINE);}
			fprintf(IW, "*** Same label for two bodies or modules.\n");
			if(verbose > 1){
			  printf("pen_quadricGeo:configure:Error: Same label for two bodies or modules.\n");
			}
		      }
		      configStatus = PEN_QUAD_GEO_SAME_LABEL_BODY; return configStatus;
		    }
		}
	    }

	  bool Eixir2 = false;          //GOTO 401
	  while(!Eixir2)
	    {
	      Eixir2 = true;
	      memset(BLINE,0,sizeof(BLINE));
          fscanf(IR,"%[^\r\n]%*[^\n]",BLINE);
          if(fgetc(IR)=='\r'){fgetc(IR);}
	      if(((BLINE[0] == 'C' || BLINE[0] == 'c') && (BLINE[1] == '\0' || BLINE[1] == CCR || BLINE[1] == CNL || BLINE[1] == ' ')) || BLINE[0] == '#')
		{
		  if(verbose > 0){	      
		    if(strlen(BLINE)>72){fprintf(IW, "%-.72s\n", BLINE);}else{fprintf(IW, "%-72s\n", BLINE);}
		  }
		  Eixir2 = false;
		  continue;
		}
	    }
	  Nombre_Elements_Escrits = sscanf(BLINE, "%8c%*c%4c", LKEYW, C4);

	  //Append end of string chars
	  LKEYW[8] = '\0';
	  C4[4] = '\0';
	  
	  if(Nombre_Elements_Escrits == 0)
	    {
	      if(verbose > 0){	      
		fprintf(IW, "%72s\n", BLINE);
		fprintf(IW, "*** Wrong input format.\n");
		if(verbose > 1){
		  printf("pen_quadricGeo:configure:Error: %72s\n", BLINE);
		  printf("pen_quadricGeo:configure:Error: Wrong input format.\n");
		}
	      }
	      configStatus = PEN_QUAD_GEO_INPUT; return configStatus;
	    }
	  if(IR == IRD)
	    {
	      strcpy(C5, C4);
	      strcat(C5, "0");
	    }
	  else if(KEEPL == 1)
	    {
	      strcpy(C5, C4);
	      strcat(C5, "0");
	    }
	  else
	    {
	      strcpy(C5, C4);
	      strcat(C5, &C1);
	    }
	  if(verbose > 0){	      
	    fprintf(IW, "# %s(%s)\n", LKEYW, C4);
	  }
	  if(strcmp(LKEYW, LMOD) != 0)
	    {
	      if(verbose > 0){	      
		fprintf(IW, "%72s\n", BLINE);
		fprintf(IW, "*** The cloned object must be a module.\n");
		if(verbose > 1){
		  printf("pen_quadricGeo:configure:Error: %72s\n", BLINE);
		  printf("pen_quadricGeo:configure:Error: The cloned object must be a module.\n");
		}
	      }
	      configStatus = PEN_QUAD_GEO_CLONED_NO_MODULE; return configStatus;
	    }
	  if(NBODYS == 0)
	    {
	      if(verbose > 0){	      
		fprintf(IW, "%72s\n", BLINE);
		fprintf(IW, "*** This module is not defined.\n");
		if(verbose > 1){
		  printf("pen_quadricGeo:configure:Error: %72s\n", BLINE);
		  printf("pen_quadricGeo:configure:Error: This module is not defined.\n");
		}
	      }
	      configStatus = PEN_QUAD_GEO_MODULE_UNDEF; return configStatus;
	    }
	  bool Eixir3 = false;  //GOTO 402
	  int KORIG;
	  for(unsigned int KB0 = 0; KB0 < NBODYS; KB0++)
	    {
	      if(strcmp(C5, ALIAB[KB0]) == 0)
		{
		  KORIG = KB0+1;
		  if(bodies[KORIG-1].KBOMO != 1)
		    {
		      if(verbose > 0){	      
			fprintf(IW, "*** The cloned object must be a module.\n");
			fprintf(IW, "*** The selected object is a body.\n");
			if(verbose > 1){
			  printf("pen_quadricGeo:configure:Error: The cloned object must be a module.\n");
			  printf("pen_quadricGeo:configure:Error: The selected object is a body.\n");			
			}
		      }
		      configStatus = PEN_QUAD_GEO_OBJECT_IS_BODY; return configStatus;
		    }
		  Eixir3 = true;
		  break;
		}
	    }
	  if(!Eixir3)
	    {
	      if(verbose > 0){	      
		fprintf(IW, "*** The label does not correspond to a module.\n");
		if(verbose > 1){
		  printf("pen_quadricGeo:configure:Error: The label does not correspond to a module.\n");
		}
	      }
	      configStatus = PEN_QUAD_GEO_LABEL_NOT_MODULE; return configStatus;	  
	    }

	  Eixir2 = false;  //GOTO 402
	  while(!Eixir2)
	    {
	      Eixir2 = true;
	      memset(BLINE,0,sizeof(BLINE));
          fscanf(IR,"%[^\r\n]%*[^\n]",BLINE);
          if(fgetc(IR)=='\r'){fgetc(IR);}
	      if(((BLINE[0] == 'C' || BLINE[0] == 'c') && (BLINE[1] == '\0' || BLINE[1] == CCR || BLINE[1] == CNL || BLINE[1] == ' ')) || BLINE[0] == '#')
		{
		  if(verbose > 0){	      
		    if(strlen(BLINE)>72){fprintf(IW, "%-.72s\n", BLINE);}else{fprintf(IW, "%-72s\n", BLINE);}
		  }
		  Eixir2 = false;
		  continue;
		}
	    }
	  sscanf(BLINE, "%8c", LKEYW);

	  //Append end of string chars
	  LKEYW[8] = '\0';
	  
	  if(strcmp(LKEYW, LONE) == 0 || strcmp(LKEYW, LNUL) == 0){}
	  else
	    {
	      if(verbose > 0){	      
		fprintf(IW, "%72s\n", BLINE);
		fprintf(IW, "*** What do you mean?\n");
		if(verbose > 1){
		  printf("pen_quadricGeo:configure:Error: %72s\n", BLINE);
		  printf("pen_quadricGeo:configure:Error: What do you mean?\n");
		}
	      }
	      configStatus = PEN_QUAD_GEO_MEAN; return configStatus;
	    }

	  //  ****  Transformation parameters.

  	  double OMEGA = 0.0;
	  double THETA = 0.0;
	  double PHI = 0.0;
	  double XSHIFT = 0.0;
	  double YSHIFT = 0.0;
	  double ZSHIFT = 0.0;
	  if(strcmp(LKEYW, LNUL) == 0){}   //GOTO 405
	  else
	    {
	      Eixir2 = false;    //GOTO 404
	      while(!Eixir2)
		{
		  Eixir2 = true;
		  memset(BLINE,0,sizeof(BLINE));
          fscanf(IR,"%[^\r\n]%*[^\n]",BLINE);
          if(fgetc(IR)=='\r'){fgetc(IR);}
		  if(((BLINE[0] == 'C' || BLINE[0] == 'c') && (BLINE[1] == '\0' || BLINE[1] == CCR || BLINE[1] == CNL || BLINE[1] == ' ')) || BLINE[0] == '#')
		    {
		      if(verbose > 0){	      
			if(strlen(BLINE)>72){fprintf(IW, "%-.72s\n", BLINE);}else{fprintf(IW, "%-72s\n", BLINE);}
		      }
		      Eixir2 = false;
		      continue;
		    }
		  double VALUE;
		  int ICHPAR;
		  char LANGLE[9];
		  char AUXSTR[9];strcpy(AUXSTR,"        ");
		  Nombre_Elements_Escrits = sscanf(BLINE, "%8c%*c%lf%*c%4d%8c", LKEYW, &VALUE, &ICHPAR, AUXSTR);

		  //Append end of string chars
		  LKEYW[8] = '\0';
		  AUXSTR[8] = '\0';
		  
		  sprintf(LANGLE,"%-8s",AUXSTR);
		  if(Nombre_Elements_Escrits == 0)
		    {
		      if(verbose > 0){	      
			fprintf(IW, "%72s\n", BLINE);
			fprintf(IW, "*** Wrong input format.\n");
			if(verbose > 1){
			  printf("pen_quadricGeo:configure:Error: %72s\n", BLINE);
			  printf("pen_quadricGeo:configure:Error: Wrong input format.\n");
			}
		      }
		      configStatus = PEN_QUAD_GEO_INPUT; return configStatus;
		    }
		  if(strcmp(LKEYW, LNUL) == 0){ break;}
		  else
  		    {

		      
		      if(strcmp(LKEYW, LOME) == 0)
			{
			  if(strcmp(LANGLE, LRAD) == 0)
			    {
			      if(verbose > 0){	      
				fprintf(IW, "# %s(%22.15E,%4d%s\n", LKEYW, VALUE, ICHPAR, LRAD);
			      }
			      OMEGA = VALUE;
			    }
			  else
			    {
			      if(verbose > 0){	      
				fprintf(IW, "# %s(%22.15E,%4d%s\n", LKEYW, VALUE, ICHPAR, LDEG);
			      }
			      OMEGA = VALUE*constants::PI/180.0;
			    }
			}
		      else if(strcmp(LKEYW, LTHE) == 0)
			{
			  if(strcmp(LANGLE, LRAD) == 0)
			    {
			      if(verbose > 0){	      
				fprintf(IW, "# %s(%22.15E,%4d%s\n", LKEYW, VALUE, ICHPAR, LRAD);
			      }
			      THETA = VALUE;
			    }
			  else
			    {
			      if(verbose > 0){	      
				fprintf(IW, "# %s(%22.15E,%4d%s\n", LKEYW, VALUE, ICHPAR, LDEG);
			      }
			      THETA = VALUE*constants::PI/180.0;
			    }
			}
		      else if(strcmp(LKEYW, LPHI) == 0)
			{
			  if(strcmp(LANGLE, LRAD) == 0)
			    {
			      if(verbose > 0){	      
				fprintf(IW, "# %s(%22.15E,%4d%s\n", LKEYW, VALUE, ICHPAR, LRAD);
			      }
			      PHI = VALUE;
			    }
			  else
			    {
			      if(verbose > 0){	      
				fprintf(IW, "# %s(%22.15E,%4d%s\n", LKEYW, VALUE, ICHPAR, LDEG);
			      }
			      PHI = VALUE*constants::PI/180.0;
			    }
			}
		      else if(strcmp(LKEYW, LXSH) == 0)
			{
			  if(verbose > 0){	      
			    fprintf(IW, "# %s(%22.15E,%4d%s\n", LKEYW, VALUE, ICHPAR, LOPEN);
			  }
			  XSHIFT = VALUE;
			}
		      else if(strcmp(LKEYW, LYSH) == 0)
			{
			  if(verbose > 0){	      
			    fprintf(IW, "# %s(%22.15E,%4d%s\n", LKEYW, VALUE, ICHPAR, LOPEN);
			  }
			  YSHIFT = VALUE;
			}
		      else if(strcmp(LKEYW, LZSH) == 0)
			{
			  if(verbose > 0){	      
			    fprintf(IW, "# %s(%22.15E,%4d%s\n", LKEYW, VALUE, ICHPAR, LOPEN);
			  }
			  ZSHIFT = VALUE;
		        }
		      else
			{
			  if(verbose > 0){	      
			    fprintf(IW, "%72s\n", BLINE);
			    fprintf(IW, "*** What do you mean?\n");
			    if(verbose > 1){
			      printf("pen_quadricGeo:configure:Error: %72s\n", BLINE);
			      printf("pen_quadricGeo:configure:Error: What do you mean?\n");
			    }
			  }
			  configStatus = PEN_QUAD_GEO_MEAN; return configStatus;
			}
		    }
		  Eixir2 = false;
		  continue;
		}
	    }     //END GOTO 404
	  
	  //  ****  Determine all the descendants of module KORIG.
  	  int ND = 1;
	  IDESC[0] = KORIG;  // Descendants of the cloned module.
	  IDONE[0] = 0;      // The descendants have not yet been identified.
  	  Eixir2 = false;    //GOTO 406
	  while(!Eixir2)
	    {
	      Eixir2 = true;
	      int NDC = ND;
	      int KDG = 0;
	      for(int I = 0; I < NDC; I++)
		{
		  if(IDONE[I] == 0)
		    {
		      KB = IDESC[I];
		      if(bodies[KB-1].KBODY[NXG-1] > 0)
			{
			  for(unsigned int J = 0; J < bodies[KB-1].KBODY[NXG-1]; J++)
			    {
			      ND = ND+1;
			      IDESC[ND-1] = bodies[KB-1].KBODY[J];  // New descendant.
			      IDONE[ND-1] = 0;
			      KDG = KDG+1;
			    }
			}
		      else if(bodies[KB-1].KDGHT[NXG-1] > 0)
			{
			  for(unsigned int J = 0; J < bodies[KB-1].KDGHT[NXG-1]; J++)
			    {
			      if(bodies[KB-1].KDGHT[J] != KB)
				{
				  ND = ND+1;
				  IDESC[ND-1] = bodies[KB-1].KDGHT[J];  // New descendant.
				  IDONE[ND-1] = 0;
				  KDG = KDG+1;
				}
			    }
			}
		      IDONE[I] = 1;  // The descendants of KB=IDESC(I) have been listed.
  		    }
		}
	      if(KDG > 0){ Eixir2 = false; continue;}    //GOTO 406
	    }

  	  int IN = 0;
	  for(unsigned int I = 0; I < NB; I++)
	    {
	      IBCL[I] = 0;  // Label of a cloned body or module.
	      IBOR[I] = 0;  // Label of the original body or module.
	    }
	  unsigned int KSD = NSURF;
	  for(unsigned int I = 0; I < NS; I++)
	    {
	      ISCL[I] = 0;  // Label of a cloned surface.
	    }
	  unsigned int KBD = NBODYS;
	  for(unsigned int KBB = 0; KBB < NBODYS; KBB++)    //GOTO 407
	    {
	      for(int ID = ND-1; ID >= 0; ID--)
		{
		  KB = IDESC[ID];
		  if(KBB+1 == KB)
		    {
		      KBD = KBD+1;
		      IBCL[KB-1] = KBD;
		      IBOR[KBD-1] = KB;
		      bodies[KBD-1].MATER = bodies[KB-1].MATER;
		      //  ****  Clone the surfaces of the original module and its descendants.
		      bodies[KBD-1].setKSURF(NXG-1,bodies[KB-1].getKSURF(NXG-1));
		      for(int KSS = bodies[KB-1].getKSURF(NXG-1)-1; KSS >= 0; KSS--)
			{
			  int KS = bodies[KB-1].getKSURF(KSS);
			  if(bodies[KB-1].getKFLAG(KSS) < 3 && ISCL[KS-1] == 0)
			    {
			      if(KSFF[KS-1] == 0)
				{
				  KSD = KSD+1;
				  if(KSD > NS)
				    {
				      if(verbose > 0){	      
					fprintf(IW, "*** The parameter NS must be increased.\n");
				      
					if(verbose > 1){
					  printf("pen_quadricGeo:configure:Error: The parameter NS must be increased.\n");
					}
				      }
				      configStatus = PEN_QUAD_GEO_NS; return configStatus;
				    }
				  std::memcpy(DEFS[KSD-1],DEFS[KS-1],100);
				  //std::strcpy(DEFS[KSD-1],DEFS[KS-1]);
				  ISCL[KS-1] = KSD;
				  strcpy(BLINE,DEFS[KSD-1]);
				  double QXX = surfaces[KS-1].AXX;
				  double QXY = surfaces[KS-1].AXY;
				  double QXZ = surfaces[KS-1].AXZ;
				  double QYY = surfaces[KS-1].AYY;
				  double QYZ = surfaces[KS-1].AYZ;
				  double QZZ = surfaces[KS-1].AZZ;
				  double QX = surfaces[KS-1].AX;
				  double QY = surfaces[KS-1].AY;
				  double QZ = surfaces[KS-1].AZ;
				  double Q0 = surfaces[KS-1].A0;
				  ROTSHF(OMEGA,THETA,PHI,XSHIFT,YSHIFT,ZSHIFT,QXX,QXY,QXZ,QYY,QYZ,QZZ,QX,QY,QZ,Q0);
				  surfaces[KSD-1].AXX = QXX;
				  surfaces[KSD-1].AXY = QXY;
				  surfaces[KSD-1].AXZ = QXZ;
				  surfaces[KSD-1].AYY = QYY;
				  surfaces[KSD-1].AYZ = QYZ;
				  surfaces[KSD-1].AZZ = QZZ;
				  surfaces[KSD-1].AX = QX;
				  surfaces[KSD-1].AY = QY;
				  surfaces[KSD-1].AZ = QZ;
				  surfaces[KSD-1].A0 = Q0;
				  KSFF[KSD-1] = 0;
				  if(verbose > 0){	      
				    fprintf(IW, "%s(%4d", LSUR, KSD);
				  
				    for(int i = 13; i < 64; i++)
				      {
					if(i==63){fprintf(IW, "%c\n", BLINE[i]);}
					else{fprintf(IW, "%c", BLINE[i]);}
				      }
				    fprintf(IW, "%s(%2d,%2d,%2d,%2d,%2d)\n", LIND, IN, IN, IN, IN, IN);
				  
				    if(fabs(surfaces[KSD-1].AXX) > 1.0E-20)
				      {
					fprintf(IW, "%s(%22.15E,%4d%s\n", LAXX, surfaces[KSD-1].AXX, IN, LOPEN);
				      }
				    if(fabs(surfaces[KSD-1].AXY) > 1.0E-20)
				      {
					fprintf(IW, "%s(%22.15E,%4d%s\n", LAXY, surfaces[KSD-1].AXY, IN, LOPEN);				      
				      }
				    if(fabs(surfaces[KSD-1].AXZ) > 1.0E-20)
				      {
					fprintf(IW, "%s(%22.15E,%4d%s\n", LAXZ, surfaces[KSD-1].AXZ, IN, LOPEN);
				      }
				    if(fabs(surfaces[KSD-1].AYY) > 1.0E-20)
				      {
					fprintf(IW, "%s(%22.15E,%4d%s\n", LAYY, surfaces[KSD-1].AYY, IN, LOPEN);
				      }
				    if(fabs(surfaces[KSD-1].AYZ) > 1.0E-20)
				      {
					fprintf(IW, "%s(%22.15E,%4d%s\n", LAYZ, surfaces[KSD-1].AYZ, IN, LOPEN);
				      }
				    if(fabs(surfaces[KSD-1].AZZ) > 1.0E-20)
				      {
					fprintf(IW, "%s(%22.15E,%4d%s\n", LAZZ, surfaces[KSD-1].AZZ, IN, LOPEN);
				      }
				    if(fabs(surfaces[KSD-1].AX) > 1.0E-20)
				      {
					fprintf(IW, "%s(%22.15E,%4d%s\n", LAX, surfaces[KSD-1].AX, IN, LOPEN);
				      }
				    if(fabs(surfaces[KSD-1].AY) > 1.0E-20)
				      {
					fprintf(IW, "%s(%22.15E,%4d%s\n", LAY, surfaces[KSD-1].AY, IN, LOPEN);
				      }
				    if(fabs(surfaces[KSD-1].AZ) > 1.0E-20)
				      {				      
					fprintf(IW, "%s(%22.15E,%4d%s\n", LAZ, surfaces[KSD-1].AZ, IN, LOPEN);
				      }
				    if(fabs(surfaces[KSD-1].A0) > 1.0E-20)
				      {	
					fprintf(IW, "%s(%22.15E,%4d%s\n", LA0, surfaces[KSD-1].A0, IN, LOPEN);
				      }
				    fprintf(IW, "0000000000000000000000000000000000000000000000000000000000000000\n");
				  }
				}
			      else
				{
				  ISCL[KS-1] = KS;
				}
			    }
			}		      
		      break;    //GOTO 407
		    }
		}
	    }    // DO 407
	  //  ****  Clone the original module and its descendants.
  	  for(KB = NBODYS+1-1; KB < KBD; KB++)
	    {
	      if(KB+1 > NB)
		{
		  if(verbose > 0){	      
		    fprintf(IW, "*** The parameter NB must be increased.\n");
		    if(verbose > 1){
		      printf("pen_quadricGeo:configure:Error: The parameter NB must be increased.\n");
		    }
		  }
		  configStatus = PEN_QUAD_GEO_NB; return configStatus;
		}
	      int KBO = IBOR[KB];
	      if(bodies[KBO-1].KMOTH > 0)
		{
		  bodies[KB].KMOTH = IBCL[bodies[KBO-1].KMOTH-1];
		}
	      else
		{
		  bodies[KB].KMOTH = 0;
		}
	      bodies[KB].KBOMO = bodies[KBO-1].KBOMO;
	      strcpy(DEFB[KB],DEFB[KBO-1]);
	      strcpy(BLINE,DEFB[KBO-1]);
	      if(bodies[KB].KBOMO == 0)
		{
		  strcpy(LKEYW,LBOD);
		  if(verbose > 0){	      
		    fprintf(IW, "%s(%4d", LBOD, KB+1);
		    for(int i = 13; i < 64; i++)
		      {
			if(i==63){fprintf(IW, "%c\n", BLINE[i]);}
			else{fprintf(IW, "%c", BLINE[i]);}
		      }
		  }
		}
	      else if(bodies[KB].KBOMO == 1)
		{
		  strcpy(LKEYW,LMOD);
		  if(verbose > 0){	      
		    fprintf(IW, "%s(%4d", LMOD, KB+1);
		    for(int i = 13; i < 64; i++)
		      {
			if(i==63){fprintf(IW, "%c\n", BLINE[i]);}
			else{fprintf(IW, "%c", BLINE[i]);}
		      }
		  }
		}
	      else
		{
		  if(verbose > 0){	      
		    fprintf(IW, "KBOMO(%4d) =%4d\n", KB+1, bodies[KB].KBOMO);
		    fprintf(IW, "*** Something wrong...\n");
		    if(verbose > 1){
		      printf("pen_quadricGeo:configure:Error: KBOMO(%4d) =%4d\n", KB+1, bodies[KB].KBOMO);
		      printf("pen_quadricGeo:configure:Error: Something wrong...\n");
		    }
		  }
		  configStatus = PEN_QUAD_GEO_WRONG; return configStatus;
		}
	      if(verbose > 0){	      
		fprintf(IW, "%s(%4d)\n", LMAT, bodies[KB].MATER);
	      }
	      int INDS;
	      for(unsigned int KS = 0; KS < bodies[KB].getKSURF(NXG-1); KS++)
		{
		  bodies[KB].setKSURF(KS,ISCL[bodies[KBO-1].getKSURF(KS)-1]);
		  bodies[KB].setKFLAG(KS,bodies[KBO-1].getKFLAG(KS));
		  
		  if(bodies[KB].getKFLAG(KS) < 3)
		    {
		      if(bodies[KB].getKFLAG(KS) == 1)
			{
			  INDS = -1;
			}
		      else
			{
			  INDS = +1;
			}
		      if(verbose > 0){	      
			fprintf(IW, "%s(%4d), SIDE POINTER=(%2d)\n", LSUR, bodies[KB].getKSURF(KS), INDS);
		      }
		    }
		}
	      if(bodies[KB].KBOMO == 0)
		{
		  bodies[KB].KBODY[NXG-1] = bodies[KBO-1].KBODY[NXG-1];
		  for(unsigned int I = 0; I < bodies[KB].KBODY[NXG-1]; I++)
		    {
		      bodies[KB].KBODY[I] = IBCL[bodies[KBO-1].KBODY[I]-1];
		      unsigned int KBB = bodies[KB].KBODY[I];
		      if(bodies[KBB-1].KBOMO == 0)
			{
			  strcpy(LKEYW,LBOD);
			  if(verbose > 0){	      
			    fprintf(IW, "%s(%4d)\n", LBOD, KBB);
			  }
			}
		      else if(KBB != KB+1)
			{
			  strcpy(LKEYW,LMOD);
			  if(verbose > 0){	      
			    fprintf(IW, "%s(%4d)\n", LMOD, KBB);
			  }
			}
		      if(KBB >= KB+1)
			{
			  if(verbose > 0){	      
			    fprintf(IW, "*** The limiting body or module is not yet defined\n");
			    if(verbose > 1){
			      printf("pen_quadricGeo:configure:Error: The limiting body or module is not yet defined\n");
			    }
			  }
			  configStatus = PEN_QUAD_GEO_LIMITING_BODY_NOT_DEF; return configStatus;
			}
		    }
		}
	      else
		{
		  bodies[KB].KDGHT[NXG-1] = bodies[KBO-1].KDGHT[NXG-1];
		  for(unsigned int I = 0; I < bodies[KB].KDGHT[NXG-1]; I++)
		    {
		      bodies[KB].KDGHT[I] = IBCL[bodies[KBO-1].KDGHT[I]-1];
		      unsigned int KBB = bodies[KB].KDGHT[I];
		      if(KBB != KB+1)
			{
			  if(bodies[KBB-1].KBOMO == 0)
			    {
			      strcpy(LKEYW,LBOD);
			      if(verbose > 0){	      
				fprintf(IW, "%s(%4d)\n", LBOD, KBB);
			      }
			    }
			  else
			    {
			      strcpy(LKEYW,LMOD);
			      if(verbose > 0){	      
				fprintf(IW, "%s(%4d)\n", LMOD, KBB);
			      }
			    }
			  if(KBB >= KB+1)
			    {
			      if(verbose > 0){	      
				fprintf(IW, "*** The limiting body or module is not yet defined\n");
				if(verbose > 1){
				  printf("pen_quadricGeo:configure:Error: The limiting body or module is not yet defined\n");
				}
			      }
			      configStatus = PEN_QUAD_GEO_LIMITING_BODY_NOT_DEF; return configStatus;
			    }
			}
		    }
		}
	      if(verbose > 0){	      
		fprintf(IW, "0000000000000000000000000000000000000000000000000000000000000000\n");
	      
		if(KB+1 == KBD)
		  {
		    fprintf(IW, "C \n");
		    fprintf(IW, "C ************  End of cloned module.\n");
		    fprintf(IW, "C <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n");
		    fprintf(IW, "C \n");
		  }
	      }
	    }
	  strcpy(ALIAB[KBD-1], C5C);
	  bodies[KBD-1].KMOTH = 0;
	  NBODYS = KBD;
	  NSURF = KSD;
	  Eixir = false;
	  continue;
	}
      else if(strcmp(LKEYW, LINC) == 0 || strcmp(LKEYW, LINA) == 0)  //GOTO 500
	{
	  //
	  //  ************  Included geometry file.
	  //
	  if(strcmp(LKEYW, LINA) == 0)
	    {
	      KEEPL = 1;  // Keep the user labels of the included file elements.
	    }
	  else
	    {
	      KEEPL = 0;
	    }
	  if(verbose > 0){	      
	    fprintf(IW, "C\n");
	    fprintf(IW, "C >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n");
	    fprintf(IW, "C ************  Included file:  \n");
	    if(KEEPL == 1)
	      {
		fprintf(IW, "C The included elements keep their user labels\n");
	      }
	    fprintf(IW, "C\n");
	    if(strlen(BLINE)>72){fprintf(IW, "%-.72s\n", BLINE);}else{fprintf(IW, "%-72s\n", BLINE);}
	  }
	  bool Eixir2 = false;    //GOTO 501
	  while(!Eixir2)
	    {
	      Eixir2 = true;
	      memset(BLINE,0,sizeof(BLINE));
          fscanf(IR,"%[^\r\n]%*[^\n]",BLINE);
          if(fgetc(IR)=='\r'){fgetc(IR);}
	      if(((BLINE[0] == 'C' || BLINE[0] == 'c') && (BLINE[1] == '\0' || BLINE[1] == CCR || BLINE[1] == CNL || BLINE[1] == ' ')) || BLINE[0] == '#')
		{
		  if(verbose > 0){	      
		    if(strlen(BLINE)>72){fprintf(IW, "%-.72s\n", BLINE);}else{fprintf(IW, "%-72s\n", BLINE);}
		  }
		  Eixir2 = false;
		  continue;
		}
	    }
	  if(BLINE[8] != '(' || BLINE[21] != ')')
	    {
	      if(verbose > 0){	      
		if(strlen(BLINE)>72){fprintf(IW, "%-.72s\n", BLINE);}else{fprintf(IW, "%-72s\n", BLINE);}
		fprintf(IW, "*** Incorrect label format.\n");
		if(verbose > 1){
		  printf("pen_quadricGeo:configure:Error: Incorrect label format.\n");
		}
	      }
	      configStatus = PEN_QUAD_GEO_LABEL_FORMAT; return configStatus;
	    }
	  int Nombre_Elements_Escrits = sscanf(BLINE, "%8c%*c%12c%8c%8c%8c%8c%8c%8c", LKEYW, GFILE, LARRAY[0], LARRAY[1], LARRAY[2], LARRAY[3], LARRAY[4], LARRAY[5]);

	  //Append end of string chars
	  LKEYW[8] = '\0';
	  GFILE[12] = '\0';
	  LARRAY[0][8] = '\0';
	  LARRAY[1][8] = '\0';
	  LARRAY[2][8] = '\0';
	  LARRAY[3][8] = '\0';
	  LARRAY[4][8] = '\0';
	  LARRAY[5][8] = '\0';
	  
	  if(Nombre_Elements_Escrits == 0)
	    {
	      if(verbose > 0){	      
		fprintf(IW, "%72s\n", BLINE);
		fprintf(IW, "*** Wrong input format.\n");
		if(verbose > 1){
		  printf("pen_quadricGeo:configure:Error: %72s\n", BLINE);
		  printf("pen_quadricGeo:configure:Error: Wrong input format.\n");
		}
	      }
	      configStatus = PEN_QUAD_GEO_INPUT; return configStatus;		  
	    }
		
	  if(strcmp(LKEYW, LFIL) != 0)
	    {
	      if(verbose > 0){	      
		fprintf(IW, "%72s\n", BLINE);
		fprintf(IW, "*** What do you mean?\n");
		if(verbose > 1){
		  printf("pen_quadricGeo:configure:Error: %72s\n", BLINE);
		  printf("pen_quadricGeo:configure:Error: What do you mean?\n");
		}
	      }
	      configStatus = PEN_QUAD_GEO_MEAN; return configStatus;
	    }
	  else
	    {
	      BLINE[strlen(BLINE)-1]='\0';
	      if(verbose > 0){	      
		if(strlen(BLINE)>72){fprintf(IW, "%-.72s\n", BLINE);}else{fprintf(IW, "%-72s\n", BLINE);}
	      }
	    }
	  Eixir2 = false;      //GOTO 502
	  while(!Eixir2)
	    {
	      Eixir2 = true;
	      memset(BLINE,0,sizeof(BLINE));
          fscanf(IR,"%[^\r\n]%*[^\n]",BLINE);
          if(fgetc(IR)=='\r'){fgetc(IR);}
	      if(((BLINE[0] == 'C' || BLINE[0] == 'c') && (BLINE[1] == '\0' || BLINE[1] == CCR || BLINE[1] == CNL || BLINE[1] == ' ')) || BLINE[0] == '#')
		{
		  if(verbose > 0){	      
		    if(strlen(BLINE)>72){fprintf(IW, "%-.72s\n", BLINE);}else{fprintf(IW, "%-72s\n", BLINE);}
		  }
		  Eixir2 = false;
		  continue;
		}
	    }
	  sscanf(BLINE, "%8c%8c%8c%8c%8c%8c%8c", LKEYW, LARRAY[0], LARRAY[1], LARRAY[2], LARRAY[3], LARRAY[4], LARRAY[5]);

	  //Append end of string chars
	  LKEYW[8] = '\0';
	  LARRAY[0][8] = '\0';
	  LARRAY[1][8] = '\0';
	  LARRAY[2][8] = '\0';
	  LARRAY[3][8] = '\0';
	  LARRAY[4][8] = '\0';
	  LARRAY[5][8] = '\0';
	  
	  if(strcmp(LKEYW, LNUL) == 0)
	    {
	      if(IR == IRI)
		{
		  if(verbose > 0){	      
		    if(strlen(BLINE)>72){fprintf(IW, "%-.72s\n", BLINE);}else{fprintf(IW, "%-72s\n", BLINE);}
		    fprintf(IW, "*** Too many include levels.\n");
		    if(verbose > 1){
		      printf("pen_quadricGeo:configure:Error: Too many include levels.\n");
		    }
		  }
		  configStatus = PEN_QUAD_GEO_LEVELS; return configStatus;
		}
	      else
		{
		  IR = IRI;
		  //  ****  The alias of elements in the included geometry are altered by
		  //  appending an extra character, to prevent having duplicated labels.
		  if(KEEPL == 0)
		    {
		      NINCL = NINCL+1;
		      if(NINCL > 35)
			{
			  if(verbose > 0){	      
			    fprintf(IW, "No. of included files =%3d\n", NINCL);
			    fprintf(IW, "*** Too many included files.\n");
			    if(verbose > 1){
			      printf("pen_quadricGeo:configure:Error: No. of included files =%3d\n", NINCL);
			      printf("pen_quadricGeo:configure:Error: Too many included files.\n");			    
			    }
			  }
			  configStatus = PEN_QUAD_GEO_LEVELS; return configStatus;
			}
		      C1 = CA[NINCL+1-1];
		    }

		  //  ****  Remove leading blanks in filename.
		  Eixir2 = false;       //GOTO 503
		  while(!Eixir2)
		    {
		      Eixir2 = true;
		      if(GFILE[0] == ' ')
			{
			  strcpy(GFILE, GFILE+1);
			  Eixir = false;
		          continue;
		        }
		    }
		  fclose(IR);
		  IR = fopen(GFILE, "r");
		  if(IR == nullptr){
		    if(verbose > 0){
		      fprintf(IW, "C Unable to open file %s. "
			      "Does the file exist?.\n", GFILE);
		      printf("pen_quadricGeo:configure:Error: Unable to open file %s. "
			     "Does the file exist?.\n", GFILE);
		    }
		    configStatus = PEN_QUAD_GEO_INPUT; return configStatus;    		    
		  }
		}
	      Eixir0 = false;
	      break;
	      //		  GOTO 1;
	    }
	}
      else if(strcmp(LKEYW, LEND) == 0)      // GOTO 600
	{
	  //	      
	  //  ************  End-line in the input file.
	  //	      
	  if(IR == IRI)
	    {
	      fclose(IR);
	      if(verbose > 0){	      
		fprintf(IW, "C \n");
		fprintf(IW, "C ************  End of included file.\n");
		fprintf(IW, "C <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\n");
		fprintf(IW, "C \n");
	      }
	      IR = IRD;
	      C1 = '0';
	      Eixir = false;      //GOTO 2
	      continue;
	    }
	  if(NBODYS == 1)
	    {
	      bodies[0].KDGHT[NXG-1] = 1;
	      bodies[0].KDGHT[0] = 1;
	    }
	      
	  for(KB = 0; KB < NBODYS; KB++)
	    {
	      strcpy(C5,ALIAB[KB]);
	      if(C5[4] == '0'){ strncpy(bodies[KB].BALIAS, C5, 4); bodies[KB].BALIAS[4]='\0';}
	    }
	  //	      
	  //  ************  Check for motherless bodies or modules.
	  //
	  int MLESS = 0;
	  int KBENC = 0;
	  for(unsigned int KBB = 0; KBB < NBODYS; KBB++)
	    {
	      if(bodies[KBB].KMOTH == 0)
		{
		  MLESS = MLESS+1;
		  KBENC = KBB+1;
		}
	    }
	  bool Eixir2 = false;    //GOTO 602
	  if(MLESS == 1)
	    {
	      if(bodies[KBENC-1].KBOMO == 1){ Eixir2 = true;} // There is a root module.
	    }

	  if(!Eixir2)
	    {
	      //		  
	      //  ************  Define the enclosure.
	      //
	      if(NBODYS > NB-1)
		{
		  if(verbose > 0){	      
		    fprintf(IW, "*** The parameter NB must be increased.\n");
		    if(verbose > 1){
		      printf("pen_quadricGeo:configure:Error: The parameter NB must be increased.\n");
		    }
		  }
		  configStatus = PEN_QUAD_GEO_NB; return configStatus;
		}
	      if(NSURF > NS-1)
		{
		  if(verbose > 0){	      
		    fprintf(IW, "*** The parameter NS must be increased.\n");
		    if(verbose > 1){
		      printf("pen_quadricGeo:configure:Error: The parameter NS must be increased.\n");
		    }
		  }
		  configStatus = PEN_QUAD_GEO_NS; return configStatus;
		}
	      //  ****  The next line serves only to avoid a warning issued by
	      //        certain compilers.
	      int NT = 0;
	      int KS = NSURF+1;
	      //  ****  The enclosure is a sphere centred at the origin and with a
	      //        radius of 1.0E7 length units.
	      surfaces[KS-1].AXX = 1.0;
	      surfaces[KS-1].AYY = 1.0;
	      surfaces[KS-1].AZZ = 1.0;
	      surfaces[KS-1].A0 = -1.0E14;
	      KB = NBODYS+1;
	      bodies[KB-1].KBOMO = 1;
	      bodies[KB-1].setKSURF(NXG-1,1);
	      bodies[KB-1].setSurf(0,KS,1);
	      bodies[KB-1].KDGHT[NXG-1] = 1;
	      bodies[KB-1].KDGHT[0] = KB;
	      unsigned int KN1,KN2,KN3,KSURF1;
	      for(unsigned int KBB = 0; KBB < NBODYS; KBB++)
		{
		  if(bodies[KBB].KMOTH == 0)
		    {
		      NT = bodies[KB-1].KDGHT[NXG-1]+1;
		      bodies[KB-1].KDGHT[NXG-1] = NT;
		      bodies[KB-1].KDGHT[NT-1] = KBB+1;
		      KN1 = bodies[KBB].getKSURF(NXG-1);
		      for(unsigned int KS1 = 0; KS1 < KN1; KS1++) //GOTO 601
			{
			  if(bodies[KBB].getKFLAG(KS1) > 3){ continue;}
			  KSURF1 = bodies[KBB].getKSURF(KS1);
			  KN2 = bodies[KB-1].getKSURF(NXG-1);
			  bool Eixir3 = false;
			  for(unsigned int KS2 = 0; KS2 < KN2; KS2++)
			    {
			      if(bodies[KB-1].getKSURF(KS2) == KSURF1){ Eixir3 = true; break;}			      
			    }
			  if(Eixir3){continue;}
			  KN2 = KN2+1;
			  if(KN2 >= NXG)
			    {
			      if(verbose > 0){	      
				fprintf(IW, "*** The parameter NXG is too small.\n");
				if(verbose > 1){
				  printf("pen_quadricGeo:configure:Error: The parameter NXG is too small.\n");
				}
			      }
			      configStatus = PEN_QUAD_GEO_NXG; return configStatus;
			    }
			  bodies[KB-1].setKSURF(NXG-1,KN2);
			  bodies[KB-1].setSurf(KN2-1,KSURF1,4);
			}
		      bodies[KBB].KMOTH = KB;
		      KN3 = KN1+1;

		      bodies[KBB].setKSURF(NXG-1,KN3);
		      bodies[KBB].setSurf(KN3-1,KS,1);
  		    }
		}
	      NSURF = KS;
	      NBODYS = KB;
	      //  ****  Sort daughters in increasing order.
	      if(NT > 1)
		{
		  unsigned int NTm1 = NT-1;
		  for(unsigned int KI = 0; KI < NTm1; KI++)
		    {
		      unsigned int KBMIN = bodies[NBODYS-1].KDGHT[KI];
		      unsigned int KMIN = KI+1;
		      for(int KJ = KI+1; KJ < NT; KJ++)
			{
			  if(bodies[NBODYS-1].KDGHT[KJ] < KBMIN)
			    {
			      KBMIN = bodies[NBODYS-1].KDGHT[KJ];
			      KMIN = KJ+1;
			    }
			}
		      if(KMIN != KI+1)
			{
			  int KSAVE = bodies[NBODYS-1].KDGHT[KI];
			  bodies[NBODYS-1].KDGHT[KI] = bodies[NBODYS-1].KDGHT[KMIN-1];
			  bodies[NBODYS-1].KDGHT[KMIN-1] = KSAVE;
			}
		    }
		}
		  
	    }     // 602 CONTNUE

	  if(verbose > 0){	      
	    fprintf(IW, "%s 0000000000000000000000000000000000000000000000000000000\n", LKEYW);
	    //
	    //  ****  Duplicated surfaces (within round-off tolerance) are removed.
	    //
	    fprintf(IW, "\n\n\n*****************************************\n");
	    fprintf(IW, "****     PENGEOM (version 2014)      ****\n");
	    fprintf(IW, "****  Constructive Quadric Geometry  ****\n");
	    fprintf(IW, "*****************************************\n\n");

	    if(NSFF > 0)
	      {
		fprintf(IW, "\nWARNING: The system contains fixed (starred) surfaces, which are\n         not affected by translations and rotations. Hence, any\n         translation or rotation that does not leave these\n         surfaces invariant will distort the system.\n\n");
	      }
	  }
	  if(NSURF < 2){}  //GOTO 704
	  else
	    {
	      int IWRITE = 0;
	      double TOL = 1.0E-14;
	      for(unsigned int KS = 0; KS < NSURF; KS++)
		{
		  KM[KS] = 0;
		}
	      for(unsigned int KS = 0; KS < NSURF-1; KS++)   //GOTO 703
		{
		  if(KM[KS] != 0){ continue;}

		  double F = (surfaces[KS].AXX > surfaces[KS].AXY ? surfaces[KS].AXX : surfaces[KS].AXY);
		  if(F < surfaces[KS].AXZ){ F = surfaces[KS].AXZ;}
		  if(F < surfaces[KS].AYY){ F = surfaces[KS].AYY;}
		  if(F < surfaces[KS].AYZ){ F = surfaces[KS].AYZ;}
		  if(F < surfaces[KS].AZZ){ F = surfaces[KS].AZZ;}
		  if(F < surfaces[KS].AX){ F = surfaces[KS].AX;}
		  if(F < surfaces[KS].AY){ F = surfaces[KS].AY;}
		  if(F < surfaces[KS].AZ){ F = surfaces[KS].AZ;}
		  if(F < surfaces[KS].A0){ F = surfaces[KS].A0;}

		  double FM = (surfaces[KS].AXX < surfaces[KS].AXY ? surfaces[KS].AXX : surfaces[KS].AXY);
		  if(FM > surfaces[KS].AXZ){ FM = surfaces[KS].AXZ;}
		  if(FM > surfaces[KS].AYY){ FM = surfaces[KS].AYY;}
		  if(FM > surfaces[KS].AYZ){ FM = surfaces[KS].AYZ;}
		  if(FM > surfaces[KS].AZZ){ FM = surfaces[KS].AZZ;}
		  if(FM > surfaces[KS].AX){ FM = surfaces[KS].AX;}
		  if(FM > surfaces[KS].AY){ FM = surfaces[KS].AY;}
		  if(FM > surfaces[KS].AZ){ FM = surfaces[KS].AZ;}
		  if(FM > surfaces[KS].A0){ FM = surfaces[KS].A0;}
			
		  if(fabs(FM) > fabs(F)){ F = FM;}
		  for(unsigned int KST = KS+1; KST < NSURF; KST++)  //GOTO 702
		    {
		      if(KM[KST] != 0){ continue;}

		      double FP = (surfaces[KST].AXX > surfaces[KST].AXY ? surfaces[KST].AXX : surfaces[KST].AXY);
		      if(FP < surfaces[KST].AXZ){ FP = surfaces[KST].AXZ;}
		      if(FP < surfaces[KST].AYY){ FP = surfaces[KST].AYY;}
		      if(FP < surfaces[KST].AYZ){ FP = surfaces[KST].AYZ;}
		      if(FP < surfaces[KST].AZZ){ FP = surfaces[KST].AZZ;}
		      if(FP < surfaces[KST].AX){ FP = surfaces[KST].AX;}
		      if(FP < surfaces[KST].AY){ FP = surfaces[KST].AY;}
		      if(FP < surfaces[KST].AZ){ FP = surfaces[KST].AZ;}
		      if(FP < surfaces[KST].A0){ FP = surfaces[KST].A0;}
			  
		      FM = (surfaces[KST].AXX < surfaces[KST].AXY ? surfaces[KST].AXX : surfaces[KST].AXY);
		      if(FM > surfaces[KST].AXZ){ FM = surfaces[KST].AXZ;}
		      if(FM > surfaces[KST].AYY){ FM = surfaces[KST].AYY;}
		      if(FM > surfaces[KST].AYZ){ FM = surfaces[KST].AYZ;}
		      if(FM > surfaces[KST].AZZ){ FM = surfaces[KST].AZZ;}
		      if(FM > surfaces[KST].AX){ FM = surfaces[KST].AX;}
		      if(FM > surfaces[KST].AY){ FM = surfaces[KST].AY;}
		      if(FM > surfaces[KST].AZ){ FM = surfaces[KST].AZ;}
		      if(FM > surfaces[KST].A0){ FM = surfaces[KST].A0;}
			  
		      if(fabs(FM) > fabs(FP)){ FP = FM;}
		      double FFP = F/FP;
		      double RFFP = 1.0/FFP;
			  
		      double TST = 0.0;
		      if(fabs(surfaces[KS].AXX) > 1.0E-16)
			{
			  double AuxDouble = fabs((surfaces[KS].AXX-surfaces[KST].AXX*FFP)/surfaces[KS].AXX);
			  if(TST < AuxDouble){ TST = AuxDouble;}
			}
		      else if(fabs(surfaces[KST].AXX) > 1.0E-16)
			{
			  double AuxDouble = fabs((surfaces[KS].AXX*RFFP-surfaces[KST].AXX)/surfaces[KST].AXX);
			  if(TST < AuxDouble){ TST = AuxDouble;}
			}
		      if(fabs(surfaces[KS].AXY) > 1.0E-16)
			{
			  double AuxDouble = fabs((surfaces[KS].AXY-surfaces[KST].AXY*FFP)/surfaces[KS].AXY);
			  if(TST < AuxDouble){ TST = AuxDouble;}
			}
		      else if(fabs(surfaces[KST].AXY) > 1.0E-16)
			{
			  double AuxDouble = fabs((surfaces[KS].AXY*RFFP-surfaces[KST].AXY)/surfaces[KST].AXY);
			  if(TST < AuxDouble){ TST = AuxDouble;}
		        }
		      if(fabs(surfaces[KS].AXZ) > 1.0E-16)
			{
			  double AuxDouble = fabs((surfaces[KS].AXZ-surfaces[KST].AXZ*FFP)/surfaces[KS].AXZ);
			  if(TST < AuxDouble){ TST = AuxDouble;}
			}
		      else if(fabs(surfaces[KST].AXZ) > 1.0E-16)
			{
			  double AuxDouble = fabs((surfaces[KS].AXZ*RFFP-surfaces[KST].AXZ)/surfaces[KST].AXZ);
			  if(TST < AuxDouble){ TST = AuxDouble;}
		        }
		      if(fabs(surfaces[KS].AYY) > 1.0E-16)
			{
			  double AuxDouble = fabs((surfaces[KS].AYY-surfaces[KST].AYY*FFP)/surfaces[KS].AYY);
			  if(TST < AuxDouble){ TST = AuxDouble;}
		        }
		      else if(fabs(surfaces[KST].AYY) > 1.0E-16)
			{
			  double AuxDouble = fabs((surfaces[KS].AYY*RFFP-surfaces[KST].AYY)/surfaces[KST].AYY);
			  if(TST < AuxDouble){ TST = AuxDouble;}
			}
		      if(fabs(surfaces[KS].AYZ) > 1.0E-16)
			{
			  double AuxDouble = fabs((surfaces[KS].AYZ-surfaces[KST].AYZ*FFP)/surfaces[KS].AYZ);
			  if(TST < AuxDouble){ TST = AuxDouble;}
			}
		      else if(fabs(surfaces[KST].AYZ) > 1.0E-16)
			{
			  double AuxDouble = fabs((surfaces[KS].AYZ*RFFP-surfaces[KST].AYZ)/surfaces[KST].AYZ);
			  if(TST < AuxDouble){ TST = AuxDouble;}
		        }
		      if(fabs(surfaces[KS].AZZ) > 1.0E-16)
			{
			  double AuxDouble = fabs((surfaces[KS].AZZ-surfaces[KST].AZZ*FFP)/surfaces[KS].AZZ);
			  if(TST < AuxDouble){ TST = AuxDouble;}			      
			}
		      else if(fabs(surfaces[KST].AZZ) > 1.0E-16)
			{
			  double AuxDouble = fabs((surfaces[KS].AZZ*RFFP-surfaces[KST].AZZ)/surfaces[KST].AZZ);
			  if(TST < AuxDouble){ TST = AuxDouble;}			      
			}
		      if(fabs(surfaces[KS].AX) > 1.0E-16)
			{
			  double AuxDouble = fabs((surfaces[KS].AX-surfaces[KST].AX*FFP)/surfaces[KS].AX);
			  if(TST < AuxDouble){ TST = AuxDouble;}
			}
		      else if(fabs(surfaces[KST].AX) > 1.0E-16)
			{
			  double AuxDouble = fabs((surfaces[KS].AX*RFFP-surfaces[KST].AX)/surfaces[KST].AX);
			  if(TST < AuxDouble){ TST = AuxDouble;}
			}
		      if(fabs(surfaces[KS].AY) > 1.0E-16)
			{
			  double AuxDouble = fabs((surfaces[KS].AY-surfaces[KST].AY*FFP)/surfaces[KS].AY);
			  if(TST < AuxDouble){ TST = AuxDouble;}
			}
		      else if(fabs(surfaces[KST].AY) > 1.0E-16)
			{
			  double AuxDouble = fabs((surfaces[KS].AY*RFFP-surfaces[KST].AY)/surfaces[KST].AY);
			  if(TST < AuxDouble){ TST = AuxDouble;}
			}
		      if(fabs(surfaces[KS].AZ) > 1.0E-16)
			{
			  double AuxDouble = fabs((surfaces[KS].AZ-surfaces[KST].AZ*FFP)/surfaces[KS].AZ);
			  if(TST < AuxDouble){ TST = AuxDouble;}
		        }
		      else if(fabs(surfaces[KST].AZ) > 1.0E-16)
			{
			  double AuxDouble = fabs((surfaces[KS].AZ*RFFP-surfaces[KST].AZ)/surfaces[KST].AZ);
			  if(TST < AuxDouble){ TST = AuxDouble;}			      
			}
		      if(fabs(surfaces[KS].A0) > 1.0E-16)
			{
			  double AuxDouble = fabs((surfaces[KS].A0-surfaces[KST].A0*FFP)/surfaces[KS].A0);
			  if(TST < AuxDouble){ TST = AuxDouble;}
		        }
		      else if(fabs(surfaces[KST].A0) > 1.0E-16)
			{
			  double AuxDouble = fabs((surfaces[KS].A0*RFFP-surfaces[KST].A0)/surfaces[KST].A0);
			  if(TST < AuxDouble){ TST = AuxDouble;}			      
			}
			  
		      if(TST < TOL)
			{
			  if(IWRITE == 0)
			    {
			      if(verbose > 0){	      
				fprintf(IW, "\n************  Removal of duplicated user-defined surfaces.\n\n");
			      }
			      IWRITE = 1;
			    }
			  if(verbose > 0){	      
			    fprintf(IW, "SURFACE (%4d) is replaced by SURFACE (%4d)\n", KST+1, KS+1);
			  }

			  //  ****  Check whether the two surface functions have the same global
			  //        sign (ISPF=0) or not (ISPF=1).
			  int ISPF = 0;
			  double TF = TOL*F;
			  if(fabs(surfaces[KS].AXX) > TF && surfaces[KS].AXX*surfaces[KST].AXX < 0.0){ ISPF = 1;}
			  if(fabs(surfaces[KS].AXY) > TF && surfaces[KS].AXY*surfaces[KST].AXY < 0.0){ ISPF = 1;}
			  if(fabs(surfaces[KS].AXZ) > TF && surfaces[KS].AXZ*surfaces[KST].AXZ < 0.0){ ISPF = 1;}
			  if(fabs(surfaces[KS].AYY) > TF && surfaces[KS].AYY*surfaces[KST].AYY < 0.0){ ISPF = 1;}
			  if(fabs(surfaces[KS].AYZ) > TF && surfaces[KS].AYZ*surfaces[KST].AYZ < 0.0){ ISPF = 1;}
			  if(fabs(surfaces[KS].AZZ) > TF && surfaces[KS].AZZ*surfaces[KST].AZZ < 0.0){ ISPF = 1;}
			  if(fabs( surfaces[KS].AX) > TF &&  surfaces[KS].AX*surfaces[KST].AX  < 0.0){ ISPF = 1;}
			  if(fabs( surfaces[KS].AY) > TF &&  surfaces[KS].AY*surfaces[KST].AY  < 0.0){ ISPF = 1;}
			  if(fabs( surfaces[KS].AZ) > TF &&  surfaces[KS].AZ*surfaces[KST].AZ  < 0.0){ ISPF = 1;}
			  if(fabs( surfaces[KS].A0) > TF &&  surfaces[KS].A0*surfaces[KST].A0  < 0.0){ ISPF = 1;}
			      
			  if(TST > 1.0E-15)
			    {
			      if(verbose > 0){	      
				fprintf(IW, "F,FP,F/FP,TST =%15.7E%15.7E%15.7E%15.7E\n", F, FP, FFP, TST);
				fprintf(IW, "AXX,AXXP,DIFF =%15.7E%15.7E%15.7E\n", surfaces[KS].AXX, surfaces[KST].AXX, surfaces[KS].AXX-surfaces[KST].AXX*FFP);
				fprintf(IW, "AXY,AXYP,DIFF =%15.7E%15.7E%15.7E\n", surfaces[KS].AXY, surfaces[KST].AXY, surfaces[KS].AXY-surfaces[KST].AXY*FFP);
				fprintf(IW, "AXZ,AXZP,DIFF =%15.7E%15.7E%15.7E\n", surfaces[KS].AXZ, surfaces[KST].AXZ, surfaces[KS].AXZ-surfaces[KST].AXZ*FFP);
				fprintf(IW, "AYY,AYYP,DIFF =%15.7E%15.7E%15.7E\n", surfaces[KS].AYY, surfaces[KST].AYY, surfaces[KS].AYY-surfaces[KST].AYY*FFP);
				fprintf(IW, "AYZ,AYZP,DIFF =%15.7E%15.7E%15.7E\n", surfaces[KS].AYZ, surfaces[KST].AYZ, surfaces[KS].AYZ-surfaces[KST].AYZ*FFP);
				fprintf(IW, "AZZ,AZZP,DIFF =%15.7E%15.7E%15.7E\n", surfaces[KS].AZZ, surfaces[KST].AZZ, surfaces[KS].AZZ-surfaces[KST].AZZ*FFP);
				fprintf(IW, "AX,AXP,DIFF   =%15.7E%15.7E%15.7E\n", surfaces[KS].AX,  surfaces[KST].AX,  surfaces[KS].AX- surfaces[KST].AX*FFP);
				fprintf(IW, "AY,AYP,DIFF   =%15.7E%15.7E%15.7E\n", surfaces[KS].AY,  surfaces[KST].AY,  surfaces[KS].AY- surfaces[KST].AY*FFP);
				fprintf(IW, "AZ,AZP,DIFF   =%15.7E%15.7E%15.7E\n", surfaces[KS].AZ,  surfaces[KST].AZ,  surfaces[KS].AZ- surfaces[KST].AZ*FFP); 	  
				fprintf(IW, "A0,A0P,DIFF   =%15.7E%15.7E%15.7E\n", surfaces[KS].A0,  surfaces[KST].A0,  surfaces[KS].A0- surfaces[KST].A0*FFP);
				fprintf(IW, " \n");
			      }
			    }

			  surfaces[KST] = surfaces[KS];

			  KM[KST] = KS+1;

			  int KBI = 1;
			  Eixir2 = false;    //GOTO 701
			  while(!Eixir2)
			    {
			      Eixir2 = true;
			      for(KB = KBI-1; KB < NBODYS; KB++)
				{
				  int KSL = 0;
				  int KSLP = 0;
				  int NSB = bodies[KB].getKSURF(NXG-1);
				  for(int K = 0; K < NSB; K++)
				    {
				      if(bodies[KB].getKSURF(K) == KS+1){ KSL = K+1;}
				      if(bodies[KB].getKSURF(K) == KST+1){ KSLP = K+1;}
				    }
				  if(KSLP > 0)
				    {
				      bodies[KB].setKSURF(KSLP-1,KS+1);
				      //  ****  If the implicit equations of surfaces KS and KST differ by a
				      //        global sign, the side pointer of KST must be reversed.
				      if(ISPF == 1)
					{
					  if(bodies[KB].getKFLAG(KSLP-1) == 1)
					    {
					      bodies[KB].setKFLAG(KSLP-1,2);
					    }
					  else if(bodies[KB].getKFLAG(KSLP-1) == 2)
					    {
					      bodies[KB].setKFLAG(KSLP-1,1);
					    }
					}
				      if(KSL > 0)
					{
					  int KFL = bodies[KB].getKFLAG(KSL-1);
					  int KFLP = bodies[KB].getKFLAG(KSLP-1);
					  if((KFL < KFLP ? KFL : KFLP) < 3)
					    {
					      bodies[KB].setKFLAG(KSL-1,(KFL < KFLP ? KFL : KFLP));

					      if(KFL != KFLP && (KFL > KFLP ? KFL : KFLP) < 3)
						{
						  if(verbose > 0){	      
						    if(bodies[KB].KBOMO == 0)
						      {
							fprintf(IW, "*** ERROR: BODY(%4d) is limited by two equivalent surfaces. Probably, \n           this body cannot be resolved because it is small and located\n           far from the origin.\n", KB+1);
							if(verbose > 1){
							  printf("pen_quadricGeo:configure:Error: BODY(%4d) is limited by two equivalent surfaces. Probably, \n           this body cannot be resolved because it is small and located\n           far from the origin.\n", KB+1);
							}
						      }
						    else
						      {
							fprintf(IW, "*** ERROR: MODULE(%4d) is limited by two equivalent surfaces. Probably, \n           this module cannot be resolved because it is small and located\n           far from the origin.\n", KB+1);
							if(verbose > 1){
							  printf("pen_quadricGeo:configure:Error: MODULE(%4d) is limited by two equivalent surfaces. Probably, \n           this module cannot be resolved because it is small and located\n           far from the origin.\n", KB+1);
							}
						      }
						  }
						  configStatus = PEN_QUAD_GEO_UNRESOLVED_BODY; return configStatus;
						}
					    }
					  else if(bodies[KB].KBOMO == 0)
					    {
					      bodies[KB].setKFLAG(KSL-1,3);
					    }
					  else if(bodies[KB].KBOMO == 1)
					    {
					      bodies[KB].setKFLAG(KSL-1,4);
					    }
					  else
					    {
					      bodies[KB].setKFLAG(KSL-1,5);
					    }
					  for(int K = KSLP-1; K < NSB-1; K++)
					    {                  
					      bodies[KB].setKSURF(K,bodies[KB].getKSURF(K+1));
					      bodies[KB].setKFLAG(K,bodies[KB].getKFLAG(K+1));
					    }

					  bodies[KB].setSurf(NSB-1,0,5);
					  bodies[KB].setKSURF(NXG-1,NSB-1);
					  
					  if(KB+1 < NBODYS)
					    {
					      KBI = KB+1;
					      Eixir2 = false;
					      break;           // GOTO 701
					    }
					}  //IF KSL
				    }  //IF KSLP
				}  //FOR KB	
			    }  //WHILE EIXIR2 701 CONTINUE
			}  //IF TST
	  	    }  //FOR KST
		}  //FOR KS
	    }  //ELSE IF NSURF<2 GOTO 704

	  //Print genealogical tree
	  if(verbose > 0){	      
	    fprintf(IW, "\n\n************  Genealogical tree. \n\n");
	  
	    for(KB = 0; KB < NBODYS; KB++)
	      {
		if(bodies[KB].KBOMO == 0)
		  {
		    fprintf(IW, "\n*** BODY   =%5d,  KMOTH =%5d,  MAT =%3d\n", KB+1, bodies[KB].KMOTH, bodies[KB].MATER);
		    if(bodies[KB].KBODY[NXG-1] > 0)
		      {
			fprintf(IW, "KBODY =");
			for(unsigned int K2 = 0; K2 < bodies[KB].KBODY[NXG-1]; K2++)
			  {
			    if(K2<15)
			      {
				fprintf(IW, "%5d", bodies[KB].KBODY[K2]);
			      }
			    else
			      {
				if(K2%15 == 0){fprintf(IW, "\n       %5d", bodies[KB].KBODY[K2]);}
				else{fprintf(IW, "%5d", bodies[KB].KBODY[K2]);}		    
			      }
			  }
			if(bodies[KB].KBODY[NXG-1]%15 == 0){fprintf(IW, "\n");}
			fprintf(IW, "\n");
		      }
		  }
		else if(bodies[KB].KBOMO == 1)
		  {
		    fprintf(IW, "\n*** MODULE =%5d,  KMOTH =%5d,  MAT =%3d\n", KB+1, bodies[KB].KMOTH, bodies[KB].MATER);
		    fprintf(IW, "KDGHT =");
		    for(unsigned int K2 = 0; K2 < bodies[KB].KDGHT[NXG-1]; K2++)
		      {
			if(K2<15)
			  {
			    fprintf(IW, "%5d", bodies[KB].KDGHT[K2]);
			  }
			else
			  {
			    if(K2%15 == 0){fprintf(IW, "\n       %5d", bodies[KB].KDGHT[K2]);}
			    else{fprintf(IW, "%5d", bodies[KB].KDGHT[K2]);}		    
			  }
		      }
		    if(bodies[KB].KDGHT[NXG-1]%15 == 0){fprintf(IW, "\n");}
		    fprintf(IW, "\n");
		  }
		else
		  {
		    fprintf(IW, "\n\n*** ERROR: the label %5d does not correspond to a body.\n", KB+1);
		    if(verbose > 1){
		      printf("pen_quadricGeo:configure:Error: the label %5d does not correspond to a body.\n", KB+1);
		    }
		    configStatus = PEN_QUAD_GEO_INCONSISTENT_BODY_LAB; return configStatus;
		  }

		fprintf(IW, "KSURF =");
		for(unsigned int KS = 0; KS < bodies[KB].getKSURF(NXG-1); KS++)
		  {
		    if(KS<15)
		      {
			fprintf(IW, "%5d", bodies[KB].getKSURF(KS));
		      }
		    else
		      {
			if(KS%15 == 0){fprintf(IW, "\n       %5d", bodies[KB].getKSURF(KS));}	
			else{fprintf(IW, "%5d", bodies[KB].getKSURF(KS));}		    
		      }
		  }
		if(bodies[KB].getKSURF(NXG-1)%15 == 0){fprintf(IW, "\n");}
		fprintf(IW, "\n");
		fprintf(IW, "KFLAG =");
		for(unsigned int KS = 0; KS < bodies[KB].getKSURF(NXG-1); KS++)
		  {
		    if(KS<15)
		      {
			fprintf(IW, "%5d", bodies[KB].getKFLAG(KS));
		      }
		    else
		      {
			if(KS%15 == 0){fprintf(IW, "\n       %5d", bodies[KB].getKFLAG(KS));}
			else{fprintf(IW, "%5d", bodies[KB].getKFLAG(KS));}		    
		      }
		  }
		if(bodies[KB].getKSURF(NXG-1)%15 == 0){fprintf(IW, "\n");}
		fprintf(IW, "\n");
	      }
	    if(MLESS == 1)
	      {
		if(bodies[KBENC-1].KBOMO == 1)
		  {
		    fprintf(IW, "\nThe module%5d is the enclosure.\n", KBENC);
		  }
	      }
	  }
	  //
	  //  ****  Surface consistency test (F. Tola).
	  //
  	  for(KB = 0; KB < NBODYS-1; KB++)
	    {
	      int KB1 = bodies[KB].KMOTH;
	      for(unsigned int I = 0; I < bodies[KB].getKSURF(NXG-1); I++)
		{
		  unsigned int KS = bodies[KB].getKSURF(I);
		  int KF = bodies[KB].getKFLAG(I);

		  for(unsigned int J = 0; J < bodies[KB1-1].getKSURF(NXG-1); J++)
	  	    {
		      if(bodies[KB1-1].getKSURF(J) == KS)
			{
			  // ---- Surface has daughter and mother at opposite sides.
			  
			  if((KF == 1 && bodies[KB1-1].getKFLAG(J) == 2) || (KF == 2 && bodies[KB1-1].getKFLAG(J) == 1))
			    {
			      if(verbose > 0){	      
				if(bodies[KB].KBOMO == 0)
				  {
				    fprintf(IW, "\n\n*** ERROR: the SURFACE (%4d), which limits BODY (%4d) and MODULE (%4d)\n           has inconsistent side pointers.\n", KS, KB+1, KB1);
				    if(verbose > 1){
				      printf("pen_quadricGeo:configure:Error: the SURFACE (%4d), which limits BODY (%4d) and MODULE (%4d)\n           has inconsistent side pointers.\n", KS, KB+1, KB1);
				    }
				  }
				else
				  {
				    fprintf(IW, "\n\n*** ERROR: the SURFACE (%4d), which limits MODULE (%4d) and MODULE (%4d)\n           has inconsistent side pointers.\n", KS, KB+1, KB1);
				    if(verbose > 1){
				      printf("pen_quadricGeo:configure:Error: the SURFACE (%4d), which limits MODULE (%4d) and MODULE (%4d)\n           has inconsistent side pointers.\n", KS, KB+1, KB1);
				    }
				  }
			      }
			      configStatus = PEN_QUAD_GEO_INCONSISTENT_SIDE; return configStatus;
			    }
			  break;
			}
		    }
		}
	    }
	  //
	  //  ****  Easiness test.
	  //
  	  unsigned int NBU = 0;
	  unsigned int NSU = 0;
	  for(KB = 0; KB < NBODYS; KB++)
	    {
	      if(bodies[KB].KBOMO == 0)
		{
		  if(NBU < bodies[KB].KBODY[NXG-1]){ NBU = bodies[KB].KBODY[NXG-1];}
		}
	      else if(bodies[KB].KBOMO == 1)
		{
		  if(NBU < bodies[KB].KDGHT[NXG-1]){ NBU = bodies[KB].KDGHT[NXG-1];}
		}
	      unsigned int NSE = 0;
	      for(unsigned int K = 0; K < bodies[KB].getKSURF(NXG-1); K++)
		{
		  if(bodies[KB].getKFLAG(K) < 5){ NSE = NSE+1;}
		}
	      if(NSU < NSE){ NSU = NSE;}
	    }
	  if(verbose > 0){	      
	    fprintf(IW, "\n\n************  Adequacy of the geometry definition.\n");
	    fprintf(IW, "\nThe largest number of bodies in a module or\n");
	    fprintf(IW, "     bodies limiting a single body is ............ %4d\n", NBU);
	    fprintf(IW, "\nThe largest number of limiting surfaces for\n");
	    fprintf(IW, "     a single body or module is .................. %4d\n\n", NSU);

	    if(NBODYS < 15 && NSURF < 15)
	      {
		fprintf(IW, "\nThe simulation of this geometry will be relatively fast,\n");
		fprintf(IW, "     no further optimization seems to be required.\n");
	      }
	    else if(NBU < 10 && NSU < 10)
	      {
		fprintf(IW, "\nThe simulation of this geometry will be relatively fast,\n");
		fprintf(IW, "      no further optimization seems to be required.\n");
	      }
	    else if(NBU < 15 && NSU < 20)
	      {
		fprintf(IW, "\nThe simulation of this geometry is expected to be slow,\n");
		fprintf(IW, "     try to split complex bodies into several modules.\n");
	      }
	    else if(NBU < 25 && NSU < 30)
	      {
		fprintf(IW, "\nThe simulation of this geometry will be very slow, you should\n");
		fprintf(IW, "     try to optimize the structure of the tree of modules.\n");
	      }
	    else
	      {
		fprintf(IW, "\nSimulating this geometry will be extremely slow.\n");
	      }
	    fprintf(IW, "\n************  The end.\n");
	    if(verbose > 1){
	      printf("pen_quadricGeo:configure:Info: Geometry loaded\n");
	    }
	  }

	  configStatus = PEN_QUAD_GEO_SUCCESS;
  	  return configStatus;	      
  	}
      else
	{
	  if(verbose > 0){	      
	    fprintf(IW, "%72s\n", BLINE);
	    fprintf(IW, "*** What do you mean?\n");
	    if(verbose > 1){
	      printf("pen_quadricGeo:configure:Error: %72s\n", BLINE);
	      printf("pen_quadricGeo:configure:Error: What do you mean?\n");
	    }
	  }
	  configStatus = PEN_QUAD_GEO_MEAN; return configStatus;
	}
    }

  //Just in case
  configStatus = PEN_QUAD_GEO_UNKNOWN_ERROR;
  return configStatus;
}

//  *********************************************************************
//                       SUBROUTINE LOCATE
//  *********************************************************************
void pen_quadricGeo::locate(pen_particleState& state) const
{
  //     This subroutine determines the body that contains the point with
  //  coordinates (X,Y,Z). The effects of numerical round-off errors are
  //  avoided by considering fuzzy surfaces, which swell or shrink slightly
  //  when the particle crosses them.
  //
  //  Input values (module TRACK_mod):
  //     X, Y, Z ... coordinates of the particle,
  //     U, V, W ... direction of movement.
  //
  //  Output values (module TRACK_mod):
  //     IBODY ..... body where the particle moves,
  //     MAT  ...... material in IBODY,
  //                    MAT=0, indicates a void region.


  //Declare array for surfaces position from particle  view (forward 1 backward 2)
  unsigned KSP[NS];

  //Init KSP to zero
  memset(KSP, 0, NSURF*(sizeof(unsigned)));
  
  unsigned KS,KF;
  double A,B,C;

  int KB0 = NBODYS;
  bool Eixir = false;
  while(!Eixir)
    {
      Eixir = true;
      for(unsigned int KSS = 0; KSS < bodies[KB0-1].getKSURF(NXG-1); KSS++)
	{
	  KS = bodies[KB0-1].getKSURF(KSS);
	  if(KSP[KS-1] != 0 || bodies[KB0-1].getKFLAG(KSS) > 4){ continue;}

	  FSURF(surfaces[KS-1],state,A,B,C);
	  double FUZZ;
	  double ABSA = fabs(A);
	  if(ABSA > 1.0E-36)
	    {
	      FUZZ = FUZZL*(B*B-4.0*A*C)/ABSA;
	    }
	  else
	    {
	      FUZZ = FUZZL*fabs(B);
	    }

	  if(C < -FUZZ)
	    {
	      KSP[KS-1] = 1;
	    }
	  else if(C > FUZZ)
	    {
	      KSP[KS-1] = 2;
	    }
	  else
	    {
	      //  ****  Point close to the surface.
	      if(B < 0.0)
		{
		  KSP[KS-1] = 1;  // Particle moving 'inwards'.
		}
	      else
		{
		  KSP[KS-1] = 2;  // Particle moving 'outwards'.
		}
	    }
	}

      //  ****  Determine the module or body that contains the point.

      for(unsigned int KBB = 0; KBB < bodies[KB0-1].KDGHT[NXG-1]; KBB++)
	{
	  int KB = bodies[KB0-1].KDGHT[KBB];
	  bool Eixir2 = false;
	  for(unsigned int KSS = 0; KSS < bodies[KB-1].getKSURF(NXG-1); KSS++)
	    {
	      bodies[KB-1].getSurf(KSS,KS,KF);
	      if(KF < 3 && KSP[KS-1] != KF){ Eixir2 = true; break;}
	    }
	  if(Eixir2){continue;}
	  if(KB == KB0)
	    {
	      state.IBODY = KB-1;  // The particle is inside the body or module KB.
	      state.MAT = bodies[KB-1].MATER;
	      return;
	    }
	  else if(bodies[KB-1].KDGHT[NXG-1] > 1)
	    {
	      KB0 = KB;    // The point is inside a submodule.
	      Eixir = false;
	      break;
	    }
	  else
	    {
	      state.IBODY = KB-1;  // The particle is inside a sister body or module.
	      state.MAT = bodies[KB-1].MATER;
	      return;
	    }
	}
    }
  state.IBODY = getBodies();  // The particle is outside the enclosure.
  state.MAT = 0;
}

//  *********************************************************************
//                       SUBROUTINE LOCATE
//  *********************************************************************
void pen_quadricGeo::step(pen_particleState& state,
			  double DS,
			  double &DSEF,
			  double &DSTOT,
			  int &NCROSS) const
{
  //     This subroutine handles the geometrical part of the track simula-
  //  tion. The particle starts from the point (X,Y,Z) and travels a length
  //  DS in the direction (U,V,W) within the material where it moves. When
  //  the track leaves the initial material, the particle is stopped just
  //  after entering the next material body (void regions with MAT=0 are
  //  crossed automatically). Furthermore, when the particle arrives from
  //  a void region, it is stopped just after entering the first material
  //  body.
  //
  //  Input values (module TRACK_mod):
  //     X, Y, Z ... coordinates of the initial point,
  //     U, V, W ... direction cosines of the displacement,
  //     IBODY ..... body where the initial point is located,
  //     MAT ....... material in body IBODY.
  //  NB: When a particle track is started, the variables IBODY and MAT
  //  must be set by calling subroutine LOCATE.
  //
  //  Input argument:
  //     DS ........ path length to travel.
  //
  //  Output arguments:
  //     DSEF....... travelled path length before leaving the initial
  //                 material or completing the jump (less than DS if the
  //                 track crosses an interface),
  //     NCROSS .... = 0 if the whole step is contained in the initial
  //                   material,
  //                 .gt.0 if the particle has crossed an interface, i.e.
  //                   if it has entered a new material.
  //
  //  Output values (module TRACK_mod):
  //     X, Y, Z ... coordinates of the final position,
  //     IBODY ..... body where the final point is located,
  //     MAT ....... material in IBODY. The value MAT=0 indicates that the
  //                 particle has escaped from the system.
  //
  //  Output values (module PENGEOM_mod):
  //     DSTOT ..... travelled path length, including path segments in void
  //                 volumes.
  //     KSLAST .... when NCROSS.ne.0, the output value of KSLAST is the
  //        label of the last surface crossed by the particle before
  //        entering a material body. KSLAST is used for rendering in 3D.
  

  unsigned KSLAST;

  //Create an array to store  <surface distance (S),surface index (IS)>
  pen_surfDS surfDS[NS2M];

  //Declare array for surfaces position from particle  view (forward 1 backward 2)
  unsigned KSP[NS];

  //Init KSP to zero
  memset(KSP, 0, NSURF*(sizeof(unsigned)));
  
  double DSRES;
  int NSCT,IERR,NERR;
  //unsigned KB1;
  const pen_quadBody* pbody;
  
  DSEF = 0.0;
  DSTOT = 0.0;
  NCROSS = 0;
  KSLAST = 0;
  int KSLST=0;

  //******************//
  // Increment IBODY  //
  // Because FORTRAN code traduction, internally,
  // PENELOPE geometry code uses 1 as first index
  ++state.IBODY;
  
  int NSC = 0;        // Number of surface crossings ahead of the particle.
  unsigned MAT0 = state.MAT;     // Initial material.
  if(state.MAT == 0)
    {
      DSRES = 1.0E35;  // In vacuum particles fly freely.
    }
  else
    {
      DSRES = DS;   // Residual path length.
    }

  if(state.IBODY > NBODYS)
    {
      //  ************  The particle enters from outside the enclosure.

      //This variable checks if the particle has entered the enclosure
      bool inEnclosure = false;
      //Calculate the number of surface crossings
      pbody = bodies + (NBODYS-1);
      STEPSI(pbody,state,KSP,surfDS,NSC,KSLAST);
      if(NSC == 0)
	{
	  //No surface has been crossed
	  double DSP=1.0E36;
	  state.MAT = 0;
	  if(KSLST > 0){KSLAST=KSLST;}
	  if(state.MAT == MAT0){ DSEF = DSEF+DSP;}
	  DSTOT = DSTOT+DSP;
	  move(DSP,state);
	  state.IBODY = NBODYS; //Return NBODYS+1-1 i.e. NBODYS
	  return;	  
	}
      NSCT = NSC;
      for(int KI = NSCT-1; KI >= 0; --KI)
	{
	  //  ****  The particle crosses a surface.
	  KSLAST = surfDS[KI].IS; //IS[KI];
	  if(KSP[KSLAST-1] == 1)
	    {
	      KSP[KSLAST-1] = 2;  //Backward to Forward
	    }
	  else
	    {
	      KSP[KSLAST-1] = 1;  //Forward to Backward
	    }
	  double DSP = surfDS[KI].S; //S[KI];
	  DSEF = DSEF+DSP;
	  DSTOT = DSTOT+DSP;
	  move(DSP,state); //Move the particle onto the surface
	  NSC = NSC-1;

	  //Update distances to surfaces
	  for(int I = 0; I < NSC; I++)
	    {
	      surfDS[I].S -= DSP;
	    }
	    

	  //Check if the particle has entered the enclosure
	  if(inBody(pbody,KSP)){

	    //  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	    //  ****  The particle enters the enclosure.
	    inEnclosure = true;
	    
	    //Go deeper or into the inner part/zone to find a material body
	    goInner(pbody,state,KSP,surfDS,NSC,KSLAST);
	    if(state.MAT != 0)
	      {
		//  ****  The particle enters a material body.
		NCROSS = 1;
		state.IBODY--;
		return;
	      }

	    //The particle remains in the void of the enclosure
	    STEPSI(pbody,state,KSP,surfDS,NSC,KSLAST);
	    break;
	  }
	}
      //  ****  At this point the program has left the DO loop.
      //  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      if(!inEnclosure)
	{
	  //The particle has not entered to the enclosure
	  double DSP=1.0E36;
	  state.IBODY = NBODYS+1;
	  state.MAT = 0;
	  if(KSLST > 0){KSLAST=KSLST;}
	  if(state.MAT == MAT0){ DSEF = DSEF+DSP;}
	  DSTOT = DSTOT+DSP;
	  move(DSP,state);
	  state.IBODY--;
	  return;
	}
    }
  else
    {
      //  ************  Surface crossings.

      //  ************  The particle is inside some body.
    
      int IBODYL = state.IBODY;
      if(LVERB){ NERR = 0;}

      //Obtain surfaces cross information
      pbody = bodies + (state.IBODY-1);
      
      STEPSI(pbody,state,KSP,surfDS,NSC,KSLAST);
      STEPLB(pbody,state,KSP,IERR);
      
      while(IERR != 0)
	{
	  //  ****  Evidence of round-off errors.
	  if(NSC > 0)
	    {
	      //  ****  When a surface is very close, we move the particle beyond it.
	      if(surfDS[NSC-1].S < 1.0E-10)
		{
		  KSLAST = surfDS[NSC-1].IS;
		  if(KSP[KSLAST-1] == 1)
		    {
		      KSP[KSLAST-1] = 2;
		    }
		  else
		    {
		      KSP[KSLAST-1] = 1;
		    }
		  double DSP = surfDS[NSC-1].S;
		  move(DSP,state);
		  if(state.MAT == MAT0)
		    {
		      DSEF = DSEF+DSP;
		      DSRES = DSRES-DSP;
		    }
		  DSTOT = DSTOT+DSP;
		  NSC = NSC-1;
		}
	    }
	  if(LVERB)
	    {
	      NERR = NERR+1;
	      if(state.MAT != 0)
		{
		  printf("*** WARNING, STEP: Accidental undershot or round-off error?\n    IBODY0 =%5d,  IBODY =%5d,  IERR =%3d,  NERR =%4d\n", int((pbody-bodies)+1), state.IBODY, IERR, NERR);		  
		  printf("(X,Y,Z) = %15.7E,%15.7E,%15.7E\n", state.X, state.Y, state.Z);
		  printf("(U,V,W) = %15.7E,%15.7E,%15.7E\n", state.U, state.V, state.W);

		  for(unsigned KSS = 0; KSS < pbody->getKSURF(NXG-1); ++KSS)
		    {
		      double A,B,C,SW;
		      unsigned KS;
		      unsigned KFLO;
		      pbody->getSurf(KSS,KS,KFLO);
		      
		      if(KFLO < 3)
			{
			  SW = 0.0;
			  for(int KI = NSC-1; KI >= 0; KI--)
			    {
			      if(KS == surfDS[KI].IS)
				{
				  SW = surfDS[KI].S;
				  break;
				}
			    }
			  FSURF(surfaces[KS-1],state,A,B,C);
			  if(KFLO == KSP[KS-1])
			    {
			      printf("Surface,FL0,FL,C,S =%5d  %3d%3d%15.7E%15.7E\n", KS, KFLO, KSP[KS-1], C, SW);
			    }
			  else
			    {
			      printf("Surface,FL0,FL,C,S =%5d  %3d%3d%15.7E%15.7E *\n", KS, KFLO, KSP[KS-1], C, SW);			      
			      KSLAST = KS;
			    }
			}
		    }
		}
	    }

	  if(state.IBODY > NBODYS){
	    //Scapes from enclosure

	    //NOTA: Açò estarà per seguretat, però no crec que puga pasar mai.
	    //L'objecte global que rodege el mon deuria tindre una superfície
	    //molt pròxima a un altra superfície d'un objecte. De fet, si aixo
	    //passa, salvat acabarà el step en el següent "if" fora del bucle,
	    //ja que state.MAT != MAT0. I dirà que no s'ha mogut (nota el DSEF=0.0),
	    //quan ell acostuma a moure la partícula 1.0e35. De moment, copie el
	    //que fa el if de baix.

	    NCROSS = 1;
	    DSEF = 0.0;
	    state.IBODY--;
	    return;	    
	  }

	  //Obtain surfaces cross information
	  pbody = bodies + (state.IBODY-1);
	  STEPSI(pbody,state,KSP,surfDS,NSC,KSLAST);
	  STEPLB(pbody,state,KSP,IERR);	  
	}

      if(bodies[state.IBODY-1].KDET != bodies[IBODYL-1].KDET || state.MAT != MAT0)
	{
	  //The particle crosses to another detector or material
	  NCROSS = 1;
	  DSEF = 0.0;
	  state.IBODY--;
	  return;
	}

      //  ****  The particle remains in the same material.

      if(state.MAT != 0 && DSRES < surfDS[NSC-1].S)
	{
	  DSEF = DSEF+DSRES;
	  DSTOT = DSTOT+DSRES;
	  move(DSRES,state);
	  state.IBODY--;
	  return;
	}
    }
  //  ************  New position.

  while(NSC > 0)
    {
      NSCT = NSC;
      unsigned MATL = state.MAT;
      unsigned IBODYL = state.IBODY;
      //Iterate over all surface intersections
      for(int KI = NSCT-1; KI >= 0; KI--)
	{
	  if(DSRES < surfDS[KI].S)
	    {
	      //  ****  The step ends within the body.
	      if(state.MAT == MAT0){ DSEF = DSEF+DSRES;}
	      DSTOT = DSTOT+DSRES;
	      //Move the particle
	      move(DSRES,state);
	      state.IBODY--;
	      return;
	    }
	  
	  //  ****  The particle crosses a surface.

	  //Update crossed surface side pointer
	  KSLAST = surfDS[KI].IS; 
	  if(KSP[KSLAST-1] == 1)
	    {
	      KSP[KSLAST-1] = 2;
	    }
	  else
	    {
	      KSP[KSLAST-1] = 1;
	    }
	  double DSP = surfDS[KI].S;

	  //Move the particle
	  move(DSP,state);

	  //Update traveled distances
	  if(state.MAT == MAT0)
	    {
	      DSEF = DSEF+DSP;
	      DSRES = DSRES-DSP;
	    }
	  DSTOT = DSTOT+DSP;

	  //Update distances to surfaces
	  NSC = NSC-1;
	  for(int I = 0; I < NSC; I++)
	    {
	      surfDS[I].S -= DSP;
	    }
	    
	  //Go inner the geometry tree to find a material body
	  if(goInner(pbody,state,KSP,surfDS,NSC,KSLAST)){
	    //  ****  The particle leaves the enclosure.
	    if(state.MAT != MATL){ ++NCROSS;}
	    DSP = 1.0E36;
	    state.MAT = 0;
	    if(KSLST > 0){KSLAST=KSLST;}
	    if(state.MAT == MAT0){ DSEF = DSEF+DSP;}
	    DSTOT = DSTOT+DSP;
	    
	    //Move the particle
	    move(DSP,state);

	    state.IBODY = NBODYS;
	    return;
	  }

	  //  ****  The particle continues flying when it enters a void region...
	  if(state.MAT == 0)
	    {
	      if(MATL != 0){KSLST=KSLAST;}
	      if(MATL == MAT0){ NCROSS = NCROSS+1;}
	      MATL = 0;
	      DSRES = 1.0E35;
	    }
	  else if(state.MAT == MATL)
	    {
	      //  ****  The particle continues flying when it enters a new body of the
	      //        same material which is not part of a different detector...      
	      if(bodies[state.IBODY-1].KDET != bodies[IBODYL-1].KDET)
		{
		  if(KSLST != 0){KSLAST=KSLST;}
		  NCROSS = NCROSS+1;
		  state.IBODY--;
		  return;		  
		}
	    }
	  else   // IF(MAT.NE.MATL) THEN
	    {
	      //  ****  ... and stops when it penetrates a new material body or a
	      //        detector.	      
	      if(MATL == 0)   // Ensure the point is inside the body.
		{
		  int ITST=2*KSP[KSLAST-1]-3;
		  double A,B,C;
		  FSURF(surfaces[KSLAST-1],state,A,B,C);
		  if(ITST*C < FUZZT)
		    {
		      double DSFZ=state.posMod()*FUZZL;  // Additional path length.
		      int ITRY=0;
		      do{
			//Move the particle
			move(DSFZ,state);
			DSTOT=DSTOT+DSFZ;
			if(ITRY < 5)
			  {
			    ++ITRY;
			    FSURF(surfaces[KSLAST-1],state,A,B,C);
			  }
			else
			  {
			    printf("*** WARNING, STEP: The particle has not entered the body after %2d trials.\n",ITRY);
			    printf("    Current conditions: IBODY =%4d, MAT =%3d, ITST*C =%13.5E.\n",state.IBODY,state.MAT,C*ITST);
			    printf("    X,Y,Z =%16.7E%16.7E%16.7E\n",state.X,state.Y,state.Z);
			    printf("    U,V,W =%16.7E%16.7E%16.7E\n",state.U,state.V,state.W);
			    break;
			  }
			    
			}while(ITST*C < FUZZT);
		    }
		}
	      if(KSLST > 0){KSLAST=KSLST;}
	      NCROSS = NCROSS+1;
	      state.IBODY--;
	      return;
	    }

	  //Get surface crosses information of current body
	  STEPSI(pbody,state,KSP,surfDS,NSC,KSLAST);
	  break;
	}
      //  ****  At this point the program has left the DO loop.
      //  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<      
    }

  if(NSC == 0)
    {
      //No surface crossings
      if(state.MAT == MAT0){ DSEF = DSEF+DSRES;}
      DSTOT = DSTOT+DSRES;
      //Move the particle
      move(DSRES,state);
      state.IBODY--;
      return;
    }  
  
  //  ************  The particle leaves the enclosure.

  double DSP = 1.0E36;
  state.IBODY = NBODYS;
  state.MAT = 0;
  if(KSLST > 0){KSLAST=KSLST;}
  if(state.MAT == MAT0){ DSEF = DSEF+DSP;}
  DSTOT = DSTOT+DSP;
  //Move the particle
  move(DSP,state);
}

//  *********************************************************************
//                       SUBROUTINE STEPSI
//  *********************************************************************
void pen_quadricGeo::STEPSI(const pen_quadBody* pbody, const pen_particleState& state, unsigned KSP[NS], pen_surfDS surfDS[NS2M], int &NSC_IO, unsigned& KSLAST_IO) const
{
  //     Calculates the intersections of the trajectory with the limiting
  //  surfaces of body KB. The intersections are added to the list and
  //  sorted in decreasing order. This subroutine works only when called
  //  from within subroutine STEP.


  // const double FUZZL = 1.0E-12;

  //Copy initial values
  int NSC = NSC_IO;
  unsigned KSLAST = KSLAST_IO;
  
  //  ************  Determine surface crossings.

  const pen_bodySurf* psurfLast = pbody->surfs+pbody->getKSURF(NXG-1);
  for(const pen_bodySurf* psurf = pbody->surfs; psurf != psurfLast; ++psurf)
    {
      //  ****  Intersections with a given surface are calculated only once.
      //        The side pointer of a surface must be changed each time the
      //        surface is crossed.
      if(psurf->KFLAG > 4){ continue;}
      const unsigned KS = psurf->KSURF;
      const unsigned KSm1 = KS-1;
      unsigned& KSP_KSm1 = KSP[KSm1];
      if(KSP_KSm1){ continue;} // if(KSP[KSm1] != 0)

      double A,B,C,ABSA;
      //FSURF(surfaces[KSm1],state,A,B,C);

      //*** FSURF inclusion ***
      //***********************
      const pen_quadSurface& surface = surfaces[KSm1];
      
      if(surface.KPLANE){ // if(surface.KPLANE != 0)
	A = ABSA = 0.0;
	B = state.U*surface.AX+state.V*surface.AY+state.W*surface.AZ;
	C = state.X*surface.AX+state.Y*surface.AY+state.Z*surface.AZ+surface.A0;
      }
      else{
	A = state.U*(surface.AXX*state.U+
		     surface.AXY*state.V+
		     surface.AXZ*state.W)+
	  state.V*(surface.AYY*state.V+
		   surface.AYZ*state.W)+
	  state.W*surface.AZZ*state.W;
      
	double XXX = surface.AXX*state.X+surface.AXY*state.Y+
	  surface.AXZ*state.Z+surface.AX;
      
	double YYY = surface.AYY*state.Y+surface.AYZ*state.Z+surface.AY;
	double ZZZ = surface.AZZ*state.Z+surface.AZ;

	B = state.U*(surface.AXX*state.X+XXX)+
	  state.V*(surface.AXY*state.X+surface.AYY*state.Y+YYY)+
	  state.W*(surface.AXZ*state.X+surface.AYZ*state.Y+surface.AZZ*state.Z+
		   ZZZ);

	C = state.X*XXX+state.Y*YYY+state.Z*ZZZ+surface.A0;
	
	ABSA = fabs(A);
      }

      //***   FSURF END    ***
      //**********************
      
      double ABSB = fabs(B);

      //  ****  Plane, a single root.

      if(ABSA < 1.0E-36)
	{
	  if(ABSB > 0.0)
	    {
	      if(C < mFUZZL)  // SP=-1
		{
		  KSP_KSm1 = 1;
		}
	      else if(C > FUZZL)  // SP=+1
		{
		  KSP_KSm1 = 2;
		}
	      else  // Point close to the surface.
		{
		  KSLAST=KS;
		  if(B < 0.0)
		    {
		      KSP_KSm1 = 1;  // Particle moving 'inwards'.
		    }
		  else
		    {
		      KSP_KSm1 = 2;  // Particle moving 'outwards'.
		    }
		  continue;
		}
	      double T1 = -C/B;
	      if(T1 > 0.0){
		surfDS[NSC++].set(T1,KS);
	      }
	    }
	  else  // Ray parallel to the plane.
	    {
	      if(C < 0.0)
		{
		  KSP_KSm1 = 1;
		}
	      else
		{
		  KSP_KSm1 = 2;
		}
	    }
	  
	  //  ****  Non-planar surface, two roots.
	}
      else
	{
	  int IAMBIG;
	  double DISCR = B*B-4.0*A*C;
	  double FUZZ = FUZZL*DISCR/ABSA;
	  if(C < -FUZZ)  // SP=-1
	    {
	      IAMBIG = 0;
	      KSP_KSm1 = 1;
	    }
	  else if(C > FUZZ)  // SP=+1
	    {
	      IAMBIG = 0;
	      KSP_KSm1 = 2;
	    }
	  else
	    {
	      IAMBIG = 1;  // Point close to the surface.
	      KSLAST=KS;
	      if(B < 0.0)
		{
		  KSP_KSm1 = 1;  // Particle moving 'inwards'.
		}
	      else
		{
		  KSP_KSm1 = 2;  // Particle moving 'outwards'.
		}
	    }

	  if(DISCR < 1.0E-36){ continue;}  // No true intersections.

	  if(IAMBIG == 0)
	    {
	      double R2A = 0.5/A;
	      double DELTA = sqrt(DISCR)*fabs(R2A);
	      double SH = -B*R2A;
	      double T1 = SH-DELTA;
	      if(T1 > 0.0){
		surfDS[NSC++].set(T1,KS);
	      }
	      double T2 = SH+DELTA;
	      if(T2 > 0.0){
		surfDS[NSC++].set(T2,KS);
	      }
	    }
	  else
	    {
	      if(B*A < 0.0){
		double R2A = 0.5/A;
		double DELTA = sqrt(DISCR)*fabs(R2A);
		double SH = -B*R2A;
		double T2 = SH+DELTA;

		surfDS[NSC++].set((T2 > 0.0 ? T2 : 0.0),KS);
	      }
	    }
	}
    }

  //  ****  Sort surface distances in decreasing order.

  if(NSC > 1)
    std::sort(surfDS,surfDS+NSC,std::greater<pen_surfDS>());
  
  /*
  if(NSC > 1)
    {
      int NSCm1 = NSC-1;
      for(int KI = 0; KI < NSCm1; ++KI)
	{
	  double SMAX = surfDS[KI].S;
	  int KMAX = KI;
	  for(int KJ = KI+1; KJ < NSC; ++KJ)
	    {
	      if(surfDS[KJ].S > SMAX)
		{
		  SMAX = surfDS[KJ].S; 
		  KMAX = KJ;
		}
	    }
	  if(KMAX != KI)
	    {
	      pen_surfDS aux = surfDS[KI];
	      surfDS[KI] = surfDS[KMAX];
	      surfDS[KMAX] = aux;
	    }
	}
    }
  */
  
  //Get final output values
  NSC_IO = NSC;
  KSLAST_IO = KSLAST;
  
}

//  *********************************************************************
//                       SUBROUTINE STEPLB
//  *********************************************************************
void pen_quadricGeo::STEPLB(const pen_quadBody* pbody, pen_particleState& state, const unsigned KSP[NS], int &IERR) const
{

  //     Helps finding the body or module that has the given side pointers
  //  for the analyzed surfaces. Subroutine STEPLB works only when invoked
  //  from within subroutine STEP. It moves through the tree of modules a
  //  single step.


  //  ****  Analyze the current body or module.
  //const pen_quadBody &body = bodies[KB-1];

  //Check if the current body is a body or a module
  if(pbody->KBOMO == 0)
    {
      //  ****  Body.

      //Get inner bodies
      int NLBOD = pbody->KBODY[NXG-1];

      for(int KBB = 0; KBB < NLBOD; KBB++)
	{
	  //Get next inner body index and the body itself
	  unsigned KBS = pbody->KBODY[KBB];
	  const pen_quadBody& bodyKBS = bodies[KBS-1];

	  //Check if the particle is in that inner body
	  if(inBody(&bodyKBS,KSP)){
	    //Is in this body, check if is a module or a final material body
	    state.IBODY = KBS;
	    if(bodyKBS.KDGHT[NXG-1] > 1)
	      {
		IERR = -1;  // The particle is inside a sister module.
	      }
	    else
	      {
		IERR = 0;   // The particle is inside a sister body.
		state.MAT = bodyKBS.MATER;
	      }
	    return;		
	  }
	}
	

      //The particle is not in any inner body nor module, check if is
      //in the current body
      if(inBody(pbody,KSP)){
	state.IBODY = int(pbody-bodies)+1;
	IERR = 0;   // The particle remains in the same body.
	state.MAT = bodies[state.IBODY-1].MATER;
	return;	
      }
    }
  else
    {
      //  ****  Module.

      //Iterate over inner daughters
      for(unsigned KBB = 0; KBB < pbody->KDGHT[NXG-1]; ++KBB)
	{
	  //Get next inner body index and object
	  unsigned KBD = pbody->KDGHT[KBB];
	  const pen_quadBody* pbodyKBD = bodies + (KBD-1); 

	  //Check if the particle is in that body
	  if(inBody(pbodyKBD,KSP)){
	    //Is inside, check if remains witihn the current module,
	    //or if is inside an other body/module
	    state.IBODY = KBD;
	    if(pbodyKBD == pbody)
	      {
		IERR = 0;  //  The particle remains within the current module.
		state.MAT = pbodyKBD->MATER;
	      }
	    else
	      {
		if(pbodyKBD->KDGHT[NXG-1] > 1)
		  {
		    IERR = -1;  // The particle is inside a submodule.
		  }
		else
		  {
		    IERR = 0;   // The particle is inside a single body.
		    state.MAT = pbodyKBD->MATER;
		  }
	      }
	    return;	    
	  }
	}
    }

  //  ****  The particle is outside the current body or module.

  IERR = 1;
  state.IBODY = pbody->KMOTH;
  if(state.IBODY == 0)
    {
      state.IBODY = NBODYS+1;
      state.MAT = 0;
    }
}


//  *********************************************************************
//                       SUBROUTINE ROTSHF
//  *********************************************************************
void ROTSHF(double OMEGA, double THETA, double PHI, double DX, double DY, double DZ, double &AXX, double &AXY, double &AXZ, double &AYY, double &AYZ, double &AZZ, double &AX, double &AY, double &AZ, double &A0)
{
  //     This subroutine rotates and shifts a quadric surface.

  //  Input parameters:
  //     OMEGA, THETA, PHI ... Euler rotation angles,
  //     DX, DY, DZ .......... components of the displacement vector,
  //     AXX, ..., A0 ........ coefficients of the initial quadric.

  //  Output parameters:
  //     AXX, ..., A0 ........ coefficients of the transformed quadric.

  double R[3][3], A2[3][3], B2[3][3], A1[3], B1[3], D1[3];
  double B0,STHETA,CTHETA,SPHI,CPHI,SOMEGA,COMEGA,A2D;

  //  ****  Initial quadric.

  B2[0][0] = AXX;
  B2[0][1] = 0.5*AXY;
  B2[0][2] = 0.5*AXZ;
  B2[1][0] = B2[0][1];
  B2[1][1] = AYY;
  B2[1][2] = 0.5*AYZ;
  B2[2][0] = B2[0][2];
  B2[2][1] = B2[1][2];
  B2[2][2] = AZZ;
  B1[0] = AX;
  B1[1] = AY;
  B1[2] = AZ;
  B0 = A0;
  D1[0] = DX;
  D1[1] = DY;
  D1[2] = DZ;

  //  ****  Rotation matrix.

  STHETA = sin(THETA);
  CTHETA = cos(THETA);
  SPHI = sin(PHI);
  CPHI = cos(PHI);
  SOMEGA = sin(OMEGA);
  COMEGA = cos(OMEGA);

  R[0][0] = CPHI*CTHETA*COMEGA-SPHI*SOMEGA;
  R[0][1] = -CPHI*CTHETA*SOMEGA-SPHI*COMEGA;
  R[0][2] = CPHI*STHETA;
  R[1][0] = SPHI*CTHETA*COMEGA+CPHI*SOMEGA;
  R[1][1] = -SPHI*CTHETA*SOMEGA+CPHI*COMEGA;
  R[1][2] = SPHI*STHETA;
  R[2][0] = -STHETA*COMEGA;
  R[2][1] = STHETA*SOMEGA;
  R[2][2] = CTHETA;

  //  ****  Rotated quadric.

  for(int I = 0; I < 3; I++)
    {
      A1[I] = 0.0;
      for(int J = 0; J < 3; J++)
	{
	  A1[I] = A1[I]+R[I][J]*B1[J];
	  A2[I][J] = 0.0;
	  for(int M = 0; M < 3; M++)
	    {
	      for(int K = 0; K < 3; K++)
		{
		  A2[I][J] = A2[I][J]+R[I][K]*B2[K][M]*R[J][M];
		}
	    }
	}
    }

  //  ****  Shifted-rotated quadric.

  for(int I = 0; I < 3; I++)
    {
      A2D = 0.0;
      for(int J = 0; J < 3; J++)
	{
	  A2D = A2D+A2[I][J]*D1[J];
	}
      B1[I] = A1[I]-2.0*A2D;
      B0 = B0+D1[I]*(A2D-A1[I]);
    }

  AXX = A2[0][0];
  AXY = A2[0][1]+A2[1][0];
  AXZ = A2[0][2]+A2[2][0];
  AYY = A2[1][1];
  AYZ = A2[1][2]+A2[2][1];
  AZZ = A2[2][2];
  AX = B1[0];
  AY = B1[1];
  AZ = B1[2];
  A0 = B0;
  if(fabs(AXX) < 1.0E-16){ AXX = 0.0;}
  if(fabs(AXY) < 1.0E-16){ AXY = 0.0;}
  if(fabs(AXZ) < 1.0E-16){ AXZ = 0.0;}
  if(fabs(AYY) < 1.0E-16){ AYY = 0.0;}
  if(fabs(AYZ) < 1.0E-16){ AYZ = 0.0;}
  if(fabs(AZZ) < 1.0E-16){ AZZ = 0.0;}
  if(fabs(AX) < 1.0E-16){ AX = 0.0;}
  if(fabs(AY) < 1.0E-16){ AY = 0.0;}
  if(fabs(AZ) < 1.0E-16){ AZ = 0.0;}
  if(fabs(A0) < 1.0E-16){ A0 = 0.0;}
}
