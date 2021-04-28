#include "pen_phaseSpaceFile.hh" 
#include "iaea_phsp.h"

//Program to convert IAEA phase space file format 
//into PenRed PSF format (binary to binary conversion)

int main(int argc, char** argv){
    
    
    if(argc < 3){
        printf("usage: %s input filename (without extension) and output file\n", argv[0]);
        return -1;
    }
    
    
  //Read IAEA phase space file and write it to PenRed format phase space file
  //Read name of input IAEA file
  char* fin;
  fin = argv[1];


  //Open output file
  FILE* fout = nullptr;
  fout = fopen(argv[2],"wb");
  if(fout == nullptr){
    printf("Error: unable to open output file '%s'.\n",argv[2]);
    return -2;
  }
  
  
  //PSF read initialization: assign an unique Id to IAEA source (psfRead)
  IAEA_I32 psfRead;
  IAEA_I32 len = 81;
  //access = 1 => opening read-only file
  const IAEA_I32 accessRead = 1;
  IAEA_I32 err;
  iaea_new_source(&psfRead,fin,&accessRead,&err,len);
  if(err<0){
      printf("Error: unable to open input IAEA PSF.\n");
      return -3;
  }
  else
  {
      printf("IAEA PSF reading initialized.\n");
      printf("\n");
  }
 
  //PSF maximum number of particles
  //If type < 0, iaeaMaxnPart is set to the total number of all particles.
  const IAEA_I32 type = -1;
  IAEA_I64 iaeaMaxnPart;
  iaea_get_max_particles(&psfRead,&type,&iaeaMaxnPart);
  if(iaeaMaxnPart<0){
      printf("Error: unable to read the total number of particles from IAEA PSF.\n");
      return -4;
  }
  
  //PSF total number of original particles
  IAEA_I64 iaeaOriginalHist;
  IAEA_I64 maxNhist=9.0e18;
  iaea_get_total_original_particles(&psfRead,&iaeaOriginalHist);
  if(iaeaOriginalHist<0){
      printf("Error: unable to read the total number of original histories from IAEA PSF.\n");
      return -5;
  }
  if(iaeaOriginalHist>maxNhist){
      printf("Error: original histories from IAEA PSF exceeds the maximum number of histories allowed.\n");
      return -6;
  }

  //PSF maximum particle energy
  IAEA_Float iaeaEmax;
  iaea_get_maximum_energy(&psfRead,&iaeaEmax);
  if(iaeaEmax < 0.0)
  {
    printf("Error: unable to get maximum energy from IAEA PSF\n");   
    return -7;
  }

  
  //Get additional floats and integers 
  IAEA_I32 nextraFloat, nextraInt;
  IAEA_I32 nextraMax=10;
  IAEA_I32 nmaxExtraLong=7;
  IAEA_I32 nmaxExtraFloat=3;
  iaea_get_extra_numbers(&psfRead,&nextraFloat,&nextraInt);
  if(nextraFloat<0 || nextraInt<0)
  {
      printf("Error: unable to read additional variables. No header file found.\n");
      return -8;
  }
  else if(nextraFloat>nmaxExtraFloat || nextraInt>nmaxExtraLong){
      printf("Error: more than %u additional float or integers variables from "
             "the IAEA PSF are read. nextraFloat = %u > nmaxExtraFloat = %u,"
             "or nextraInt = %u > nmaxExtraLong = %u.\n",
             nextraMax, nextraFloat, nmaxExtraFloat, nextraInt, nmaxExtraLong);
      printf("Check the PSF original format of the Monte Carlo code "
              "used to create this IAEA PSF.\n");
      return -9;
  }
  else
  {
        printf("Number of extra variables read: %u floats and %u integers.\n",nextraFloat, nextraInt);
        if(nextraInt>0)
            printf("Extra integers at PSF:\n");
  }
  
  //Get the type of all additional variables
  //Declare variables to save the type and variables to save the values of the additional variables
  IAEA_I32 *extraIntTypes,*extraFloatTypes, *extraInt;
  IAEA_Float *extraFloat;
  
  //Get memory
  extraFloatTypes = (IAEA_I32*)malloc(nextraFloat*sizeof(IAEA_I32));
  extraFloat      = (IAEA_Float*)malloc(nextraFloat*sizeof(IAEA_Float));
  extraIntTypes   = (IAEA_I32*)malloc(nextraInt*sizeof(IAEA_I32));
  extraInt        = (IAEA_I32*)malloc(nextraInt*sizeof(IAEA_I32));
  
  iaea_get_type_extra_variables(&psfRead,&err,extraIntTypes,extraFloatTypes);
  
  bool extraVariables[7]={false};
  unsigned long extraPosition[7] = {0};
  //Int extraIntTypes where ILB values are stored if exists
  for(IAEA_I32 i=0; i<nextraInt; ++i){
      if(extraIntTypes[i]==0){
        printf("User defined generic integer type.\n");
      }
      else if(extraIntTypes[i]==1){
          //Incremental history number (EGS, PENELOPE)
          //=0 nonprimary event
          //>0 primary particles
          printf("Incremental history number, "
                "indicating primary or secondary particles.\n");
          extraVariables[0]= true;
          extraPosition[0]=i;
      }
      else if(extraIntTypes[i]==2){
          //LATCH, EGS variable, unused by PenRed
          extraVariables[1]= true;
          extraPosition[1]=i;
      }
      else if(extraIntTypes[i]==3){
          //ILB5
          printf("ILB(5), from PENELOPE simulation.\n");
          extraVariables[2]= true;
          extraPosition[2]=i;
      }
      else if(extraIntTypes[i]==4){
          //ILB4
          printf("ILB(4), from PENELOPE simulation.\n");
          extraVariables[3]= true;
          extraPosition[3]=i;
      }
      else if(extraIntTypes[i]==5){
          //ILB3
          printf("ILB(3), from PENELOPE simulation.\n");
          extraVariables[4]= true;
          extraPosition[4]=i;
      }
      else if(extraIntTypes[i]==6){
          //ILB2
          printf("ILB(2), from PENELOPE simulation.\n");
          extraVariables[5]= true;
          extraPosition[5]=i;
      }
      else if(extraIntTypes[i]==7){
          //ILB1
          printf("ILB(1), from PENELOPE simulation.\n");
          extraVariables[6]= true;
          extraPosition[6]=i;
      }
      else{
          printf("PSF contains undefined variable.\n");
      }
      
  }
  printf("\n");
  
  //No prints are shown for FloatTypes. 
  //These types are from EGS simulation and unused by PenRed
  
  
  //Check file size and byte order. It must be equal to the checksum header value
  printf("IAEA PSF size information:\n");
  iaea_check_file_size_byte_order(&psfRead,&err);
  if(err==-1){
      printf("Error: unable to check file size, header does not exist.\n");
      return -10;
  }
  else if(err==-2){
      printf("Error: unable to check file size, function fseek fails for some reason.\n");
      return -11;
  }
  else if(err==-3){
      printf("Error: file size mismatch. PSF size and checksum of the header differs.\n");
      return -12;
  }
  else if(err==-4){
      printf("Error: byte order mismatch. PSF byte order and the byte order of the machine being run on differs.\n");
      return -13;
  }
  else if(err==-5){
      printf("Error: both file size and byte order mismatch.\n");
      printf("PSF size and checksum of the header differs.\n");
      printf("PSF byte order and the byte order of the machine being run on differs.\n");
      return -14;
  }
      
      

  //Create phase space file writer
  pen_psfwriter psfWriter;
  
  //Create particle state and kpar
  pen_particleState state;
  pen_KPAR kpar;
  
  //Declare hist increment and total number of histories
  unsigned long long dhist;
  unsigned long long nhist = 1; //Ensure non zero history number
  
  //Create chunk counter
  size_t chunk = 0;
  
  //Create particle counter
  unsigned long part[3];
  part[0]=0;
  part[1]=0;
  part[2]=0;
  
  printf("\n");
  printf("Converting IAEA PSF into PenRed PSF format.\n");
  printf("\n");
  
  //Loop to read particles from IAEA PSF 
  for(IAEA_I32 partReaded=0; partReaded<iaeaMaxnPart; ++partReaded)
  {
    //Get IAEA particle state  
    IAEA_Float iaeaE, iaeaWght,iaeaX,iaeaY,iaeaZ,iaeaU,iaeaV,iaeaW;
    IAEA_I32 iaeaKpar;
        
    iaea_get_particle(&psfRead,&err,&iaeaKpar,&iaeaE,&iaeaWght,&iaeaX,&iaeaY,&iaeaZ,&iaeaU,&iaeaV,&iaeaW,extraFloat,extraInt);
        
    if(err==-1){
        printf("Error: unable to read particle number %u.\n", partReaded+1);
        return -15;
    }
    else if(err==-2){
        printf("Error: end of file of the PSF source reached attempting to read particle number %u.\n", partReaded+1);
        return -16;
    }  
        
    //Store state data
    state.E = iaeaE*1.0E6; //Energy must be converted from MeV to eV
    state.X = iaeaX;
    state.Y = iaeaY;
    state.Z = iaeaZ;
    state.U = iaeaU;
    state.V = iaeaV;
    state.W = iaeaW;
    state.WGHT = iaeaWght;
    //Particle age is not stored in IAEA format
    state.PAGE = 0.0; 

    //Check if extra variables are present at IAEA PSF. First extra variables stores incremental history.
    if(extraVariables[0]){
        dhist = (unsigned long)extraInt[extraPosition[0]];
    }
    else{
        //History does not increment
        dhist=0; 
        if(err>0){
            //Set to 1. Each particle is a new history
            dhist=1; 
        }
    }
        
    nhist +=dhist;  
    
    //kpar number must be converted form IAEA format to PenRed format
    if(iaeaKpar==1){
        kpar  = PEN_PHOTON; //photons
        part[0]++;
    }
    else if(iaeaKpar==2){
        kpar  = PEN_ELECTRON; //electrons
        part[1]++;
    }
    else if(iaeaKpar==3){
        kpar  = PEN_POSITRON; //positrons
        part[2]++;
    }
    else{
        printf("Error: iaeaKpar = %u for particle number %u, invalid for PenRed simulation. \n", iaeaKpar, partReaded+1);
        printf("iaeaKapar = 4 for neutrons, iaeaKpar = 5 for protons. These particles are not allowed in PenRed simulations.\n");
        return -17;
    }
            
    //This value is not stored in IAEA format
    state.LAGE = 0.0; 

        
    //Check if ILB values are stored at IAEA PSF format
    if(extraVariables[6]){
        state.ILB[0] = extraInt[extraPosition[6]];
    }
    else{
        //If not information are stored at IAEA PSF, 
        //PenRed assumes primary particle
        state.ILB[0] = 1;
    }
    
    //If not information are stored at IAEA PSF, 
    //The other ILB values are assumed null
    if(extraVariables[5]){
        state.ILB[1] = extraInt[extraPosition[5]];
    }
    else{
        state.ILB[1] = 0;
    }
    
    if(extraVariables[4]){
        state.ILB[2] = extraInt[extraPosition[4]];
    }
    else{
        state.ILB[2] = 0;
    }
        
    if(extraVariables[3]){
        state.ILB[3] = extraInt[extraPosition[3]];
    }
    else{
        state.ILB[3] = 0;
    }
        
    if(extraVariables[2]){
        state.ILB[4] = extraInt[extraPosition[2]];
    }
    else{
        state.ILB[4] = 0;
    }
        
    //Note: extraVariables[1] is not used here, it corresponds to EGS LATCH varible which is unused by PenRed.
    
    
    //Write particles at PSF PenRed file
    if(psfWriter.store(nhist,kpar,state) == -1)  //Buffer is full
    {
        //printf("Chunk %lu full (%lu).\n", chunk, psfWriter.stored());
        //printf("Concluded particles (%lu).\n",psfWriter.partConcluded()); 
        //printf("Last history particles (%lu).\n",psfWriter.partLastHist());
        
        //Buffer is full, write articles info into output file
        psfWriter.write(fout,1);
        //Clear buffer
        psfWriter.clearBuffer();
        //Save last particle state that could not be written because of full buffer
        psfWriter.store(nhist,kpar,state);
        
        chunk++;
    }
  }
  
  //Write extra particles read while the last chunk is not full
  //printf("Remaining particles: %lu\n", psfWriter.stored());  
  chunk++;
  //Write remaining states
  psfWriter.write(fout,1,true);
  psfWriter.clear();
  //Close output file
  fclose(fout);
  
  
  
  long int totalPart =  part[0]+part[1]+part[2];
  printf("IAEA PSF read and converted to PenRed format.\n");
  printf("\n");
  printf("Maximum energy from IAEA PSF: %6.2E eV\n",iaeaEmax*1.0E06);   
  printf("Photons read from PSF: %lu\n",part[0]);
  printf("Electrons read from PSF: %lu\n",part[1]);
  printf("Positrons read from PSF: %lu\n",part[2]);
  printf("Total particles read from PSF: %lu\n",totalPart);
  printf("Total particles declared in IAEA PSF header: %llu\n",iaeaMaxnPart);
  printf("Total original histories declared in IAEA PSF header: %llu\n",iaeaOriginalHist);
  printf("Total number of chunks writted at penRed format PSF: %lu\n",chunk);
  
  if(iaeaMaxnPart!=totalPart)
  {
      printf("Error: number of total particles declared in IAEA PSF header (%llu)"
           " and number of total particles read from PSF (%lu), differs.\n",iaeaMaxnPart,totalPart);
      printf("Conversion terminated with errors\n");
      return -18;
  }

  
  iaea_destroy_source(&psfRead,&err);
  if(err<0)
  {
    printf("Error: unable to close IAEA PSF.\n");
    return -19;
  }

  printf("\n");
  printf("Conversion terminated successfully.\n");  
}
  
  
  
  
  
