#include "pen_phaseSpaceFile.hh" 
#include "iaea_phsp.h"

//Program to convert PenRed PSF format 
//into IAEA PSF format (binary to binary conversion)

int main(int argc, char** argv){
    
    
    if(argc < 3){
        printf("usage: %s input filename and output file (without extension)\n", argv[0]);
        return -1;
    }
    
    
    //Open input PSF pneRed file
    FILE* fin = nullptr;
    fin = fopen(argv[1],"rb");
    if(fin == nullptr){
        printf("Error: unable to open input file '%s'.\n",argv[1]);
        return -2;
    }
    
    //Read name of output IAEA file
    char* fout;
    fout = argv[2];
    
    
    
    //PSF write initialization: assign an unique Id to IAEA source (psfWrite)
    IAEA_I32 psfWrite;
    IAEA_I32 len = 81;
    //access = 2 => opening file for writing
    const IAEA_I32 accessWrite = 2;
    IAEA_I32 err;
    iaea_new_source(&psfWrite,fout,&accessWrite,&err,len);
    if(err<0){
        printf("Error: unable to open output IAEA PSF.\n");
        return -3;
    }
    else
    {
        printf("IAEA PSF writing initialized.\n");
        printf("\n");
    }
    
    
    //Set additional floats and integers 
    IAEA_I32 nextraFloat = 0; //These additional variables are from EGS
                              //simulation, set to 0 by penRed
    IAEA_I32 nextraInt   = 6; //These additional variables corresponds to
                              //incremental history number and ILBs, 
                              //set to corresponding values by penRed
                              
    iaea_set_extra_numbers(&psfWrite,&nextraFloat, &nextraInt);

    
    //Set a type type of the extra long variable corresponding to the "index" number
    IAEA_I32 extraIntTypes[nextraInt],type[nextraInt],extraInt[nextraInt];
    IAEA_Float extraFloat[1];
    
    type[0] = 1;  //Incremental history number
    //type[1] = 0;  //User defined generic type
    ///type[2] = 0;
    //type[3] = 0;
    //type[4] = 0;
   
    
    //Using IAEA routines that recognize ILB1 to ILB4 for PENELOPE codes
    type[1] = 7;  //ILB1
    type[2] = 6;  //ILB2
    type[3] = 5;  //ILB3
    type[4] = 4;  //ILB4
    

    type[5] = 3;  //ILB5

    
    for(IAEA_I32 index = 0; index<nextraInt; ++index)
    {
        extraIntTypes[index] = type[index];
        
        iaea_set_type_extralong_variable(&psfWrite,&index,&extraIntTypes[index]);
        if(extraIntTypes[index]==-1)
        {
           printf("Error: unable to set additional variable type.\n");
           return -5;
        }
        else if(extraIntTypes[index]==-2 || extraIntTypes[index]==-3)
        {
           printf("Error: additional variable type out of range.\n");
           return -6;
        }
    }
    
    
    
    //Create phase space file reader
    pen_psfreader psfReader;
    
    //Create particle state and kpar
    pen_particleState state;
    //pen_KPAR kpar;
    unsigned kpar;
    
    //Declare hist increment and total number of histories
    unsigned long dhist = 0;
    unsigned long nhist = 0;
    
    //Create chunk counter
    unsigned nchunks = 0;
    
    //Create particle counter
    unsigned long part[3];
    part[0]=0;
    part[1]=0;
    part[2]=0;
    
    printf("\n");
    printf("Converting PenRed PSF format into IAEA format.\n");
    printf("\n");
    
    //Read input file until the end
    while(psfReader.read(fin,1) == PEN_PSF_SUCCESS)
    {
        nchunks++;
        while(psfReader.get(dhist,kpar,state) > 0)
        {
            //Set IAEA particle state  
            IAEA_Float iaeaE, iaeaWght,iaeaX,iaeaY,iaeaZ,iaeaU,iaeaV,iaeaW;
            IAEA_I32 iaeaKpar, iaean_stat;
            
            //Convert kpar from penRed format to IAEA
            if(kpar == PEN_PHOTON){
                iaeaKpar=1;//photons
                part[0]++;
            }
            else if(kpar == PEN_ELECTRON){
                iaeaKpar=2; //electrons
                part[1]++;
            }
            else if(kpar == PEN_POSITRON){
                iaeaKpar=3; //positrons
                part[2]++;
            }
            
            //Set state data
            iaeaE = state.E*1.0E-06; //Energy must be converted from eV to MeV
            iaeaX = state.X;
            iaeaY = state.Y;
            iaeaZ = state.Z;
            iaeaU = state.U;
            iaeaV = state.V;
            iaeaW = state.W;
            iaeaWght = state.WGHT;
            
            extraInt[0] = dhist;
            extraInt[1] = state.ILB[0];
            extraInt[2] = state.ILB[1];
            extraInt[3] = state.ILB[2];
            extraInt[4] = state.ILB[3];
            extraInt[5] = state.ILB[4];
            
            extraFloat[0] = 0.0; //No float additional variables for PenRed simulations
            
            //History increment; 0 if same history and greater if else
            iaean_stat = dhist;
            
            nhist += dhist;
             
            //Write particle in IAEA PSF format
            iaea_write_particle(&psfWrite,&iaean_stat,&iaeaKpar,&iaeaE,&iaeaWght,&iaeaX,&iaeaY,&iaeaZ,&iaeaU,&iaeaV,&iaeaW,extraFloat,extraInt);
            
            if(iaean_stat==-1){
                printf("Error: unable to write particle into IAEA format.\n");
                return -7;
            }
        }
    }
        
    
    //Write header IAEA file
    IAEA_I64 iaeaOriginalHist = nhist;
    iaea_set_total_original_particles(&psfWrite,&iaeaOriginalHist);
    
    if(iaeaOriginalHist<0){
        printf("Error: unable to write total number of original histories to IAEA header file.\n");
        return -8;
    }
    
    //Update header IAEA file
    iaea_update_header(&psfWrite,&err);
    if(err<0)
    {
        printf("Error: unable to update IAEA header file.\n");
        return -9;   
    }
    
    
    //Close psfReader
    iaea_destroy_source(&psfWrite,&err);
    if(err<0)
    {
        printf("Error: unable to close IAEA PSF\n");
        return -10;
    }
    
 fclose(fin);
 
 
  long int totalPart =  part[0]+part[1]+part[2];
  printf("PenRed PSF read and converted to IAEA format.\n");
  printf("\n");   
  printf("Photons read from PSF: %lu\n",part[0]);
  printf("Electrons read from PSF: %lu\n",part[1]);
  printf("Positrons read from PSF: %lu\n",part[2]);
  printf("Total particles read from PSF: %lu\n",totalPart);
  printf("Total original histories declared in IAEA PSF header: %llu\n",iaeaOriginalHist);
  printf("Total number of chunks read from penRed format PSF: %u\n",nchunks);
 
  printf("\n");
  printf("Conversion terminated successfully.\n");  
 
}
  
