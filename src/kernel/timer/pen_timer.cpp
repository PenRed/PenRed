
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


#include "pen_timer.hh"

//  *********************************************************************
//                       SUBROUTINE TIMER
//  *********************************************************************
double pen_timer::timer()
{
  //  This subroutine returns the execution time in seconds
  //  elapsed since the last call to TIME0. The output value.

  time_t time_end = time(NULL);

  return difftime(time_end,time_beginning);
}

//  *********************************************************************
//                       FUNCTION CPUTIM
//  *********************************************************************
double CPUtime()
{
  //  This function returns the CPU (user) time used since the start of the
  //  calling program, in seconds.

  return(((double)clock())/CLOCKS_PER_SEC);
}
//  *********************************************************************
//                       SUBROUTINE PDATET
//  *********************************************************************
void PDATET(char* DATE23)
{
  //  Delivers a 23-character string with the date and time.
  //  (ddth Mon yyyy. hh:mm:ss)
  
  char LORD[3];
  char MONTH[12][4];
//  char zone[5];
//  char date[8];
//  char c_time[10];
//  int values(8),IORD;
//  int  IORD;

  tm* values = 0;

  strcpy(MONTH[0] ,"Jan");
  strcpy(MONTH[1] ,"Feb");
  strcpy(MONTH[2] ,"Mar");
  strcpy(MONTH[3] ,"Apr");
  strcpy(MONTH[4] ,"May");
  strcpy(MONTH[5] ,"Jun");
  strcpy(MONTH[6] ,"Jul");
  strcpy(MONTH[7] ,"Aug");
  strcpy(MONTH[8] ,"Sep");
  strcpy(MONTH[9] ,"Oct");
  strcpy(MONTH[10],"Nov");
  strcpy(MONTH[11],"Dec");
  
  time_t actualTime = time(NULL);
  values = localtime(&actualTime);  

  if(values->tm_mday == 0){
    strcpy(LORD,"st");
  }
  else if(values->tm_mday == 1){
    strcpy(LORD,"nd");
  }
  else if(values->tm_mday == 2){
    strcpy(LORD,"rd");
  }
  else{
    strcpy(LORD,"th");
  }

  //  (ddth Mon yyyy. hh:mm:ss)  
  sprintf(DATE23,"%2d%2s %3s%5d. %02d:%02d:%02d\n",values->tm_mday,LORD,MONTH[values->tm_mon],1900+values->tm_year,values->tm_hour,values->tm_min,values->tm_sec);
  
  return;
}
