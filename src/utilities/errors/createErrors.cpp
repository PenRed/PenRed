
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


#include <cstdio>  //FILENAME_MAX
#include <ctype.h>

#if defined _MSC_VER
//  Microsoft Visual C++
#include <filesystem>   
namespace fs = std::filesystem;

#include <cstring>

int main()
{
	//Get actual dir
	std::error_code ec;
	const std::filesystem::path& dirpath = fs::current_path(ec);
	if (ec.value() != 0)
	{
		printf("Error: Could not get actual dir path.\n");
		return -2;
	}

	if ((fs::exists(dirpath)) == true)
	{

		//Open output files
		FILE* fenum = 0;
		fenum = fopen("eenum.h", "w");
		FILE* fmess = 0;
		fmess = fopen("emess.h", "w");

		if (fenum == 0 || fmess == 0)
		{
			printf("Error: Could not open output files 'eenum.h' and 'emess.h'.\n");
			return -3;
		}

		fprintf(fenum, "#ifndef __PEN_ERROR_ENUM__\n");
		fprintf(fenum, "#define __PEN_ERROR_ENUM__\n");

		fprintf(fenum, "enum pen_errCode{\n");
		fprintf(fenum, "      PEN_SUCCESS = 0");

		fprintf(fmess, "#ifndef __PEN_ERROR_MESS__\n");
		fprintf(fmess, "#define __PEN_ERROR_MESS__\n");

		fprintf(fmess, "const char* __pen_errorMessages[] = { \n");
		fprintf(fmess, "              \"Success\"");


		for (auto& strDir : fs::directory_iterator(dirpath))
		{

			//Skip hide files
			if (strDir.path().filename().string()[0] == '.') {
				printf("Skipping hide file: '%s'.\n", strDir.path().filename().string().c_str());
				continue;
			}
			//Avoid files with no .err extension
			if (strDir.path().extension().string() != ".err")
			{
				printf("Skipping non '.err' file: '%s'.\n", strDir.path().extension().string().c_str());
				continue;
			}

			printf("Adding errors list from '%s'.\n", strDir.path().filename().string().c_str());
			FILE* fin = 0;
			fin = fopen(strDir.path().filename().string().c_str(), "r");
			char line[30000];
			while (fgets(line, sizeof(line), fin))
			{
				//Check for \n character
				char* psn = strchr(line, '\n');
				if (psn == NULL)
				{
					//discard remaining line
					fscanf(fin, "%*[^\n]");
					getc(fin);
				}
				else
				{
					//Remove \n char
					*psn = '\0';
				}

				//get enum string
				char* pline = line;
				//skip whitespaces
				while (isspace(*pline)) pline++;
				char* penum = pline;
				//Search next whitespace
				while (!isspace(*pline)) pline++;
				//Print enum string
				*pline = '\0';
				fprintf(fenum, ",\n");
				fprintf(fenum, "      ERR_%s", penum);
				*pline = ' ';
				//skip whitespaces
				while (isspace(*pline)) pline++;
				//Print error descriptiom
				fprintf(fmess, ",\n");
				fprintf(fmess, "              \"%s\"", pline);
			}

		}

		// Add extra element to enumeration to count number of errors
		fprintf(fenum, ",\n");
		fprintf(fenum, "      ERR_LAST___");

		fprintf(fmess, ",\n");
		fprintf(fmess, "              \"Last error.\"");

		fprintf(fenum, "};\n");
		fprintf(fmess, "};\n");
		//closedir(dir);
		fprintf(fenum, "#endif\n");
		fprintf(fmess, "#endif\n");
		return 0;
	}
	else
	{
		printf("Error: Colud not open directory '%s'\n", dirpath.string().c_str());
		return -1;
	}
}

#else
#include <dirent.h>

#ifdef _WIN32
    #define getcwd _getcwd
#else
    #include <unistd.h>
#endif
#include <cstring>

int main()
{
  //Get actual dir
  char dirpath[FILENAME_MAX];
  if(getcwd(dirpath,FILENAME_MAX) == 0)
    {
      printf("Error: Could not get actual dir path.\n");
      return -2;
    }
  
  //Read all files in "dirpath"
  DIR* dir;
  struct dirent* strDir;
  if((dir = opendir(dirpath)) != NULL)
    {

      //Open output files
      FILE* fenum = 0;
      fenum = fopen("eenum.h","w");
      FILE* fmess = 0;
      fmess = fopen("emess.h","w");
      
      if(fenum == 0 || fmess == 0)
	{
	  printf("Error: Could not open output files 'eenum.h' and 'emess.h'.\n");
	  return -3;
	}

      fprintf(fenum,"#ifndef __PEN_ERROR_ENUM__\n");
      fprintf(fenum,"#define __PEN_ERROR_ENUM__\n");

      fprintf(fenum,"enum pen_errCode{\n");
      fprintf(fenum,"      PEN_SUCCESS = 0");
      
      fprintf(fmess,"#ifndef __PEN_ERROR_MESS__\n");
      fprintf(fmess,"#define __PEN_ERROR_MESS__\n");
      
      fprintf(fmess,"const char* __pen_errorMessages[] = { \n");
      fprintf(fmess,"              \"Success\"");

      
      while((strDir = readdir(dir)) != NULL)
	{

	  //Skip hide files
	  if(strDir->d_name[0] == '.'){
	    printf("Skipping hide file: '%s'.\n",strDir->d_name);
	    continue;
	  }
	  //Avoid files with no .err extension
	  int len = strlen(strDir->d_name);
	  if(strcmp(strDir->d_name+len-4,".err") != 0){
	    printf("Skipping non '.err' file: '%s'.\n",strDir->d_name);
	    continue;
	  }
	  
	  printf("Adding errors list from '%s'.\n",strDir->d_name);
	  FILE* fin = 0;
	  fin = fopen(strDir->d_name,"r");
	  char line[30000];
	  while(fgets(line,sizeof(line),fin))
	    {
	      //Check for \n character
	      char* psn = strchr(line, '\n');
	      if(psn == NULL)
		{
		  //discard remaining line
		  fscanf(fin,"%*[^\n]");
		  getc(fin);
		}
	      else
		{
		  //Remove \n char
		  *psn = '\0';
		}

	      //get enum string
	      char* pline = line;
	      //skip whitespaces
	      while(isspace(*pline)) pline++;
	      char* penum = pline;
	      //Search next whitespace
	      while(!isspace(*pline)) pline++;
	      //Print enum string
	      *pline = '\0';
	      fprintf(fenum,",\n");
	      fprintf(fenum,"      ERR_%s",penum);
	      *pline = ' ';
	      //skip whitespaces
	      while(isspace(*pline)) pline++;	      
	      //Print error descriptiom
	      fprintf(fmess,",\n");
	      fprintf(fmess,"              \"%s\"",pline);	      
	    }
	  
	}

      // Add extra element to enumeration to count number of errors
      fprintf(fenum,",\n");
      fprintf(fenum,"      ERR_LAST___");

      fprintf(fmess,",\n");
      fprintf(fmess,"              \"Last error.\"");	      
      
      fprintf(fenum,"};\n");
      fprintf(fmess,"};\n");
      closedir(dir);
      fprintf(fenum,"#endif\n");
      fprintf(fmess,"#endif\n");
      return 0;
    }
  else
    {
      printf("Error: Colud not open directory '%s'\n",dirpath);
      return -1;
    }
}
#endif

