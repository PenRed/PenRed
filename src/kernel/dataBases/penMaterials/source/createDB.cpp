 
#include <fstream>


int main(int argc, char** argv){

  if(argc < 4){
    printf("usage: %s path/to/file/with/db/files/list "
	   "path/to/db/dir path/to/output/dir\n", argv[0]);
    return 1;
  }

  //Try to open file with the list of files in the DB
  std::ifstream fList(argv[1]);
  if(!fList.is_open()){
    printf("Error: Unable to open data base file list (file '%s')\n",
	   argv[1]);
    return -1;
  }

  //Open the header file to be generated
  std::string headerFilename(argv[3]);
  headerFilename.append("/database.hh");
  std::ofstream fheader(headerFilename);
  fheader << "#ifndef __PEN_MAT_DATABASE__" << std::endl;
  fheader << "#define __PEN_MAT_DATABASE__" << std::endl;

  fheader << "\nnamespace penred{\n";
  fheader << "   namespace penMatDB{" << std::endl;
    
  //Read and process all specified files
  std::string line;
  while(std::getline(fList, line)){
    //Try to open the specified file
    std::string pathToFile(argv[2]);
    pathToFile.append("/").append(line);
    std::ifstream fin(pathToFile);
    if(!fin.is_open()){
      printf("Unable to open input file '%s'\n", pathToFile.c_str());
    }

    //Open the output file
    std::string pathToOut(argv[3]);
    pathToOut.append("/").append(line);
    std::ofstream fout(pathToOut);
    if(!fout.is_open()){
      printf("Unable to open output file '%s'\n", pathToOut.c_str());
    }

    //Both files opened, process it
    fout << "R\"***(";

    std::string line2;
    while(std::getline(fin,line2)){
      fout << line2 << std::endl;
    }
    
    fout << ")***\"" << std::endl;

    fin.close();
    fout.close();

    //Create the variable name with no points
    std::string varName = line;
    std::string::size_type pointPos;
    pointPos = varName.find(".") ;
    while(pointPos != std::string::npos){
      varName.replace(pointPos,1,"_");
      pointPos = varName.find(".") ;
    }
    
    //Add the line to the header file
    fheader << "      const char* const " << varName << " = {\n";
    fheader << "                  #include \"" << line << "\"\n";
    fheader << "      };" << std::endl;
  }

  //Write function to convert filename to data pointer
  fList.close();
  fList.open(argv[1]);

  fheader << "// To covert to 'constexpr' when updating standard version" << std::endl;
  fheader << "\n      inline const char* readDataBaseFile(const std::string& filename){\n";
  while(std::getline(fList, line)){
    //Create the variable name with no points
    std::string varName = line;
    std::string::size_type pointPos;
    pointPos = varName.find(".") ;
    while(pointPos != std::string::npos){
      varName.replace(pointPos,1,"_");
      pointPos = varName.find(".") ;
    }
    
    fheader << "            if(filename.compare(\"" << line << "\") == 0){\n";
    fheader << "              return " << varName << ";\n";
    fheader << "            }" << std::endl;    
  }

  fheader << "            return nullptr;\n       }" << std::endl;

  fheader << "   };\n";
  fheader << "};\n";
  fheader << "#endif" << std::endl;
  fheader.close();

  fList.close();
  
  return 0;
}
