 
#include <fstream>
#include <string>

#include "dataBasesCommon.hh"

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

  fheader << "\n#include <array>" << std::endl;
  fheader << "\n#include <string>" << std::endl;

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
      return -2;
    }

    //Create string literal files
    unsigned nSubFiles;
    std::string pathToOut(argv[3]);
    pathToOut.append("/").append(line);
    int err = penred::dataBases::createStringLiteralFiles(fin, pathToOut, nSubFiles);
    if(err != penred::dataBases::errors::SUCCESS){
      printf("Unable to create string literal files for '%s'\n"
	     "%s\n",
	     pathToFile.c_str(), penred::dataBases::errorMessage(err));
      return -3;
    }

    //Close database input file
    fin.close();

    //Create the variable name with no points
    std::string varName = penred::dataBases::toVariableName(line);

    //Define a variable corresponding to a each subfile in the header file
    fheader << "// " << varName << std::endl;
    for (unsigned i = 0; i < nSubFiles; ++i) {
        fheader << "      extern const char* const " << varName << "_" << i << ";" << std::endl;
    }
    fheader << std::endl;
    //Define also the array containing all element data in the header file
    fheader << "      extern const std::array<const char* const, " << nSubFiles << "> " << varName << ";\n" << std::endl;

    //Create the source file for this database element
    std::string srcFileName(argv[3]);
    srcFileName += "/" + varName + ".cpp";
    std::ofstream fsrc(srcFileName);
    fsrc << "#include \"database.hh\"" << std::endl;
    for (unsigned i = 0; i < nSubFiles; ++i) {
        fsrc << "      const char* const penred::penMatDB::" << varName << "_" << i << " = {\n";
        fsrc << "                  #include \"" << line << "_" << i << "\"\n";
        fsrc << "      };" << std::endl;
    }
    fsrc << std::endl;

    //Create an array with all subfiles (src)
    fsrc << "      const std::array<const char* const, " << nSubFiles << "> penred::penMatDB::" << varName << " = {\n"
        << "                                                 penred::penMatDB::" << varName << "_0";
    for (unsigned i = 1; i < nSubFiles; ++i) {
        fsrc << ",\n"
            << "                                                 penred::penMatDB::" << varName << "_" << i;
    }
    fsrc << "\n                                                 };\n" << std::endl;
    fsrc.close();
  }

  //Write function to convert filename to data pointer
  fList.close();
  fList.open(argv[1]);

  fheader << "\n      inline const char* readDataBaseFile(const std::string& filename, size_t subFile){\n";
  while(std::getline(fList, line)){
    //Create the variable name with no points
    std::string varName = penred::dataBases::toVariableName(line);
    
    fheader << "            if(filename.compare(\"" << line << "\") == 0){\n";
    fheader << "              return " << varName << ".size() > subFile ? " << varName << "[subFile] : nullptr;\n";
    fheader << "            }" << std::endl;    
  }

  fheader << "            return nullptr;\n       }" << std::endl;

  fheader << "   }\n";
  fheader << "}\n";
  fheader << "#endif" << std::endl;
  fheader.close();

  fList.close();
  
  return 0;
}
