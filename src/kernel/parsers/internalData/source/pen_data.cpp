
//
//
//    Copyright (C) 2019-2020 Universitat de València - UV
//    Copyright (C) 2019-2020 Universitat Politècnica de València - UPV
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


#include "pen_data.hh"

const char* pen_parserError(const int err){
  switch(err){
  case INTDATA_SUCCESS: return "Success";break;
  case INTDATA_INVALID_KEY: return "Invalid key";break;
  case INTDATA_NOT_A_SCALAR: return "Scalar expected";break;
  case INTDATA_NOT_A_ARRAY: return "Array expected";break;
  case INTDATA_NOT_A_STRING: return "String expected";break;
  case INTDATA_NOT_A_SECTION: return "Section expected";break;
  case INTDATA_KEY_IS_NOT_PATH: return "Key is not a path";break;
  case INTDATA_NOT_A_CHAR: return "Character expected";break;
  case INTDATA_NOT_A_INT: return "Integer expected";break;
  case INTDATA_NOT_A_DOUBLE: return "Double expected";break;
  case INTDATA_NOT_A_BOOL: return "Boolean expected";break;
  case INTDATA_OUT_OF_RANGE: return "Out of range";break;
  case INTDATA_BAD_ALLOCATION: return "Bad allocation";break;
  case INTDATA_USED_KEY: return "Used key";break;
  case INTDATA_UNUSED_KEY: return "Unused key";break;
  case INTDATA_EMPTY_KEY: return "Empty key";break;
  case INTDATA_NOT_A_ELEMENT: return "Is not an element";break;
  case INTDATA_INVALID_PREFIX: return "Invalid prefix";break;
  case INTDATA_KEY_NO_EXIST: return "Key doesn't exists";break;
  case INTDATA_PARSE_STRING_EMPTY: return "Empty string";break;
  case INTDATA_PARSE_STRING_PARTIAL_UNUSED:
    return "Parsed string partially unused";break;
  case INTDATA_PARSE_STRING_INVALID_CHARACTER:
    return "Invalid parsed character";break;
  case INTDATA_PARSE_STRING_INVALID_DATA: return "Invalid parsed data";break;
  case INTDATA_PARSE_STRING_INVALID_ARRAY: return "Invalid parsed array";break;
  case INTDATA_PARSE_STRING_INVALID_STRING: return "Invalid parsed string";break;
  case INTDATA_PARSE_STRING_INVALID_ELEMENT: return "Invalid parsed element";break;
  case INTDATA_PARSE_STRING_INVALID_VALUE: return "Invalid parsed value";break;
  case INTDATA_NULL_STRING: return "Null string";break;
  case INTDATA_OPEN_FILE: return "Null file pointer";break;
  case INTDATA_UNKNOWN_ERROR: return "Unknown error";break;
  default: return "Non tabulated error";
  }
}


pen_parserData::pen_parserData(){
  tag = CHAR;
  c = '\0';
}
pen_parserData::pen_parserData(const char val){
  c = val;
  tag = CHAR;
}
pen_parserData::pen_parserData(const int val){
  i = val;
  tag = INT;
}
pen_parserData::pen_parserData(const double val){
  d = val;
  tag = DOUBLE;
}
pen_parserData::pen_parserData(const bool val){
  b = val;
  tag = BOOL;
}
pen_parserData::pen_parserData(const pen_parserData& val){
  tag = val.readTag();
  switch(tag){
  case CHAR: c = val.c; break;
  case INT: i = val.i; break;
  case DOUBLE: d = val.d; break;
  case BOOL: b = val.b; break;
  }

}

pen_parserData& pen_parserData::operator=(const char val){
  c = val;
  tag = CHAR;  
  return *this;
}
pen_parserData& pen_parserData::operator=(const int val){
  i = val;
  tag = INT;
  return *this;
}
pen_parserData& pen_parserData::operator=(const double val){
  d = val;
  tag = DOUBLE;  
  return *this;
}
pen_parserData& pen_parserData::operator=(const bool val){
  b = val;
  tag = BOOL;  
  return *this;
}

int pen_parserData::read(char& vget) const{
  if(tag == CHAR){
    vget = c;
    return INTDATA_SUCCESS;
  }
  return INTDATA_NOT_A_CHAR;
}
int pen_parserData::read(int& vget) const{
  if(tag == INT){
    vget = i;
    return INTDATA_SUCCESS;
  }
  return INTDATA_NOT_A_INT;  
}
int pen_parserData::read(double& vget) const{
  if(tag == DOUBLE){
    vget = d;
    return INTDATA_SUCCESS;
  }
  else if(tag == INT){
    vget = double(i);
    return INTDATA_SUCCESS;
  }
  return INTDATA_NOT_A_DOUBLE;  
}
int pen_parserData::read(bool& vget) const{
  if(tag == BOOL){
    vget = b;
    return INTDATA_SUCCESS;
  }
  return INTDATA_NOT_A_BOOL;
}
int pen_parserData::read(pen_parserData& vget) const{
    vget = *this;
    return INTDATA_SUCCESS;
}

//Stringify function
void pen_parserData::stringify(std::string& strout) const{
  char buffer[20];
  switch(tag){
  case CHAR: {
    buffer[0] = '\'';
    buffer[1] = c;
    buffer[2] = '\'';
    buffer[3] = '\0';
    break;
  }
  case INT: snprintf(buffer,sizeof(buffer),"%d",i); break;
  case DOUBLE: snprintf(buffer,sizeof(buffer),"%14.5E",d); break;
  case BOOL: snprintf(buffer,sizeof(buffer),"%s", b ? "true" : "false"); break;
  }
  strout.assign(buffer);
}

//Parse function
int pen_parserData::parse(const std::string& strIn){

  using std::stoi;
  using std::stod;
  using std::isspace;
    
  //Clear string
  std::string strClear = trim(strIn);
  
  //Check string size
  size_t strlength = strClear.length();
  if(strlength == 0)
    return INTDATA_PARSE_STRING_EMPTY;

  //Check if is a character
  if(strClear[0] == '\''){
    if(strlength != 3) //Chars have following format 'c' or
      return INTDATA_PARSE_STRING_INVALID_CHARACTER;
      
    if(strClear[2] == '\''){
      *this = strClear[1];
      return INTDATA_SUCCESS;
    }
    else{
      return INTDATA_PARSE_STRING_INVALID_CHARACTER;
    }
  }
  //Check if is a bool
  if(strClear.compare("True") == 0 || strClear.compare("true") == 0){
    *this = bool(true);
    return INTDATA_SUCCESS;
  }
  if(strClear.compare("False") == 0 || strClear.compare("false") == 0){
    *this = bool(false);
    return INTDATA_SUCCESS;
  }

  //Check if is a number
  bool isInt = true;
  bool isDoub = true;

  char* end_l;
  char* end_d;

  int intVar = 0;
  long int longVar = strtol(strClear.c_str(), &end_l, 0);
  size_t long_str_rep_length = strlength - strlen(end_l);

  if (longVar == 0)
  {
      if (long_str_rep_length == 0)
      {
          //Is not an integer
          isInt = false;
      }else{ //Is an integer and is zero
          intVar = 0;
      }
  }
  else if (longVar == LONG_MIN || longVar == LONG_MAX)
  {
      //Is not an integer
      isInt = false;
  }
  else if(longVar > INT_MAX || longVar < INT_MIN){
      isInt = false; //Can't fit in an integer
  }
  else
  {   //Is an integer
      intVar = static_cast<int>(longVar);
  }

  double doubVar = strtod(strClear.c_str(), &end_d);
  size_t double_str_rep_length = strlength - strlen(end_d);

  if (doubVar == 0.0)
  {
      if (double_str_rep_length == 0)
      {
          //Is not a double
          isDoub = false;
      } 
      // else{ Is a double and is zero} 
  }
  else if (doubVar == HUGE_VAL)
  {
      //Is not a double
      isDoub = false;
  }
  // else{  Is a double }

  if (!isDoub && !isInt) {
      //Is a string
      return INTDATA_PARSE_STRING_INVALID_DATA;
  }

  //Is a integer or a double, check read positions
  if (long_str_rep_length == double_str_rep_length)
  {
      //Number not contain '.', 'e' etc, is a integer
      *this = intVar;
      if (long_str_rep_length < strlength)
          return INTDATA_PARSE_STRING_PARTIAL_UNUSED;
  }
  else {
      //Is a double
      *this = doubVar;
      if (double_str_rep_length < strlength)
          return INTDATA_PARSE_STRING_PARTIAL_UNUSED;
  }
  return INTDATA_SUCCESS;
}

// pen_parserArray
///////////////////////////

pen_parserArray::pen_parserArray(){}
pen_parserArray::pen_parserArray(const pen_parserArray& c) : vect(c.vect){}

//Append element
void pen_parserArray::append(const pen_parserData val){
  vect.push_back(val);
}
void pen_parserArray::append(const char val){
  const pen_parserData data = val;
  vect.push_back(data);
}
void pen_parserArray::append(const int val){
  const pen_parserData data = val;
  vect.push_back(data);
}
void pen_parserArray::append(const double val){
  const pen_parserData data = val;
  vect.push_back(data);
}
void pen_parserArray::append(const bool val){
  const pen_parserData data = val;
  vect.push_back(data);
}

//Remove elements
int pen_parserArray::remove(const unsigned index){

  if(index >= vect.size())
    return INTDATA_OUT_OF_RANGE;
  
  vect.erase(vect.begin()+index);
  return INTDATA_SUCCESS;
}
int pen_parserArray::remove(const unsigned begin, const unsigned end){

  const unsigned vecSize = vect.size();

  if(begin >= vecSize)
    return INTDATA_OUT_OF_RANGE;
  
  unsigned limit = end;
  if(limit > vecSize)
    limit = vecSize;
  
  vect.erase(vect.begin()+begin, vect.begin()+limit);
  return INTDATA_SUCCESS;
}

//Access elements
pen_parserData& pen_parserArray::operator[](const unsigned index){
  unsigned length = vect.size();
  if(index >= length){
    char error[400];
    sprintf(error,"pen_parserArray: index %d out of range (%d)",index, length);
    throw std::out_of_range (error);
  }
  return vect[index];
}

//Assignment operator
pen_parserArray& pen_parserArray::operator=(const pen_parserArray& c){
  vect = c.vect;
  return *this;
}

//Assignment functions
int pen_parserArray::assign(const pen_parserElement& c){
  if(c.readTag() == STRING){
    return INTDATA_NOT_A_ARRAY;
  }
  assign(c.array);
  return INTDATA_SUCCESS;
}

//Stringify function
void pen_parserArray::stringify(std::string& strout) const{
  strout.assign("[");
  std::string straux;

  for(unsigned i = 0; i < size(); i++){
    vect[i].stringify(straux);
    if(i > 0)
      strout += ',';
    strout += straux;
  }
    
  strout += ']';
}

//Parse function
int pen_parserArray::parse(const std::string& strIn){

  //Clear whitecharacters
  std::string strClear = trim(strIn);
  
  //Check size
  size_t strlength = strClear.length();
  if(strlength == 0)
    return INTDATA_PARSE_STRING_EMPTY;

  //Clear array
  clear();
  
  //Check if first character is '[' and last ']'
  if(strClear[0] != '[' || strClear[strlength-1] != ']')
    return INTDATA_PARSE_STRING_INVALID_ARRAY;

  //Split elements separated with non character ','
  size_t posLastCom = 0;
  size_t firstPos = 1;
  size_t posCom = 0;
  while(posCom != std::string::npos){

    //Get position of next ','
    posCom = strClear.find(',',posLastCom+1);
    if(posCom != std::string::npos){
      posLastCom = posCom;
      //Check if this ',' is a data character
      if( strClear[posCom-1] == '\'' && strClear[posCom+1] == '\'' ){
	//Is not a element separator
	posLastCom = posCom;
	continue;
      }
    }
    else{
      //There are not more ','
      posLastCom = strlength-1; //Avoid last character ']'
    }

    //Is a element separator
    //Parse data string
    pen_parserData auxData;
    int err = auxData.parse(strClear.substr(firstPos,posLastCom-firstPos));
    if(err != INTDATA_SUCCESS){
      return INTDATA_PARSE_STRING_INVALID_ARRAY;
    }
    //Append data to array
    append(auxData);
    //Store new first position
    firstPos = posLastCom+1;    
  }
  return INTDATA_SUCCESS;
}

// pen_parserElement
//////////////////////////

//Constructors
pen_parserElement::pen_parserElement(){
  tag = SCALAR;
  array.append('\0');
}

pen_parserElement::pen_parserElement(const pen_parserElement& c) : tag(c.tag),
								   array(c.array),
								   str(c.str)
{}

pen_parserElement::pen_parserElement(const pen_parserData& c) : tag(SCALAR){
  array.append(c);
}

pen_parserElement::pen_parserElement(const pen_parserArray& c) : tag(ARRAY),
								 array(c)
{}

pen_parserElement::pen_parserElement(const std::string& c) : tag(STRING),
							     str(c)
{}

pen_parserElement::pen_parserElement(const char* c) : tag(STRING),
							     str(c)
{
}

pen_parserElement::pen_parserElement(const char c) : tag(SCALAR){
  array.append(c);
}
pen_parserElement::pen_parserElement(const int c) : tag(SCALAR){
  array.append(c);
}
pen_parserElement::pen_parserElement(const double c) : tag(SCALAR){
  array.append(c);
}
pen_parserElement::pen_parserElement(const bool c) : tag(SCALAR){
  array.append(c);
}

//Remove data
int pen_parserElement::remove(const unsigned index){
  if(tag != ARRAY)
    return INTDATA_NOT_A_ARRAY;

  return array.remove(index);
}
int pen_parserElement::remove(const unsigned begin, const unsigned end){
  if(tag != ARRAY)
    return INTDATA_NOT_A_ARRAY;

  return array.remove(begin,end);
}

//Assignment operators
pen_parserElement& pen_parserElement::operator=(const pen_parserElement& c){
  tag = c.tag;
  array = c.array;
  str.assign(c.str);
  return *this;
}
pen_parserElement& pen_parserElement::operator=(const pen_parserArray& c){
  if(tag != ARRAY){
    clear();
    tag = ARRAY;
  }
  array = c;
  return *this;  
}
pen_parserElement& pen_parserElement::operator=(const pen_parserData& c){
  if(tag != SCALAR){
    clear();
    tag = SCALAR;
    array.append(c);
  }
  else{
    array[0] = c;
  }
  return *this;
}
pen_parserElement& pen_parserElement::operator=(const std::string& c){
  if(tag != STRING){
    clear();
    tag = STRING;
  }
  str.assign(c);
  return *this;
}
pen_parserElement& pen_parserElement::operator=(const char* c){
  if(tag != STRING){
    clear();
    tag = STRING;
  }
  str.assign(c);
  return *this;
}
pen_parserElement& pen_parserElement::operator=(const char c){
  *this = pen_parserData(c);
  return *this;
}
pen_parserElement& pen_parserElement::operator=(const int c){
  *this = pen_parserData(c);
  return *this;
}
pen_parserElement& pen_parserElement::operator=(const double c){
  *this = pen_parserData(c);
  return *this;
}
pen_parserElement& pen_parserElement::operator=(const bool c){
  *this = pen_parserData(c);
  return *this;
}

//Access elements
pen_parserData& pen_parserElement::operator[](const unsigned index){
  if(tag == ARRAY || tag == SCALAR){
    return array[index];
  }
  char error[400];
  sprintf(error,"pen_parserElement: invalid operator [] access to string type element.");
  throw std::out_of_range (error);
  
}

//Stringify function
void pen_parserElement::stringify(std::string& strout) const{
  if(tag == SCALAR){
    array.at(0).stringify(strout);
  }
  else if(tag == ARRAY){
    array.stringify(strout);
  }
  else if(tag == STRING){
    strout.assign(str);
  }
}

//Parse function
int pen_parserElement::parse(const std::string& strIn){

  //Clear element
  clear();

  //Clear string
  std::string strClear = trim(strIn);
  size_t strlength = strClear.length();
  
  //Check size
  if(strlength == 0)
    return INTDATA_PARSE_STRING_EMPTY;
  
  //Try to parse element as a array
  int errArray = array.parse(strClear);
  if(errArray == INTDATA_SUCCESS){
    tag = ARRAY;
    return INTDATA_SUCCESS;
  }
  
  //Clear array
  array.clear();
  tag = SCALAR;
  array.append('\0');
  
  //Try to parse element as data
  pen_parserData auxData;
  int errScalar = auxData.parse(strClear);
  if(errScalar == INTDATA_SUCCESS){
    *this = auxData;
    return INTDATA_SUCCESS;
  }

  //Check if is a string
  if(strClear[0] == '"'){
    if(strClear[strlength-1] != '"'){
      return INTDATA_PARSE_STRING_INVALID_STRING;
    }
    *this = strClear.substr(1,strlength-2); //Avoid limiting '"'
    return INTDATA_SUCCESS;
  }
  return INTDATA_PARSE_STRING_INVALID_ELEMENT;
}

// pen_parserSection
//////////////////////////

//Constructors
pen_parserSection::pen_parserSection(){
}

pen_parserSection::pen_parserSection(const pen_parserSection& c){
  //Copy elements from "c"
  elements = c.elements;
}

//Key format check
int pen_parserSection::checkKey(const char* key) const{
  std::string aux(key);
  return checkKey(aux);
}
int pen_parserSection::checkKey(std::string& key) const{

  //Clear key white-space characters
  key = trim(key);
  
  //Check key size
  size_t len = key.size();
  if(len == 0)
    return INTDATA_EMPTY_KEY;

  //Check if "key" has a '/' as last char
  char last = key.back();
  size_t limit = len;
  if(last == '/'){
    limit--; //Avoid last '/'
  }

  //Check if first character is a '/'
  if(key[0] == '/'){
    key = key.substr(1,limit-1);
  }
  else{
    key = key.substr(0,limit);
  }

  //Check key size after "substr" call
  len = key.size();
  if(len == 0)
    return INTDATA_EMPTY_KEY;

  //Check new extremes
  if(key[0] == '/' || key[limit-1] == '/' || isspace(key[0]) || isspace(key[limit-1])){
    return INTDATA_INVALID_KEY;
  }
  
  //Check if before and after each '/' there is, at last, one character
  size_t posb = 0; // position 0 can't be a '/' (checked before)
  do{
    posb = key.find('/',posb+1);
    
    if(posb != std::string::npos){
      //Is not the last name

      if(key[posb-1] == '/' || key[posb+1] == '/' || isspace(key[posb-1]) || isspace(key[posb+1])){
	return INTDATA_INVALID_KEY;
      }
      //This subpath can't be a element , check it
      std::string name = key.substr(0,posb);
      if(isElement(name))
	return INTDATA_NOT_A_SECTION;      
    }
  }while(posb != std::string::npos);

  return INTDATA_SUCCESS;
}

//Find key functions
bool pen_parserSection::isSection(const std::string& key) const{

  //First check if is an element
  if(isElement(key)){
    return false;
  }
  
  //Get iterator to lower and upper bound of key
  elementMap::const_iterator itlow, itup;

  itlow = elements.lower_bound(key);
  if(itlow == elements.end())
    return false;

  //Is not an element, so, the introduced key and key pointed by
  //"itlow" must be different. The two possibilities are:
  //
  //  -Introduced key is a substring of itlow pointed key. Return true
  //  -Introduced key is completly different of itlow pointed key. Return false
  //

  //So, find key in iterator key
  if(itlow->first.find(key) == 0){
    //Coincidence at first position, is a section
    return true;
  }

  //Is not a section
  return false;
}
bool pen_parserSection::isElement(const std::string& key) const{
  if(elements.count(key) > 0)
    return true;
  return false;
}
int pen_parserSection::trusted_sectionRange(const std::string& key, elementMap::const_iterator& itlow, elementMap::const_iterator& itup) const{

  //Check if is a element
  if(isElement(key)){
    itlow = elements.end();
    itup = elements.end();
    return INTDATA_NOT_A_SECTION;
  }
  
  itlow = elements.lower_bound(key);  
  if(itlow == elements.end()){
    itup = elements.end();
    return INTDATA_NOT_A_SECTION;
  }
  
  //find key as substring in "itlow" key
  if(itlow->first.find(key) != 0){
    //key is not contained in iterator pointed key.
    //This section doesn't exists
    itlow = elements.end();
    itup = elements.end();
    return INTDATA_KEY_NO_EXIST;
  }

  itup = itlow;
  itup++;
  while(itup != elements.end()){
    if(itup->first.find(key) != 0){
      break;
    }
    itup++;
  }
  return INTDATA_SUCCESS;
}


//Create array
int pen_parserSection::createArray(std::string& key, const bool overwrite){
  pen_parserArray dummy;
  return set(key,dummy,overwrite);
}

//Remove elements
int pen_parserSection::trusted_remove(const std::string& key){
  elements.erase(key);

  //If is a section remove all inner elements and sections
  elementMap::const_iterator itlow, itup;
  int err = trusted_sectionRange(key,itlow,itup);
  if(err != INTDATA_SUCCESS){
    return err;
  }
  elements.erase(itlow,itup);
  return INTDATA_SUCCESS;
}
int pen_parserSection::remove(const char* key){
  std::string aux(key);
  return remove(aux);
}
int pen_parserSection::remove(std::string& key){

  //Check key
  int err = checkKey(key);
  if(err != INTDATA_SUCCESS)
    return err;

  return trusted_remove(key);
  
}
int pen_parserSection::remove(std::string& arraykey, const unsigned index){

  //Check key
  int err = checkKey(arraykey);
  if(err != INTDATA_SUCCESS)
    return err;
  
  //Check if exists
  if(!isElement(arraykey)){
    return INTDATA_NOT_A_ELEMENT;
  }

  //Try to remove specified index
  return elements[arraykey].remove(index);
}
int pen_parserSection::remove(std::string& arraykey, const unsigned begin, const unsigned end){

  //Check key
  int err = checkKey(arraykey);
  if(err != INTDATA_SUCCESS)
    return err;  

  //Check if exists
  if(!isElement(arraykey)){
    return INTDATA_NOT_A_ELEMENT;
  }

  //Try to remove specified range
  return elements[arraykey].remove(begin,end);
}

//Assign from input section
void pen_parserSection::assign(const pen_parserSection& c){
  elements = c.elements;
}

//Read subsection
int pen_parserSection::readSubsection(std::string& key, pen_parserSection& secOut, const bool removeKey) const{
  
  //Check key
  int err = checkKey(key);
  if(err != INTDATA_SUCCESS){
    return err;
  }  

  //Get section iterator begin and end
  elementMap::const_iterator itlow, itup;
  err = trusted_sectionRange(key,itlow,itup);
  if(err != INTDATA_SUCCESS){
    return err;
  }

  //Clear output section
  secOut.clear();

  if(removeKey){
    //Insert elements removing subsection key
    size_t keyLen = key.length()+1; //+1 to avoid '/'
    elementMap::const_iterator it;
    for(it = itlow; it != itup; it++){
      std::string aux(it->first,keyLen);
      secOut.elements[aux] = it->second;
    }
  }
  else{
    //Insert iterator range
    secOut.elements.insert(itlow,itup);
  }
  return INTDATA_SUCCESS;
}

//Get section

//Clear function
void pen_parserSection::clear(){
  elements.clear();
}

//Stringify function
void pen_parserSection::stringify(std::string& strout) const{

  //Prepare strings
  strout.clear();
  std::string aux;
  //Iterate all elements
  elementMap::const_iterator it;
  for(it = elements.begin(); it != elements.end(); it++){
    it->second.stringify(aux);
    strout += it->first + ": " + aux + '\n';
  }
}

//Parse function
int pen_parserSection::parse(const std::string& strIn, const bool overwrite){

  //Clear string
  std::string strClear = trim(strIn);
  size_t strlength = strClear.length();
  
  //Check size
  if(strlength == 0)
    return INTDATA_PARSE_STRING_EMPTY;

  //Get value substring
  size_t valueBeg = strlength-1;
  if(strClear[strlength-1] == '"'){
    //Value is a string, get beginning
    while(valueBeg != std::string::npos){

      //Get previous '"'
      valueBeg = strClear.rfind('"',valueBeg-1);

      if(valueBeg != std::string::npos){
	//Check if is a slashed "
	if(strClear[valueBeg-1] != '\\'){
	  //Is the string beginning
	  break;
	}
      }
    }
  }
  else if(strClear[strlength-1] == ']'){
    //Value is an array
    while(valueBeg != std::string::npos){

      //Get previous '['
      valueBeg = strClear.rfind('[',valueBeg-1);

      if(valueBeg != std::string::npos){
	//Check if is a character stored in the array
	if(strClear[valueBeg-1] != '\'' && strClear[valueBeg+1] != '\''){
	  //Is the array beginning
	  break;
	}
      }
    }
  }
  else{
    //Value is regular data, search for space or tabulation
    size_t spacepos = strClear.rfind(' ');
    size_t tabpos = strClear.rfind('\t');
    if(tabpos == std::string::npos){
      valueBeg = spacepos;
    }
    else if(spacepos == std::string::npos){
      valueBeg = tabpos;
    }
    else if(tabpos > spacepos){
      valueBeg = tabpos;
    }
    else{
      valueBeg = spacepos;
    }
  }

  //Check retrieved beginning position
  if(valueBeg == std::string::npos){
    return INTDATA_PARSE_STRING_INVALID_VALUE;
  }

  //Get key substring
  std::string key = strClear.substr(0,valueBeg);
  std::string value = strClear.substr(valueBeg);

  //Parse ant set element
  return parse(key,value,overwrite);
}

int pen_parserSection::parse(std::string& key, std::string& value, const bool overwrite){

  //Pase element value
  pen_parserElement auxElement;
  int err = auxElement.parse(value);
  if(err != INTDATA_SUCCESS)
    return err;

  //Try to set element
  return set(key,auxElement,overwrite);
}

//Parse function


//Assigment operator
pen_parserSection& pen_parserSection::operator=(const pen_parserSection& c){

  //Copy maps
  elements = c.elements;
  return *this;
}

//ls function
int pen_parserSection::ls(std::string& key, std::vector<std::string>& vect) const{

  //Check if key is empty
  if(key.length() == 0){
    //Take all names at root section
    std::string previous("//");
    elementMap::const_iterator it;
    for(it = elements.begin(); it != elements.end(); it++){
      //Find next '/'
      size_t barp = it->first.find('/');
      std::string aux(it->first,0,barp);

      //Check if this name has already added to array
      if(aux.compare(previous) != 0){
	//Is a new name
	vect.push_back(aux);
	previous.assign(aux);
      }    
    }
    return 0;
  }
  
  //Check key
  int err = checkKey(key);
  if(err != INTDATA_SUCCESS){
    return err;
  }

  //Check if is a single element
  if(isElement(key)){
    return INTDATA_NOT_A_SECTION;
  }

  //Get section iterator begin and end
  elementMap::const_iterator itlow, itup;
  err = trusted_sectionRange(key,itlow,itup);
  if(err != INTDATA_SUCCESS){
    return err;
  }

  vect.clear();
  
  //Take names contained in this section
  std::string previous("//");
  size_t keyLen = key.length()+1; //+1 to avoid '/'
  elementMap::const_iterator it;
  for(it = itlow; it != itup; it++){
    //Find next '/'
    size_t barp = it->first.find('/',keyLen);
    std::string aux(it->first,keyLen,barp-keyLen);

    //Check if this name has already added to array
    if(aux.compare(previous) != 0){
      //Is a new name
      vect.push_back(aux);
      previous.assign(aux);
    }
    
  }
  return 0;
}

//Destructor
pen_parserSection::~pen_parserSection(){
  clear();
}


//Auxiliar functions
int parseFile(const char* filename,
	      pen_parserSection& section,
	      std::string& errorString,
	      long unsigned& errorLine){

  if(filename == nullptr){
    return INTDATA_NULL_STRING;
  }

  FILE* fin = nullptr;
  fin = fopen(filename,"r");

  if(fin == nullptr){
    return INTDATA_OPEN_FILE;
  }

  //Clear output section
  section.clear();

  //Read and parse file lines
  char line[50000];
  long unsigned lineNum = 0;
  unsigned long read;
  while(pen_getLine(fin,50000,line,read) == 0){
    lineNum += read;
    int err = section.parse(line);
    if(err != INTDATA_SUCCESS)
      if(err != INTDATA_PARSE_STRING_EMPTY){
	errorLine = lineNum;
	errorString.assign(line);
	return err;
      }
  }

  return INTDATA_SUCCESS;
}

int pen_getLine(FILE* fin,
		const unsigned size,
		char* line,
		unsigned long& nlines){

  //'line' will be filled with the next non comment line.
  //The number of read lines will be stored at 'nlines'
  //The function will return 0 on success or a
  //negative value on error.

  nlines = 0;
  
  //Check if input file is a null pointer
  if(fin == nullptr || line == nullptr){
    return -1;
  }

  while(fgets(line,size,fin) != nullptr){
    nlines++;
    //Get pointer to \n if exists
    char* pn = strchr(line,'\n');

    if(pn == nullptr){
      //Read remaining line
      char dummy[500];
      while(pn == nullptr){
        if(fgets(dummy,500,fin) == nullptr)
            break;
        pn = strchr(dummy,'\n');
      }
    }

    //Check if is a comment
    char firstc = ' ';
    sscanf(line," %c",&firstc);

    if(firstc != '#'){
      return 0;
    }
  }

  return -2;
}

