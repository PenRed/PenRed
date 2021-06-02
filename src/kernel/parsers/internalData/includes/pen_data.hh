
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

 
#ifndef __PENRED_PARSER_CLASSES__
#define __PENRED_PARSER_CLASSES__

#include <map>
#include <cstring>
#include <vector>
#include <exception>
#include <stdexcept>
#include <string>
#include <cctype>
#include <climits>
#include <cmath>

enum pen_parserErrors{INTDATA_SUCCESS = 0,
		      INTDATA_FALSE = 0,
		      INTDATA_TRUE = 1,
		      INTDATA_INVALID_KEY = -1,
		      INTDATA_NOT_A_SCALAR = -2,
		      INTDATA_NOT_A_ARRAY = -3,
		      INTDATA_NOT_A_STRING = -4,
		      INTDATA_NOT_A_SECTION = -5,
		      INTDATA_KEY_IS_NOT_PATH = -6,
		      INTDATA_NOT_A_CHAR = -7,
		      INTDATA_NOT_A_INT = -8,
		      INTDATA_NOT_A_DOUBLE = -9,
		      INTDATA_NOT_A_BOOL = -10,
		      INTDATA_OUT_OF_RANGE = -11,
		      INTDATA_BAD_ALLOCATION = -12,
		      INTDATA_USED_KEY = -13,
		      INTDATA_UNUSED_KEY = -14,
		      INTDATA_EMPTY_KEY = -15,
		      INTDATA_NOT_A_ELEMENT = -16,
		      INTDATA_INVALID_PREFIX = -17,
		      INTDATA_KEY_NO_EXIST = -18,
		      INTDATA_PARSE_STRING_EMPTY = -19,
		      INTDATA_PARSE_STRING_PARTIAL_UNUSED = -20,
		      INTDATA_PARSE_STRING_INVALID_CHARACTER = -21,
		      INTDATA_PARSE_STRING_INVALID_DATA = -22,
		      INTDATA_PARSE_STRING_INVALID_ARRAY = -23,
		      INTDATA_PARSE_STRING_INVALID_STRING = -24,
		      INTDATA_PARSE_STRING_INVALID_ELEMENT = -25,
		      INTDATA_PARSE_STRING_INVALID_VALUE = -26,
		      INTDATA_NULL_STRING = -27,
		      INTDATA_OPEN_FILE = -28,
		      INTDATA_UNKNOWN_ERROR = -99
};

const char* pen_parserError(const int err);

enum pen_parserDataTypes {CHAR, INT, DOUBLE, BOOL};
enum pen_parserElementTypes {SCALAR,ARRAY,STRING};

struct pen_parserData;
struct pen_parserArray;
struct pen_parserElement;
struct pen_parserSection;

struct pen_parserData{
private:
  enum pen_parserDataTypes tag;
  union{
    char c;
    int i;
    double d;
    bool b;
  };
public:

  pen_parserData();
  pen_parserData(const char val);
  pen_parserData(const int val);
  pen_parserData(const double val);
  pen_parserData(const bool val);
  pen_parserData(const pen_parserData& val);
  
  pen_parserData& operator=(const char);
  pen_parserData& operator=(const int);
  pen_parserData& operator=(const double);
  pen_parserData& operator=(const bool);

  int read(char&) const;
  int read(int&) const;
  int read(double&) const;
  int read(bool&) const;
  int read(pen_parserData&) const;

  //Stringify function
  inline std::string stringify() const{
    std::string aux;
    stringify(aux);
    return aux;
  }
  void stringify(std::string&) const;

  //Parse function
  inline int parse(const char* data){
    std::string aux(data);
    return parse(aux);
  }
  int parse(const std::string&);
  
  inline pen_parserDataTypes readTag() const{
    return tag;
  }
  
};

struct pen_parserArray{
private:
  std::vector<pen_parserData> vect;

  //Read functions template
  template<class T>
  int readT(T& value, const unsigned index) const{
    if(index >= vect.size())
      return INTDATA_OUT_OF_RANGE;

    return vect[index].read(value);
  }
  
  //Assignment functions template
  template<class T>
  int setT(const T value, const unsigned index){
    if(index >= vect.size())
      return INTDATA_OUT_OF_RANGE;

    vect[index] = value;
    
    return INTDATA_SUCCESS;
  }
  
public:
  pen_parserArray();
  pen_parserArray(const pen_parserArray&);
  //Append element
  void append(const pen_parserData);
  void append(const char);
  void append(const int);
  void append(const double);
  void append(const bool);

  //Remove elements
  int remove(const unsigned);
  int remove(const unsigned, const unsigned);

  //Read functions
  inline int read(char& value, const unsigned index) const{
    return readT(value,index);
  }
  inline int read(int& value, const unsigned index) const{
    return readT(value,index);
  }
  inline int read(double& value, const unsigned index) const{
    return readT(value,index);
  }
  inline int read(bool& value, const unsigned index) const{
    return readT(value,index);
  }
  inline int read(pen_parserData& value, const unsigned index) const{
    return readT(value,index);
  }

  //Assignment functions
  inline int set(const char value, const unsigned index){
    return setT(value,index);
  }
  inline int set(const int value, const unsigned index){
    return setT(value,index);
  }
  inline int set(const double value, const unsigned index){
    return setT(value,index);
  }
  inline int set(const bool value, const unsigned index){
    return setT(value,index);
  }
  inline int set(const pen_parserData value, const unsigned index){
    return setT(value,index);
  }
  
  //Access elements
  pen_parserData& operator[](const unsigned);

  //Get function
  const pen_parserData& at(const unsigned index) const{
    return vect.at(index);
  }
  
  //Assignment operators  
  pen_parserArray& operator=(const pen_parserArray&);

  //Assignment functions
  inline int assign(const pen_parserArray& c){
    vect = c.vect;
    return INTDATA_SUCCESS;
  }
  int assign(const pen_parserElement& c);
  
  inline void clear(){vect.clear();}
  inline unsigned size() const {return vect.size();}

  //Stringify function
  inline std::string stringify() const{
    std::string aux;
    stringify(aux);
    return aux;
  }  
  void stringify(std::string&) const;

  //Parse function
  inline int parse(const char* data){
    std::string aux(data);
    return parse(aux);
  }  
  int parse(const std::string&);
  
};

struct pen_parserElement{
  friend class pen_parserArray;
private:
  enum pen_parserElementTypes tag;
  pen_parserArray array;
  std::string str;

  //Append template functions
  template<class T>
  int appendT(const T data){
    if(tag == SCALAR){
      tag = ARRAY;
    } else if(tag == STRING){
      return INTDATA_NOT_A_ARRAY;
    }

    array.append(data);

    return INTDATA_SUCCESS;
  }

  //Read functions
  template<class T>
  int readT(T& value, const unsigned index = 0) const{
    if(tag == STRING)
      return INTDATA_NOT_A_ARRAY;
    
    return array.read(value,index);
  }

  //Assignment functions for arrays
  template<class T>
  int setT(const T value, const unsigned index = 0){
    if(tag == STRING)
      return INTDATA_NOT_A_ARRAY;

    return array.set(value,index);
  }
  
  
public:

  //Constructors
  pen_parserElement();
  pen_parserElement(const pen_parserElement&);
  pen_parserElement(const pen_parserArray&);
  pen_parserElement(const pen_parserData&);
  pen_parserElement(const std::string&);
  pen_parserElement(const char*);
  pen_parserElement(const char);
  pen_parserElement(const int);
  pen_parserElement(const double);
  pen_parserElement(const bool);

  //Append data
  inline int append(const char data){return appendT(data);}
  inline int append(const int data){return appendT(data);}
  inline int append(const double data){return appendT(data);}
  inline int append(const bool data){return appendT(data);}
  inline int append(const pen_parserData data){return appendT(data);}

  //Remove data
  int remove(const unsigned);
  int remove(const unsigned, const unsigned);

  //Assignment operators
  pen_parserElement& operator=(const pen_parserElement&);
  pen_parserElement& operator=(const pen_parserArray&);
  pen_parserElement& operator=(const pen_parserData&);
  pen_parserElement& operator=(const std::string&);
  pen_parserElement& operator=(const char*);
  pen_parserElement& operator=(const char);
  pen_parserElement& operator=(const int);
  pen_parserElement& operator=(const double);
  pen_parserElement& operator=(const bool);

  //Read functions
  inline int read(char& value, const unsigned index = 0) const {
    return readT(value,index);
  }
  inline int read(int& value, const unsigned index = 0) const {
    return readT(value,index);
  }
  inline int read(double& value, const unsigned index = 0) const {
    return readT(value,index);
  }
  inline int read(bool& value, const unsigned index = 0) const {
    return readT(value,index);
  }
  inline int read(pen_parserData& value, const unsigned index = 0) const {
    return readT(value,index);
  }
  inline int read(std::string& value, const unsigned = 0) const{
    if(tag != STRING)
      return INTDATA_NOT_A_STRING;

    value.assign(str);
    return INTDATA_SUCCESS;
  }
  inline int read(pen_parserArray& value, const unsigned = 0) const{
    if(tag != ARRAY)
      return INTDATA_NOT_A_ARRAY;

    value = array;
    return INTDATA_SUCCESS;
  }

  
  //Assignment functions for arrays
  inline int set(const char value, const unsigned index = 0){
    return setT(value,index);
  }
  inline int set(const int value, const unsigned index = 0){
    return setT(value,index);
  }
  inline int set(const double value, const unsigned index = 0){
    return setT(value,index);
  }
  inline int set(const bool value, const unsigned index = 0){
    return setT(value,index);
  }
  inline int set(const pen_parserData value, const unsigned index = 0){
    return setT(value,index);
  }
  
  //Access array elements
  pen_parserData& operator[](const unsigned);
  
  inline void clear(){
    array.clear();
    str.clear();
    tag = SCALAR;
    array.append('\0');
  }
  inline unsigned size() const {
    switch(tag){
    case ARRAY: return array.size();
    case SCALAR: return 1;
    case STRING: return str.length();
    }
  }
  
  inline pen_parserElementTypes readTag() const{
    return tag;
  }

  //Stringify function
  inline std::string stringify() const{
    std::string aux;
    stringify(aux);
    return aux;
  }  
  void stringify(std::string&) const;

  //Parse function
  inline int parse(const char* data){
    std::string aux(data);
    return parse(aux);
  }  
  int parse(const std::string&);
  
};

struct pen_parserSection{

  typedef std::map<std::string,pen_parserElement> elementMap;
  
private:
  elementMap elements;

  int trusted_remove(const std::string&);
  int trusted_sectionRange(const std::string&, elementMap::const_iterator&, elementMap::const_iterator&) const;

  //Functions to append elements to an existing array
  template<class T>
  int appendT(std::string& key, const T data){

    //Check key
    int err = checkKey(key);
    if(err != INTDATA_SUCCESS)
      return err;
    
    //Check if exists
    if(!isElement(key)){
      return INTDATA_NOT_A_ELEMENT;
    }

    //Try to append data
    return elements[key].append(data);
  }

  //Functions to set elements
  template<class T>
  int setT(std::string& key, const T data, const bool overwrite = false){

    //Check key
    int err = checkKey(key);
    if(err != INTDATA_SUCCESS)
      return err;

    //Check if the key can be overwritten
    if(overwrite){
      //Overwrite enabled
      //Try to remove previous key
      trusted_remove(key);
    }
    else{
      //Overwrite disabled
      //Check if this key already exists
      if(exists(key)){
	return INTDATA_USED_KEY;
      }
    }
    
    //Set element
    elements[key] = data;

    //Return success
    return INTDATA_SUCCESS;    
  }

  //Set elements in array
  template<class T>
  int setPosT(std::string& key, const unsigned index, const T data){
    
    //Check if is a valid key
    int err = checkKey(key);
    if(err != INTDATA_SUCCESS){
      return err;
    }

    //Check if this key is an existing array
    if(!isElement(key)){
      return INTDATA_NOT_A_ELEMENT;
    }
    
    return elements[key].set(data,index);
  }

  //Read elements
  template<class T>
  int readT(std::string& key, T& data, const unsigned index = 0) const{

    //Check if is a valid key
    int err = checkKey(key);
    if(err != INTDATA_SUCCESS){
      return err;
    }

    //Check if this key is an existing element
    if(!isElement(key)){
      return INTDATA_NOT_A_ELEMENT;
    }

    //Try to get the value
    return elements.at(key).read(data,index);
  }

  
public:

  pen_parserSection();
  pen_parserSection(const pen_parserSection&);

  //Check functions
  int checkKey(const char*) const;  
  int checkKey(std::string&) const;
  //Find key functions
  bool isSection(const std::string&) const;
  bool isElement(const std::string&) const;
  inline bool exists(const std::string& key) const{
    return (isElement(key) || isSection(key));  
  }
  
  //Create array
  int createArray(std::string&, const bool);

  //Functions to append elements to an existing array
  inline int append(std::string& key, const char data){
    return appendT(key,data);
  }
  inline int append(std::string& key, const int data){
    return appendT(key,data);
  }
  inline int append(std::string& key, const double data){
    return appendT(key,data);
  }
  inline int append(std::string& key, const bool data){
    return appendT(key,data);
  }
  inline int append(std::string& key, const pen_parserData data){
    return appendT(key,data);
  }

  inline int append(const char* key, const char data){
    std::string aux(key);
    return appendT(aux,data);
  }
  inline int append(const char* key, const int data){
    std::string aux(key);
    return appendT(aux,data);
  }
  inline int append(const char* key, const double data){
    std::string aux(key);
    return appendT(aux,data);
  }
  inline int append(const char* key, const bool data){
    std::string aux(key);
    return appendT(aux,data);
  }
  inline int append(const char* key, const pen_parserData data){
    std::string aux(key);
    return appendT(aux,data);
  }
  
  //Remove elements
  int remove(const char*);
  int remove(std::string&);
  int remove(std::string&, const unsigned);
  int remove(std::string&, const unsigned, const unsigned);

  //Assign from input section
  void assign(const pen_parserSection&);
  
  //Set elements

  inline int set(std::string& key, const char data, const bool overwrite = false){
    return setT(key,data,overwrite);
  }
  inline int set(std::string& key, const int data, const bool overwrite = false){
    return setT(key,data,overwrite);
  }  
  inline int set(std::string& key, const double data, const bool overwrite = false){
    return setT(key,data,overwrite);
  }
  inline int set(std::string& key, const bool data, const bool overwrite = false){
    return setT(key,data,overwrite);
  }
  inline int set(std::string& key, const pen_parserData data, const bool overwrite = false){
    return setT(key,data,overwrite);
  }
  inline int set(std::string& key, const pen_parserElement& data, const bool overwrite = false){
    return setT(key,data,overwrite);
  }  
  inline int set(std::string& key, const pen_parserArray& data, const bool overwrite = false){
    return setT(key,data,overwrite);
  }
  inline int set(std::string& key, const std::string& data, const bool overwrite = false){
    return setT(key,data,overwrite);
  }  

  inline int set(const char* key, const char data, const bool overwrite = false){
    std::string aux(key);
    return setT(aux,data,overwrite);
  }
  inline int set(const char* key, const int data, const bool overwrite = false){
    std::string aux(key);
    return setT(aux,data,overwrite);
  }  
  inline int set(const char* key, const double data, const bool overwrite = false){
    std::string aux(key);
    return setT(aux,data,overwrite);
  }
  inline int set(const char* key, const bool data, const bool overwrite = false){
    std::string aux(key);
    return setT(aux,data,overwrite);
  }
  inline int set(const char* key, const pen_parserData data, const bool overwrite = false){
    std::string aux(key);
    return setT(aux,data,overwrite);
  }
  inline int set(const char* key, const pen_parserElement& data, const bool overwrite = false){
    std::string aux(key);
    return setT(aux,data,overwrite);
  }  
  inline int set(const char* key, const pen_parserArray& data, const bool overwrite = false){
    std::string aux(key);
    return setT(aux,data,overwrite);
  }
  inline int set(const char* key, const std::string& data, const bool overwrite = false){
    std::string aux(key);
    return setT(aux,data,overwrite);
  }    
    
  
  //Set elements in array
  inline int setPos(std::string& key, const unsigned index, const char data){
    return setPosT(key,index,data);
  }
  inline int setPos(std::string& key, const unsigned index, const int data){
    return setPosT(key,index,data);
  }
  inline int setPos(std::string& key, const unsigned index, const double data){
    return setPosT(key,index,data);
  }
  inline int setPos(std::string& key, const unsigned index, const bool data){
    return setPosT(key,index,data);
  }
  inline int setPos(std::string& key, const unsigned index, const pen_parserData data){
    return setPosT(key,index,data);
  }

  inline int setPos(const char* key, const unsigned index, const char data){
    std::string aux(key);
    return setPosT(aux,index,data);
  }
  inline int setPos(const char* key, const unsigned index, const int data){
    std::string aux(key);
    return setPosT(aux,index,data);
  }
  inline int setPos(const char* key, const unsigned index, const double data){
    std::string aux(key);
    return setPosT(aux,index,data);
  }
  inline int setPos(const char* key, const unsigned index, const bool data){
    std::string aux(key);
    return setPosT(aux,index,data);
  }
  inline int setPos(const char* key, const unsigned index, const pen_parserData data){
    std::string aux(key);
    return setPosT(aux,index,data);
  }
  
  //Read elements
  inline int read(std::string& key, char& data, const unsigned index = 0) const{
    return readT(key,data,index);
  }
  inline int read(std::string& key, int& data, const unsigned index = 0) const{
    return readT(key,data,index);
  }
  inline int read(std::string& key, double& data, const unsigned index = 0) const{
    return readT(key,data,index);
  }
  inline int read(std::string& key, bool& data, const unsigned index = 0) const{
    return readT(key,data,index);
  }
  inline int read(std::string& key, pen_parserData& data, const unsigned index = 0) const{
    return readT(key,data,index);
  }
  inline int read(std::string& key, std::string& data, const unsigned index = 0) const{
    return readT(key,data,index);
  }  
  inline int read(std::string& key, pen_parserArray& data, const unsigned = 0) const{
    //Check if is a valid key
    int err = checkKey(key);
    if(err != INTDATA_SUCCESS){
      return err;
    }

    //Check if this key is an existing element
    if(!isElement(key)){
      return INTDATA_NOT_A_ELEMENT;
    }

    //Try to copy element to array
    return data.assign(elements.at(key));
  }
  inline int read(std::string& key, pen_parserSection& data, const unsigned = 0) const{
    return readSubsection(key,data);
  }  
  
  inline int read(const char* key, char& data, const unsigned index = 0) const{
    std::string aux(key);
    return readT(aux,data,index);
  }
  inline int read(const char* key, int& data, const unsigned index = 0) const{
    std::string aux(key);
    return readT(aux,data,index);
  }
  inline int read(const char* key, double& data, const unsigned index = 0) const{
    std::string aux(key);
    return readT(aux,data,index);
  }
  inline int read(const char* key, bool& data, const unsigned index = 0) const{
    std::string aux(key);
    return readT(aux,data,index);
  }
  inline int read(const char* key, pen_parserData& data, const unsigned index = 0) const{
    std::string aux(key);
    return readT(aux,data,index);
  }
  inline int read(const char* key, std::string& data, const unsigned index = 0) const{
    std::string aux(key);
    return readT(aux,data,index);
  }    
  inline int read(const char* key, pen_parserArray& data, const unsigned = 0) const{
    std::string aux(key);
    return read(aux,data);
  }
  inline int read(const char* key, pen_parserSection& data, const unsigned = 0) const{
    std::string aux(key);
    return readSubsection(aux,data);
  }
  
  inline int readSubsection(const char* key, pen_parserSection& secOut, const bool removeKey = true) const{
    std::string aux(key);
    return readSubsection(aux,secOut,removeKey);
  }
  int readSubsection(std::string&, pen_parserSection&, const bool = true) const;
  
  void clear();
  
  pen_parserSection& operator=(const pen_parserSection&);

  //Stringify function
  inline std::string stringify() const{
    std::string aux;
    stringify(aux);
    return aux;
  }  
  void stringify(std::string&) const;

  //Parse function
  inline int parse(const char* data, const bool overwrite = false){
    std::string aux(data);
    return parse(aux,overwrite);
  }
  int parse(const std::string&, const bool overwrite = false);
  int parse(std::string&, std::string&, const bool overwrite = false);

  //ls function
  inline int ls(std::vector<std::string>& vect) const{
    
    std::string aux("");
    return ls(aux,vect);
  }  
  inline int ls(const char* key, std::vector<std::string>& vect) const{

    std::string aux(key);
    return ls(aux,vect);
  }
  int ls(std::string&, std::vector<std::string>&) const;
  
  ~pen_parserSection();
};

inline std::string triml(const std::string& strin){
  const std::string delimiters = " \n\r\t\f\v";
  size_t first = strin.find_first_not_of(delimiters);
  return (first == std::string::npos) ? "" : strin.substr(first);
}

inline std::string trimr(const std::string& strin){
  const std::string delimiters = " \n\r\t\f\v";
  size_t last = strin.find_last_not_of(delimiters);
  return (last == std::string::npos) ? "" : strin.substr(0,last+1);
}

inline std::string trim(const std::string& strin){
  return trimr(triml(strin));
}

//Auxiliar functions
int parseFile(const char* filename,
	      pen_parserSection& section,
	      std::string& errorString,
	      long unsigned& errorLine);

int pen_getLine(FILE* fin,
		const unsigned size,
		char* line,
		unsigned long& nlines);

#endif
