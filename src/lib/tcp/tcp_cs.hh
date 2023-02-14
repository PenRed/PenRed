
//
//
//    Copyright (C) 2020 Universitat de València - UV
//    Copyright (C) 2020 Universitat Politècnica de València - UPV
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
//        vicent.gimenez.alventosa@gmail.com (Vicent Giménez Alventosa)
//        vicente.gimenez@uv.es (Vicente Giménez Gómez)
//    
//

#ifndef __PEN_TCP_CS__
#define __PEN_TCP_CS__

#include <cstdio>
#include "asio.hpp"

#ifdef _PEN_USE_SSL_
#include "asio/ssl.hpp"
#endif

enum PEN_TCP_ERR{
		 PEN_TCP_SUCCESS,
		 PEN_TCP_BAD_SSL_CONTEXT_OPTIONS,
		 PEN_TCP_BAD_SSL_VERIFY_MODE,
		 PEN_TCP_CALLBACK_SET_FAIL,
		 PEN_TCP_INVALID_CA_CERT,
		 PEN_TCP_INVALID_CERT,
		 PEN_TCP_INVALID_KEY,
		 PEN_TCP_INVALID_DH_PARAM,
		 PEN_TCP_ERROR_ON_BIND,
		 PEN_TCP_WRITE_FAIL,
		 PEN_TCP_READ_FAIL,
		 PEN_TCP_ACCEPT_CONNECTION_FAIL,
		 PEN_TCP_HANDSHAKE_FAILED,
		 PEN_TCP_SSL_STREAM_READ_FAILED,
		 PEN_TCP_EXCEPTION_CATCH,
		 PEN_TCP_SOCKET_READ_FAILED,
		 PEN_TCP_ACCEPTOR_OPEN_FAIL,
		 PEN_TCP_ACCEPTOR_LISTEN_FAIL,
		 PEN_TCP_SOCKET_ALREADY_ACTIVE,
		 PEN_TCP_INACTIVE_SOCKET,
		 PEN_TCP_NULL_POINTER,
		 PEN_TCP_CONNECTION_FAIL,
		 PEN_TCP_UNABLE_TO_RESOLVE_HOSTNAME,
		 PEN_TCP_UNABLE_TO_CREATE_REQUEST,
		 PEN_TCP_SSL_SUPPORT_DISABLED,
};

namespace pen_tcp{

  //TCP message length
  static const size_t messageLen = 512;

  //Fill message function to fit the expected message length
  inline static void fillmessage(char msg[messageLen],
				 long int offset = -1){
    size_t len;
    if(offset <= 0)
      len = strlen(msg);
    else
      len = offset;
    if(len < messageLen){
      size_t tofill = messageLen-len;
      memset(&msg[len],' ',tofill);
    }
  }
  
  std::string trim(const std::string& str,
		   const std::string& whitespace = " \t\n");
  
  static inline bool filterLog(FILE* flog,
			       const unsigned verbose,
			       const unsigned required) {
    if(flog != nullptr && verbose > required )
      return true;
    return false;
  }
  
  //*********************** SSL *************************//
#ifdef _PEN_USE_SSL_

  void run(asio::io_context& io_context,
	   std::chrono::seconds timeout,
	   asio::ssl::stream<asio::ip::tcp::socket>*& stream,
	   FILE* flog = nullptr, const unsigned verbose = 1);
  
  int listen(asio::io_context& io_context,
	     asio::ssl::context& ssl_context,
	     asio::ip::tcp::acceptor& acceptor,
	     asio::ssl::stream<asio::ip::tcp::socket>*& stream,
	     std::chrono::seconds timeout = std::chrono::seconds(5),
	     FILE* flog = nullptr, const unsigned verbose = 1);
  
  int connect(asio::io_context& io_context,
	      asio::ssl::context& ssl_context,
	      const std::vector<asio::ip::tcp::endpoint>& endpoints,
	      asio::ssl::stream<asio::ip::tcp::socket>*& stream,
	      std::chrono::seconds timeout,
	      FILE* flog = nullptr, const unsigned verbose = 1);

  int receive(std::string& message,
	      asio::io_context& io_context,
	      asio::ssl::stream<asio::ip::tcp::socket>*& stream,
	      std::chrono::seconds timeout,
	      FILE* flog = nullptr, const unsigned verbose = 1);
  
  int write(char message[messageLen], long int len,
	    asio::io_context& io_context,
	    std::chrono::seconds timeout,
	    asio::ssl::stream<asio::ip::tcp::socket>*& stream,
	    FILE* flog = nullptr, const unsigned verbose = 1);
  
  void close(asio::io_context& io_context,
	   asio::ssl::stream<asio::ip::tcp::socket>*& stream);
#endif
  //*********************** SSL END *************************//
  
void run(asio::io_context& io_context,
	 std::chrono::seconds timeout,
	 asio::ip::tcp::socket*& socket,
	 FILE* flog = nullptr, const unsigned verbose = 1);

int listen(asio::io_context& io_context,
	     asio::ip::tcp::acceptor& acceptor,
	     asio::ip::tcp::socket*& socket,
	     FILE* flog = nullptr, const unsigned verbose = 1);  
  
  int connect(asio::io_context& io_context,
	      const std::vector<asio::ip::tcp::endpoint>& endpoints,
	      asio::ip::tcp::socket*& socket,
	      std::chrono::seconds timeout,
	      FILE* flog = nullptr, const unsigned verbose = 1);
  
  int receive(std::string& message,
	      asio::io_context& io_context,
	      asio::ip::tcp::socket*& socket,
	      std::chrono::seconds timeout,
	      FILE* flog = nullptr, const unsigned verbose = 1);

  int write(char message[messageLen], long int len,
	    asio::io_context& io_context,
	    std::chrono::seconds timeout,
	    asio::ip::tcp::socket*& socket,
	    FILE* flog = nullptr, const unsigned verbose = 1);
  
  void close(asio::ip::tcp::socket*& socket);

  
  class server_seq{

  private:
    
    FILE* flog;

    asio::io_context io_context;
    
#ifdef _PEN_USE_SSL_
    asio::ssl::context ssl_context;
#endif
    asio::ip::tcp::acceptor acceptor;
    
    std::chrono::seconds timeout;

#ifdef _PEN_USE_SSL_
    asio::ssl::stream<asio::ip::tcp::socket>* stream;
#endif
    asio::ip::tcp::socket* socket;
    unsigned verbose;

    std::chrono::steady_clock::time_point tinit;
    bool sslEnabled;

  public:

    inline bool usingSSL() const {return sslEnabled;}
    
    server_seq(const unsigned verb = 1);

    inline std::chrono::seconds::rep timeStamp(const std::chrono::steady_clock::time_point time) const{
      return std::chrono::duration_cast
	<std::chrono::seconds>(time-tinit).count();
    }
    
    inline std::chrono::seconds::rep timeStamp() const{
      return timeStamp(std::chrono::steady_clock::now());
    }

    inline void setInit(const std::chrono::steady_clock::time_point& t){
      tinit = t;
    }
    
    inline bool filterLog(const unsigned required) const {
      if(flog != nullptr && verbose > required )
	return true;
      return false;
    }

    inline void setLogFile(FILE* newflog){flog = newflog;}

    inline void setVerbose(const unsigned verb){verbose = verb;}

#ifdef _PEN_USE_SSL_
    int init(const bool local,
	     const unsigned short port,
	     const char* CAfilename,
	     const char* certFilename,
	     const char* keyFilename,
	     const char* dhFilename,
	     const char* password);
#endif
    int init(const bool local,
	     const unsigned short port);
    
    inline int listen(){
      if(filterLog(1)){
	fprintf(flog,"%07ld s - listening for incomming connection\n",
		static_cast<long int>(timeStamp()));
	fflush(flog);
      }

#ifdef _PEN_USE_SSL_
      if(sslEnabled)
	return pen_tcp::listen(io_context,ssl_context,acceptor,
			       stream,timeout,flog,verbose);
#endif
      return pen_tcp::listen(io_context,acceptor,socket,flog,verbose);
    }

    inline int receive(std::string& message){
      if(filterLog(1)){
	fprintf(flog,"%07ld s - Receiving message\n",
		static_cast<long int>(timeStamp()));
	fflush(flog);
      }
      
#ifdef _PEN_USE_SSL_
      if(sslEnabled)
	return pen_tcp::receive(message,io_context,
				stream,timeout,flog,verbose);
#endif
      return pen_tcp::receive(message,io_context,
			      socket,timeout,flog,verbose);
    }
        
    inline int write(char message[messageLen], long int len = -1){
      if(filterLog(1)){
	fprintf(flog,"%07ld s - Writing message\n",
		static_cast<long int>(timeStamp()));
	fflush(flog);
      }
      
#ifdef _PEN_USE_SSL_
      if(sslEnabled)
	return pen_tcp::write(message,len,io_context,
			      timeout,stream,flog,verbose);
#endif
      return pen_tcp::write(message,len,io_context,
			    timeout,socket,flog,verbose);
    }

    inline void close(){
      if(filterLog(1)){
	fprintf(flog,"%07ld s - Close connection\n",
		static_cast<long int>(timeStamp()));
	fflush(flog);
      }
      
#ifdef _PEN_USE_SSL_
      if(sslEnabled){
	pen_tcp::close(io_context,stream);
	return;
      }
#endif
      pen_tcp::close(socket);
         
    }

    void clear();
    
    ~server_seq(){
      clear();
    }
  };

  class client{

  private:

    FILE* flog;
    
    asio::io_context io_context;
#ifdef _PEN_USE_SSL_
    asio::ssl::context ssl_context;
#endif
    std::chrono::seconds timeout;

#ifdef _PEN_USE_SSL_
    asio::ssl::stream<asio::ip::tcp::socket>* stream;
#endif
    asio::ip::tcp::socket* socket;

    std::vector<asio::ip::tcp::endpoint> serverEndpoints;
    
    unsigned verbose;

    bool sslEnabled;
    
  public:

    inline bool usingSSL() const {return sslEnabled;}
    
    client(const unsigned verb = 1);

    int setHost(const char* host,
		const char* port);
    
#ifdef _PEN_USE_SSL_
    int initSSL(const char* CAfilename,
		const char* certFilename,
		const char* keyFilename,
		const char* password,
		const char* hostname);
#endif
    
    inline bool filterLog(const unsigned required) const {
      if(flog != nullptr && verbose > required )
	return true;
      return false;
    }

    inline void setLogFile(FILE* newflog){flog = newflog;}

    inline void setVerbose(const unsigned verb){verbose = verb;}
    
    inline int connect(){
      if(filterLog(3)){
	fprintf(flog,"Connecting\n");
	fflush(flog);
      }

#ifdef _PEN_USE_SSL_
      if(sslEnabled)
	return pen_tcp::connect(io_context,ssl_context,serverEndpoints,
			      stream,timeout,flog,verbose);
#endif
      return pen_tcp::connect(io_context,serverEndpoints,
			      socket,timeout,flog,verbose);
    }

    inline int receive(std::string& message){
      if(filterLog(3)){
	fprintf(flog,"Receiving message\n");
	fflush(flog);
      }

#ifdef _PEN_USE_SSL_
      if(sslEnabled)
	return pen_tcp::receive(message,io_context,
				stream,timeout,flog,verbose);
#endif
      return pen_tcp::receive(message,io_context,
			      socket,timeout,flog,verbose);
    }
        
    inline int write(char message[messageLen], long int len = -1){
      if(filterLog(3)){
	fprintf(flog,"Writing message\n");
	fflush(flog);
      }
      
#ifdef _PEN_USE_SSL_
      if(sslEnabled)
	return pen_tcp::write(message,len,io_context,
			      timeout,stream,flog,verbose);
#endif
      return pen_tcp::write(message,len,io_context,
			    timeout,socket,flog,verbose);
    }

    inline void close(){
      if(filterLog(3)){
	fprintf(flog,"Close connection\n");
	fflush(flog);
      }
      
#ifdef _PEN_USE_SSL_
      if(sslEnabled){
	pen_tcp::close(io_context,stream);
	return;
      }
#endif    
      pen_tcp::close(socket);
            
    }
    
    void clear();
    
    ~client(){
      clear();
    }
  };
  
}

#endif
