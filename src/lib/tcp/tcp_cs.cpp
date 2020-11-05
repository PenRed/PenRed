
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

#include "tcp_cs.hh"

  //*********************** SSL *************************//
#ifdef _PEN_USE_SSL_

pen_tcp::server_seq::server_seq(const unsigned verb) :
  ssl_context(asio::ssl::context::sslv23),
  acceptor(io_context),
  stream(nullptr), socket(nullptr),
  verbose(verb), sslEnabled(false){
      
  tinit = std::chrono::steady_clock::now();
  timeout = std::chrono::seconds(5);      
}

int pen_tcp::server_seq::init(const bool local,
			      const unsigned short port,
			      const char* CAfilename,
			      const char* certFilename,
			      const char* keyFilename,
			      const char* dhFilename,
			      const char* password){

  asio::error_code ec;

  //Close acceptor, just in case
  acceptor.close(ec);
  
  if(filterLog(1)){
    fprintf(flog,"%07ld s - Starting server monitor listening at port %u\n",
	    static_cast<long int>(timeStamp()),port);
    fflush(flog);
  }

  sslEnabled = false;
  if(certFilename != nullptr && keyFilename != nullptr &&
     dhFilename != nullptr   && CAfilename != nullptr)
    sslEnabled = true;
    
  asio::ip::tcp::endpoint endpoint;
  if(local)
    endpoint = asio::ip::tcp::endpoint(asio::ip::address::from_string("127.0.0.1"),port);
  else
    endpoint = asio::ip::tcp::endpoint(asio::ip::tcp::v4(),port);
  
  if(sslEnabled){

    if(filterLog(1)){
      fprintf(flog,"%07ld s - Enabled SSL communication\n"
	      "        CA Cert: %s\n"
	      "     Server Key: %s\n"
	      "    Server Cert: %s\n"
	      "             dh: %s\n",
	      static_cast<long int>(timeStamp()),
	      CAfilename,keyFilename,certFilename,dhFilename);
      fflush(flog);
    }
    
    std::string spassword;
    if(password != nullptr)
      spassword.assign(password);
    
    
    ssl_context.set_options(asio::ssl::context::default_workarounds |
			    asio::ssl::context::no_sslv2            |
			    asio::ssl::context::single_dh_use,
			    ec);

    if(ec){
      if(filterLog(0)){
	fprintf(flog,"%07ld s - Error: Unable to set SSL context options\n"
		"                   Error: %s\n"
		"              Error code: %d\n",
		static_cast<long int>(timeStamp()),
		ec.message().c_str(),ec.value());
	fflush(flog);
      }
      return PEN_TCP_BAD_SSL_CONTEXT_OPTIONS;
    }

    //Require peer certificate
    ssl_context.set_verify_mode(asio::ssl::verify_peer |
				asio::ssl::verify_fail_if_no_peer_cert,
				ec);
    
    if(ec){
      if(filterLog(0)){
	fprintf(flog,"%07ld s - Error: Unable to set SSL verify mode\n"
		"                   Error: %s\n"
		"              Error code: %d\n",
		static_cast<long int>(timeStamp()),
		ec.message().c_str(),ec.value());
	fflush(flog);
      }
      return PEN_TCP_BAD_SSL_VERIFY_MODE;
    }

    if(password != nullptr){
      ssl_context.set_password_callback([spassword](std::size_t max_length,
						    asio::ssl::context::password_purpose /*purpose*/) -> std::string{
					  if(spassword.length() > max_length){
					    printf("Warning: provided password"
						   "length exceeds the limit.\n"
						   "   length: %lu\n"
						   "      max: %lu\n",
						   static_cast<unsigned long>
						   (spassword.length()),
						   static_cast<unsigned long>
						   (max_length));
					    return spassword.substr(0,max_length);
					  }
					  else
					    return spassword;
					},ec);

      if(ec){
	if(filterLog(0)){
	  fprintf(flog,"%07ld s - Error: Unable to set password callback\n"
		  "                   Error: %s\n"
		  "              Error code: %d\n",
		  static_cast<long int>(timeStamp()),
		  ec.message().c_str(),ec.value());
	  fflush(flog);
	}
	return PEN_TCP_CALLBACK_SET_FAIL;
      }
      
    }

    //Load CA certificate
    ssl_context.load_verify_file(CAfilename,ec);
    if(ec){
      if(filterLog(0)){
	fprintf(flog,"%07ld s - Error loading CA certificate\n"
		"                   Error: %s\n"
		"              Error code: %d\n",
		static_cast<long int>(timeStamp()),
		ec.message().c_str(),ec.value());
	fflush(flog);
      }
      return PEN_TCP_INVALID_CA_CERT;
    }

    //Add server certificate file
    ssl_context.use_certificate_chain_file(certFilename,ec);
    if(ec){
      if(filterLog(0)){
	fprintf(flog,"%07ld s - Error loading SSL certificate chain file\n"
		"                   Error: %s\n"
		"              Error code: %d\n",
		static_cast<long int>(timeStamp()),
		ec.message().c_str(),ec.value());
	fflush(flog);
      }
      return PEN_TCP_INVALID_CERT;
    }

    //Add server key file
    ssl_context.use_private_key_file(keyFilename,asio::ssl::context::pem,ec);
    if(ec){
      if(filterLog(0)){
	fprintf(flog,"%07ld s - Error loading SSL key file\n"
		"                   Error: %s\n"
		"              Error code: %d\n",
		static_cast<long int>(timeStamp()),
		ec.message().c_str(),ec.value());
	fflush(flog);
      }
      return PEN_TCP_INVALID_KEY;
    }
    
    ssl_context.use_tmp_dh_file(dhFilename,ec);
    if(ec){
      if(filterLog(0)){
	fprintf(flog,"%07ld s - Error loading Diffie-Hellman parameters file\n"
		"                   Error: %s\n"
		"              Error code: %d\n",
		static_cast<long int>(timeStamp()),
		ec.message().c_str(),ec.value());
	fflush(flog);
      }
      return PEN_TCP_INVALID_DH_PARAM;
    }

  }

  //Configure the acceptor

  //Open acceptor using tcp protocol and IPv4 
  acceptor.open(asio::ip::tcp::v4(), ec);
  if(ec){
    if(filterLog(0)){
      fprintf(flog,"%07ld s - Error at acceptor open\n"
	      "                   Error: %s\n"
	      "              Error code: %d\n",
	      static_cast<long int>(timeStamp()),
	      ec.message().c_str(),ec.value());
      fflush(flog);
    }
    return PEN_TCP_ACCEPTOR_OPEN_FAIL;
  }

  //Bind the acceptor
  acceptor.bind(endpoint,ec);
  if(ec){
    if(filterLog(0)){
      fprintf(flog,"%07ld s - Error binding to specified local endpoint\n"
	      "                   Error: %s\n"
	      "              Error code: %d\n",
	      static_cast<long int>(timeStamp()),
	      ec.message().c_str(),ec.value());
      fflush(flog);
    }
    return PEN_TCP_ERROR_ON_BIND;    
  }

  acceptor.listen(asio::socket_base::max_listen_connections, ec);
  if(ec){
    if(filterLog(0)){
      fprintf(flog,"%07ld s - Error on acceptor listen\n"
	      "                   Error: %s\n"
	      "              Error code: %d\n",
	      static_cast<long int>(timeStamp()),
	      ec.message().c_str(),ec.value());
      fflush(flog);
    }
    return PEN_TCP_ACCEPTOR_LISTEN_FAIL;    
  }

  return PEN_TCP_SUCCESS;
}

void pen_tcp::run(asio::io_context& io_context,
		  std::chrono::seconds timeout,
		  asio::ssl::stream<asio::ip::tcp::socket>*& stream,
		  FILE* flog, const unsigned verbose){
  //Restart context
  io_context.restart();
  
  //Run context with the timeout
  io_context.run_for(timeout);
  
  //Check if context has been stopped
  if(!io_context.stopped()){
    //Context has not been stoped, the timeout has reached
    if(filterLog(flog,verbose,0)){
      fprintf(flog,"Timeout on SSL action\n");
      fflush(flog);
    }

    if(stream != nullptr){
      std::error_code ec;
      //Close the connection
      //stream->next_layer().cancel();
      stream->lowest_layer().close(ec); //close TCP
      if(ec){
	if(filterLog(flog,verbose,0)){
	  fprintf(flog,"Unable to shutdown TCP\n"
		  "              Error: %s\n"
		  "         Error code: %d\n",
		  ec.message().c_str(),ec.value());
	  fflush(flog);
	}
      }
    }

    //Run context to finish the operation
    io_context.run();      
  }
}

int pen_tcp::write(char message[messageLen], long int len,
		   asio::io_context& io_context,
		   std::chrono::seconds timeout,
		   asio::ssl::stream<asio::ip::tcp::socket>*& stream,
		   FILE* flog, const unsigned verbose){
  try{

    if(stream == nullptr)
      return PEN_TCP_INACTIVE_SOCKET;
    
    fillmessage(message,len);
    
    //Add the write handler
    asio::error_code ec;
      
    asio::async_write(*stream,asio::buffer(message,messageLen),
		      [&ec](const std::error_code& error,
			    std::size_t /*result_n*/){
			ec = error;
		      });
    run(io_context,timeout,stream,flog,verbose);

    if(ec){
      if(filterLog(flog,verbose,0)){
	fprintf(flog,"Error: Unable to send message\n"
		"              Error: %s\n"
		"         Error code: %d\n",
		ec.message().c_str(),ec.value());
	fflush(flog);
      }
      return PEN_TCP_WRITE_FAIL;
    }
    else if(filterLog(flog,verbose,3)){
      fprintf(flog,"Message sent\n");
      fflush(flog);
    }    
    return PEN_TCP_SUCCESS;    
  }
  catch (std::exception& e){
    if(filterLog(flog,verbose,0)){
      fprintf(flog,"Exception on message send: %s\n",
	      e.what());
      fflush(flog);
    }
    return PEN_TCP_EXCEPTION_CATCH;
  }

}

int pen_tcp::listen(asio::io_context& io_context,
		    asio::ssl::context& ssl_context,
		    asio::ip::tcp::acceptor& acceptor,
		    asio::ssl::stream<asio::ip::tcp::socket>*& stream,
		    std::chrono::seconds timeout,
		    FILE* flog, const unsigned verbose){

  try{

    if(stream != nullptr)
      return PEN_TCP_SOCKET_ALREADY_ACTIVE;

    //Create a new stream
    stream = new asio::ssl::stream<asio::ip::tcp::socket>(io_context,ssl_context);
    
    asio::error_code ec;
    //Accept next communication
    acceptor.accept(stream->lowest_layer(),ec);    
    if(ec){
      if(filterLog(flog,verbose,0)){
	fprintf(flog,"Connection accept failed\n"
		"              Error: %s\n"
		"         Error code: %d\n",
		ec.message().c_str(),ec.value());
	fflush(flog);
      }
      return PEN_TCP_ACCEPT_CONNECTION_FAIL;
    }

    if(filterLog(flog,verbose,1)){
      fprintf(flog,"Incoming communication from %s\n",
	      stream->lowest_layer().remote_endpoint().address().to_string().c_str());
      fflush(flog);
    }

    //Perform the handshake
    stream->async_handshake(asio::ssl::stream_base::server,
			    [&ec](const asio::error_code& error){
				ec = error;
			    });

    //Start timeout
    run(io_context,timeout,stream,flog,verbose);

    if(ec){
      if(filterLog(flog,verbose,0)){
	fprintf(flog,"Handshake failed\n"
		"              Error: %s\n"
		"         Error code: %d\n",
		ec.message().c_str(),ec.value());
	fflush(flog);
      }
      return PEN_TCP_HANDSHAKE_FAILED;
    }
    
    //Successfully handshake.    
    return PEN_TCP_SUCCESS;
    
  }catch (std::exception& e){
    if(filterLog(flog,verbose,0)){
      fprintf(flog,"Exception on SSL stream receive: %s\n",
	      e.what());
      fflush(flog);
    }
    return PEN_TCP_EXCEPTION_CATCH;
  }
}

int pen_tcp::connect(asio::io_context& io_context,
		     asio::ssl::context& ssl_context,
		     const std::vector<asio::ip::tcp::endpoint>& endpoints,
		     asio::ssl::stream<asio::ip::tcp::socket>*& stream,
		     std::chrono::seconds timeout,
		     FILE* flog, const unsigned verbose){
  
  if(stream != nullptr)
    return PEN_TCP_SOCKET_ALREADY_ACTIVE;
  
  //Create a new stream
  stream = new asio::ssl::stream<asio::ip::tcp::socket>(io_context,ssl_context);
  
  std::error_code ec;
  asio::async_connect(stream->lowest_layer(), endpoints,
		      [&ec](const std::error_code& error,
			    const asio::ip::tcp::endpoint& /*endpoint*/){
			ec = error;
		      });

  run(io_context,timeout,stream,flog,verbose);
  
  if(ec){
    if(filterLog(flog,verbose,0)){
      fprintf(flog,"Connection failed\n"
	      "              Error: %s\n"
	      "         Error code: %d\n",
	      ec.message().c_str(),ec.value());
      fflush(flog);
    }
    return PEN_TCP_CONNECTION_FAIL;
  }

  //Perform the handshake
  ec.clear();
  stream->async_handshake(asio::ssl::stream_base::client,
			  [&ec](const asio::error_code& error){
			    ec = error;
			  });

  //Start timeout
  run(io_context,timeout,stream,flog,verbose);

  if(ec){
    if(filterLog(flog,verbose,0)){
      fprintf(flog,"Handshake failed\n"
	      "              Error: %s\n"
	      "         Error code: %d\n",
	      ec.message().c_str(),ec.value());
      fflush(flog);
    }
    return PEN_TCP_HANDSHAKE_FAILED;
  }

  return PEN_TCP_SUCCESS;
}

int pen_tcp::receive(std::string& message,
		     asio::io_context& io_context,
		     asio::ssl::stream<asio::ip::tcp::socket>*& stream,
		     std::chrono::seconds timeout,
		     FILE* flog, const unsigned verbose){

  if(stream == nullptr){
    if(filterLog(flog,verbose,0)){
      fprintf(flog,"Invalid stream (null pointer)\n");
      fflush(flog);
    }    
    return PEN_TCP_INACTIVE_SOCKET;
  }
  
  asio::error_code ec;
  //Read from stream
  message.clear();
  size_t read = 0;
  for(;;){    
    char buf[messageLen+1];
    size_t len;
    stream->async_read_some(asio::buffer(buf),
			    [&ec,&len](const asio::error_code& error,
				       std::size_t bytes_transferred){
			      ec = error;
			      len = bytes_transferred;
			    });

    //Start timeout
    run(io_context,timeout,stream,flog,verbose);

    if(ec == asio::error::eof){
      break;
    }
    else if(ec){
      if(filterLog(flog,verbose,0)){
	fprintf(flog,"Stream read failed\n"
		"             Error: %s\n"
		"        Error code: %d\n",
		ec.message().c_str(),ec.value());
	fflush(flog);
      }
      return PEN_TCP_SSL_STREAM_READ_FAILED;
    }
    read += len;
    message.append(buf,len);
    if(read >= messageLen){
      message = message.substr(0,messageLen);
      break;
    }
  }
  //trim message
  message = trim(message);
  return PEN_TCP_SUCCESS;
}

void pen_tcp::close(asio::io_context& io_context,
		    asio::ssl::stream<asio::ip::tcp::socket>*& stream){
  try{
    if(stream != nullptr){
      //stream->next_layer().cancel();
      stream->async_shutdown([](const asio::error_code& /*error*/){});
      run(io_context,std::chrono::seconds(5),stream);
      stream->lowest_layer().close();
      delete stream;
      stream = nullptr;
    }
  }catch (std::exception& /*e*/){
    if(stream != nullptr){
      delete stream;
      stream = nullptr;
    }
  }  
}

  //*********************** SSL END *************************//
#else
pen_tcp::server_seq::server_seq(const unsigned verb) :
  acceptor(io_context), socket(nullptr),
  verbose(verb), sslEnabled(false){
      
  tinit = std::chrono::steady_clock::now();
  timeout = std::chrono::seconds(5);      
}
#endif

int pen_tcp::server_seq::init(const bool local,
			      const unsigned short port){

  asio::error_code ec;

  //Close acceptor, just in case
  acceptor.close(ec);
  
  if(filterLog(1)){
    fprintf(flog,"%07ld s - Starting server monitor listening at port %u\n",
	    static_cast<long int>(timeStamp()),port);
    fflush(flog);
  }

  sslEnabled = false;
    
  asio::ip::tcp::endpoint endpoint;
  if(local)
    endpoint = asio::ip::tcp::endpoint(asio::ip::address::from_string("127.0.0.1"),port);
  else
    endpoint = asio::ip::tcp::endpoint(asio::ip::tcp::v4(),port);
  
  //Configure the acceptor

  //Open acceptor using tcp protocol and IPv4 
  acceptor.open(asio::ip::tcp::v4(), ec);
  if(ec){
    if(filterLog(0)){
      fprintf(flog,"%07ld s - Error at acceptor open\n"
	      "                   Error: %s\n"
	      "              Error code: %d\n",
	      static_cast<long int>(timeStamp()),
	      ec.message().c_str(),ec.value());
      fflush(flog);
    }
    return PEN_TCP_ACCEPTOR_OPEN_FAIL;
  }

  //Bind the acceptor
  acceptor.bind(endpoint,ec);
  if(ec){
    if(filterLog(0)){
      fprintf(flog,"%07ld s - Error binding to specified local endpoint\n"
	      "                   Error: %s\n"
	      "              Error code: %d\n",
	      static_cast<long int>(timeStamp()),
	      ec.message().c_str(),ec.value());
      fflush(flog);
    }
    return PEN_TCP_ERROR_ON_BIND;    
  }

  acceptor.listen(asio::socket_base::max_listen_connections, ec);
  if(ec){
    if(filterLog(0)){
      fprintf(flog,"%07ld s - Error on acceptor listen\n"
	      "                   Error: %s\n"
	      "              Error code: %d\n",
	      static_cast<long int>(timeStamp()),
	      ec.message().c_str(),ec.value());
      fflush(flog);
    }
    return PEN_TCP_ACCEPTOR_LISTEN_FAIL;    
  }

  return PEN_TCP_SUCCESS;
}

void pen_tcp::run(asio::io_context& io_context,
		  std::chrono::seconds timeout,
		  asio::ip::tcp::socket*& socket,
		  FILE* flog, const unsigned verbose){
  //Restart context
  io_context.restart();
  
  //Run context with the timeout
  io_context.run_for(timeout);
  
  //Check if context has been stopped
  if(!io_context.stopped()){
    //Context has not been stoped, the timeout has reached
    if(filterLog(flog,verbose,0)){
      fprintf(flog,"Timeout on raw TCP action\n");
      fflush(flog);
    }

    if(socket != nullptr){
      std::error_code ec;
      //Close the connection
      socket->close(ec); //close TCP
      if(ec){
	if(filterLog(flog,verbose,0)){
	  fprintf(flog,"Unable to shutdown TCP"
		  "              Error: %s\n"
		  "         Error code: %d\n",
		  ec.message().c_str(),ec.value());
	  fflush(flog);
	}
      }
    }

    //Run context to finish the operation
    io_context.run();      
  }
}

int pen_tcp::write(char message[messageLen], long int len,
		   asio::io_context& io_context,
		   std::chrono::seconds timeout,
		   asio::ip::tcp::socket*& socket,
		   FILE* flog, const unsigned verbose){
  try{

    if(socket == nullptr)
      return PEN_TCP_INACTIVE_SOCKET;
    
    fillmessage(message,len);
    
    //Add the write handler
    asio::error_code ec;
      
    asio::async_write(*socket,asio::buffer(message,messageLen),
		      [&ec](const std::error_code& error,
			    std::size_t /*result_n*/){
			ec = error;
		      });
    run(io_context,timeout,socket,flog,verbose);
    

    if(ec){
      if(filterLog(flog,verbose,0)){
	fprintf(flog,"Error: Unable to send message\n"
		"              Error: %s\n"
		"         Error code: %d\n",
		ec.message().c_str(),ec.value());
	fflush(flog);
      }
      return PEN_TCP_WRITE_FAIL;
    }
    else if(filterLog(flog,verbose,3)){
      fprintf(flog,"Message sent\n");
      fflush(flog);
    }    
    return PEN_TCP_SUCCESS;    
  }
  catch (std::exception& e){
    if(filterLog(flog,verbose,0)){
      fprintf(flog,"Exception on message send: %s\n",
	      e.what());
      fflush(flog);
    }
    return PEN_TCP_EXCEPTION_CATCH;
  }
}

int pen_tcp::listen(asio::io_context& io_context,
		    asio::ip::tcp::acceptor& acceptor,
		    asio::ip::tcp::socket*& socket,
		    FILE* flog, const unsigned verbose){

  try{

    if(socket != nullptr)
      return PEN_TCP_SOCKET_ALREADY_ACTIVE;
    
    //Create a new socket
    socket = new asio::ip::tcp::socket(io_context);
    
    asio::error_code ec;
  
    acceptor.accept(*socket,ec);
    if(ec){
      if(filterLog(flog,verbose,0)){
	fprintf(flog,"Connection accept failed\n"
		"              Error: %s\n"
		"         Error code: %d\n",
		ec.message().c_str(),ec.value());
	fflush(flog);
      }
      return PEN_TCP_ACCEPT_CONNECTION_FAIL;
    }

    if(filterLog(flog,verbose,1)){
      fprintf(flog,"Incoming communication from %s\n",
	      socket->remote_endpoint().address().to_string().c_str());
      fflush(flog);
    }
    
    return PEN_TCP_SUCCESS;
  }catch (std::exception& e){
    if(filterLog(flog,verbose,0)){
      fprintf(flog,"Exception on socket receive: %s\n",
	      e.what());
      fflush(flog);
    }
    return PEN_TCP_EXCEPTION_CATCH;
  }
}

int pen_tcp::connect(asio::io_context& io_context,
		     const std::vector<asio::ip::tcp::endpoint>& endpoints,
		     asio::ip::tcp::socket*& socket,
		     std::chrono::seconds timeout,
		     FILE* flog, const unsigned verbose){

  if(socket != nullptr)
    return PEN_TCP_SOCKET_ALREADY_ACTIVE;

  //Create a new socket
  socket = new asio::ip::tcp::socket(io_context);
  
  std::error_code ec;
  asio::async_connect(*socket, endpoints,
		      [&ec](const std::error_code& error,
			    const asio::ip::tcp::endpoint& /*endpoint*/){
			ec = error;
		      });

  run(io_context,timeout,socket,flog,verbose);
  
  if(ec){
    if(filterLog(flog,verbose,0)){
      fprintf(flog,"Connection failed\n"
	      "              Error: %s\n"
	      "         Error code: %d\n",
	      ec.message().c_str(),ec.value());
      fflush(flog);
    }
    return PEN_TCP_CONNECTION_FAIL;
  }

  return PEN_TCP_SUCCESS;
}

int pen_tcp::receive(std::string& message,
		     asio::io_context& io_context,
		     asio::ip::tcp::socket*& socket,
		     std::chrono::seconds timeout,
		     FILE* flog, const unsigned verbose){

  if(socket == nullptr){
    if(filterLog(flog,verbose,0)){
      fprintf(flog,"Invalid socket (null pointer)\n");
      fflush(flog);
    }
    return PEN_TCP_INACTIVE_SOCKET;
  }
  asio::error_code ec;
  //Read the message
  message.clear();
  size_t read = 0;
  for(;;){
    char buf[messageLen+1];
    size_t len;
    socket->async_read_some(asio::buffer(buf),
			    [&ec,&len](const asio::error_code& error,
				       std::size_t bytes_transferred){
			      ec = error;
			      len = bytes_transferred;
			    });
      
    //Start timeout
    run(io_context,timeout,socket,flog,verbose);
      
    if(ec == asio::error::eof){
      break;
    }
    else if(ec){
      if(filterLog(flog,verbose,0)){
	fprintf(flog,"Socket read failed\n"
		"               Error: %s\n"
		"          Error code: %d\n",
		ec.message().c_str(),ec.value());
	fflush(flog);
      }
      return PEN_TCP_SOCKET_READ_FAILED;
    }
    read += len;
    message.append(buf,len);
    if(read >= messageLen){
      message = message.substr(0,messageLen);
      break;
    }
  }

  //trim message
  message = trim(message);

  return PEN_TCP_SUCCESS;
}

void pen_tcp::close(asio::ip::tcp::socket*& socket){
  try{
    if(socket != nullptr){
      socket->close();
      delete socket;
      socket = nullptr;
    }
  }catch (std::exception& /*e*/){
    if(socket != nullptr){
      delete socket;
      socket = nullptr;
    }
  }
}

std::string pen_tcp::trim(const std::string& str,
			  const std::string& whitespace){
    const auto strBegin = str.find_first_not_of(whitespace);
    if(strBegin == std::string::npos)
      return ""; // no content

    const auto strEnd = str.find_last_not_of(whitespace);
    const auto strRange = strEnd - strBegin + 1;

    return str.substr(strBegin, strRange);
}

void pen_tcp::server_seq::clear(){
#ifdef _PEN_USE_SSL_
  pen_tcp::close(io_context,stream);
#endif
  pen_tcp::close(socket);
  asio::error_code ec;
  acceptor.close(ec);
  flog = nullptr;
}

  //*********************** SSL *************************//
#ifdef _PEN_USE_SSL_
pen_tcp::client::client(const unsigned verb) :
  ssl_context(asio::ssl::context::sslv23),
  stream(nullptr), socket(nullptr),
  verbose(verb), sslEnabled(false){

  timeout = std::chrono::seconds(20);      
}

int pen_tcp::client::initSSL(const char* CAfilename,
			     const char* certFilename,
			     const char* keyFilename,
			     const char* password,
			     const char* hostname){

  asio::error_code ec;
  
  sslEnabled = true;
  
  if(CAfilename == nullptr ||
     certFilename == nullptr ||
     keyFilename == nullptr){

    if(filterLog(0)){
      fprintf(flog,"Client SSL configuration error: "
	      "Null pointer provided as input\n");
      fflush(flog);
    }
    return PEN_TCP_NULL_POINTER;
  }

  if(filterLog(1)){
    fprintf(flog,"clien SSL enabled\n");
    fflush(flog);
  }

  //Configure SSL
  std::string spassword;
  if(password != nullptr)
    spassword.assign(password);
      
  ssl_context.set_options(asio::ssl::context::default_workarounds |
			  asio::ssl::context::no_sslv2,
			  ec);      
  if(ec){
    if(filterLog(0)){
      fprintf(flog,"Error: Unable to set SSL context options\n"
	      "                 Error: %s\n"
	      "            Error code: %d\n",
	      ec.message().c_str(),ec.value());
      fflush(flog);
    }
    return PEN_TCP_BAD_SSL_CONTEXT_OPTIONS;
  }

  ssl_context.set_verify_mode(asio::ssl::verify_peer |
			      asio::ssl::verify_fail_if_no_peer_cert,
			      ec);

  if(ec){
    if(filterLog(0)){
      fprintf(flog,"Error: Unable to set SSL context verify mode\n"
	      "                 Error: %s\n"
	      "            Error code: %d\n",
	      ec.message().c_str(),ec.value());
      fflush(flog);
    }
    return PEN_TCP_BAD_SSL_CONTEXT_OPTIONS;
  }

  if(password != nullptr){
    ssl_context.set_password_callback([spassword](std::size_t max_length,
						  asio::ssl::context::password_purpose /*purpose*/) -> std::string{
					if(spassword.length() > max_length){
					  printf("Warning: provided password"
						 "length exceeds the limit.\n"
						 "   length: %lu\n"
						 "      max: %lu\n",
						 static_cast<unsigned long>
						 (spassword.length()),
						 static_cast<unsigned long>
						 (max_length));
					  return spassword.substr(0,max_length);
					}
					else
					  return spassword;
				      },ec);

    if(ec){
      if(filterLog(0)){
	fprintf(flog,"Error: Unable to set password callback\n"
		"                   Error: %s\n"
		"              Error code: %d\n",
		ec.message().c_str(),ec.value());
	fflush(flog);
      }
      return PEN_TCP_BAD_SSL_CONTEXT_OPTIONS;
    }
  }

  //Load CA certificate
  ssl_context.load_verify_file(CAfilename,ec);
  if(ec){
    if(filterLog(0)){
      fprintf(flog,"Error loading CA certificate\n"
	      "              Error: %s\n"
	      "         Error code: %d\n",
	      ec.message().c_str(),ec.value());
      fflush(flog);
    }
    return PEN_TCP_INVALID_CA_CERT;
  }

  //Load client certificate file
  ssl_context.use_certificate_chain_file(certFilename,ec);
  if(ec){
    if(filterLog(0)){
      fprintf(flog,"Error loading SSL certificate chain file\n"
	      "              Error: %s\n"
	      "         Error code: %d\n",
	      ec.message().c_str(),ec.value());
      fflush(flog);
    }
    return PEN_TCP_INVALID_CERT;
  }

  //Load client key file
  ssl_context.use_private_key_file(keyFilename,
				   asio::ssl::context::pem,
				   ec);
  if(ec){
    if(filterLog(0)){
      fprintf(flog,"Error loading SSL key file\n"
	      "              Error: %s\n"
	      "         Error code: %d\n",
	      ec.message().c_str(),ec.value());
      fflush(flog);
    }
    return PEN_TCP_INVALID_KEY;
  }      

  //Set verify callback
  std::string shostname;
  if(hostname != nullptr)
    shostname.assign(hostname);
  ssl_context.set_verify_callback([this,shostname](bool preverified,
						   asio::ssl::verify_context& ctx){
				    char subject_name[256];
				    X509* cert = X509_STORE_CTX_get_current_cert(ctx.native_handle());
				    X509_NAME_oneline(X509_get_subject_name(cert), subject_name, 256);
				    if(filterLog(3)){
				      fprintf(flog,"Server name: %s\n",
					      subject_name);
				      fflush(flog);
				    }
				    //Check if coincide with expected hostname
				    if(shostname.length() > 0){
				      if(shostname.compare(subject_name) != 0){
					if(filterLog(0)){
					  fprintf(flog,"Warning: Unexpected "
						  "hostname (%s), close "
						  "connection\n",
						  subject_name);
					  fflush(flog);
					}
					preverified = false;
				      }
				    }
				    return preverified;					
				  },ec);
  if(ec){
    if(filterLog(0)){
      fprintf(flog,"Error creating verify callback\n"
	      "              Error: %s\n"
	      "         Error code: %d\n",
	      ec.message().c_str(),ec.value());
      fflush(flog);
    }
    return PEN_TCP_CALLBACK_SET_FAIL;
  }

  return PEN_TCP_SUCCESS;
}
#else
  //*********************** SSL END *************************//
pen_tcp::client::client(const unsigned verb) :
  socket(nullptr), verbose(verb), sslEnabled(false){

  timeout = std::chrono::seconds(20);      
}
#endif

int pen_tcp::client::setHost(const char* host,
			     const char* port){
  if(host == nullptr){
    if(filterLog(0)){
      fprintf(flog,"client: init: Error: Provided null pointer as hostname\n");
      fflush(flog);
    }
    return PEN_TCP_NULL_POINTER;
  }

  if(port == nullptr){
    if(filterLog(0)){
      fprintf(flog,"client: init: Error: Provided null pointer as port\n");
      fflush(flog);
    }
    return PEN_TCP_NULL_POINTER;
  }

  asio::error_code ec;
  
  if(filterLog(1)){
    fprintf(flog,"Configuring host located at %s:%s\n",
	    host,port);
    fflush(flog);
  }

  //Clear previous endpoints
  serverEndpoints.clear();
#ifndef _PEN_BUILD_STATIC_
  //Static build is dissable

  //Resolve the hostname to a list of endpoints
  asio::ip::tcp::resolver resolver(io_context);
  asio::ip::tcp::resolver::results_type list = resolver.resolve(host,port,ec);
  //Copy endpoints
  for(const auto& endpoint : list)
    serverEndpoints.push_back(endpoint);
#else
  //Static build is neabled. As "resolve" function calls "getaddrinfo",
  //it can't be statically linked if we are using glibc. To avoid that
  //call we will require the host IP instead of a hostname
  auto adress = asio::ip::make_address(host,ec);
  if(!ec){
    serverEndpoints.push_back(asio::ip::tcp::endpoint(adress,std::atoi(port)));
  }
#endif
  if(ec){
    if(filterLog(0)){
      fprintf(flog,"Error resolving hostname\n"
	      "           Error: %s\n"
	      "           Error code: %d\n",
	      ec.message().c_str(),ec.value());
      fflush(flog);
    }
    return PEN_TCP_UNABLE_TO_RESOLVE_HOSTNAME;
  }
  
  return PEN_TCP_SUCCESS;
}

void pen_tcp::client::clear(){
#ifdef _PEN_USE_SSL_
  pen_tcp::close(io_context,stream);
#endif
  pen_tcp::close(socket);
  flog = nullptr;
}
