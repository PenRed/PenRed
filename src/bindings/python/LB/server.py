import sys
from RUPER_LB import Server

if len(sys.argv) < 6:
	print("usage: %s ID procedence port nworkers nIter CA-cert-file "
           "server-cert-file server-key-file dh-params-file key-password(optional)\n\n"
	   "procedence must be set to 0 to accept only connections "
	   "from localhost. Otherwise set it to 1\n\n"
       "NOTE: CA-cert-file, server-cert-file, server-key-file and dh-params-file "
       "are optional parameters. If ALL of them are provided, the server will "
       "enable SSL comunications. Otherwise the comunications will be performed "
       "with no security.\n",sys.argv[0])
	sys.exit()

if sys.argv[2] == 0:
	local = True
else:
	local = False

port = int(sys.argv[3])
if port < 0:
	print("Invalid port: %d" % port)
	sys.exit(1)

nw = int(sys.argv[4])
if nw < 1:
	print("Invalid number of workers: %d" % nw)
	sys.exit(2)

nIter = int(sys.argv[5])
if nIter < 1:
	print("Invalid number of iterations: %d" %nIter)
	sys.exit(3)

if len(sys.argv) >= 10:
	caCertfile = str(sys.argv[6])
	serverCertfile = str(sys.argv[7])
	serverKeyfile = str(sys.argv[8])
	DHparamfile = str(sys.argv[9])
	if len(sys.argv) > 10:
		keyPassword = str(sys.argv[10])
	else:
		keyPassword = None
else:
	caCertfile = None
	serverCertfile = None
	serverKeyfile = None
	DHparamfile = None
	keyPassword = None

#Create server
server = Server()
server.setID(str(sys.argv[1]))
server.init(nw,nIter,"serverLog.txt",3)

#Start monitor
server.monitor(local,port,caCertfile,serverCertfile,serverKeyfile,DHparamfile,keyPassword,3)
