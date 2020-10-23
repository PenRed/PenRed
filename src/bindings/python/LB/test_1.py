import sys
from RUPER_LB import Server

#Number of workers
nw = 2
#Number of iterations
nIter = 10000000000
#Log filename
logfilename = 'logs.txt'
#Threshold
threshold = 100

#Create server
server = Server()

#Init server
err = server.init(nw,nIter,logfilename,2)
if err != 0:
	print("Unable to init server. Error code: %d" % err)
	sys.exit()
#Set threshold
server.setThreshold(threshold)

#Start one worker
err,assignw0 = server.receiveStart(0,10,2)
if err != 0:
	print("Error on worker 0 start: %d" % err)
	sys.exit()
else:
	print("Assigned to worker 0: %d" % assignw0)

#Simulate server elapsed time
server.moveInit(-25)

#Perform a report
err, assignw0 = server.receiveReport(0,1000,30,2)
if err != 0:
	print("Error on worker report: %d" % err)
	sys.exit()

#Elapse time
server.moveInit(-30)

#Perform a second report
err, assignw0 = server.receiveReport(0,30000,55,2)
if err != 0:
	print("Error on worker second report: %d" % err)
	sys.exit()

#Start worker 1
err,assignw1 = server.receiveStart(1,5,2)
if err != 0:
        print("Error on worker 1 start: %d" % err)
        sys.exit()
else:
        print("Assigned to worker 1: %d" % assignw1)


#Perform a report on both workers
server.moveInit(-30)
server.receiveReport(1,40000,30,2)
server.receiveReport(0,45000,80,2)

server.moveInit(-30)
err, assignw0 = server.receiveReport(1,1000000,70,2)
err, assignw1 = server.receiveReport(0,1300000,130,2)

print("Report assigns: %d %d" % (assignw0,assignw1))

#Send a finish petition on worker 0

server.moveInit(-40)
server.receiveFinish(0,2000000,160,2)
server.receiveReport(1,2400000,165,2)

#Finish the second worker
server.moveInit(-40)
server.receiveFinish(0,3000000,185,2)
server.receiveFinish(1,5000000,190,2)

#Start a new worker
server.moveInit(-5)
server.receiveStart(4,10,2)

server.moveInit(-30)
server.receiveReport(4,300000,35,2)

#Start worker 3 and 5
server.receiveStart(3,50,2)
server.receiveStart(5,10,2)

server.receiveReport(3,20000000,30,2)

server.moveInit(-30)
server.receiveReport(4,7000000,65,2)
server.receiveReport(5,4000000,30,2)

#Perform some reports
server.moveInit(-60)
server.receiveReport(4,10000000,120,2)
server.receiveReport(5,10500000,95,2)

server.moveInit(-60)
server.receiveReport(4,17000000,185,2)
server.receiveReport(5,18000000,160,2)

server.moveInit(-70)
server.receiveReport(4,24000000,250,2)
server.receiveReport(5,23000000,230,2)

server.moveInit(-70)
server.receiveReport(4,32000000,310,2)
server.receiveReport(5,34000000,295,2)

server.moveInit(-70)
server.receiveReport(4,37000000,380,2)
server.receiveReport(5,40000000,370,2)

server.moveInit(-70)
server.receiveReport(4,46000000,460,2)
server.receiveReport(5,50000000,450,2)

server.moveInit(-70)
server.receiveReport(4,56000000,550,2)
server.receiveReport(5,60000000,540,2)
