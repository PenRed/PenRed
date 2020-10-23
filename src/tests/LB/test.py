from subprocess import Popen, PIPE, STDOUT
import sys
import time

print(sys.argv)

p = Popen(['./stdinServer'] + sys.argv[1:], stdout=PIPE, stdin=PIPE)

#Start worker 1
p.stdin.write(b'2 1 10\n') 
p.stdin.flush()

#Start worker 0
p.stdin.write(b'2 0 15\n') 
p.stdin.flush()

#Elapse time
time.sleep(5)

#Report worker 0
p.stdin.write(b'1 0 1000 16\n')
p.stdin.flush()

#Elapse time
time.sleep(3)

#Report both workers
p.stdin.write(b'1 0 2400 21\n')
p.stdin.flush()
p.stdin.write(b'1 1 1000 20\n')
p.stdin.flush()


#Elapse time
time.sleep(5)

#Report both workers
p.stdin.write(b'1 0 3600 25\n')
p.stdin.flush()
p.stdin.write(b'1 1 2000 24\n')
p.stdin.flush()

#Elapse time
time.sleep(7)

#Start a new worker
p.stdin.write(b'2 2 4\n')
p.stdin.flush()

#Report both workers
p.stdin.write(b'1 0 5600 32\n')
p.stdin.flush()
p.stdin.write(b'1 1 3000 30\n') 
p.stdin.flush()

#Elapse time
time.sleep(5)

#Report workers 0 2
p.stdin.write(b'1 0 6600 37\n')
p.stdin.flush()
p.stdin.write(b'1 2 1000 10\n') 
p.stdin.flush()

#Elapse time
time.sleep(5)

#Report workers 0 2
p.stdin.write(b'1 0 7500 41\n')
p.stdin.flush()
p.stdin.write(b'1 2 1800 16\n') 
p.stdin.flush()

#Elapse time
time.sleep(10)

#Report workers 0 2
p.stdin.write(b'1 0 9600 52\n')
p.stdin.flush()
p.stdin.write(b'1 2 2700 25\n') 
p.stdin.flush()

#Elapse time
time.sleep(10)

#Report workers 2
p.stdin.write(b'1 2 3600 35\n') 
p.stdin.flush()

#Elapse time
time.sleep(10)

#Finish worker 0
p.stdin.write(b'3 0 12400 72\n')
p.stdin.flush()

#Report workers 2
p.stdin.write(b'1 2 4800 46\n') 
p.stdin.flush()

#Dump server state
p.stdin.write(b'5 dump.dat\n') 
p.stdin.flush()

#Elapse time
time.sleep(8)

#Report worker 2
p.stdin.write(b'1 2 5300 55\n')
p.stdin.flush()

#Elapse time
time.sleep(5)

#Start a new worker
p.stdin.write(b'2 6 3\n')
p.stdin.flush()

#Report worker 2
p.stdin.write(b'1 2 6000 60\n')
p.stdin.flush()

#Dump server state
p.stdin.write(b'5 finalDump.dat\n') 
p.stdin.flush()

#Close server
p.stdin.write(b'0\n')
p.stdin.flush()
stdout_data = p.stdout.read()
print(stdout_data.decode('utf-8'))
