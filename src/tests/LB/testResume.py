from subprocess import Popen, PIPE, STDOUT
import sys
import time

print(sys.argv)

p = Popen(['./stdinServer'] + sys.argv[1:], stdout=PIPE, stdin=PIPE)

#Load dumped server state
p.stdin.write(b'4 dump.dat\n') 
p.stdin.flush()

#Dump server state
p.stdin.write(b'5 dump2.dat\n') 
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
p.stdin.write(b'5 finalDump2.dat\n') 
p.stdin.flush()

#Close server
p.stdin.write(b'0\n')
p.stdin.flush()
stdout_data = p.stdout.read()
print(stdout_data.decode('utf-8'))
