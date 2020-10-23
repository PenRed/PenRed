 
import socket
import sys
import time

if len(sys.argv) < 3:
    print("usage: python %s address port message" % sys.argv[0])
    sys.exit(0)

print("Send '%s' to %s:%d" % (sys.argv[3],sys.argv[1],int(sys.argv[2])))

sock = socket.socket(socket.AF_INET,socket.SOCK_STREAM)
server_address = (sys.argv[1],int(sys.argv[2]))
sock.connect(server_address)

length = len(sys.argv[3])
if length >= 512:
    message = sys.argv[3][0:512]
else:
    remaining = 512-length
    message = sys.argv[3] + " " * remaining

toSend = message.encode('utf-8')
sock.sendall(toSend)
print("Waiting for response")

data = sock.recv(512)
response = data.decode("utf-8")
print(response)

sock.close()
