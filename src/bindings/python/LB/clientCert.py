 
import socket
import sys
import ssl
import time

if len(sys.argv) < 3:
    print("usage: python %s address port message" % sys.argv[0])
    sys.exit(0)

print("Send '%s' to %s:%d" % (sys.argv[3],sys.argv[1],int(sys.argv[2])))

context = ssl.SSLContext(ssl.PROTOCOL_TLS_CLIENT)
context.load_verify_locations('CA/ca1.cert')
context.load_cert_chain("client.pem",password="1234")

sock = socket.socket(socket.AF_INET,socket.SOCK_STREAM)
sslSock = context.wrap_socket(sock,server_hostname="penred-server")
server_address = (sys.argv[1],int(sys.argv[2]))
sslSock.connect(server_address)

length = len(sys.argv[3])
if length >= 512:
    message = sys.argv[3][0:512]
else:
    remaining = 512-length
    message = sys.argv[3] + " " * remaining

toSend = message.encode('utf-8')
sslSock.sendall(toSend)
print("Waiting for response")

data = sslSock.recv(512)
response = data.decode("utf-8")
print(response)

sslSock.close()
