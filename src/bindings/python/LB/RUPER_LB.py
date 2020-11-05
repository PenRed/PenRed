import ctypes

lib = ctypes.cdll.LoadLibrary('./libRUPER_LB_C.so')

class Server():
    def __init__(self):

        #Constructor
        #lib.taskServer_new.argtypes = [ctypes.c_void_p]
        lib.taskServer_new.restype  = ctypes.c_void_p

        #Set init function
        lib.taskServer_init.argtypes = [ctypes.c_void_p,    #server
                                        ctypes.c_ulong,     #nw
                                        ctypes.c_ulonglong, #nIter
                                        ctypes.c_char_p,    #logfilename
                                        ctypes.c_uint]      #verbose
        lib.taskServer_init.restype  =  ctypes.c_int
        
        #Set ETA function
        lib.taskServer_ETA.argtypes = [ctypes.c_void_p]    #server
        lib.taskServer_ETA.restype  =  ctypes.c_ulonglong
        
        #Set threshold function
        lib.taskServer_setThreshold.argtypes = [ctypes.c_void_p,    #server
                                                ctypes.c_ulonglong] #threshold
        lib.taskServer_setThreshold.restype  =  ctypes.c_void_p
        
        #Set log file function
        lib.taskServer_setLogFile.argtypes = [ctypes.c_void_p,  #server
                                              ctypes.c_char_p]  #log filename
        lib.taskServer_setLogFile.restype  =  ctypes.c_int
        
        #Set receive report function
        lib.taskServer_receiveReport.argtypes = [ctypes.c_void_p,    #server
                                                 ctypes.c_ulong,     #iw
                                                 ctypes.c_ulonglong, #nIter
                                                 ctypes.c_long,      #elapsed
                                                 ctypes.POINTER(ctypes.c_ulonglong), #newAssign
                                                 ctypes.c_uint]      #verbose
        lib.taskServer_receiveReport.restype  =  ctypes.c_int

        #Set receive start function
        lib.taskServer_receiveStart.argtypes = [ctypes.c_void_p,    #server
                                                ctypes.c_ulong,     #iw
                                                ctypes.c_long,      #elapsed
                                                ctypes.POINTER(ctypes.c_ulonglong), #assigned                                                
                                                ctypes.c_uint]      #verbose
        lib.taskServer_receiveStart.restype  =  ctypes.c_int        

        #Set receive finish function
        lib.taskServer_receiveFinish.argtypes = [ctypes.c_void_p,    #server
                                                 ctypes.c_ulong,     #iw
                                                 ctypes.c_ulonglong, #nIter
                                                 ctypes.c_long,      #elapsed
                                                 ctypes.c_uint]      #verbose
        lib.taskServer_receiveFinish.restype  =  ctypes.c_int

        #Set monitor function
        lib.taskServer_monitor.argtypes = [ctypes.c_void_p,  #server
                                           ctypes.c_bool,    #local
                                           ctypes.c_ushort,  #port
                                           ctypes.c_char_p,  #CAfilename
                                           ctypes.c_char_p,  #certfilename
                                           ctypes.c_char_p,  #keyfilename
                                           ctypes.c_char_p,  #dhfilename
                                           ctypes.c_char_p,  #password
                                           ctypes.c_uint]    #verbose
        lib.taskServer_monitor.restype  =  ctypes.c_void_p

        #Set speed function
        lib.taskServer_speed.argtypes = [ctypes.c_void_p, #server
                                         ctypes.POINTER(ctypes.c_float)]  #globspeed
        lib.taskServer_speed.restype  =  ctypes.c_ulonglong

        #Set print report function
        lib.taskServer_printReport.argtypes = [ctypes.c_void_p,  #server
                                               ctypes.c_char_p]  #filename
        lib.taskServer_printReport.restype  =  ctypes.c_int

        #Set set init function
        lib.taskServer_setInit.argtypes = [ctypes.c_void_p]  #server
        lib.taskServer_setInit.restype  =  ctypes.c_void_p

        #Set set init function
        lib.taskServer_moveInit.argtypes = [ctypes.c_void_p,  #server
                                            ctypes.c_long]    #elapsed
        lib.taskServer_moveInit.restype  =  ctypes.c_void_p        
        
        #Set log file function
        lib.taskServer_setID.argtypes = [ctypes.c_void_p,  #server
                                         ctypes.c_char_p]  #newID
        lib.taskServer_setID.restype  =  ctypes.c_int
        
        self.obj = lib.taskServer_new()

    def init(self,nw,nIter,logFilename,verbose):
        logFilenameC = logFilename.encode('utf-8')
        return lib.taskServer_init(self.obj,nw,nIter,
                                   logFilenameC,verbose)
        
    def ETA(self):
        return lib.taskServer_ETA(self.obj)

    def setThreshold(self,t):
        lib.taskServer_setThreshold(self.obj,t)

    def setLogFile(self,filename):
        logFilenameC = filename.encode('utf-8')
        return lib.taskServer_setLogFile(self.obj,logFilenameC)

    def receiveReport(self,iw,nIter,elapsed,verbose):

        auxAssign = ctypes.c_ulonglong(0)
        ret = lib.taskServer_receiveReport(self.obj,iw,nIter,elapsed,
                                           ctypes.byref(auxAssign),verbose)
        return [ret,auxAssign.value]
        
    def receiveStart(self,iw,elapsed,verbose):
        
        auxAssign = ctypes.c_ulonglong(0)
        ret = lib.taskServer_receiveStart(self.obj,iw,elapsed,
                                          ctypes.byref(auxAssign),verbose)
        return [ret,auxAssign.value]

    def receiveFinish(self,iw,nIter,elapsed,verbose):
        return lib.taskServer_receiveFinish(self.obj,iw,nIter,elapsed,verbose)

    def monitor(self,local,port,CAfilename,certFilename,keyFilename,
                dhFilename,password,verbose):
        if CAfilename is None:
            CAfilenameC = None
        else:
            CAfilenameC = CAfilename.encode('utf-8')
            
        if certFilename is None:
            certFilenameC = None
        else:
            certFilenameC = certFilename.encode('utf-8')
            
        if keyFilename is None:
            keyFilenameC = None
        else:
            keyFilenameC = keyFilename.encode('utf-8')
            
        if dhFilename is None:
            dhFilenameC = None
        else:
            dhFilenameC = dhFilename.encode('utf-8')
            
        if password is None:
            passwordC = None
        else:
            passwordC = password.encode('utf-8')
            
        lib.taskServer_monitor(self.obj,local,port,CAfilenameC,certFilenameC,
                           keyFilenameC,dhFilenameC,passwordC,verbose)

    def speed(self):
        auxSpeed = ctypes.c_float(0.0)
        done = lib.taskServer_speed(self.obj,ctypes.byref(auxSpeed))
        return [done,auxSpeed.value]

    def printReport(self,filename):
        if filename is None:
            filenameC = None
        else:
            filenameC = filename.encode('utf-8')

        return lib.taskServer_printReport(self.obj,filenameC)
        
    def setInit(self):
        lib.taskServer_setInit(self.obj)
        
    def moveInit(self,elapsed):
        lib.taskServer_moveInit(self.obj,elapsed)
        
    def setID(self,newID):
        newIDC = newID.encode('utf-8')
        lib.taskServer_setID(self.obj,newIDC)
