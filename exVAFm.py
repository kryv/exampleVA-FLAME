import sys
import math
import numpy 
import time
import cPickle
import re
import fnmatch

home_dir = './build'
python_dir  = home_dir + '/python'
file_name = python_dir + '/flame/test/to_strl.lat'

sys.path.append(python_dir)

from flame import Machine

PS_X = 0; PS_PX = 1; PS_Y = 2; PS_PY = 3; PS_S = 4; PS_PS = 5
MeVtoeV = 1e6

class VAF:
    """
    Contains:
    getelem(int)
    getindex(str)
    getindexu(str)
    setelem(int,str,float)
    prt_lat_prms()
    prt_lat()
    tcs()
    tcs(int,int)
    looptcs()
    SaveBeam(S)
    """
    def __init__(self, fname=file_name):
        self.name = fname
        self.itr = 0

        with open(self.name, 'rb') as inf:

            self.M = Machine(inf)

            self.lat = self.M.conf()['elements']

            #Flag for limited longitudinal run
            self.clng = 0

            S=self.M.allocState({})
            self.M.propagate(S,0,1)

            self.refIonZ = S.ref_IonZ
            self.IonZ_io = S.IonZ# list of real IonZ

            self.refIonEk = S.ref_IonEk

            self.BC0 = S.moment0  
            self.ENV0 = S.moment1

            # memory space for all beam data
            # self.LD = [[0.0]*8 for i in range(len(self.M))]
            self.LD = numpy.zeros((len(self.M),8))

            # store element position data
            self.bpmls=[]
            self.corls=[]
            self.cavls=[]
            self.solls=[]
            self.cavpara=[]
            self.solpara=[]

            for i in range(len(self.M)):
                elem = self.M.conf()['elements'][i]
                if elem['type'] == 'bpm' : self.bpmls.append(i)
                elif elem['type'] == 'orbtrim' : self.corls.append(i)
                elif elem['type'] == 'rfcavity' :
                    self.cavls.append(i)
                    self.cavpara.append([elem['scl_fac'],elem['phi']])
                elif elem['type'] == 'solenoid' :
                    self.solls.append(i)
                    self.solpara.append([elem['B']])

            # output data for plotting
            with open('plot.info','w') as f:
                iteminfo=[len(self.M),self.bpmls,self.corls,self.cavls,self.solls]
                cPickle.dump(iteminfo,f)


    def getelem(self,num):
        """
        Get parameter of lattice element
        getelem(index of lattice element)
        """
        print self.M.conf()['elements'][num]


    def getvalue(self,num,name):
        """
        Get parameter of lattice element
        getelem(index of lattice element, parameter name)
        """
        print self.M.conf(num)[naem]


    def getindex(self,name,searchby='name'):
        """
        Get index list of lattice elements by python style regular expression
        getindex(name or type,
                 searchby = 'name')        
        """
        name = name.replace(':','_').lower()
        pat = re.compile(name)
        result = []

        for (i,elem) in enumerate(self.lat):
            if pat.search(elem[searchby]):
                result.append(i)
        return result
    
    def getindexu(self,name,searchby='name'):
        """
        Get index list of lattice elements by unix style regular expression
        getindex(name or type,
                 searchby = 'name')        
        """
        name = name.replace(':','_').lower()
        result = []

        for (i,elem) in enumerate(self.lat):
            if fnmatch.fnmatch(elem[searchby],name):
                result.append(i)
        return result


    def setelem(self,num,name,val):
        """
        Set parameter of lattice element
        setelem(position number of lattice element, 
                name of parameter
                value of parameter)
        """
        #D = self.M.conf()['elements'][num]
        #D = self.lat[num]
        #D[name] = float(val)
        self.M.reconfigure(num,{name:float(val)})


    """
    # Developing
    def genRandomNoise(self):
        for (i,cv) in enumerate(self.cavls):
            D = self.M.conf()['elements'][cv]
            for (j,st) in enumerate(['scl_fac','phi']):
                defv = self.cavpara[i][j]
                D[st] = defv + numpy.random.normal(0.0,abs(defv)*1e-3)
                self.M.reconfigure(cv, D)
            D.clear()

        for (i,cv) in enumerate(self.solls):
            D = self.M.conf()['elements'][cv]
            for (j,st) in enumerate(['B']):
                defv = self.solpara[i][j]
                D[st] = defv + numpy.random.normal(0.0,abs(defv)*1e-3)
                self.M.reconfigure(cv, D)
            D.clear()
    
    def loopRN(self):
        while self.itr < 100:
            self.genRandomNoise()
            self.itr += 1

    """

    
    def prt_lat_prms(self):
        """Print initial beam parameters in lattice file"""
        print '\nIon Charge States = ', self.M.conf()['IonChargeStates']
        print 'IonEs [MeV]       = ', self.M.conf()['IonEs']/1e6
        print 'IonEk [MeV]       = ', self.M.conf()['IonEk']/1e6
        print '\nBaryCenter 0:\n', self.M.conf()['BaryCenter0']
        print '\nBaryCenter 1:\n', self.M.conf()['BaryCenter1']
        print '\nBeam Envelope 0:\n', self.M.conf()['S0']
        print '\nBeam Envelope 1:\n', self.M.conf()['S1']

    def prt_lat(self):
        """Print all lattice elements"""
        for (i,elem) in enumerate(self.lat):
            print(i, elem['name'], elem['type'])


    def tcs(self,lpf=0, opf=1):
        """
        Main function of one through calculation
        Calculation data is stored at LD.
        LD[i] contains:
        [pos, xcen, ycen, zrad, xrms, yrms, ref_phis, ref_IonEk]
        
        """
        S = self.M.allocState({})
        self.M.propagate(S, 0, 1)

        # set initial beam data
        S.ref_IonZ   = self.refIonZ
        S.IonZ       = self.IonZ_io

        S.moment0    = self.BC0
        S.moment1    = self.ENV0

        S.ref_IonEk  = self.refIonEk
        
        S.phis       = S.moment0[PS_S,:]
        S.IonEk      = S.moment0[PS_PS,:]*MeVtoeV + S.ref_IonEk

        S.clng = self.clng

        fin = len(self.M)

        
        # store initial beam data
        self.LD[0][0] = S.pos

        #Mean data
        self.LD[0][1] = S.moment0_env[0]
        self.LD[0][2] = S.moment0_env[2]
        self.LD[0][3] = S.moment0_env[4]
        self.LD[0][4] = S.moment0_rms[0]
        self.LD[0][5] = S.moment0_rms[2]
        self.LD[0][6] = S.ref_phis
        self.LD[0][7] = S.ref_IonEk
        

        # propagate step by step and store beam data
        for i in range(1,len(self.M)):
            self.M.propagate(S, i, 1)
            
            
            self.LD[i][0] = S.pos
            #Mean data
            self.LD[i][1] = S.moment0_env[0]
            self.LD[i][2] = S.moment0_env[2]
            self.LD[i][3] = S.moment0_env[4]
            self.LD[i][4] = S.moment0_rms[0]
            self.LD[i][5] = S.moment0_rms[2]
            self.LD[i][6] = S.ref_phis
            self.LD[i][7] = S.ref_IonEk
            

        #output data for plotting
        if opf: numpy.savetxt('ldata.txt',self.LD)

        if not lpf: return S

        
    def tcs_seg(self,start,end,SD=None,allcs=0,opf=1):
        """
        Main function of segmented section calculation
        (S, RD) = tcs_seg(start_index, end_index, SD=beam_data)
        beam_data can be generated by SaveBeam(S)
        RD[i] contains:
        [pos, xcen, ycen, zrad, xrms, yrms, ref_phis, ref_IonEk]
        i: index of the element

        For all charge outputs: allcs=1,
        (S, RD, AD) = tcs_seg(start_index, end_index, SD=beam_data, allcs=1)
        AD[k][i] contains:
        [pos, xcen, ycen, zrad, xrms, yrms]
        k: index of the charge state
        """
        S = self.M.allocState({})
        self.M.propagate(S, 0, 1)
        
        if SD == None :
            # set initial beam data
            S.ref_IonZ   = self.refIonZ
            S.IonZ       = self.IonZ_io

            S.moment0    = self.BC0
            S.moment1    = self.ENV0

            S.ref_IonEk  = self.refIonEk
            
            S.phis       = S.moment0[PS_S,:]
            S.IonEk      = S.moment0[PS_PS,:]*MeVtoeV + S.ref_IonEk
        
        else :
            # set initial beam data
            S.ref_IonZ    = SD.refIonZ
            S.IonZ        = SD.IonZ_io

            S.moment0     = SD.BC0
            S.moment1     = SD.ENV0

            S.ref_phis    = SD.refphis
            S.phis        = SD.phis_io

            S.ref_IonEk   = SD.refIonEk
            S.IonEk       = SD.BC0[PS_PS,:]*MeVtoeV + SD.refIonEk
            S.pos         = SD.pos

        phis_ini      = S.ref_phis

        S.clng = self.clng

        fin = end - start + 1
        RD = numpy.zeros((fin,8))

        #if allcs: AD = [[[0.0]*6 for i in range(fin)] for j in range(len(S.IonZ))]
        if allcs: AD = numpy.zeros((len(S.IonZ),fin,8))

        # store initial beam data
        RD[0][0] = S.pos
        RD[0][1] = S.moment0_env[0]
        RD[0][2] = S.moment0_env[2]
        RD[0][3] = S.moment0_env[4]
        RD[0][4] = S.moment0_rms[0]
        RD[0][5] = S.moment0_rms[2]
        RD[0][6] = S.ref_phis - phis_ini
        RD[0][7] = S.ref_IonEk

        if allcs:
            for k in range(len(S.IonZ)):
                AD[k][0][0] = S.pos
                AD[k][0][1] = S.moment0[0,k]
                AD[k][0][2] = S.moment0[2,k]
                AD[k][0][3] = S.moment0[4,k]
                AD[k][0][4] = numpy.sqrt(S.moment1[0,0,k])
                AD[k][0][5] = numpy.sqrt(S.moment1[2,2,k])
           
        # propagate step by step and store beam data
        for (j,i) in enumerate(range(start,end)):
            self.M.propagate(S, i+1, 1)

            RD[j+1][0] = S.pos
            RD[j+1][1] = S.moment0_env[0]
            RD[j+1][2] = S.moment0_env[2]
            RD[j+1][3] = S.moment0_env[4]
            RD[j+1][4] = S.moment0_rms[0]
            RD[j+1][5] = S.moment0_rms[2]
            RD[j+1][6] = S.ref_phis - phis_ini
            RD[j+1][7] = S.ref_IonEk

            if allcs:
                for k in range(len(S.IonZ)):
                    AD[k][j+1][0] = S.pos
                    AD[k][j+1][1] = S.moment0[0,k]
                    AD[k][j+1][2] = S.moment0[2,k]
                    AD[k][j+1][3] = S.moment0[4,k]
                    AD[k][j+1][4] = numpy.sqrt(S.moment1[0,0,k])
                    AD[k][j+1][5] = numpy.sqrt(S.moment1[2,2,k])

        if opf: numpy.savetxt('ldata.txt',RD)

        if allcs: return (S,RD,AD)
        else : return (S,RD)        
        
        
    def looptcs(self):
        """
        Main loop for Virtual Accelerator
        This function should be called as background job.     
        Loop is stopped  by setting itr = 1.
        """ 
        while self.itr < 1: 
            #self.genRandomNoise() #developing
            self.tcs(lpf=1)
            #self.itr +=1 

    class SaveBeam:
        def __init__(self,S_in):
            self.refIonZ  = S_in.ref_IonZ
            self.IonZ_io  = S_in.IonZ

            self.refphis  = S_in.ref_phis
            self.phis_io  = S_in.phis
            
            self.BC0      = S_in.moment0
            self.ENV0     = S_in.moment1
            
            self.refIonEk = S_in.ref_IonEk
            self.pos      = S_in.pos



class SaveBeam:
    def __init__(self,S_in):
        self.refIonZ  = S_in.ref_IonZ
        self.IonZ_io  = S_in.IonZ

        self.refphis  = S_in.ref_phis
        self.phis_io  = S_in.phis
        
        self.BC0      = S_in.moment0
        self.ENV0     = S_in.moment1
        
        self.refIonEk = S_in.ref_IonEk
        self.pos      = S_in.pos
