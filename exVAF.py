import sys
import math
import numpy 
import time
import cPickle
import re
import fnmatch

home_dir = './build'
python_dir  = home_dir + '/python'
file_name = python_dir + '/uscsi/test/to_strl.lat'

sys.path.append(python_dir)

from uscsi import Machine

PS_X = 0; PS_PX = 1; PS_Y = 2; PS_PY = 3; PS_S = 4; PS_PS = 5
MeVtoeV = 1e6

class VAF:
    """
    Contains:
    getelem(num)
    setelem(num,name,value)
    prt_lat_prms()
    prt_lat()
    tcs()
    looptcs()
    """
    def __init__(self, fname=file_name):
        self.name = fname
        self.itr = 0

        with open(self.name, 'rb') as inf:

            self.M = Machine(inf.read())

            self.lat = self.M.conf()['elements']

            self.refIonZ = self.M.conf()['IonChargeStates'][0]
            self.realIonZ = self.M.conf()['IonChargeStates'][0]

            self.refIonEk = self.M.conf()['IonEk']

            self.BC0 = self.M.conf()['BaryCenter0']
            self.ENV0 = numpy.split(self.M.conf()['S0'],7)

            # memory space for all beam data
            self.ldata = [[0.0]*8 for i in range(len(self.M))]

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
        getelem(position number of lattice element)
        """
        print self.M.conf()['elements'][num]


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
        D = self.lat[num]
        D[name] = float(val)
        self.M.reconfigure(num, D)


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


    def tcs(self,lpf=0):
        """
        Main function of one through calculation
        Calculation data is stored at ldata.
        ldata[i] contains:
        [pos, xcen, ycen, zrad, xrms, yrms, ref_phis, ref_IonEk]
        
        """
        S = self.M.allocState({})
        self.M.propagate(S, 0, 1)

        # set initial beam data
        S.ref_IonZ    = self.refIonZ
        S.real_IonZ   = self.realIonZ

        S.moment0[:]  = self.BC0
        S.state[:]    = self.ENV0

        S.ref_IonEk   = self.refIonEk
        S.real_IonEk  = self.refIonEk

        #S.ref_IonW    = self.ref_IonEk + self.ref_IonEk
        #S.real_IonW   = self.ref_IonW

        #S.real_gamma  = S.real_IonW/S.real_IonEs;
        #S.real_beta   = math.sqrt(1e0-1e0/S.real_gamma**2.0);
        #S.real_bg     = S.real_beta*S.real_gamma;

        S.real_phis   = S.moment0[PS_S];
        S.real_IonEk  += S.moment0[PS_PS]*MeVtoeV;

        fin = len(self.M)

        # store initial beam data
        self.ldata[0][0] = S.pos
        self.ldata[0][1] = S.moment0[0]
        self.ldata[0][2] = S.moment0[2]
        self.ldata[0][3] = S.moment0[4]
        self.ldata[0][4] = numpy.sqrt(S.state[0][0])
        self.ldata[0][5] = numpy.sqrt(S.state[2][2])
        self.ldata[0][6] = S.ref_phis
        self.ldata[0][7] = S.ref_IonEk

        # propagate step by step and store beam data
        for i in range(1,len(self.M)):
            self.M.propagate(S, i, 1)
            
            self.ldata[i][0] = S.pos
            self.ldata[i][1] = S.moment0[0]
            self.ldata[i][2] = S.moment0[2]
            self.ldata[i][3] = S.moment0[4]
            self.ldata[i][4] = numpy.sqrt(S.state[0][0])
            self.ldata[i][5] = numpy.sqrt(S.state[2][2])
            self.ldata[i][6] = S.ref_phis
            self.ldata[i][7] = S.ref_IonEk
            
        #output data for plotting
        numpy.savetxt('ldata.txt',self.ldata)

        if not lpf: return S

        
    def tcs_seg(self,start,end,SD=None):
        """
        Main function of segmented section calculation
        (S, rdata) = tcs_seg(start_index, end_index, SD=beam_data)
        beam_data can be generated by SaveBeam(S)
        rdata[i] contains:
        [pos, xcen, ycen, zrad, xrms, yrms, ref_phis, ref_IonEk]
        
        """
        S = self.M.allocState({})
        self.M.propagate(S, 0, 1)
        
        if SD == None :
            # set initial beam data
            S.ref_IonZ    = self.refIonZ
            S.real_IonZ   = self.realIonZ

            S.moment0[:]  = self.BC0
            S.state[:]    = self.ENV0

            S.ref_IonEk   = self.refIonEk
            S.real_IonEk  = self.refIonEk
        
        else :
            # set initial beam data
            S.ref_IonZ    = SD.refIonZ
            S.real_IonZ   = SD.realIonZ

            S.moment0[:]  = SD.BC0
            S.state[:]    = SD.ENV0

            S.ref_IonEk   = SD.refIonEk
            S.real_IonEk  = SD.refIonEk
            S.pos         = SD.pos

        S.real_phis   = S.moment0[PS_S];
        S.real_IonEk  += S.moment0[PS_PS]*MeVtoeV;
        
        fin = end - start + 1
        rdata = [[0.0]*8 for i in range(fin)]

        # store initial beam data
        rdata[0][0] = S.pos
        rdata[0][1] = S.moment0[0]
        rdata[0][2] = S.moment0[2]
        rdata[0][3] = S.moment0[4]
        rdata[0][4] = numpy.sqrt(S.state[0][0])
        rdata[0][5] = numpy.sqrt(S.state[2][2])
        rdata[0][6] = S.ref_phis
        rdata[0][7] = S.ref_IonEk

        
        # propagate step by step and store beam data
        for (j,i) in enumerate(range(start,end)):
            self.M.propagate(S, i, 1)

            rdata[j+1][0] = S.pos
            rdata[j+1][1] = S.moment0[0]
            rdata[j+1][2] = S.moment0[2]
            rdata[j+1][3] = S.moment0[4]
            rdata[j+1][4] = numpy.sqrt(S.state[0][0])
            rdata[j+1][5] = numpy.sqrt(S.state[2][2])
            rdata[j+1][6] = S.ref_phis
            rdata[j+1][7] = S.ref_IonEk
        
        return (S,rdata)        
        
        
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
        self.realIonZ = S_in.real_IonZ
        
        self.BC0      = S_in.moment0
        self.ENV0     = S_in.state
        
        self.refIonEk = S_in.ref_IonEk
        self.pos      = S_in.pos
