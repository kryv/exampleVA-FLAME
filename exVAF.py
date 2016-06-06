import sys
import math
import numpy 
import time
import cPickle
import copy

home_dir = './build'
python_dir  = home_dir + '/python/'
file_name = './to_strl.lat'

sys.path.append(python_dir)

from uscsi import Machine

PS_X = 0; PS_PX = 1; PS_Y = 2; PS_PY = 3; PS_S = 4; PS_PS = 5

MeVtoeV = 1e6


class VAF:
    def __init__(self, fname=file_name):
        self.name = fname

        self.itr = 0

        with open(self.name, 'rb') as inf:

            self.M = Machine(inf.read())

            self.lat = self.M.conf()['elements']

            self.refIonZ = self.M.conf()['IonChargeStates'][0]
            self.realIonZ = self.M.conf()['IonChargeStates'][0]
            self.BC0 = self.M.conf()['BaryCenter0']
            self.ENV0 = numpy.split(self.M.conf()['S0'],7)

            #self.PRNseed = self.M.conf()['PRN_seed']
            #self.PRNvalue = self.M.conf()['PRN_value']

            self.ldata = [[0.0]*6 for i in range(len(self.M))]

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

            with open('plot.info','w') as f:
                iteminfo=[len(self.M),self.bpmls,self.corls,self.cavls,self.solls]
                cPickle.dump(iteminfo,f)




    def getelem(self,num):
        print self.M.conf()['elements'][num]


    def setelem(self,num,name,val):
        #D = self.M.conf()['elements'][num]
        D = self.lat[num]
        D[name] = val
        self.M.reconfigure(num, D)


    """
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
        print '\nIon Charge States = ', self.M.conf()['IonChargeStates']
        print 'IonEs [MeV]       = ', self.M.conf()['IonEs']/1e6
        print 'IonEk [MeV]       = ', self.M.conf()['IonEk']/1e6
        print '\nBaryCenter 1:\n', self.M.conf()['BaryCenter0']
        print '\nBaryCenter 2:\n', self.M.conf()['BaryCenter1']
        print '\nBeam Envelope 1:\n', self.M.conf()['S0']
        print '\nBeam Envelope 2:\n', self.M.conf()['S1']

    def prt_lat(self):
        for (i,elem) in enumerate(self.M.conf()['elements']):
            print(i, elem['name'], elem['type'])



    def tcs(self,lpf=0):
        #global M, refIonZ, realIonZ, BC0, ENV0, ldata
        
        S = self.M.allocState({})
        self.M.propagate(S, 0, 1)

        #S.error_value      = self.errorv

        S.ref_IonZ    = self.refIonZ
        S.real_IonZ   = self.realIonZ

        S.moment0[:]  = self.BC0
        S.state[:]    = self.ENV0

        #S.PRN_seed     = self.PRNseed
        #S.PRN_value    = self.PRNvalue

        S.real_gamma  = S.real_IonW/S.real_IonEs;
        S.real_beta   = math.sqrt(1e0-1e0/S.real_gamma**2.0);
        S.real_bg     = S.real_beta*S.real_gamma;

        S.real_phis   = S.moment0[PS_S];
        S.real_IonEk  += S.moment0[PS_PS]*MeVtoeV;

        fin = len(self.M)

        self.ldata[0][0] = S.pos
        self.ldata[0][1] = S.moment0[0]
        self.ldata[0][2] = S.moment0[2]
        self.ldata[0][3] = S.moment0[4]
        self.ldata[0][4] = numpy.sqrt(S.state[0][0])
        self.ldata[0][5] = numpy.sqrt(S.state[2][2])

        for i in range(1,len(self.M)):
            #print i
            self.M.propagate(S, i, 1)
            
            self.ldata[i][0] = S.pos
            self.ldata[i][1] = S.moment0[0]
            self.ldata[i][2] = S.moment0[2]
            self.ldata[i][3] = S.moment0[4]
            self.ldata[i][4] = numpy.sqrt(S.state[0][0])
            self.ldata[i][5] = numpy.sqrt(S.state[2][2])

        numpy.savetxt('ldata.txt',self.ldata)

        if not lpf: return S

    def looptcs(self):
        while self.itr < 1:
            #self.genRandomNoise()
            self.tcs(lpf=1)
            #self.itr +=1


