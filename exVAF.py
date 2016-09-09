import sys
import math
import numpy as np
import time
import cPickle
import re
import fnmatch
import collections

#home_dir = '/home/fukushim/software/flame_main_local/build'
home_dir = '/home/fukushim/software/flame_local/build'
python_dir  = home_dir + '/python'
#file_name = python_dir + '/flame/test/to_strl.lat'
file_name = '/home/fukushim/software/flame_main_local/lat/to_strl.lat'

sys.path.append(python_dir)

from flame import Machine

PS_X = 0; PS_PX = 1; PS_Y = 2; PS_PY = 3; PS_S = 4; PS_PS = 5
MeVtoeV = 1e6

cen = ['x','xp','y','yp', 'z', 'zp']
rms = ['xrms', 'xprms', 'yrms', 'yprms','zrms', 'zprms']

class VAF:
    """
VAF
===

Tool box for FLAME python interface

You can start VAF by
>>> va = vaf.VAF(fname=${lattice_file_path})

Initial beam parameters 
(default is the same as lattice file)
-------------------------------------
BC0: 7*n array (n: number of charge state)
    Barycenter vectors for each charge state.
    contains:[x,  px,  y,  py,  z,   pz,    1]
    unit:    [mm, rad, mm, rad, rad, MeV/u, 1]

ENV0: 7*7*n array (n: number of charge state)
    Envelope matrixes for each charge state.

IonZ: n array (n: number of charge state)
    Charge to mass ratio for each charge state.

RefIonEk: float
    Reference kinetic energy at starting point [eV/u]



This tool box contains
(Check each documentation for detail.)
--------------------------------------
    prt_beam_prms()
    prt_lat_list()
    getelem(int)
    getvalue(int,str)
    getindex(str)
    getindexu(str)
    setvalue(int,str,float)
    tcs()
    tcs_seg()
    SaveBeam(S)

On developping
--------------
    genRN()
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
            self.IonZ = S.IonZ # list of real IonZ

            self.refIonEk = S.ref_IonEk

            self.BC0 = S.moment0  
            self.ENV0 = S.moment1

            # memory space for all beam data
            # self.LD = [[0.0]*8 for i in range(len(self.M))]
            self.LD = np.zeros((len(self.M),9))

            self.LD2 = np.zeros((len(self.M),7))

            # store element position data
            self.bpmls=[]
            self.corls=[]
            self.cavls=[]
            self.solls=[]
            self.cavpara=[]
            self.solpara=[]

            for i in range(len(self.M)):
                #elem = self.M.conf()['elements'][i]
                elem = self.lat[i]
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
getelem(a)

    Get base parameter of lattice element.
    
    Parameters
    ----------
    a : int
        Input index of the element.
    
    Return
    ------
    out: dict
        Python dict style lattice element.

    Examples
    --------
    >>> print va.getelem(8)
    OrderedDict([('B', 5.34), ('L', 0.1), ('aper', 0.02), ('name', 'ls1_ca01_sol1_d1131_1'), ('type', 'solenoid')])
        """
        #return self.M.conf()['elements'][num]
        return self.lat[num]


    def getvalue(self,num,name):
        """
getvalue(a, b) 
    
    Get latest parameter of lattice element.
    
    Parameters
    ----------
    a: int
        Input index of the element
    b: str
        Input name of the parameter

    Return
    ------
    out: depends on parameter name
        Value of the parameter.

    Example
    -------
    >>> print va.getvalue(8,'B')
    5.34
        """
        return self.M.conf(num)[name]



    def getindex(self,name,searchfrom='name'):
        """
getindex(a, searchby='name')
    
    Get index list of lattice elements by python style regular expression.

    Parameters
    ----------
    a: str
        Keyword for search.
    searchfrom: str
        Search from 'name' or 'type'.
    
    Return
    ------
    out: list
        Search result.
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
getindex(a, searchby='name')
    
    Get index list of lattice elements by unix style regular expression.

    Parameters
    ----------
    a: str
        Keyword for search.
    searchfrom: str
        Search from 'name' or 'type'.
        
    Return
    ------
    out: list
        Search result.
        """
        name = name.replace(':','_').lower()
        result = []

        for (i,elem) in enumerate(self.lat):
            if fnmatch.fnmatch(elem[searchby],name):
                result.append(i)
        return result


    def setvalue(self,num,name,val):
        """
setvalue(a, b, c)

    Set parameter of lattice element.
    
    Parameters
    ----------
    a: int
        Input index of the element.
    b: str
        Input name of the parameter.
    c: float 
        Set value to the parameter.
        """
        self.M.reconfigure(num,{name:float(val)})


    
    # Developing
    def genRN(self,scl=0.0005,seed=None):
        np.random.seed(seed)
        for (i,elem) in enumerate(self.lat):
            rnds = np.random.randn(5)*scl #rms 0.5 mm and rms 0.5 mrad
            etype = elem['type']
            if etype == 'rfcavity' or etype == 'solenoid' or etype == 'quadrupole':
               self.M.reconfigure(i,{'dx':float(rnds[0])})
               self.M.reconfigure(i,{'dy':float(rnds[1])})
               #self.M.reconfigure(i,{'pitch':float(rnds[2])})
               #self.M.reconfigure(i,{'yaw':float(rnds[3])})
               #self.M.reconfigure(i,{'tilt':float(rnds[4])})

    
    def prt_beam_prms(self):
        """
Print initial beam parameters in lattice file
        """
        print '\nIon Charge States = ', self.M.conf()['IonChargeStates']
        print 'IonEs [MeV]       = ', self.M.conf()['IonEs']/1e6
        print 'IonEk [MeV]       = ', self.M.conf()['IonEk']/1e6
        print '\nBaryCenter 0:\n', self.M.conf()['BaryCenter0']
        print '\nBaryCenter 1:\n', self.M.conf()['BaryCenter1']
        print '\nBeam Envelope 0:\n', self.M.conf()['S0']
        print '\nBeam Envelope 1:\n', self.M.conf()['S1']

    def prt_lat_list(self):
        """
Print all lattice elements with index.
    returns (index, name, type)
        """
        for (i,elem) in enumerate(self.lat):
            print(i, elem['name'], elem['type'])


    def tcs(self,lpf=0, opf=1):
        """
tcs(lpf=0, opf=1)
    
    Main function for one through calculation.
    Calculation data is stored as VAF.LD.

        VAF.LD[i]: 
        ----------
        i: int
            index of the element
        
        contains:
            [pos, xcen, ycen, zcen, xrms, yrms, zrms, ref_phis, ref_IonEk]
        units:
            [m,   mm,   mm,   mm,   mm,   mm,   mm,   rad,      eV/u]
            
    Parameters
    ----------
    lpf: bool
        flag for loop run.
    opf: bool
        flag for file type output.

    Returns
    -------
    S: object
        Beam parameters at the ending point.
        """
        S = self.M.allocState({})
        self.M.propagate(S, 0, 1)

        # set initial beam data
        S.ref_IonZ   = self.refIonZ
        S.IonZ       = self.IonZ

        S.moment0    = self.BC0
        S.moment1    = self.ENV0

        S.ref_IonEk  = self.refIonEk
        
        S.phis       = S.moment0[PS_S,:]
        S.IonEk      = S.moment0[PS_PS,:]*MeVtoeV + S.ref_IonEk

        #S.clng = self.clng

        fin = len(self.M)

        
        # store initial beam data
        self.LD[0][0] = S.pos

        #Mean data
        self.LD[0][1] = S.moment0_env[0]
        self.LD[0][2] = S.moment0_env[2]
        self.LD[0][3] = S.moment0_env[4]
        self.LD[0][4] = S.moment0_rms[0]
        self.LD[0][5] = S.moment0_rms[2]
        self.LD[0][6] = S.moment0_rms[4]
        self.LD[0][7] = S.ref_phis
        self.LD[0][8] = S.ref_IonEk

        # store initial beam data
        self.LD2[0][0] = S.pos
        #Mean data
        self.LD2[0][1] = S.moment0_env[1]
        self.LD2[0][2] = S.moment0_env[3]
        self.LD2[0][3] = S.moment0_env[5]
        self.LD2[0][4] = S.moment0_rms[1]
        self.LD2[0][5] = S.moment0_rms[3]
        self.LD2[0][6] = S.moment0_rms[5]


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
            self.LD[i][6] = S.moment0_rms[4]
            self.LD[i][7] = S.ref_phis
            self.LD[i][8] = S.ref_IonEk

            self.LD2[i][0] = S.pos
            #Mean data
            self.LD2[i][1] = S.moment0_env[1]
            self.LD2[i][2] = S.moment0_env[3]
            self.LD2[i][3] = S.moment0_env[5]
            self.LD2[i][4] = S.moment0_rms[1]
            self.LD2[i][5] = S.moment0_rms[3]
            self.LD2[i][6] = S.moment0_rms[5]

        #output data for plotting
        if opf: np.savetxt('ldata.txt',self.LD)

        if not lpf: return S


    def tcs2(self):
        """
tcs(lpf=0, opf=1)
    
    Main function for one through calculation.
    Calculation data is stored as VAF.LD.

        VAF.LD[i]: 
        ----------
        i: int
            index of the element
        
        contains:
            [pos, xcen, ycen, zcen, xrms, yrms, zrms, ref_phis, ref_IonEk]
        units:
            [m,   mm,   mm,   mm,   mm,   mm,   mm,   rad,      eV/u]
            
    Parameters
    ----------
    lpf: bool
        flag for loop run.
    opf: bool
        flag for file type output.

    Returns
    -------
    S: object
        Beam parameters at the ending point.
        """
        S = self.M.allocState({})
        self.M.propagate(S, 0, 1)

        # set initial beam data
        S.ref_IonZ   = self.refIonZ
        S.IonZ       = self.IonZ

        S.moment0    = self.BC0
        S.moment1    = self.ENV0

        S.ref_IonEk  = self.refIonEk
        
        S.phis       = S.moment0[PS_S,:]
        S.IonEk      = S.moment0[PS_PS,:]*MeVtoeV + S.ref_IonEk

        ### Reconstract data
        U = collections.OrderedDict()
        T = collections.OrderedDict()
        
        U,T = self.save(U, T, S)


        #S.clng = self.clng

        fin = len(self.M)

        H = self.M.propagate(S, 1, fin, observe=range(fin))
        
        ### Reconstract data
        U = collections.OrderedDict()
        T = collections.OrderedDict()

        for i in range(fin-1):
            U,T = self.save(U, T, H[i][1])

        return U,T

        
    def tcss(self,start=None,end=None,S_in=None,allcs=0,opf=1):

        if start == None : start = 1
        if end == None : end = len(self.M) -1
        
        if S_in == None :

            S = self.M.allocState({})
            self.M.propagate(S, 0, 1)
            # set initial beam data
            S.ref_IonZ   = self.refIonZ
            S.IonZ       = self.IonZ

            S.moment0    = self.BC0
            S.moment1    = self.ENV0

            S.ref_IonEk  = self.refIonEk
            
            S.phis       = S.moment0[PS_S,:]
            S.IonEk      = S.moment0[PS_PS,:]*MeVtoeV + S.ref_IonEk
        
        else :
            # set initial beam data
            S = S_in.clone()      

        #S.clng = self.clng

        H_ini = (0, S.clone())

        fin = end - start + 1

        H = self.M.propagate(S, start, fin, observe=range(len(self.M)))

        return [H_ini] + H
        #return H
        
        
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


    def save(self,u, t, S, label = 'total'):
        for (k,j) in enumerate(S.IonZ):
            cs = str(int(round(j*238.0)))
            try:
                u[cs].pos.append(S.pos)
            except:
                u[cs] = self.SPI()
                u[cs].pos.append(S.pos)

            for i,(n,m) in enumerate(zip(cen,rms)):
                getattr(u[cs], n).append(S.moment0[i,k])
                getattr(u[cs], m).append(np.sqrt(S.moment1[i,i,k]))
                
        try:
            t[label].pos.append(S.pos)
        except:
            t[label] = self.SPI()
            t[label].pos.append(S.pos)

        for i,(n,m) in enumerate(zip(cen,rms)):
            getattr(t[label], n).append(S.moment0_env[i])
            getattr(t[label], m).append(S.moment0_rms[i])

        t[label].ek.append(S.ref_IonEk)
        t[label].phis.append(S.ref_phis)   

        return u,t 



    class SaveBeam:
        """
SaveBeam(S_in)
    Store beam parameters for VAF.tcs_seg().
    
    Parameters
    ----------
    S_in: object
        Beam parameters made by VAF.tcs() or VAF.tcs_seg().
    
    Returns
    -------
    out: object
        Beam parameters object for VAF.tcs_seg().
        """
        def __init__(self,S_in):
            self.refIonZ  = S_in.ref_IonZ
            self.IonZ_io  = S_in.IonZ

            self.refphis  = S_in.ref_phis
            self.phis_io  = S_in.phis
            
            self.BC0      = S_in.moment0
            self.ENV0     = S_in.moment1
            
            self.refIonEk = S_in.ref_IonEk
            self.pos      = S_in.pos

    class SPI:
        def __init__(self):
            self.pos = []
            self.x = []
            self.xp = []
            self.y = []
            self.yp = []
            self.z = []
            self.zp = []
            self.xrms = []
            self.xprms = []
            self.yrms = []
            self.yprms = []            
            self.zrms = []
            self.zprms = []
            self.ek = []
            self.phis = []

        def n(self, name):
            self.x_n = np.array(self.x)


    #---------------------------------------------------------

    def rms2(self,xy,S,cs=0):
        if xy == 'x': i=0
        elif xy == 'px': i=1
        elif xy == 'y': i=2
        elif xy == 'py': i=3
        elif xy == 'z': i=4
        elif xy == 'pz': i=5
        return S.moment1[i,i,cs]

    def rms(self,xy,S,cs=0):
        return np.sqrt(self.rms2(xy,S,cs))

    def emit(self,xy,S,cs=0):
        if xy == 'x': i=0
        elif xy == 'y': i=2
        else : raise NameError('First argument must be \'x\' or \'y\'')
        xx = S.moment1[i,i,cs]*1e-6
        xxp = S.moment1[i,i+1,cs]*1e-3
        xpxp = S.moment1[i+1,i+1,cs]
        return np.sqrt(np.abs(xx*xpxp - xxp*xxp))

    def nemit(self,xy,S,cs=0):
        bg = S.beta[cs]*S.gamma[cs]
        return bg*self.emit(xy,S,cs)

    def bfunc(self,xy,S,cs=0):
        return self.rms2(xy,S,cs)*1e-6/self.emit(xy,S,cs)

    def bfuncp(self,xy,S,cs=0):
        if xy == 'x': i = 0
        elif xy == 'y': i = 2
        return S.moment1[i,i+1,cs]*1e-3/self.emit(xy,S,cs)


    def hbfunc(self,xy,f,cs=0):
        ans = np.zeros(len(f))
        for i in range(len(f)):
            ans[i] = self.bfunc(xy,f[i][1],cs)
        return ans

    def hbfuncp(self,xy,f,cs=0):
        ans = np.zeros(len(f))
        for i in range(len(f)):
            ans[i] = self.bfuncp(xy,f[i][1],cs)
        return ans


    def hpos(self,f):
        ans = np.zeros(len(f))
        for i in range(len(f)):
            ans[i] = f[i][1].pos
        return ans

    def hcen(self,xy,f,cs=0):
        if xy == 'x': i=0
        elif xy == 'px': i=1
        elif xy == 'y': i=2
        elif xy == 'py': i=3
        elif xy == 'z': i=4
        elif xy == 'pz': i=5

        ans = np.zeros(len(f))
        if cs == -1:
            for k in range(len(f)):
                ans[k] = f[k][1].moment0_env[i]
        else:
            for k in range(len(f)):
                ans[k] = f[k][1].moment0[i,cs]

        return ans


    def hrms(self,xy,f,cs=0):
        if xy == 'x': i=0
        elif xy == 'px': i=1
        elif xy == 'y': i=2
        elif xy == 'py': i=3
        elif xy == 'z': i=4
        elif xy == 'pz': i=5

        ans = np.zeros(len(f))
        if cs == -1:
            for k in range(len(f)):
                ans[k] = f[k][1].moment0_rms[i]
        else:
            for k in range(len(f)):
                ans[k] = self.rms(xy,f[k][1],cs)

        return ans

    def dpp(self,S,cs=0):
        gm = S.gamma[cs]
        iw = S.IonEk[cs]
        dw = S.moment0[5][cs]*1e6
        return gm/(1+gm)*dw/iw

    #---------------------------------------------------------


    def hdpp(self,f,cs=0):
        ans = np.zeros(len(f))
        for k in range(len(f)):
            ans[k] = self.dpp(f[k][1],cs=cs)
        return ans

    def hdsp(self,xy,f,cs=0):
        if xy == 'x' : i = 0
        elif xy == 'y' : i = 2
        ans = np.zeros(len(f))
        for k in range(len(f)):
            gm = f[k][1].gamma[cs]
            iw = f[k][1].IonEk[cs]
            xw = f[k][1].moment1[i,5,cs]*1e-3
            ww = f[k][1].moment1[5,5,cs]*1e6*gm/(1+gm)/iw
            ans[k] = xw/ww          
        return ans



    def getek(self,f,ls):
        ans = np.zeros(len(ls))
        for (i,j) in enumerate(ls):
            ans[i] = f[j][1].ref_IonEk
        return ans

    def getid0(self,id,f,ls):
        ans = np.zeros(len(ls))
        for (i,j) in enumerate(ls):
            ans[i] = f[j][1].moment0_env[id]
        return ans

    def getid00(self,cs,id,f,ls):
        ans = np.zeros(len(ls))
        for (i,j) in enumerate(ls):
            ans[i] = f[j][1].moment0[id,cs]
        return ans

    def getid1(self,id,f,ls):
        ans = np.zeros(len(ls))
        for (i,j) in enumerate(ls):
            ans[i] = f[j][1].moment0_rms[id]
        return ans


    def getid01(self,cs,id,f,ls):
        ans = np.zeros(len(ls))
        for (i,j) in enumerate(ls):
            ans[i] = np.sqrt(f[j][1].moment1[id,id,cs])
        return ans

