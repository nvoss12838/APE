import numpy as np
import scipy as sp
import matplotlib.pylab as plt
'''
objects for storing pest related control file
'''

class controldata(object):
  '''
  object containing paramaters in the control data section of the pestcontrol file
  
  Attributes:
    RSTFLE : "restart or norestart" Instructs PEST whether to write restart data
    PESTMODE : "esimation","prediction","regularisation","pareto" PEST mode of operation
    NPAR : Number of Paramaters
    NOBS : Number of Observations
    NPARGP : Number of Paramapter Groups
    NPRIOR : Number of Prior Information
    NOBSG : Number of Observation Groups
    MAXCOMPDIM : Number of Elements in the compresed Jacobian matrix
    NTPLFLE : Number of Template Files
    NINSFLE : Number of Instruction Files
    PRECIS : "single or double" Precision
    DPOINT : "point or nopoint" Whether to use Decimal Points
    NUMCOM : Number of command lines used to run model
    JACFILE : 0,1,-1 whether the model uses an external derivitives file
    MESSFILE : 0,1 whether PEST writes a PEST-to-Model message file
    RLAMBDA1 : initial Marquart lambda
    RLAMFAC : dictates Marqaurt lambda adjustment process
    PHIRATSUF : fractional objective function sufficient for end of current iterations
    PHIREDLAM : termination criterion for Marquardt lamda search
    NUMLAM : maximum number of marquart lambda test
    JACUPDATE : activation of Broyden;s Jacobian Update procedure
    LAMFORGIVE : treat model run failure during the lamda search as high objective function
    RELPARMAX : paramater relative change limit
    FACPARMAX : paramater factor change limit
    FACORIG : minimum fraction of original paramearcheter value in evaulating relative change
    IBOUNDSTICK : instruct PEST not to compute derivatives for parameter at its bounds
    UPVECBEND : instructs PEST to bend parameter upgrade vector if parameter hits bounds
    ABSPARMAX : paramater absolute change limit
    PHIREDSWH : sets objective function change for introduction of central derivatives
    NOPTMAX : -2,-1,0,or any number greater then zero number of optimisation iterations
    PHIREDSTP : relative objective function reduction triggering termination
    NPHISTP : number of successive iterations over which PHIREDSTP applies
    NPHINORED : number of iterations since last drop in objective function to trigger termination
    RELPARSTP : maximum relative paramater change triggering termination
    NRELPAR : number of successive iteration over which RELPARSTP applies
    PHISTOPTHRESH : objective function threshold triggering termination
    LASTRUN : 0, 1 instructs PEST to undertake final model run with best paramaters 
    PHIABANDON : objective function value at which to abandon optimisation 
    ICOV : 0,1 record covariance matrix 
    ICOR : 0,1 record correlation matrix
    IEIG : 0,1 record eigenvectiors in matrix file 
    IRES : 0,1 record resolution data
    JCOSAVE : "jcosave or nojcosave" save best Jacobian fileas JCO file 
    VERBOSEREC : "verboserec or noverboserec" omit obeservation data from rec file
    JCOSAVEITN : "jcosaveitn" or nojcosaveitn" write current jacobian matrix to iteration specific JCO file
    REISAVEITN : "reisaveitn or noreisaveitn" Store residuals to iteration specific residuals file
    PARSAVEITN : "parsaveitn or noparsaveitin" Store iteration specific paramater value files
    '''
    def __init__(self,rstfile,pestmode,npar,nobs,npargp,nprior,nobsgp,maxcompdim,ntplfle,ninsfle,precis,dpoint,numcom,jacfile,messfile \
      ,obsreref,rlambda1,rlamfac,phiratsuf,phiredlam,numlam,jacupdate,lamforgive,derforgive,relparmax,facparmax,facorig,iboundstick,upvecbound,absparmax \
	,phiredswh,noptmax,phiredstp,nphistp,nphinored,relparstp,nrelpar,phistropthresh,lastrun,phiaboandon,icov,icor,ieig,ires,jcosave,verboserec \
	  ,jcosaveitn,reisaveitn,parsaveitn,parsaverun):
	self.rstfile = rstfile
	self.pestmode = pestmode
	self.npar = npar
	self.nobs = nobs
	self.npargp = npargp
	self.nprior = nprior
	self.nobsgp = nobsgp
	self.maxcompdim = maxcompdim
	self.ntplfle = ntplfle
	self.ninsfle = ninsfle
	self.precis = precis
	self.dpoint = dpoint
	self.numcom = numcom
	self.jacfile = jacfile
	self.messfile = messfile
	self.obsreref = obsreref
	self.rlambda1 = rlambda1
	self.rlamfac = rlamfac
	self.phiratsuf = phiratsuf
	self.phiredlam = phiredlam
	self.numlam = numlam
	self.jacupdate = jacupdate
	self.lamforgive = lamforgive
	self.derforgive = deforgive
	self.relparmax = relparmax
	self.facparmax = facparmax
	self.facorig = facorig
	self.iboundstick = iboundstick
	self.upvecbound = upvecbound
	self.absparmax = absparmax
	self.phiredswh = phiredswh
	self.noptmax = noptmax
	self.phiredstp = phiredstp
	self.nphistp = nphistp
	self.nphinored = nphinored
	self.relparstp = relparstp
	self.nrelpar = nrelpar
	self.phistropthresh = phistropthresh
	self.lastrun = lastrun
	self.phiabandon = phiabandon
	self.icov = icov
	self.icor = icor
	self.ieig = ieig
	self.ires = ires
	self.jcosave = jcosave
	self.verboserec = verboserec
	self.jcosaveitn = jcosaveitn
	self.reisaveitn = reisaveitn
	self.parsaveitn = parsaveitn
	self.parsaverun = rarsaverun
	return
      
class svd(object):
  '''
  Create a singular value decomposition
  
  Attributes:
    SVDMODE : 0,1 activates truncated singular value decomposition fir solution of the inverse problem
    MAXSING : number of singular values at which truncation occurs
    EIGTHRESH : <1 eiqenvalue ratio threshold for truncation
    EIGWRITE : 0,1 determines content of scd output file
  '''
  def __init__(self,svdmode,maxsing,eigthresh,eigwrite):
      self.svdmode = svdmode
      self.maxsing = maxsing
      self.eigthresh = eigthresh
      self.eigwrite = eigwrite
      return

class pargroup(object):
  '''
  create a pargroup line
  
  Attributes:
    PARGPNAME  : name of the paramater group name
    INCTYPE : "relative", "absolute" or "rel_to_max" method by which parameter increments are calculate
    DERINC : absolutre or relative parameter increment 
    FORCEN : "switch,always_2,always_3,always_3,always_5,switch_5" determines whether central derivatives calculation is undertaken
    DERINCMUL : derivative increment multiplier when undertaking centrl derivatives calculation
    DERMTHD : "parabolic,outside_pts,best_fit,minvar,maxprec" method of central derivatives calculation
  '''
  def __init__(self,pargpname,inctyp,derinc,forcen,derincmul, dermthd):
      self.pargpname = pargpname
      self.inctyp = inctyp
      self.derinc = derinc
      self.forcen = forcen
      self.derincmul = derincmul
      self.dermthd = dermthd
      return
class paramater(object):
  '''
  create a paramater line for the PCF
  
  Attributes:
    PARNME : Parameter name
    PARTRANS : "log,none,fixed,tied" paramater transformation
    PARCHGLIM :  "relative,factor,absolute(N) type of paramater change limit
    PARVAL1 : inital paramater value
    PARLBND : parameter lower bound
    PARUBND : paramater upper bound
    PARGP : paramater group name
    SCALE : multiplication factor for paramater
    OFFSET : number to addd to paramater
    DERCOM : model command line used in computing paramater increments
  '''
  def __init__(self,parnme,partrans,parchglim,parval1,parlbnd,parupbnd,pargp,scale,offset,dercom):
      self.parnme = parnme
      self.partrans = partrans
      self.parchglim = parchglim
      self.parval1 = parval1
      self.parlbnd = parlbnd
      self.parupbnd = parupbnd
      self.pargp = pargp
      self.scale = scale
      self.offset = offset
      self.dercom = dercom
      return
  
class obsgroup(object):
  '''
  create the obsgroup line
  
  Attributes:
  OBGNME  : Name of observation group
  '''
  def __init__(seld,obgnme):
    self.obgnme = obgnme 
    return
  
  
class observation(object):
   '''
   create the observation data object
    
   Attributes:
      OBSNME : observation name
      OBSVAL : observation value
      WEIGHT : observation weight
      OBGNME : group that observation belongs
  '''
  def __init__(self,obsnme,obsval,weight,obgnme):
    self.obsnme = obsnme
    self.obsval = obsval
    self.weight = weight
    self.obgnme = obgnme

class commandline(object):
  '''
  object containing method for running the model command line
  
  Attributes:
    COMMAND : Text which specifies how model is run
  '''
  def __init__(self,command):
    self.command = command
    return

class modelinout(object):
  '''
  contain the information on how to read and write model input file
  
  Attributes :
      TEMPFLE : Template File
      INFLE : Input File
      INSFLE : Instruction file for reading model output
      OUTFLE : Output File
  '''
  def __init__(self,tempfle,infle,insfle,outfle):
    self.tempfle = tempfle
    self.infle = infle
    self.insfle = insfle
    self.outfle = outfle
    return

clas prior(object):
  '''
  container for a prior information line
  
  Attributes:
    PILBL : Prior information label
    EQ : PFAC*PARNME + PFIC * log(PARNME) = PIVAL Text representing a prior information equation
    WEIGHT : weight of the prior information
    OBGPNME: Observation group that prior information belongs to 
 '''
 def __init__(self,piblbl,eq,weight,obgpnme):
    self.piblbl = piblbl
    self.eq = eq
    self.weight = weight
    self.obgpnme = obgpnme 
    return

class regularisation(object):
  '''
  Regularisation section of the PEST control file
  
  Attritutes:
    PHIMLIM : target measurement objective function
    PHIMACCEPT : acceptable measurement objective function
    WFINIT : inital regulristion weight factor
    WFMIN : minimum regularisation weight factor
    WFMAX : maximum regularisation weight factor
    WFFAC : regularisation weight factor adjjustment factor
    WFTOL : covergence criterion for regulatrisation weight factor
    IRGADJ : 0,1,2,3,4,5  instructs PEST to perform iner-regularisation group weight factor adjustments
   '''
   def __init__(self,phimlim,phimaccept,wfinit,wfmin,wfmax,wffac,wftol,irgadj):
      self.phimlim = phimlim
      self.phimaccept  = phumaccept
      self.wfinit = wfinit
      self.wfmin = wfmin 
      self.wfmax = wfmax
      self.wffac = wffac
      self.wftol = wftol
      self.irgadj = irgadj
      return
    
def createPCF(filName,controlObj,svdObj,pargpObj,parDataObjList,obsgpObj,obsObjList,commandObj,modinoutObj,piObj,regularObj):
    '''
    create a PCF by compiling the proper objects and writing them to file
    '''
    f = open(fileName, 'w')
    f.write('pcf\n')
    f.write('* control data\n')
    
