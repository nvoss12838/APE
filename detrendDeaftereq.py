import numpy as np
import matplotlib.pylab as plt
from scipy.optimize import leastsq
import os 

def readData(filename):
  '''
  read in data in Decater format 
  date location uncertainty?
  '''
  #read in filename splitting into useful python arrays 
  date,loc,uncert = np.genfromtxt(filename,dtype=float,usecols = (0,1,2), unpack = True)
  dateNum = [float(d) for d in date]
  loc = [float(l) for l in loc]
  uncert = [float(x) for x in uncert]
  
  return dateNum,loc,uncert

def heavy(t,t0):
   if t<t0:
	 x=0
   else:
	 x = 1
   return x 
def heavyInv(t,t0):
   if t<t0:
	 x=1
   else:
	 x = 0
   return x 
def trend(t,m,b):
    #fit trend : y = mx+b 
    y = m*t + b #m = trend and b is the bias 
    return y
def period(t,m,b,A1,A2,A3,A4):
    #semi-annual and annual periods trend should be removed first
    y = A1*np.cos(np.radians(2.0*np.pi*t)) + A2*np.sin(np.radians(2.0*pi*t))+A3*np.cos(np.radians(4.0*np.pi*t)) + A4*np.sin(np.radians(4.0*pi*t))
    return y 
def jump(t,t0,G):
    #fit for jump (either antenna change or earthquake
    y = G*heavy(t,t0)
    return y 

def exponential(t,C,eqt,tau):
    #fit the postsiesmic deformation
    #we could also do a logLarithmic (aught to implement in future)
    y = C*(1.0-np.exp(-(t-eqt)/(tau/365.0)))
    return y 

def fun(t,B,C,tau1,D,tau2,E,tau3):#t02,tau3,U):#t,G,E,
    y = []
    #fit arctan as approximation for sse from Holtkamp and Brudzinski 2008
    #t01 is the median time of the sse
    #tau is the period over which the sse takes place.
    eqt1 = 2012.6817
    eqt2 = 2012.8159
    t0 = 2014.12
    dur = 0.3
    #tau1 = 7.0
    #tau2 = 70.0
    #tau3 = 420.0
    for i in range(len(t)):
      #y.append(jump(t[i],eqt1,G)+jump(t[i],eqt2,E) + trend(t[i],m,b)*heavyInv(t[i],eqt1))#+sse(t[i],t02,tau3,U))
	  y.append(jump(t1[i],eqt2,B)+exponential(t1[i],C,eqt1,tau1)*heavy(t1[i],eqt1)+exponential(t1[i],D,eqt2,tau1)*heavy(t1[i],eqt1)+exponential(t1[i],E,eqt1,tau3)*heavy(t1[i],eqt1))#*heavyInv(t[i],2013.081))
    return y

def residuals(par):
	eqt1 = 2012.6817
	eqt2 = 2012.8159
	#tau1 = 7.0
	#tau2 = 70.0
	#tau3 = 420.0
	res = x-fun(t1,par[0],par[1],par[2],par[3],par[4],par[5],par[6])
	return res
def resTrend(par):
	eqt1 = 2012.6817
	eqt2 = 2012.8159
	t0 = 2014.12
	dur=0.3
	#tau1 = 7.0
	#tau2 = 70.0
	#tau3 = 420.0
	y = []
	for i in range(len(t1)):
	  mean = np.mean(x)
	 # print mean
      #y.append(jump(t[i],eqt1,G)+jump(t[i],eqt2,E) + trend(t[i],m,b)*heavyInv(t[i],eqt1))#+sse(t[i],t02,tau3,U))
	  y.append(x[i]-jump(t1[i],eqt2,par[0])-exponential(t1[i],par[1],eqt1,par[2])*heavy(t1[i],eqt1)-exponential(t1[i],par[3],eqt2,par[4])*heavy(t1[i],eqt1)-exponential(t1[i],par[5],eqt1,par[6])*heavy(t1[i],eqt1))#+sse(t1[i],dur,t0,par[7]))
	return y

def sse(t,t0,tau,U):
    #fit arctan as approximation for sse from Holtkamp and Brudzinski 2008
    #t01 is the median time of the sse
    #tau is the period over which the sse takes place. 
    y = 0.5*U*((np.tanh(t-t0)/(tau)-1))
    return y 

def getFit(station):
  global t1
  global x
  date,locat,uncert = readData(station)
  #dateOrd = [ dt.date.toordinal(i) for i in date]
  #x = [(i - np.mean(x))*-1.0 for i in x]
  #these are our initial guesses 
  #t,G,C,tau,D,tau2,E,t02,tau3,U
  A = 0.0  
  B = -.5 
  C = -3.0  
  D = 0.0 # 
  E = 2.5
  m = -0.3
  b = 0.09
  tau1 = 0.1
  tau2 = 1.5
  tau3 = 38.0
  t0 = 2014.15
  dur = 0.26
  U = 3.0
  t1= np.array(date)
  x = np.array(locat)
  time,loc =[],[]
  for i in range(len(t1)):
  	  if t1[i]<2013.5:
  	  	time.append(t1[i])
  	  	loc.append(x[i])
  t1 = np.array(time)
  x = np.array(loc)

#G,E,C,tau,D,tau2,t02,tau3,U
  fitParams, fitCov = leastsq(residuals,[B,C,tau1,D,tau2,E,tau3])#,D,E,tau1,tau2,tau3])#t02,tau3,U,G,E])
  print fitParams
  t1 = np.array(date)
  plt.plot(t1,fun(t1,fitParams[0],fitParams[1],fitParams[2],fitParams[3],fitParams[4],fitParams[5],fitParams[6]))
  plt.scatter(t1,locat)
  plt.show()
  x = np.array(locat)
  cleaned = resTrend(fitParams)
  return t1,cleaned,uncert

def plotData(date,loc,uncert,fil,dir):
	os.chdir(dir)
	fig = plt.figure()
	plt.scatter(date,loc)
	uncert = [u*3 for u in uncert]
	plt.errorbar(date, loc, uncert, ls='none')
	plt.savefig(fil + '.png')
	plt.close(fig)
	os.chdir('..')
	return

def writeData(t,x,u,dir,filname):
	#print os.getcwd()
	os.chdir(dir)
	#print os.getcwd()
	#f = open(filname, 'w')
	with open(filname, 'w') as f:
		for i in range(len(x)):
			f.write(str(t[i]) + ' ' + str(x[i]) + ' ' + str(u[i]) + '\n')
	#f.close()
	os.chdir('..')
	return