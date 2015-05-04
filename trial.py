from pcf import *

control = controldata()
svd = svd()

b = pargroup('b','absolute',0.01,0.0,'switch',2.0,'parabolic')
a1 = pargroup('a1','relative',0.01,0.0,'switch',2.0,'parabolic')
a2 = pargroup('a2','relative',0.0,0.0,'switch',2.0,'parabolic')
a3 = pargroup('a3','relative',0.01,0.0,'switch',2.0,'parabolic')
a4 = pargroup('a4','relative',0.01,0.0,'switch',2.0,'parabolic')
g = pargroup('g','relative',0.01,0.0,'switch' ,2.0, 'parabolic')
h = pargroup('h','absolute',0.01,0.0,'switch',2.0,'parabolic')
c = pargroup('c','relative',0.01,0.0,'switch',2.0,'parabolic')
d = pargroup('d','relative',0.01,0.0,'switch',2.0,'parabolic')
e = pargroup('e','relative',0.01,0.0,'switch',2.0,'parabolic')
t04 = pargroup('t04','relative',0.0,0.0,'switch',2.0,'parabolic')
t4 = pargroup('t4','relative',0.01,0.0,'switch',2.0,'parabolic')
u4 = pargroup('u4','relative',0.01,0.0,'switch',2.0,'parabolic')
t01 = pargroup('t01','relative',0.01,0.0,'switch' ,2.0, 'parabolic')
t1= pargroup('t1','relative',0.01,0.0,'switch',2.0,'parabolic')
u1 = pargroup('u1','relative',0.01,0.0,'switch',2.0,'parabolic')
m = pargroup('m','absolute',0.01,0.0,'switch' ,2.0, 'parabolic')

pargp = [b,a1,a2,a3,a4,g,h,c,d,e,t04,t4,u4,t01,t1,u1,m]

pb = paramater('none','relative',0.31,-1.000000E+10,1.0,'b',1.0000,0.0000,1,0)
pa1 = paramater('none','relative',0.121970,-1.000000E+10,1.000000E+10,'a1',1.0000,0.0000,1,0)
pa2 = paramater('none','relative',3.449000E-02,-1.000000E+10,1.000000E+10,'a2',1.0000,0.0000,1,0)
pa3 = paramater('none','relative',-1.640000E-02,-1.000000E+10,1.000000E+10,'a3',1.0000,0.0000,1,0)
pa4 = paramater('none','relative',2.038000E-02,-1.000000E+10,1.000000E+10,'a4',1.0000,0.0000,1,0)
pg = paramater('none','relative',-19.9829,-1.000000E+10,1.000000E+10,'g',1.0000,0.0000,1,0)
ph = paramater('none','relative',-2.71254,-1.000000E+10,0.0,'h',1.0000,0.0000,1,0)
pc = paramater('none','relative',-5.29094,-1.000000E+10,1.000000E+10,'c',1.0000,0.0000,1,0)
pd = paramater('none','relative',-7.14580,-1.000000E+10,1.000000E+10,'d',1.0000,0.0000,1,0)
pe = paramater('none','relative',-17.2090,-1.000000E+10,1.000000E+10,'e',1.0000,0.0000,1,0)
pt04 = paramater('none','relative',2014.14,2014.05,2014.4 ,'t04',1.0000,0.0000,1,0)
pt4 = paramater('none','relative',30.0,0.01,80.0,'t4',1.0000,0.0000,1,0)
pu4 = paramater('none','relative',-1.15363,-1.000000E+10,1.000000E+10,'u4',1.0000,0.0000,1,0)
pt01 = paramater('none','relative',2011.49,2011.45,2014.55,'t01',1.0000,0.0000,1,0)
pt1 = paramater('none','relative',20.0,10.0,30.0,'t1',1.0000 ,0.0000 ,1,0)
pu1 = paramater('none','relative',-0.24,-0.6,-0.0,'u1',1.0000,0.0000,1,0)
pm  = paramater('none','relative',0.0,-1.000000E+10,1.000000E+10, 'm' , 1.0000,0.0000,1,0)

pars = [pb,pa1,pa2,pa3,pa4,pg,ph,pc,pd,pe,pt04,pt4,pu4,pt01,pt1,pu1,pm]

obgrp = [obsgroup('obsgroup')]

obs = []
import detrendDeaftereq as cr
t,l,u = cr.readData('../Geodesy/CR_2014/cr_2014long/IND1.lon')
#f.write('l1 w !obs1! \n')
for i in range(len(t)):
	obs.append(observation(' obs' + str(i+1),l[i],np.min(u)/u[i],'obsgroup'))

command = [commandline('python pestTimeSeries2.py cr_2014long/PNE2.rad paramsFile.txt')]
inout = modelinout('paramsFile.tp1','paramsFile.txt','model.ins','model.txt')
pi = [prior('','','','')]

reg = regularisation('','','','','','','','')

createPCF('trial.pst',control,svd,pargp,pars,obgrp,obs,command,inout,pi,reg)