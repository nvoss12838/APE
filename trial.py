from pcf import *

control = controldata()
svd = svd()

m = pargroup('m','absolute',0.01,0.0,'switch' ,2.0, 'parabolic')
b = pargroup('b','absolute',0.01,0.0,'switch',2.0,'parabolic')
a1 = pargroup('a1','relative',0.01,0.0,'switch',2.0,'parabolic')
a2 = pargroup('a2','relative',0.0,0.0,'switch',2.0,'parabolic')
a3 = pargroup('a3','relative',0.01,0.0,'switch',2.0,'parabolic')
a4 = pargroup('a4','relative',0.01,0.0,'switch',2.0,'parabolic')
g = pargroup('g','relative',0.01,0.0,'switch' ,2.0, 'parabolic')
h = pargroup('h','absolute',0.01,0.0,'switch',2.0,'parabolic')
c = pargroup('c','relative',0.01,0.0,'switch',2.0,'parabolic')
cTau = pargroup('cTau','relative',0.01,0.0,'switch' ,2.0, 'parabolic')
d = pargroup('d','relative',0.01,0.0,'switch',2.0,'parabolic')
dTau = pargroup('dTau','relative',0.01,0.0,'switch' ,2.0, 'parabolic')
e = pargroup('e','relative',0.01,0.0,'switch',2.0,'parabolic')
eTau = pargroup('eTau','relative',0.01,0.0,'switch' ,2.0, 'parabolic')
t01 = pargroup('t01','relative',0.01,0.0,'switch' ,2.0, 'parabolic')
t1= pargroup('t1','relative',0.01,0.0,'switch',2.0,'parabolic')
u1 = pargroup('u1','relative',0.01,0.0,'switch',2.0,'parabolic')
t02 = pargroup('t02','relative',0.01,0.0,'switch' ,2.0, 'parabolic')
t2= pargroup('t2','relative',0.01,0.0,'switch',2.0,'parabolic')
u2 = pargroup('u2','relative',0.01,0.0,'switch',2.0,'parabolic')
t03 = pargroup('t03','relative',0.01,0.0,'switch' ,2.0, 'parabolic')
t3= pargroup('t3','relative',0.01,0.0,'switch',2.0,'parabolic')
u3 = pargroup('u3','relative',0.01,0.0,'switch',2.0,'parabolic')
t04 = pargroup('t04','relative',0.0,0.0,'switch',2.0,'parabolic')
t4 = pargroup('t4','relative',0.01,0.0,'switch',2.0,'parabolic')
u4 = pargroup('u4','relative',0.01,0.0,'switch',2.0,'parabolic')
t05 = pargroup('t05','relative',0.0,0.0,'switch',2.0,'parabolic')
t5 = pargroup('t5','relative',0.01,0.0,'switch',2.0,'parabolic')
u5 = pargroup('u5','relative',0.01,0.0,'switch',2.0,'parabolic')
t06 = pargroup('t06','relative',0.0,0.0,'switch',2.0,'parabolic')
t6 = pargroup('t6','relative',0.01,0.0,'switch',2.0,'parabolic')
u6 = pargroup('u6','relative',0.01,0.0,'switch',2.0,'parabolic')
t07 = pargroup('t07','relative',0.0,0.0,'switch',2.0,'parabolic')
t7 = pargroup('t7','relative',0.01,0.0,'switch',2.0,'parabolic')
u7 = pargroup('u7','relative',0.01,0.0,'switch',2.0,'parabolic')
t08 = pargroup('t08','relative',0.0,0.0,'switch',2.0,'parabolic')
t8 = pargroup('t8','relative',0.01,0.0,'switch',2.0,'parabolic')
u8 = pargroup('u8','relative',0.01,0.0,'switch',2.0,'parabolic')
t09 = pargroup('t08','relative',0.0,0.0,'switch',2.0,'parabolic')
t9 = pargroup('t8','relative',0.01,0.0,'switch',2.0,'parabolic')
u9 = pargroup('u8','relative',0.01,0.0,'switch',2.0,'parabolic')

pargp = [m,b,a1,a2,a3,a4,g,h,c,cTau,d,dTau,e,eTau,t01,t1,u1,t02,t2,u2,t03,t3,u3,t04,t4,u4,t05,t5,u5,t06,t6,u6,t07,t7,u7,t08,t8,u8,t09,t9,u9]

pm  = paramater('m','none','relative',0.0,-1.000000E+10,1.000000E+10, 'm' , 1.0000,0.0000,1)
pb = paramater('b','none','relative',0.31,-1.000000E+10,1.0,'b',1.0000,0.0000,1)
pa1 = paramater('a1','none','relative',0.121970,-1.000000E+10,1.000000E+10,'a1',1.0000,0.0000,1)
pa2 = paramater('a2','none','relative',3.449000E-02,-1.000000E+10,1.000000E+10,'a2',1.0000,0.0000,1)
pa3 = paramater('a3','none','relative',-1.640000E-02,-1.000000E+10,1.000000E+10,'a3',1.0000,0.0000,1)
pa4 = paramater('a4','none','relative',2.038000E-02,-1.000000E+10,1.000000E+10,'a4',1.0000,0.0000,1)
pg = paramater('g','none','relative',-19.9829,-1.000000E+10,1.000000E+10,'g',1.0000,0.0000,1)
ph = paramater('h','none','relative',-2.71254,-1.000000E+10,0.0,'h',1.0000,0.0000,1)
pc = paramater('c','none','relative',-5.29094,-1.000000E+10,1.000000E+10,'c',1.0000,0.0000,1)
cTau = paramater('cTau','log','relative',-19.9829,-1.000000E+10,1.000000E+10,'cTau',1.0000,0.0000,1)
pd = paramater('d','none','relative',-7.14580,-1.000000E+10,1.000000E+10,'d',1.0000,0.0000,1)
dTau = paramater('dTau','log','relative',-19.9829,-1.000000E+10,1.000000E+10,'dTau',1.0000,0.0000,1)
pe = paramater('e','none','relative',-17.2090,-1.000000E+10,1.000000E+10,'e',1.0000,0.0000,1)
eTau = paramater('eTau','log','relative',-19.9829,-1.000000E+10,1.000000E+10,'eTau',1.0000,0.0000,1)
pt01 = paramater('t01','log','relative',2011.49,2011.45,2014.55,'t01',1.0000,0.0000,1)
pt1 = paramater('t1','log','relative',20.0,10.0,30.0,'t1',1.0000 ,0.0000 ,1)
pu1 = paramater('u1','none','relative',-0.24,-0.6,-0.0,'u1',1.0000,0.0000,1)
pt02 = paramater('t02','log','relative',2011.49,2011.45,2014.55,'t02',1.0000,0.0000,1)
pt2 = paramater('t2','log','relative',20.0,10.0,30.0,'t2',1.0000 ,0.0000 ,1)
pu2 = paramater('u2','none','relative',-0.24,-0.6,-0.0,'u2',1.0000,0.0000,1)
pt03 = paramater('t03','log','relative',2011.49,2011.45,2014.55,'t03',1.0000,0.0000,1)
pt3 = paramater('t3','log','relative',20.0,10.0,30.0,'t3',1.0000 ,0.0000 ,1)
pu3 = paramater('u3','none','relative',-0.24,-0.6,-0.0,'u3',1.0000,0.0000,1)
pt04 = paramater('t04','log','relative',2014.14,2014.05,2014.4 ,'t04',1.0000,0.0000,1)
pt4 = paramater('t4','log','relative',30.0,0.01,80.0,'t4',1.0000,0.0000,1)
pu4 = paramater('u4','none','relative',-1.15363,-1.000000E+10,1.000000E+10,'u4',1.0000,0.0000,1)
pt05 = paramater('t05','log','relative',2011.49,2011.45,2014.55,'t05',1.0000,0.0000,1)
pt5 = paramater('t5','log','relative',20.0,10.0,30.0,'t5',1.0000 ,0.0000 ,1)
pu5 = paramater('u5','none','relative',-0.24,-0.6,-0.0,'u5',1.0000,0.0000,1)
pt06 = paramater('t06','log','relative',2011.49,2011.45,2014.55,'t06',1.0000,0.0000,1)
pt6 = paramater('t6','log','relative',20.0,10.0,30.0,'t6',1.0000 ,0.0000 ,1)
pu6 = paramater('u6','none','relative',-0.24,-0.6,-0.0,'u6',1.0000,0.0000,1)
pt07 = paramater('t07','log','relative',2011.49,2011.45,2014.55,'t07',1.0000,0.0000,1)
pt7 = paramater('t7','log','relative',20.0,10.0,30.0,'t7',1.0000 ,0.0000 ,1)
pu7 = paramater('u7','none','relative',-0.24,-0.6,-0.0,'u7',1.0000,0.0000,1)
pt08 = paramater('t08','log','relative',2011.49,2011.45,2014.55,'t08',1.0000,0.0000,1)
pt8 = paramater('t8','log','relative',20.0,10.0,30.0,'t8',1.0000 ,0.0000 ,1)
pu8 = paramater('u8','none','relative',-0.24,-0.6,-0.0,'u8',1.0000,0.0000,1)
pt09 = paramater('t09','log','relative',2011.49,2011.45,2014.55,'t09',1.0000,0.0000,1)
pt9 = paramater('t9','log','relative',20.0,10.0,30.0,'t9',1.0000 ,0.0000 ,1)
pu9 = paramater('u9','none','relative',-0.24,-0.6,-0.0,'u9',1.0000,0.0000,1)


pars = [pm,pb,pa1,pa2,pa3,pa4,pg,ph,pc,pd,pe,pt01,pt1,pu1,,pt02,pt2,pu2,pt03,pt3,pu3,pt04,pt4,pu4,pt05,pt5,pu5,pt06,pt6,pu6,pt07,pt7,pu7,pt08,pt8,pu8,pt09,pt9,pu9]

obgrp = [obsgroup('obsgroup')]

obs = []
import detrendDeaftereq as cr
t,l,u = cr.readData('../gpsFIT/Data/bija/BIJA.lon')
#f.write('l1 w !obs1! \n')
for i in range(len(t)):
	obs.append(observation(' obs' + str(i+1),l[i],np.min(u)/u[i],'obsgroup'))

command = [commandline('python pestTimeSeries2.py cr_2014long/PNE2.rad paramsFile.txt')]
inout = modelinout('paramsFile.tp1','paramsFile.txt','model.ins','model.txt')
pi = [prior('','','','')]

reg = regularisation('','','','','','','','')
control.npar = len(pars)
control.nobs = len(obs)
control.npargp = len(pargp)
control.nprior = len(pi)
createPCF('../gpsFIT/Data/bija/bijalon.pst',control,svd,pargp,pars,obgrp,obs,command,inout,pi,reg)