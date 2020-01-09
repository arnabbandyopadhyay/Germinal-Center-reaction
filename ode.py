#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 09:41:47 2019

@author: abp19
"""

import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from lmfit import minimize, Parameters, Parameter, report_fit
#from scipy.optimize import minimize
# function that returns dy/dt
def model(dydt,t):
    bn, bgc, bm, bpc, uc, ic, v, nk = dydt
    if t>168+72:
        k=0
    else:
        k=1
    if t<168:
        kk=0
    else:
        kk=1
    s=(bm)*k*kk
    kmb=1
    km = 0.01
    kdm=0.1
    kv=10
    kmp=0.005
    p=0.05
    kv=0.001
    knk=0.1
    kdv=0.01
    kdv2=0.0001
    kdv3=0.1
    kin=0.003
    kdi=0.1
    kdnk=0.01
    eq = [s-p*bn, 2*p*bn-km*bgc-kmp*bgc, kmb*v/(kv+v) +km*bgc-kdm*bm-kdv2*bm*ic, kmp*bgc-kdv3*bpc*ic, -kin*uc*v, kin*uc*v -kdi*ic*nk-kdv2*bm*ic-kdv3*bpc*ic, kv*ic-kdv*v, knk-kdnk*nk-kdi*nk*ic]
    return eq

# initial condition
bn0=1
bgc0=0
bm0=0
bpc0=0
uc0=1000
ic0=0
v0=1
nk0=100
y0 = [bn0,bgc0,bm0,bpc0,uc0,ic0,v0,nk0]

# time points
t = np.linspace(0,500,1000)

# solve ODE
y = odeint(model,y0,t)

t=t/24

# plot results
plt.plot(t,y[:,0],'b',label='naive')
plt.plot(t,y[:,1],'g',label='Bgc')
plt.plot(t,y[:,2],'r',label='Bm')
plt.plot(t,y[:,3],'k',label='Bpc')

plt.xlabel('time')
plt.ylabel('y(t)')
plt.legend()
plt.show()


# plot results
plt.plot(t,y[:,4],'b',label='Uc')
plt.plot(t,y[:,5],'g',label='Ic')
plt.xlabel('time')
plt.ylabel('y(t)')
plt.legend()
plt.show()

# plot results
plt.plot(t,y[:,6],'b',label='V')
plt.plot(t,y[:,7],'g',label='NK')
plt.xlabel('time')
plt.ylabel('y(t)')
plt.legend()
plt.show()





def model(dydt,t):
    bn, bgc, bmigm, bmigg, bpc, uc, ic, v, nk = dydt

    if t<300:
        kk=1
    else:
        kk=0
        
    s=1*(bmigm)*kk
      
    kmb=0.01;kmb2=0.001;km = 0.0005;kdm=0.01;kmp=0.0005;p=0.05;
    kvm=800;kv=0.01;knk=0.1;kdv=0.01;kdv2=0.0001;kdv3=0.1;kin=0.0005;kdi=0.01;kdnk=0.01;
    
    dbn = s - p*bn
    dbgc = 2*p*bn - km*bgc - kmp*bgc
    dbmigm = kmb2*bgc + kmb*v/(kvm+v)  - kdm*bmigm 
    dbmigg = km*bgc - kdv2*bmigg*ic
    dbpc = kmp*bgc - kdv3*bpc*ic
    duc = - kin*uc*v
    dic = kin*uc*v - kdi*ic*nk - kdv2*bmigg*ic - kdv3*bpc*ic
    dv = kv*ic - kdv*v
    dnk = knk - kdnk*nk - kdi*nk*ic
    
    eq = [dbn, dbgc, dbmigm, dbmigg, dbpc, duc, dic, dv, dnk]
    return eq

# initial condition
bn0=1
bgc0=0
bmigm0=0
bmigg0=0
bpc0=0
uc0=1000
ic0=0
v0=1
nk0=100
y0 = [bn0,bgc0,bmigm0,bmigg0,bpc0,uc0,ic0,v0,nk0]

# time points
t = np.linspace(0,1000,5000)

# solve ODE
y = odeint(model,y0,t)

t=t/24

def plott(n,pos,col,label):
    for i in range(n):
        plt.plot(t,y[:,pos[i]],col[i],label=label[i])
        plt.xlabel('time')
        plt.ylabel('y(t)')
        plt.legend()


# plot results
plt.plot(t,y[:,0],'b',label='naive')
plt.plot(t,y[:,1],'g',label='Bgc')

plt.xlabel('time')
plt.ylabel('y(t)')
plt.legend()
plt.show()

plt.plot(t,y[:,2],'r',label='Bmigm')
plt.plot(t,y[:,3],'m',label='Bmigg')
plt.plot(t,y[:,4],'k',label='Bpc')

plt.xlabel('time')
plt.ylabel('y(t)')
plt.legend()
plt.show()


# plot results
plt.plot(t,y[:,5],'b',label='Uc')
plt.plot(t,y[:,6],'g',label='Ic')
plt.xlabel('time')
plt.ylabel('y(t)')
plt.legend()
plt.show()

# plot results
plt.plot(t,y[:,7],'b',label='V')
plt.plot(t,y[:,8],'g',label='NK')
plt.xlabel('time')
plt.ylabel('y(t)')
plt.legend()
plt.show()

n=3
pos=[2,3,4]
col=['r','g','b','m','k']
label=['Bmigm','Bmigg','Bpc']


        

fig = plt.figure()
plott(n,pos,col,label)
fig = plt.figure()
plott(n,pos,col,label)







def modelN(dydt,t,params):
    bn, bgc, bmigm, bmigg, bpigm, bpigg, uc, ic, v, nk = dydt

    if t<300:
        kk=1
    else:
        kk=0
        
    s=1*(bmigm)*kk
    p=0.2;

    kmigm = params[0]; kmigg = params[1]; kpigm = params[2]; kpigg = params[3]; kdm=params[4];
    
#    kswim = 0.01; kswimm = 1; kswig = 0.01; kswigm = 1;
    kswim = params[5]; kswimm = params[6]; kswig =params[7]; kswigm = params[8];
#    kcigm=0.01; kcigg=0.01;
    kcigm=params[9];kcigg=params[10];
#    kmb=0.01;kmb2=0.001;kmp=0.0005;
#    kin=0.0005;kdi=0.01; kv=0.01;kdv=0.01; knk=0.1;kdnk=0.01;
    kin=params[11];kdi=params[12];kv=params[13];kdv=params[14];knk=params[15];kdnk=params[16];
    
    dbn = s - p*bn
    dbgc = 2*p*bn - kmigm*bgc - kmigg*bgc - kpigm*bgc - kpigg*bgc
    dbmigm = kmigm*bgc  - kdm*bmigm - kswim*v*bn/(kswimm+v) 
    dbmigg = kmigg*bgc - kdm*bmigg - kswig*v*bmigg/(kswigm+v)
    dbpigm = kpigm*bgc + kswim * bn*v/(kswimm+v) - kdm*bpigm  - kcigm*bpigm*ic
    dbpigg = kpigg*bgc - kdm*bpigg + kswig*v*bmigg/(kswigm+v) - kcigg*bpigg*ic
    duc = - kin*uc*v
    dic = kin*uc*v - kdi*ic*nk - kcigm*bpigm*ic - kcigg*bpigg*ic
    dv = kv*ic - kdv*v
    dnk = knk - kdnk*nk - kdi*nk*ic
    
    eq = [dbn, dbgc, dbmigm, dbmigg, dbpigm, dbpigg, duc, dic, dv, dnk]
    return eq


# initial condition
bn0=1
bgc0=0
bmigm0=0
bmigg0=0
bpigm0=0
bpigg0=0
uc0=1000
ic0=0
v0=1
nk0=100
y0 = [bn0,bgc0,bmigm0,bmigg0,bpigm0,bpigg0,uc0,ic0,v0,nk0]

##good param set 
#p=0.2;s=1;
#kmigm = 0.001; kmigg = 0.0005; kpigm = 0.0003; kpigg = 0.0001; kdm=0.01;
#kswim = 0.01; kswimm = 10; kswig =0.01; kswigm = 10;
#kcigm=kcigg=0.00001;
#kin=0.0005;kdi=0.01; kv=0.01;kdv=0.01; knk=0.1;kdnk=0.01;





kmigm = 0.001; kmigg = 0.0005; kpigm = 0.0003; kpigg = 0.0001; kdm=0.01;
    
#    kswim = 0.01; kswimm = 1; kswig = 0.01; kswigm = 1;
kswim = 0.01; kswimm = 10; kswig =0.01; kswigm = 10;
#    kcigm=0.01; kcigg=0.01;
kcigm=0.00001;kcigg=0.00005
#    kmb=0.01;kmb2=0.001;kmp=0.0005;
kin=0.0005;kdi=0.01; kv=0.01;kdv=0.01; knk=0.1;kdnk=0.01;
#kin=kdi=kv=kdv=knk=kdnk=0.0;

params=[kmigm, kmigg,kpigm, kpigg, kdm, kswim, kswimm, kswig, kswigm, kcigm, kcigg,kin, kdi, kv, kdv, knk, kdnk]




# time points
t = np.linspace(0,1000,5000)

# solve ODE
y = odeint(modelN,y0,t,args=(params,))

t=t/24
n=2;pos=[0,1]
col=['r','g','b','m','k']
label=['naive','Bgc']
fig = plt.figure()
plott(n,pos,col,label)

n=4;pos=[2,3,4,5]
label=['Bmigm','Bmigg','Bpigm','Bpigg']
fig = plt.figure()
plott(n,pos,col,label)

n=2;pos=[6,7]
label=['Uc','Ic']
fig = plt.figure()
plott(n,pos,col,label)   

n=2;pos=[8,9]
label=['V','Nk']
fig = plt.figure()
plott(n,pos,col,label)         

def modelN2(y, t, paras):
    """
    Your system of differential equations
    """

    bn = y[0]
    bgc = y[1]
    bmigm = y[2]
    bmigg = y[3]
    bpigm = y[4]
    bpigg = y[5]
    uc = y[6]
    ic = y[7]
    v = y[8]
    nk = y[9]
    
    if t<300:
        kk=1
    else:
        kk=0
        
    s=1*(bmigm)*kk
    pp=0.2;
        
#    kmigm = 0.001; kmigg = 0.0005; kpigm = 0.0003; kpigg = 0.0001; 
    kdm=0.01;
    
#    kswim = 0.01; kswimm = 1; kswig = 0.01; kswigm = 1;
#    kswim = 0.01; kswimm = 10; kswig =0.01; kswigm = 10;
#    kcigm=0.01; kcigg=0.01;
    kcigm=0.00001;kcigg=0.00005
#    kmb=0.01;kmb2=0.001;kmp=0.0005;
    kin=0.0005;kdi=0.01; kv=0.01;kdv=0.01; knk=0.1;kdnk=0.01;

    try:
        kmigm = paras['kmigm'].value
        kmigg = paras['kmigg'].value
        kpigm = paras['kpigm'].value
        kpigg = paras['kpigg'].value
        kswim = paras['kswim'].value
        kswimm = paras['kswimm'].value
        kswig = paras['kswig'].value
        kswigm = paras['kswigm'].value

    except KeyError:
        kmigm, kmigg, kpigm, kpigg, kswim, kswimm, kswig, kswigm = paras

#    kmigm, kmigg, kpigm, kpigg = paras
    # the model equations
    dbn = s - pp*bn
    dbgc = 2*pp*bn - kmigm*bgc - kmigg*bgc - kpigm*bgc - kpigg*bgc
    dbmigm = kmigm*bgc  - kdm*bmigm - kswim*v*bn/(kswimm+v) 
    dbmigg = kmigg*bgc - kdm*bmigg - kswig*v*bmigg/(kswigm+v)
    dbpigm = kpigm*bgc + kswim * bn*v/(kswimm+v) - kdm*bpigm  - kcigm*bpigm*ic
    dbpigg = kpigg*bgc - kdm*bpigg + kswig*v*bmigg/(kswigm+v) - kcigg*bpigg*ic
    duc = - kin*uc*v
    dic = kin*uc*v - kdi*ic*nk - kcigm*bpigm*ic - kcigg*bpigg*ic
    dv = kv*ic - kdv*v
    dnk = knk - kdnk*nk - kdi*nk*ic
    
    return [dbn, dbgc, dbmigm, dbmigg, dbpigm, dbpigg, duc, dic, dv, dnk]





def g(t, x0, paras):
    """
    Solution to the ODE x'(t) = f(t,x,k) with initial condition x(0) = x0
    """
    x = odeint(modelN2, x0, t, args=(paras,))
    return x


def residual(paras, t, data):

    """
    compute the residual between actual data and fitted data
    """

#    x0 = paras['bn0'].value, paras['bgc0'].value, paras['bm0'].value, paras['bpc0'].value, paras['uc0'].value, paras['ic0'].value,paras['v0'].value, paras['nk0'].value
    model = g(t, y0, paras)

    # you only have data for one of your variables
    x2_model = model[:, 2]
    return (x2_model - data).ravel()


# measured data
t_measured = np.array([1, 4, 7, 10, 13, 16, 19, 22, 25, 28, 31, 34, 37])
t_measured = t_measured*24
#x2_measured = np.array([0.000, 0.416, 0.489, 0.595, 0.506, 0.493, 0.458, 0.394, 0.335, 0.309])

igmm = np.array([1, 3.69, 8.98, 12.88, 16.27, 19.37, 22.96, 26.46, 29.55, 32.15, 34.45, 36.05, 37.46])
iggm = np.array([0.40,1.80,5.00,8.59,12.39, 16.78, 22.56, 27.75, 33.14, 37.73, 41.63, 45.12, 48.72])
iggpc = np.array([0.20, 0.71,1.81, 3.31, 4.82, 6.12, 8.12, 10.82, 14.11, 18.21, 22.4, 26.89, 31.48])

#data=[igmm,iggm,iggpc]

plt.figure()
plt.scatter(t_measured, igmm, marker='o', color='b', label='measured data', s=75)

# set parameters including bounds; you can also fix parameters (use vary=False)
params = Parameters()
params.add('bn0', value=bn0, vary=False)
params.add('bgc0', value=bgc0, vary=False)
params.add('bmigm0', value=bmigm0, vary=False)
params.add('bmigg0', value=bmigg0, vary=False)
params.add('bpigm0', value=bpigm0, vary=False)
params.add('bpigg0', value=bpigg0, vary=False)
params.add('uc0', value=uc0, vary=False)
params.add('ic0', value=ic0, vary=False)
params.add('v0', value=v0, vary=False)
params.add('nk0', value=nk0, vary=False)
params.add('kmigm', value=0.01, min=0.00001, max=1000)
params.add('kmigg', value=0.02, min=0.00001, max=1000)
params.add('kpigm', value=0.03, min=0.00001, max=1000)
params.add('kpigg', value=0.03, min=0.00001, max=1000)
params.add('kswim', value=0.01, min=0.00001, max=1000)
params.add('kswimm', value=0.02, min=0.00001, max=1000)
params.add('kswig', value=0.03, min=0.00001, max=1000)
params.add('kswigm', value=0.03, min=0.00001, max=1000)

# fit model
result = minimize(residual, params, args=(t_measured, igmm), method='leastsq')  # leastsq nelder
# check results of the fit
data_fitted = g(np.linspace(0., 1000., 1000), y0, result.params)

# plot fitted data
plt.plot(np.linspace(0., 1000., 1000), data_fitted[:, 1], '-', linewidth=2, color='red', label='fitted data')
plt.legend()
plt.xlim([0, max(t_measured)])
plt.ylim([0, 50])
#plt.ylim([0, 1.1 * max(data_fitted[:, 1])])
# display fitted statistics
report_fit(result)

plt.show()

##1
#def f(y, t, paras):
#    """
#    Your system of differential equations
#    """
#
#    bn = y[0]
#    bgc = y[1]
#    bm = y[2]
#    bpc = y[3]
#    uc = y[4]
#    ic = y[5]
#    v = y[6]
#    nk = y[7]
#    
#    if t>72:
#        k=0
#    else:
#        k=1
#        
#    kv=0.001
#    knk=0.1
#    kdv=0.01
#    kin=0.003
#    kdi=0.1
#
#    try:
#        ky = paras['ky'].value
#        km = paras['km'].value
#        kmp = paras['kmp'].value
#        p = paras['p'].value
#
#    except KeyError:
#        ky, km, kmp, p = paras
#    # the model equations
#    bnp = ky*k-p*bn
#    bgcp = 2*p*bn-km*bgc -kmp*bgc
#    bmp = km*bgc
#    bpcp = kmp*bgc
#    ucp = -kin*uc*v
#    icp = kin*uc*v -kdi*ic*nk
#    vp = kv*ic-kdv*v
#    nkp = knk-kdi*nk*ic
#    
#    return [bnp, bgcp, bmp, bpcp, ucp, icp, vp, nkp]
#
#def g(t, x0, paras):
#    """
#    Solution to the ODE x'(t) = f(t,x,k) with initial condition x(0) = x0
#    """
#    x = odeint(f, x0, t, args=(paras,))
#    return x
#
#
#def residual(paras, t, data):
#
#    """
#    compute the residual between actual data and fitted data
#    """
#
#    x0 = paras['bn0'].value, paras['bgc0'].value, paras['bm0'].value, paras['bpc0'].value, paras['uc0'].value, paras['ic0'].value,paras['v0'].value, paras['nk0'].value
#    model = g(t, x0, paras)
#
#    # you only have data for one of your variables
#    x2_model = model[:, 1]
#    return (x2_model - data).ravel()
#
## initial condition
#bn0=1
#bgc0=0
#bm0=0
#bpc0=0
#uc0=1000
#ic0=0
#v0=1
#nk0=100
#y0 = [bn0,bgc0,bm0,bpc0,uc0,ic0,v0,nk0]
#
#
#
## measured data
#t_measured = np.array([1, 4, 7, 10, 13, 16, 19, 22, 25, 28, 31, 34, 37])
##x2_measured = np.array([0.000, 0.416, 0.489, 0.595, 0.506, 0.493, 0.458, 0.394, 0.335, 0.309])
#
#igmm = np.array([1, 3.69, 8.98, 12.88, 16.27, 19.37, 22.96, 26.46, 29.55, 32.15, 34.45, 36.05, 37.46])
#iggm = np.array([0.40,1.80,5.00,8.59,12.39, 16.78, 22.56, 27.75, 33.14, 37.73, 41.63, 45.12, 48.72])
#iggpc = np.array([0.20, 0.71,1.81, 3.31, 4.82, 6.12, 8.12, 10.82, 14.11, 18.21, 22.4, 26.89, 31.48])
#
##data=[igmm,iggm,iggpc]
#
#
#plt.figure()
#plt.scatter(t_measured, igmm, marker='o', color='b', label='measured data', s=75)
#
## set parameters including bounds; you can also fix parameters (use vary=False)
#params = Parameters()
#params.add('bn0', value=bn0, vary=False)
#params.add('bgc0', value=bgc0, vary=False)
#params.add('bm0', value=bm0, vary=False)
#params.add('bpc0', value=bpc0, vary=False)
#params.add('uc0', value=uc0, vary=False)
#params.add('ic0', value=ic0, vary=False)
#params.add('v0', value=v0, vary=False)
#params.add('nk0', value=nk0, vary=False)
#params.add('ky', value=ky, min=0.001, max=1000)
##params.add('ky', value=10, vary=False)
#params.add('km', value=km, min=0.001, max=1000)
#params.add('kmp', value=0.2, min=0.001, max=1000)
#params.add('p', value=0.3, min=0.001, max=1000)
#
## fit model
#result = minimize(residual, params, args=(t_measured, igmm), method='leastsq')  # leastsq nelder
## check results of the fit
#data_fitted = g(np.linspace(0., 50., 100), y0, result.params)
#
## plot fitted data
#plt.plot(np.linspace(0., 50., 100), data_fitted[:, 1], '-', linewidth=2, color='red', label='fitted data')
#plt.legend()
#plt.xlim([0, max(t_measured)])
#plt.ylim([0, 50])
##plt.ylim([0, 1.1 * max(data_fitted[:, 1])])
## display fitted statistics
#report_fit(result)
#
#plt.show()





###initialize the data
#x_data = np.linspace(0,9,10)
#y_data = np.array([0.000,0.416,0.489,0.595,0.506,0.493,0.458,0.394,0.335,0.309])
#
#
#def f(y, t, k): 
#    """define the ODE system in terms of 
#        dependent variable y,
#        independent variable t, and
#        optinal parmaeters, in this case a single variable k """
#    return (-k[0]*y[0],
#          k[0]*y[0]-k[1]*y[1],
#          k[1]*y[1])
#
#def my_ls_func(x,teta):
#    """definition of function for LS fit
#        x gives evaluation points,
#        teta is an array of parameters to be varied for fit"""
#    # create an alias to f which passes the optional params    
#    f2 = lambda y,t: f(y, t, teta)
#    # calculate ode solution, retuen values for each entry of "x"
#    r = integrate.odeint(f2,y0,x)
#    #in this case, we only need one of the dependent variable values
#    return r[:,1]
#
#def f_resid(p):
#    """ function to pass to optimize.leastsq
#        The routine will square and sum the values returned by 
#        this function""" 
#    return y_data-my_ls_func(x_data,p)
##solve the system - the solution is in variable c
#guess = [0.2,0.3] #initial guess for params
#y0 = [1,0,0] #inital conditions for ODEs
#(c,kvg) = optimize.leastsq(f_resid, guess) #get params
#
#print("parameter values are ",c)
#
## fit ODE results to interpolating spline just for fun
#xeval=np.linspace(min(x_data), max(x_data),30) 
#gls = interpolate.UnivariateSpline(xeval, my_ls_func(xeval,c), k=3, s=0)
#
##pick a few more points for a very smooth curve, then plot 
##   data and curve fit
#xeval=np.linspace(min(x_data), max(x_data),200)
##Plot of the data as red dots and fit as blue line
#pp.plot(x_data, y_data,'.r',xeval,gls(xeval),'-b')
#pp.xlabel('xlabel',{"fontsize":16})
#pp.ylabel("ylabel",{"fontsize":16})
#pp.legend(('data','fit'),loc=0)
#pp.show()


































##initialize the data
#
#2
#import pylab as pp
#import numpy as np
#from scipy import integrate, interpolate
#from scipy import optimize
#from scipy.optimize import least_squares
#
##x_data = np.linspace(0,9,10)
##y_data = np.array([0.000,0.416,0.489,0.595,0.506,0.493,0.458,0.394,0.335,0.309])
## measured data
#x_data = np.array([1, 4, 7, 10, 13, 16, 19, 22, 25, 28, 31, 34, 37])
##x2_measured = np.array([0.000, 0.416, 0.489, 0.595, 0.506, 0.493, 0.458, 0.394, 0.335, 0.309])
#
#y_data = np.array([1, 3.69, 8.98, 12.88, 16.27, 19.37, 22.96, 26.46, 29.55, 32.15, 34.45, 36.05, 37.46])
#y_data = np.array([0.40,1.80,5.00,8.59,12.39, 16.78, 22.56, 27.75, 33.14, 37.73, 41.63, 45.12, 48.72])
#
#def f(y, t, kk): 
#    bn, bgc, bm, bpc, uc, ic, v, nk = y
#    km,kmp,ky,p=kk
#    if t>72:
#        k=0
#    else:
#        k=1
##    s=100*k
##    p=0.05
#    kv=0.001
#    knk=0.1
#    kdv=0.01
#    kin=0.003
#    kdi=0.1
#    eq = [kk[2]*k-kk[3]*bn, 2*kk[3]*bn-kk[0]*bgc-kk[1]*bgc, kk[0]*bgc, kk[1]*bgc, -kin*uc*v, kin*uc*v -kdi*ic*nk, kv*ic-kdv*v, knk-kdi*nk*ic]
#    return eq
#
#
#def my_ls_func(x,teta):
#    """definition of function for LS fit
#        x gives evaluation points,
#        teta is an array of parameters to be varied for fit"""
#    # create an alias to f which passes the optional params    
#    f2 = lambda y,t: f(y, t, teta)
#    # calculate ode solution, retuen values for each entry of "x"
#    r = integrate.odeint(f2,y0,x)
#    #in this case, we only need one of the dependent variable values
#    return r[:,2]
#
#def f_resid(p):
#    """ function to pass to optimize.leastsq
#        The routine will square and sum the values returned by 
#        this function""" 
#    return y_data-my_ls_func(x_data,p)
##solve the system - the solution is in variable c
#guess = [0.2,0.3,0.1,0.1] #initial guess for params
#y0 = [1,0,0,0,1000,0,1,100]
#(c,kvg) = optimize.leastsq(f_resid, guess) #get params
#
#
#print("parameter values are ",c)
#
## fit ODE results to interpolating spline just for fun
#xeval=np.linspace(min(x_data), max(x_data),30) 
#gls = interpolate.UnivariateSpline(xeval, my_ls_func(xeval,c), k=3, s=0)
#
##pick a few more points for a very smooth curve, then plot 
##   data and curve fit
#xeval=np.linspace(min(x_data), max(x_data),200)
##Plot of the data as red dots and fit as blue line
#pp.plot(x_data, y_data,'.r',xeval,gls(xeval),'-b')
#pp.xlabel('xlabel',{"fontsize":16})
#pp.ylabel("ylabel",{"fontsize":16})
#pp.legend(('data','fit'),loc=0)
#pp.show()
#
#
#opt= least_squares(f_resid,(0.2,0.3,0.1,0.1),bounds=([0,1000]))
#gls = interpolate.UnivariateSpline(xeval, my_ls_func(xeval,opt.x), k=3, s=0)
#print("parameter values are ",opt.x)
##pick a few more points for a very smooth curve, then plot 
##   data and curve fit
#xeval=np.linspace(min(x_data), max(x_data),200)
##Plot of the data as red dots and fit as blue line
#pp.plot(x_data, y_data,'.r',xeval,gls(xeval),'-b')
#pp.xlabel('xlabel',{"fontsize":16})
#pp.ylabel("ylabel",{"fontsize":16})
#pp.legend(('data','fit'),loc=0)
#pp.show()










#import numpy as np
#import pandas as pd
#import matplotlib.pyplot as plt
##matplotlib inline
#from scipy.integrate import odeint
#from scipy.optimize import minimize
#from scipy.interpolate import interp1d
#
## initial guesses
#x0 = np.zeros(4)
#x0[0] = 0.8 # Kp
#x0[1] = 0.2 # Kd
#x0[2] = 150.0 # taup
#x0[3] = 10.0 # thetap
#
## Import CSV data file
## try to read local data file first
#try:
#    filename = 'data.csv'
#    data = pd.read_csv(filename)
#except:
#    filename = 'http://apmonitor.com/pdc/uploads/Main/tclab_data2.txt'
#    data = pd.read_csv(filename)
#Q1_0 = data['Q1'].values[0]
#Q2_0 = data['Q2'].values[0]
#T1_0 = data['T1'].values[0]
#T2_0 = data['T2'].values[0]
#t = data['Time'].values - data['Time'].values[0]
#Q1 = data['Q1'].values
#Q2 = data['Q2'].values
#T1 = data['T1'].values
#T2 = data['T2'].values
#
## specify number of steps
#ns = len(t)
#delta_t = t[1]-t[0]
## create linear interpolation of the u data versus time
#Qf1 = interp1d(t,Q1)
#Qf2 = interp1d(t,Q2)
#
## define first-order plus dead-time approximation    
#def fopdt(T,t,Qf1,Qf2,Kp,Kd,taup,thetap):
#    #  T      = states
#    #  t      = time
#    #  Qf1    = input linear function (for time shift)
#    #  Qf2    = input linear function (for time shift)
#    #  Kp     = model gain
#    #  Kd     = disturbance gain
#    #  taup   = model time constant
#    #  thetap = model time constant
#    # time-shift Q
#    try:
#        if (t-thetap) <= 0:
#            Qm1 = Qf1(0.0)
#            Qm2 = Qf2(0.0)
#        else:
#            Qm1 = Qf1(t-thetap)
#            Qm2 = Qf2(t-thetap)
#    except:
#        Qm1 = Q1_0
#        Qm2 = Q2_0
#    # calculate derivative
#    dT1dt = (-(T[0]-T1_0) + Kp*(Qm1-Q1_0) + Kd*(T[1]-T[0]))/taup
#    dT2dt = (-(T[1]-T2_0) + (Kp/2.0)*(Qm2-Q2_0) + Kd*(T[0]-T[1]))/taup
#    return [dT1dt,dT2dt]
#
## simulate FOPDT model
#def sim_model(x):
#    # input arguments
#    Kp,Kd,taup,thetap = x
#    # storage for model values
#    T1p = np.ones(ns) * T1_0
#    T2p = np.ones(ns) * T2_0
#    # loop through time steps    
#    for i in range(0,ns-1):
#        ts = [t[i],t[i+1]]
#        T = odeint(fopdt,[T1p[i],T2p[i]],ts,args=(Qf1,Qf2,Kp,Kd,taup,thetap))
#        T1p[i+1] = T[-1,0]
#        T2p[i+1] = T[-1,1]
#    return T1p,T2p
#
## define objective
#def objective(x):
#    # simulate model
#    T1p,T2p = sim_model(x)
#    # return objective
#    return sum(np.abs(T1p-T1)+np.abs(T2p-T2))
#
## show initial objective
#print('Initial SSE Objective: ' + str(objective(x0)))
#print('Optimizing Values...')
#
## optimize without parameter constraints
##solution = minimize(objective,x0)
#
## optimize with bounds on variables
#bnds = ((0.4, 1.5), (0.1, 0.5), (50.0, 200.0), (0.0, 30.0))
#solution = minimize(objective,x0,bounds=bnds,method='SLSQP')
#
## show final objective
#x = solution.x
#iae = objective(x)
#Kp,Kd,taup,thetap = x
#print('Final SSE Objective: ' + str(objective(x)))
#print('Kp: ' + str(Kp))
#print('Kd: ' + str(Kd))
#print('taup: ' + str(taup))
#print('thetap: ' + str(thetap))
## save fopdt.txt file
#fid = open('fopdt.txt','w')
#fid.write(str(Kp)+'\n')
#fid.write(str(Kd)+'\n')
#fid.write(str(taup)+'\n')
#fid.write(str(thetap)+'\n')
#fid.write(str(T1_0)+'\n')
#fid.write(str(T2_0)+'\n')
#fid.close()
#
## calculate model with updated parameters
#T1p,T2p = sim_model(x)
#
#plt.figure(1,figsize=(15,7))
#plt.subplot(2,1,1)
#plt.plot(t,T1,'r.',linewidth=2,label='Temperature 1 (meas)')
#plt.plot(t,T2,'b.',linewidth=2,label='Temperature 2 (meas)')
#plt.plot(t,T1p,'r--',linewidth=2,label='Temperature 1 (pred)')
#plt.plot(t,T2p,'b--',linewidth=2,label='Temperature 2 (pred)')
#plt.ylabel(r'T $(^oC)$')
#plt.text(200,20,'Integral Abs Error: ' + str(np.round(iae,2)))
#plt.text(400,35,r'$K_p$: ' + str(np.round(Kp,2)))  
#plt.text(400,30,r'$K_d$: ' + str(np.round(Kd,2)))  
#plt.text(400,25,r'$\tau_p$: ' + str(np.round(taup,1)) + ' sec')  
#plt.text(400,20,r'$\theta_p$: ' + str(np.round(thetap,1)) + ' sec')  
#plt.legend(loc=2)
#plt.subplot(2,1,2)
#plt.plot(t,Q1,'b--',linewidth=2,label=r'Heater 1 ($Q_1$)')
#plt.plot(t,Q2,'r:',linewidth=2,label=r'Heater 2 ($Q_2$)')
#plt.legend(loc='best')
#plt.xlabel('time (sec)')
#plt.show()

## example 2

#import pylab as pp
#import numpy as np
#from scipy import integrate, interpolate
#from scipy import optimize
#
###initialize the data
#x_data = np.linspace(0,9,10)
#y_data = np.array([0.000,0.416,0.489,0.595,0.506,0.493,0.458,0.394,0.335,0.309])
#
#
#def f(y, t, k): 
#    """define the ODE system in terms of 
#        dependent variable y,
#        independent variable t, and
#        optinal parmaeters, in this case a single variable k """
#    return (-k[0]*y[0],
#          k[0]*y[0]-k[1]*y[1],
#          k[1]*y[1])
#
#def my_ls_func(x,teta):
#    """definition of function for LS fit
#        x gives evaluation points,
#        teta is an array of parameters to be varied for fit"""
#    # create an alias to f which passes the optional params    
#    f2 = lambda y,t: f(y, t, teta)
#    # calculate ode solution, retuen values for each entry of "x"
#    r = integrate.odeint(f2,y0,x)
#    #in this case, we only need one of the dependent variable values
#    return r[:,1]
#
#def f_resid(p):
#    """ function to pass to optimize.leastsq
#        The routine will square and sum the values returned by 
#        this function""" 
#    return y_data-my_ls_func(x_data,p)
##solve the system - the solution is in variable c
#guess = [0.2,0.3] #initial guess for params
#y0 = [1,0,0] #inital conditions for ODEs
#(c,kvg) = optimize.leastsq(f_resid, guess) #get params
#
#print("parameter values are ",c)
#
## fit ODE results to interpolating spline just for fun
#xeval=np.linspace(min(x_data), max(x_data),30) 
#gls = interpolate.UnivariateSpline(xeval, my_ls_func(xeval,c), k=3, s=0)
#
##pick a few more points for a very smooth curve, then plot 
##   data and curve fit
#xeval=np.linspace(min(x_data), max(x_data),200)
##Plot of the data as red dots and fit as blue line
#pp.plot(x_data, y_data,'.r',xeval,gls(xeval),'-b')
#pp.xlabel('xlabel',{"fontsize":16})
#pp.ylabel("ylabel",{"fontsize":16})
#pp.legend(('data','fit'),loc=0)
#pp.show()



## example 3

# cleaned up a bit to get my head around it - thanks for sharing 
#import pylab as pp
#import numpy as np
#from scipy import integrate, optimize
#
#class Parameterize_ODE():
#    def __init__(self):
#        self.X = np.linspace(0,9,10)
#        self.y = np.array([0.000,0.416,0.489,0.595,0.506,0.493,0.458,0.394,0.335,0.309])
#        self.y0 = [1,0,0] # inital conditions ODEs
#    def ode(self, y, X, p):
#        return (-p[0]*y[0],
#                 p[0]*y[0]-p[1]*y[1],
#                           p[1]*y[1])
#    def model(self, X, p):
#        return integrate.odeint(self.ode, self.y0, X, args=(p,))
#    def f_resid(self, p):
#        return self.y - self.model(self.X, p)[:,1]
#    def optim(self, p_quess):
#        return optimize.leastsq(self.f_resid, p_guess) # fit params
#
#po = Parameterize_ODE(); p_guess = [0.2, 0.3] 
#c, kvg = po.optim(p_guess)
#
## --- show ---
#print("parameter values are ", c, kvg)
#x = np.linspace(min(po.X), max(po.X), 2000)
#pp.plot(po.X, po.y,'.r',x, po.model(x, c)[:,1],'-b')
#pp.xlabel('X',{"fontsize":16}); pp.ylabel("y",{"fontsize":16}); pp.legend(('data','fit'),loc=0); pp.show()






## example 4
#from lmfit import minimize, Parameters, Parameter, report_fit
#from scipy.integrate import odeint
#
#def f(xs, t, ps):
#    """Lotka-Volterra predator-prey model."""
#    try:
#        a = ps['a'].value
#        b = ps['b'].value
#        c = ps['c'].value
#        d = ps['d'].value
#    except:
#        a, b, c, d = ps
#
#    x, y = xs
#    return [a*x - b*x*y, c*x*y - d*y]
#
#def g(t, x0, ps):
#    """
#    Solution to the ODE x'(t) = f(t,x,k) with initial condition x(0) = x0
#    """
#    x = odeint(f, x0, t, args=(ps,))
#    return x
#
#def residual(ps, ts, data):
#    x0 = ps['x0'].value, ps['y0'].value
#    model = g(ts, x0, ps)
#    return (model - data).ravel()
#
#t = np.linspace(0, 10, 100)
#x0 = np.array([1,1])
#
#a, b, c, d = 3,1,1,1
#true_params = np.array((a, b, c, d))
#data = g(t, x0, true_params)
#data += np.random.normal(size=data.shape)
#
## set parameters incluing bounds
#params = Parameters()
#params.add('x0', value= float(data[0, 0]), min=0, max=10)
#params.add('y0', value=float(data[0, 1]), min=0, max=10)
#params.add('a', value=2.0, min=0, max=10)
#params.add('b', value=1.0, min=0, max=10)
#params.add('c', value=1.0, min=0, max=10)
#params.add('d', value=1.0, min=0, max=10)
#
## fit model and find predicted values
#result = minimize(residual, params, args=(t, data), method='leastsq')
#final = data + result.residual.reshape(data.shape)
#
## plot data and fitted curves
#plt.plot(t, data, 'o')
#plt.plot(t, final, '-', linewidth=2);
#
## display fitted statistics
#report_fit(result)













































