
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import trapz,cumtrapz
import re
import matplotlib as mpl
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.optimize import curve_fit, minimize
from scipy import interpolate
from SharedFunctions import return2d,ImportGridue,plotWALL
import numpy.ma as ma

#define plot parameters
plt.rcParams["font.family"] = "serif"
params = {'legend.fontsize': 'medium',
         'axes.labelsize': 'medium',
         'axes.titlesize':'medium',
         'xtick.labelsize':'medium',
         'ytick.labelsize':'medium',
         }
plt.rcParams.update(params)

# define exponential decay of heat flux
def EICH(x,LAMBDA,peakq):
    if x>=0:
        q = peakq*np.exp(-x/LAMBDA)
    else:
        q = 0
    return q

def GAUSS(x,S):
    y = (1/(S*np.sqrt(np.pi)))*np.exp(-(x/S)**2)
    
    return y
xtest = np.linspace(-100,100,10000)
print(np.trapz(GAUSS(xtest,0.1),xtest))


folder0 = "patchData"
dataue = ImportGridue(folder0+"\\gridue")

def expFunc(x,A,B):
    return A*np.exp((-x)/B)


# Single Null and XPT balance files
# file0 = "balFiles//balance0.nc"
# file1 = "balFiles//balance.nc"
file0 = "balFiles//balanceSN.nc"
file1 = "balFiles//balanceXPT.nc"
#unpack two balance files
quantities2d0,rootgrp0 = return2d(file0)
quantities2d1,rootgrp1 = return2d(file1)

#define midplane x point indices
midplaneix0 = np.argmax(quantities2d0["r"][-1])
midplaneix1 = np.argmax(quantities2d1["r"][-1])
Xpointix0=91
Xpointix1=dataue["ix_cut2"]+1

JSEP = rootgrp1["jsep"][0]+2

jsep0 = rootgrp0["jsep"][0]+2
#define R-Rsep for each radial surface
Rrsep0 = quantities2d0["r"]-quantities2d0["r"][jsep0,midplaneix0]
Rrsep1 = quantities2d1["r"]-quantities2d1["r"][JSEP,midplaneix1]

#calculate connection lengths for each radial surface
LArray0 = []
for i in range(jsep0,len(quantities2d0["qpar"])-1):
    Lpar = trapz(np.abs(quantities2d0["sdiff"][i,midplaneix0:-1]))
    LArray0.append(Lpar)

#calculate connection lengths for each radial surface on XPT
LArray1 = []
qdsArray1 = []
Tt1 = []
Nt = []
for i in range(dataue["iyseparatrix2"]+1,dataue["iyseparatrix3"]+1):
    Lpar = trapz(np.abs(quantities2d1["sdiff"][i,midplaneix1:dataue['ix_cut3']+1]))
    Lpar = Lpar+trapz(np.abs(quantities2d1["sdiff"][i,dataue['ix_cut4']+1:]))
    qds = trapz(np.abs(quantities2d1["sdiff"][i,midplaneix1:dataue['ix_cut3']+1]*quantities2d1["cond"][i,midplaneix1:dataue['ix_cut3']+1]))
    qds =  qds+trapz(np.abs(quantities2d1["sdiff"][i,dataue['ix_cut4']+1:]*quantities2d1["cond"][i,dataue['ix_cut4']+1:]))
    qdsArray1.append(qds)
    LArray1.append(Lpar)
    Tt1.append(quantities2d1["te"][i][-1])
    Nt.append(quantities2d1["ne"][i][-2])
    plt.plot(quantities2d1["r"][i,midplaneix1:dataue['ix_cut3']+1],quantities2d1["z"][i,midplaneix1:dataue['ix_cut3']+1],color="C0",linewidth=0.2)
    plt.plot(quantities2d1["r"][i,dataue['ix_cut4']+1:],quantities2d1["z"][i,dataue['ix_cut4']+1:],color="C0",linewidth=0.2)
for i in range(dataue["iyseparatrix3"]+1,len(quantities2d1["qpar"])-1):
    Lpar = trapz(np.abs(quantities2d1["sdiff"][i,midplaneix1:dataue["ix_plate2"]+1]))
    LArray1.append(Lpar)
    qds = trapz(np.abs(quantities2d1["sdiff"][i,midplaneix1:dataue["ix_plate2"]+1]*quantities2d1["cond"][i,midplaneix1:dataue["ix_plate2"]+1]))
    qdsArray1.append(qds)
    Tt1.append(quantities2d1["te"][i,dataue["ix_plate2"]+1])
    Nt.append(quantities2d1["ne"][i][dataue["ix_plate2"]])
    plt.plot(quantities2d1["r"][i,midplaneix1:dataue["ix_plate2"]+1],quantities2d1["z"][i,midplaneix1:dataue["ix_plate2"]+1],color="C2",linewidth=0.2)

plt.savefig("1.png",dpi=1000,bbox_inches='tight')
plt.show()


x0 = Rrsep0[jsep0:-1,midplaneix0]
x0 = x0*1000

x1 = Rrsep1[JSEP:-1,midplaneix1]
x1 = x1*1000


y0 = (quantities2d0["qpar"])[jsep0:-1,Xpointix0]
popt, pcov = curve_fit(expFunc, x0, y0,p0=[np.amax(y0),4])
y1 = (quantities2d1["qpar"])[JSEP:-1,Xpointix1]
popt1, pcov1 = curve_fit(expFunc, x1, y1,p0=[np.amax(y0),1])

# popt1, pcov = curve_fit(expFunc, x, y1,p0=[np.amax(y1),0.002])
# # plt.plot(x,y0,marker="o",linestyle="")
print(trapz(y0,x0))
print(trapz(y1,x1))
fig, axs = plt.subplots(1, 1)
# twin1 = axs.twinx()
# twin1.plot(x0,LArray0,marker='o',color="#5F4690")
axs.plot(x0,y0,marker="o",color="#EDAD08",label="SN")
# axs.plot(x1,y1,marker="o",color="#EDAD08",alpha=0)
# axs.plot(x0,expFunc(x0, *popt),label="fit with "+r"$\lambda_{q}$"+"="+str(np.round(popt[1],1))+"mm",color="black")
# axs.tick_params(axis='y', labelcolor="#EDAD08")
axs.set_ylabel("q"+r'$_{||,u}$'+" (MWm"+r"$^{-2}$"+")")
# twin1.tick_params(axis='y', labelcolor="#5F4690")
# twin1.set_ylabel("L (m)")
axs.set_xlabel(r"$r_{u} (mm)$")
axs.legend()
axs.set_xlim([-0.07,2.1])
# twin1.set_ylim([0,140])
# axs.set_ylim([-0.05,90])
# plt.savefig("images//SNupstreamqpar.png",dpi=1000,bbox_inches='tight')
# plt.show()

# fig, axs = plt.subplots(1, 1)
# twin1 = axs.twinx()
# twin1.plot(x1,LArray1,marker='o',color="#5F4690")
axs.plot(x1,y1,marker="o",color="#5F4690",label="XPT")
# axs.plot(x0,y0,marker="o",color="#EDAD08",alpha=0)
# axs.plot(x1,expFunc(x1, *popt1),label="fit with "+r"$\lambda_{q}$"+"="+str(np.round(popt1[1],1))+"mm",color="black")
# axs.tick_params(axis='y', labelcolor="#EDAD08")
axs.set_ylabel("q"+r'$_{||,u}$'+" (MWm"+r"$^{-2}$"+")")
# twin1.tick_params(axis='y', labelcolor="#5F4690")
# twin1.set_ylabel("L (m)")
axs.set_xlabel("r"+r"$_{u}$"+" (mm)")
axs.legend()
axs.set_xlim([-0.07,2.1])
# twin1.set_ylim([0,140])
# axs.set_ylim([-0.05,1.15])
plt.savefig("images//XPTupstreamqpar.png",dpi=1000,bbox_inches='tight')
plt.show()


x0 = Rrsep0[1:-1,midplaneix0]
x0 = x0*1000
print("extent is")
print(x0[0])
print(x0[-1])
x1 = Rrsep1[1:-1,midplaneix1]
x1 = x1*1000
print(x1[0])
print(x1[-1])
y0 = (quantities2d0["ne"])[1:-1,midplaneix0]
y1 = (quantities2d1["ne"])[1:-1,midplaneix1]


plt.plot(x0,y0,marker="o",label="SN",color="#38A6A5")
plt.plot(x1,y1,marker="^",label="XPT",color="#CC503E")
plt.ylabel("n"+r'$_{e,u}$'+" (m"+r"$^{-3}$"+")")
plt.xlabel("r"+r"$_{u}$"+" (mm)")
plt.legend()
plt.xlim([-0.5,2.1])
plt.savefig("images//upstreamdens.png",dpi=1000,bbox_inches='tight')
plt.show()

y0 = (quantities2d0["te"])[1:-1,midplaneix0]
y1 = (quantities2d1["te"])[1:-1,midplaneix1]


plt.plot(x0,y0,marker="o",label="SN",color="#38A6A5")
plt.plot(x1,y1,marker="^",label="XPT",color="#CC503E")
# plt.plot(x,expFunc(x, *popt1),label="fit with "+r"$\lambda_{n}$"+"="+str(np.round(popt1[1],0))+"mm")
plt.ylabel("T"+r'$_{e,u}$'+" (eV)")
plt.xlabel("r"+r"$_{u}$"+" (mm)")
plt.ylim([30,130])
plt.xlim([-0.5,2.1])
plt.legend()
plt.savefig("images//upstreamTemp.png",dpi=1000,bbox_inches='tight')
plt.show()

echarge = 1.60E-19
x0 = Rrsep0[jsep0:-1,midplaneix0]
x0 = x0*1000
y0 = (echarge*quantities2d0["te"]*quantities2d0["ne"])[jsep0:-1,-2]
x1 = Rrsep1[JSEP:-1,midplaneix1]
x1 = x1*1000
y1 = (echarge*quantities2d1["te"]*quantities2d1["ne"])[JSEP:-1,-2]
y12 = (echarge*quantities2d1["te"]*quantities2d1["ne"])[JSEP:-1,dataue["ix_plate2"]]
plt.plot(x0,y0,marker="o",label="SN",color="#38A6A5")
plt.plot(x1,y1,marker="o",label="XPT Target 1",color="#CC503E")
plt.plot(x1,y12,marker="^",label="XPT Target 2",color="#CC503E")
plt.plot([Rrsep1[dataue["iyseparatrix2"]+1,midplaneix1]*1000,Rrsep1[dataue["iyseparatrix2"]+1,midplaneix1]*1000],[0.5*np.min(y1),2*np.max(y1)],label="Primary Separatrix",linestyle="--",color="#3E434A")
plt.plot([Rrsep1[dataue["iyseparatrix4"]+1,midplaneix1]*1000,Rrsep1[dataue["iyseparatrix4"]+1,midplaneix1]*1000],[0.5*np.min(y1),2*np.max(y1)],label="Secondary Separatrix",linestyle="-.",color="#3E434A")
plt.ylim([np.min(y1)*0.8,np.max(y1)*1.1])
plt.ylabel("P"+r'$_{e,t}$'+" (Pa)"  )
plt.xlabel("r"+r"$_{u}$"+" (mm)")
plt.legend()
plt.savefig("images//targetPressure.png",dpi=1000,bbox_inches='tight')
plt.show()

tdiff1actual0 = quantities2d0["te"][JSEP:-1,Xpointix0]-quantities2d0["te"][JSEP:-1,-2]

tdiff1actual1 = quantities2d1["te"][JSEP:-1,Xpointix1]-quantities2d1["te"][JSEP:-1,dataue["ix_plate2"]]
x0 = Rrsep0[jsep0:-1,midplaneix0]*1000
y0 = (quantities2d0["te"])[jsep0:-1,-2]
x1 = Rrsep1[JSEP:-1,midplaneix1]*1000
y1 = (quantities2d1["te"])[JSEP:-1,-2]
y12 = (quantities2d1["te"])[JSEP:-1,dataue["ix_plate2"]]
plt.plot(x0,y0,marker="o",label="SN",color="#38A6A5")

plt.plot([Rrsep1[dataue["iyseparatrix2"]+1,midplaneix1]*1000,Rrsep1[dataue["iyseparatrix2"]+1,midplaneix1]*1000],[0.5*np.min(y1),2*np.max(y1)],label="Primary Separatrix",linestyle="--",color="#3E434A")
plt.plot([Rrsep1[dataue["iyseparatrix4"]+1,midplaneix1]*1000,Rrsep1[dataue["iyseparatrix4"]+1,midplaneix1]*1000],[0.5*np.min(y1),2*np.max(y1)],label="Secondary Separatrix",linestyle="-.",color="#3E434A")
plt.ylim([0,np.max(y0)*1.1])
plt.plot(x1,y1,marker="o",label="XPT Target 1",color="#CC503E")
plt.plot(x1,y12,marker="^",label="XPT Target 2",color="#CC503E")
plt.ylabel("T"+r'$_{e,t}$'+" (eV)"  )
plt.xlabel("r"+r"$_{u}$"+" (mm)")
plt.legend()
plt.savefig("images//targetTemp.png",dpi=1000,bbox_inches='tight')
plt.show()

def returnEICH(LQ,peakq,S):
    x = np.linspace(-200,200,20000)
    y0fun = []
    y1fun = []
    for i in range(len(x)):
        y0fun.append(EICH(x[i],LAMBDA=LQ,peakq=peakq))
        y1fun.append(GAUSS(x[i],S))
    # plt.plot(x,y0)
    # plt.plot(x,y1)
    x = np.array(x)
    yfinal = np.array(np.convolve(y0fun,y1fun,mode='same'))
    interpfunc = interpolate.interp1d(x,yfinal)
    return interpfunc

def FIT(params,xfit,ydata):
    [LQ,peakq,S] = params


    xfit = np.array(xfit)

    interpfunc = returnEICH(LQ,peakq,S)

    Yout = interpfunc(xfit)
    return np.mean((Yout-np.array(ydata))**2)#/(np.mean(Yout+ydata))

# plot target heat flux density and diffusion model fits.
x0 = ma.getdata(Rrsep0[jsep0:-1,midplaneix0]*1000)
y0 = (quantities2d0["qpar"])[jsep0:-1,-2]
xplot = np.linspace(0,2,1000)
x1 = ma.getdata(Rrsep1[JSEP:-1,midplaneix1]*1000)
y1 = (quantities2d1["qpar"])[JSEP:-1,-2]
y12 = (quantities2d1["qpar"])[JSEP:-1,dataue["ix_plate2"]]
bounds = ((0.01,0.8),(0.1,150),(0.01,1))
dataout = minimize(FIT,x0=[0.5,25,0.1],args=(x0,y0),method='TNC',bounds=bounds,options={'maxiter':5000,'gtol':1E-9})
paramsknown = dataout['x']
interpfunc0 = returnEICH(paramsknown[0],paramsknown[1],paramsknown[2])
plt.plot(xplot,interpfunc0(xplot),color="#EDAD08")
plt.plot(x0,y0,marker="o",label="SN",color="#EDAD08",linestyle="")
bounds = ((0.01,0.8),(1,150),(0.01,1))
dataout = minimize(FIT,x0=[0.5,25,0.1],args=(x1,y1),method='TNC',bounds=bounds,options={'maxiter':5000,'gtol':1E-9})
paramsknown = dataout['x']
print(paramsknown)
interpfunc0 = returnEICH(paramsknown[0],paramsknown[1],paramsknown[2])
plt.plot(xplot,interpfunc0(xplot),color="#5F4690")
plt.plot(x1,y1,marker="o",label="XPT Target 1",color="#5F4690",linestyle='')
bounds = ((0.01,0.8),(0.1,10),(0.01,1))
x1fit = x1-Rrsep1[dataue["iyseparatrix4"]+1,midplaneix1]*1000
dataout = minimize(FIT,x0=[0.5,9,0.1],args=(x1fit,y12),method='TNC',bounds=bounds,options={'maxiter':1000})
paramsknown = dataout['x']
print(dataout)
interpfunc0 = returnEICH(paramsknown[0],paramsknown[1],paramsknown[2])
plt.plot(xplot,interpfunc0(xplot-Rrsep1[dataue["iyseparatrix4"]+1,midplaneix1]*1000),color="#5F4690")
plt.plot(x1,y12,marker="^",label="XPT Target 2",color="#5F4690",linestyle='')
plt.plot([0],[0],color="black",label="diffusion model")
plt.plot([Rrsep1[dataue["iyseparatrix2"]+1,midplaneix1]*1000,Rrsep1[dataue["iyseparatrix2"]+1,midplaneix1]*1000],[0.5*np.min(y1),2*np.max(y1)],label="Primary Separatrix",linestyle="--",color="#3E434A")
plt.plot([Rrsep1[dataue["iyseparatrix4"]+1,midplaneix1]*1000,Rrsep1[dataue["iyseparatrix4"]+1,midplaneix1]*1000],[0.5*np.min(y1),2*np.max(y1)],label="Secondary Separatrix",linestyle="-.",color="#3E434A")
plt.ylim([0,np.max(y0)*1.1])
plt.ylabel("q"+r'$_{||,t}$'+" (MWm"+r"$^{-2}$"+")")
plt.xlabel("r"+r"$_{u}$"+" (mm)")
plt.legend()
plt.savefig("images//targetQpar.png",dpi=1000,bbox_inches='tight')
plt.show()


# plot target density as a function of radial separation mapped to midplane.
x0 = Rrsep0[jsep0:-1,midplaneix0]*1000
y0 = (quantities2d0["ne"])[jsep0:-1,-2]
x1 = Rrsep1[JSEP:-1,midplaneix1]*1000
y1 = (quantities2d1["ne"])[JSEP:-1,-2]
y12 = (quantities2d1["ne"])[JSEP:-1,dataue["ix_plate2"]]
plt.plot(x0,y0,marker="o",label="SN",color="#38A6A5")
plt.plot(x1,y1,marker="o",label="XPT Target 1",color="#5F4690")
plt.plot(x1,y12,marker="^",label="XPT Target 2",color="#5F4690")
plt.plot([Rrsep1[dataue["iyseparatrix2"]+1,midplaneix1]*1000,Rrsep1[dataue["iyseparatrix2"]+1,midplaneix1]*1000],[0.5*np.min(y1),2*np.max(y1)],label="Primary Separatrix",linestyle="--",color="#3E434A")
plt.plot([Rrsep1[dataue["iyseparatrix4"]+1,midplaneix1]*1000,Rrsep1[dataue["iyseparatrix4"]+1,midplaneix1]*1000],[0.5*np.min(y1),2*np.max(y1)],label="Secondary Separatrix",linestyle="-.",color="#3E434A")
plt.ylim([0,np.max(y1)*1.1])
plt.ylabel("n"+r'$_{e,t}$'+" (m"+r"$^{-3}$"+")")
plt.xlabel("r"+r"$_{u}$"+" (mm)")
plt.legend()
plt.savefig("images//targetNe.png",dpi=1000,bbox_inches='tight')
plt.show()


# perform analysis of target temperature
tdiff0 = (LArray0*(quantities2d0["qpar"][jsep0:-1,Xpointix0]))
tdiff1 = (LArray1*(quantities2d1["qpar"][JSEP:-1,Xpointix1]))
q0 = (quantities2d0["qpar"][jsep0:-1,Xpointix0])
q1 = (quantities2d1["qpar"][JSEP:-1,Xpointix1])
qint1 = tdiff1
TT1 = quantities2d0["te"][JSEP:-1,-2]
TU1 = quantities2d1["te"][JSEP:-1,midplaneix1]#**(7/2)
tdiff1actual0 = quantities2d0["te"][jsep0:-1,Xpointix0]**(7/2)-quantities2d0["te"][jsep0:-1,-2]**(7/2)
tdiff1actual1 = quantities2d1["te"][JSEP:-1,midplaneix1]**(7/2)-np.array(Tt1)**(7/2)
x0 = Rrsep0[jsep0:-1,midplaneix0]*1000
y0 = (quantities2d0["te"])[jsep0:-4,-2]
x1 = Rrsep1[JSEP:-1,midplaneix1]*1000
y1 = (quantities2d1["te"])[JSEP:-1,-2]
y12 = (quantities2d1["te"])[JSEP:-1,dataue["ix_plate2"]]
kappa = tdiff1actual1/qdsArray1
print(kappa)
TTCALC1 = (TU1**(7/2)-kappa*(q1*np.mean(LArray1)))**(2/7)
TTCALC2 = (TU1**(7/2)-kappa*(q1*LArray1))**(2/7)
plt.plot(x1,LArray1/LArray1[0],marker="^",label="L",linestyle="--", color = "#5C8001")
plt.plot(x1,tdiff1/tdiff1[0],marker="o",label="q"+r'$_{u}$'"L",linestyle="--", color = "#CC503E")
plt.plot(x1,tdiff1actual1/tdiff1actual1[0],marker="o",label="T"+r"$_{u}^{7/2}$"+"-T"+r"$_{t}^{7/2}$", color = "#38A6A5")
nt = quantities2d1["ne"][JSEP:-1,-2]/(quantities2d1["ne"][JSEP:-1,-2][0])
nu = quantities2d1["ne"][JSEP:-1,Xpointix1]/(quantities2d1["ne"][JSEP:-1,Xpointix1][0])
qu = quantities2d1["qpar"][JSEP:-1,Xpointix1]/(quantities2d1["qpar"][JSEP:-1,Xpointix1][0])
Fmom = quantities2d1["ne"][JSEP:-1,-2]*quantities2d1["te"][JSEP:-1,-2]/(quantities2d1["ne"][JSEP:-1,Xpointix1]*quantities2d1["te"][JSEP:-1,Xpointix1])
plt.xlabel("r"+r"$_{u}$"+" (mm)")
plt.legend()
plt.savefig("images//targetTempAnalysis1.png",dpi=1000,bbox_inches='tight')
plt.show()
