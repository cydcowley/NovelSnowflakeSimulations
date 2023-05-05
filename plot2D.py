from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import os
import matplotlib.image as m
import imageio as io
from scipy import interpolate
from PIL import Image
import glob
import re
from scipy.integrate import trapz,cumtrapz
import sys
import matplotlib as mpl
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.optimize import curve_fit
plt.rcParams["font.family"] = "serif"
params = {'legend.fontsize': 'medium',
         'axes.labelsize': 'medium',
         'axes.titlesize':'medium',
         'xtick.labelsize':'medium',
         'ytick.labelsize':'medium',
        #  'figure.figsize': (4,3.2),
         }
plt.rcParams.update(params)

def ImportGridue(fname: str = 'gridue') -> dict:
        """
        Import UEDGE grid file as dictionary.

        Parameters
        ----------
        fname : str, optional
            Path/file name to gridue formatted file.

        Returns
        -------
            A dict containing header and body information from the gridue file.

        """
        try:
            f = open(fname, mode='r')
            Values = []
            for i in range(5):
                Values.append([int(x) for x in next(f).split()])
            gridtype = "dn"
            HeaderItems = 0
            
            if gridtype == "dn":
                HeaderItems = ['nxm', 'nym','iyseparatrix1', 'iyseparatrix2',
                'ix_plate1', 'ix_cut1', '_FILLER_', 'ix_cut2', 'ix_plate2',
                'iyseparatrix3', 'iyseparatrix4',
                'ix_plate3', 'ix_cut3', '_FILLER_', 'ix_cut4', 'ix_plate4']
            
            else:
                HeaderItems = ['nxm', 'nym', 'ixpt1', 'ixpt2', 'iyseptrx1']
            flat_list = [item for sublist in Values for item in sublist]
            print(flat_list)
            Values = flat_list
            gridue_settings = dict(zip(HeaderItems, Values))
            print("settings are",gridue_settings)
            next(f)
            BodyItems = ['rm', 'zm', 'psi', 'br', 'bz', 'bpol', 'bphi', 'b']
            Str = {i: [] for i in BodyItems}
            k = iter(Str.keys())
            Key = next(k)
            for line in f:
                if line == 'iogridue\n':
                    continue
                if line == '\n':
                    try:
                        Key = next(k)
                    except:
                        continue

                else:
                    Str[Key].append(line)
            f.close()
            nx = gridue_settings['nxm'] + 2
            ny = gridue_settings['nym'] + 2
            for k, v in Str.items():
                L = (''.join(v).replace('\n', '').replace('D', 'e')).split()
                _l = iter(L)
                vv = next(_l)

                data_ = np.zeros((nx, ny, 5))
                for n in range(5):
                    for j in range(ny):
                        for i in range(nx):

                            data_[i][j][n] = float(vv)

                            try:
                                vv = next(_l)
                            except:
                                continue
                gridue_settings[k] = data_
            return gridue_settings
        except Exception as e:
            print(repr(e))

def plotWALL(filewalls,axis):
    datawall = open(filewalls)
    tline = datawall.readlines(1)
    linenum = 0
    while True:
        tline = datawall.readlines(1)
        if  '*** 3b. Data for additional surfaces' in tline[0]:
            break
    ns = int(datawall.readlines(1)[0])
    coords = []
    ic = 1
    for i in range(ns):
        displayname = datawall.readlines(1)
        rlbnd = int(datawall.readlines(1)[0][1])
        tmp = datawall.readlines(1)
        iliin = int(tmp[0][1:6])

        if iliin <0:
            tmp = datawall.readlines(1)
            continue
        if rlbnd == 2:
            check = datawall.readlines(1)[0]
            co = re.findall("\D+\d+\.\d+\D+\d+",check)
            for i in range(len(co)):
                co[i] = float(co[i])/100
            coords.append(co)
            stype = datawall.readlines(1)
            axis.plot([co[0],co[3]],[co[1],co[4]],color="black",linewidth=0.5)

folder0 = "patchData"
dataue = ImportGridue(folder0+"\\gridue")
def expFunc(x,A,B):
    return A*np.exp((-x)/B)

def plot2d(file,quantity,mode):
    rootgrp =Dataset(file, "r", format="NETCDF4")
    # print(rootgrp)
    bb = rootgrp['bb']
    dv = np.array(rootgrp['vol'])
    quantities2d = {}
    #grid coordinates:
    quantities2d["r"] = (rootgrp['crx'][0]+rootgrp['crx'][1]+rootgrp['crx'][2]+rootgrp['crx'][3])/4

    quantities2d["z"] = (rootgrp['cry'][0]+rootgrp['cry'][1]+rootgrp['cry'][2]+rootgrp['cry'][3])/4

    hx = rootgrp['hx']

    #total magnetic field:
    quantities2d["TotalField"] = np.array(bb[3])
    quantities2d["Bpol"] = np.array(bb[0])
    #length of grid
    s = np.array(hx)*np.abs(np.array(bb[3])/np.array(bb[0]))
    quantities2d["sdiff"] = s
    quantities2d["sdiffpol"] = np.array(hx)
    # Parallel area:
    quantities2d["Area"] = np.array(dv)/s
    #Grid volume
    quantities2d["V"] = dv
    #specific flux ring to focus on
    # print(rootgrp['jsep'][0])
    sep = rootgrp['jsep'][0]
    ring = sep+5
    quantities2d["ring"] = ring
    #electron density (m^{-3})
    quantities2d["ne"] = rootgrp["ne"]
    #ion density (m^{-3})
    quantities2d["ni"] = rootgrp["na"][1]
    #Conductive electron heat flux (Wm^{-2}):
    fhe_cond = rootgrp['fhe_cond'][0]/np.abs(quantities2d["Area"])
    quantities2d["cond"] = np.abs(fhe_cond)
    #electron temperature (eV):
    te = np.array(rootgrp["te"])
    quantities2d["te"] = te/(1.60*10**(-19))
    
    #ion temperature (eV)
    ti = np.array(rootgrp["ti"])
    quantities2d["ti"] = ti/(1.60*10**(-19))
    #artificial impurity radiation (W):
    #imprad = np.sum(rootgrp["b2stel_she_bal"][2:],axis=0)
    imprad = rootgrp["b2stel_she_bal"][1]
    quantities2d["imprad"] = imprad
    #flow velocity
    vfluid = rootgrp["ua"][1]
    quantities2d["vfluid"] = vfluid
    # dab2 = rootgrp['dab2']
    # quantities2d["n0"] = dab2[0]
    quantities2d["qpar"] = rootgrp['fhe_cond'][0]+rootgrp['fhe_32'][0]+rootgrp['fhe_52'][0]+rootgrp['fhe_thermj'][0]+rootgrp['fhe_dia'][0]+rootgrp['fhe_ecrb'][0]
    quantities2d["qpar"] =quantities2d["qpar"] +rootgrp['fhe_strange'][0]+rootgrp['fhe_pschused'][0]
    quantities2d["qpar"] = np.abs(quantities2d["qpar"])/quantities2d["Area"]

    #recombination particle source 
    # R = rootgrp['crx'][1]
    # # print(R[-1])
    # quantities2d["fna"] = rootgrp["fna_tot"][0][1]#/quantities2d["Area"]
    # Z = rootgrp['cry'][1]
    # fnx = rootgrp["fna_pll"][0,0,:,1]
    # fnx2 = rootgrp["fna_pll"][0,1,:,4]

    # midplaneix = np.argmax(quantities2d["r"][-1])
    # # for j in range(midplaneix,len(rootgrp['crx'][0][0])):#len(R[0])-1):
    # #     cutR = []
    # #     cutZ = []
    # #     cutiX = []
    # #     cutiY = []
    # #     for i in range(0,len(rootgrp['crx'][0])):
    # Tu = []
    # for i in range(14,21):
    #     heatflux = np.append(np.abs(quantities2d["qpar"][i,midplaneix:111]),np.abs(quantities2d["qpar"][i,133:]))
    #     Spar = np.append(np.abs(quantities2d["sdiff"][i,midplaneix:111]),np.abs(quantities2d["sdiff"][i,133:]))
    #     Spar = cumtrapz(Spar,initial=0)
    #     Tu.append(trapz(heatflux,Spar))
    # # plt.plot(np.abs(quantities2d["qpar"][18,midplaneix:111]))
    # # plt.plot(np.abs(quantities2d["qpar"][18,133:]))
    #     plt.plot(Spar,heatflux)
    # plt.savefig("heat flux profile.png",dpi=1000)
    # plt.show()

    # print("mid is",midplaneix)
    # Lpar = trapz(np.abs(quantities2d["sdiff"][10,midplaneix:]))
    # print("Lpar is",Lpar)
    # Lpar = trapz(np.abs(quantities2d["sdiff"][10,midplaneix:dataue['ix_cut3']+1]))
    # Lpar = Lpar+trapz(np.abs(quantities2d["sdiff"][10,dataue['ix_cut4']+1:]))
    
    # print("Lpar is",Lpar)
    # # plt.plot(quantities2d["te"][9:,57]/np.amax(quantities2d["te"][9:,57]))
    # # plt.plot(quantities2d["ne"][9:,57]/np.amax(quantities2d["ne"][9:,57]))

    # Rrsep = quantities2d["r"]-0.5*(rootgrp['crx'][2])[8,midplaneix]-0.5*(rootgrp['crx'][3])[8,midplaneix]
    # # plt.plot(Rrsep[14:19,midplaneix]*1000,quantities2d["te"][14:19,midplaneix]**(7/2)-quantities2d["te"][14:19,-3]**(7/2))
    # # plt.plot(quantities2d["qpar"][:,72])
    # # plt.plot(Rrsep[14:,midplaneix]*1000,(quantities2d["ne"]*quantities2d["te"]*quantities2d["vfluid"])[14:,120],marker="o")
    # # plt.plot(Rrsep[14:21,midplaneix]*1000,Tu)
    # plt.plot(Rrsep[10:-1,65],quantities2d["qpar"][10:-1,87],marker="o")
    # # plt.show()
    # # plt.plot((rootgrp['crx'][2])[8,:],(rootgrp['cry'][2])[8,:])
    # # plt.show()
    # # plt.plot(Rrsep[8:,midplaneix]*1000,np.abs(quantities2d["qpar"][8:,midplaneix+13]))
    # # plt.plot(Rrsep[18,midplaneix]*1000,np.abs(quantities2d["qpar"][18,midplaneix+13]),marker="o")
    # # plt.ylim([0.78E10,2.22E10])
    # # axs.set_ylim([-1.8,2.5])
    # plt.ylabel("elec. density at midplane")
    # plt.xlabel("R-Rsep (mm)")
    # plt.savefig("upstream.png")
    # plt.show()

    # plt.plot(Rrsep[10:-1,77],quantities2d["ne"][10:-1,77],marker="o")

    # plt.ylabel("elec. density at midplane")
    # plt.xlabel("R-Rsep (mm)")
    # plt.savefig("upstream.png")
    # plt.show()
    fig, axs = plt.subplots(1, 1)
    # axs.set_xlim([1.3,2.5])
    # axs.set_ylim([-1.8,2.5])
    axs.set_aspect('equal')
    prev = 0

    plotquantity = quantities2d[quantity]
    print(dataue['ix_cut3'])
    cutix0 =0#dataue['ix_cut2']+1
    # cutix0 =73
    # norm = mpl.colors.Normalize(vmin=10, vmax=80)
    norm = 0
    cmap = 0
    # if mode=="rectangular":
    #     norm = mpl.colors.Normalize(vmin=0, vmax=len(rootgrp['crx'][0]))
    #     cmap = cm.plasma
    # else:
    #     norm = mpl.colors.Normalize(vmin=np.amin(plotquantity[1:-1,cutix0:-1]), vmax=np.amax(plotquantity[1:-1,cutix0:-1]))
    #     cmap = cm.plasma
    norm = mpl.colors.Normalize(vmin=np.amin(plotquantity[1:-1,cutix0:-1]), vmax=np.amax(plotquantity[1:-1,cutix0:-1]))
    cmap = cm.plasma
    m = cm.ScalarMappable(norm=norm, cmap=cmap)


    # print(m.to_rgba(x))
    # plotWALL("C:\\Users\\cyd cowley\\Desktop\\PhD\\collaboratory\\Analysis\\input.dat",axs)
    Settings = []
    print(len(rootgrp['crx'][0][0]))
    # print(rootgrp["rightix"][0])
    # plt.plot(quantities2d["r"][10,midplaneix:dataue['ix_cut3']+1],quantities2d["z"][10,midplaneix:dataue['ix_cut3']+1])
    # plt.plot(quantities2d["r"][10,dataue['ix_cut4']+1:],quantities2d["z"][10,dataue['ix_cut4']+1:])
    
    for j in range(cutix0,len(rootgrp['crx'][0][0])):#len(R[0])-1):
        cutR = []
        cutZ = []
        cutiX = []
        cutiY = []
       
        for i in range(0,len(rootgrp['crx'][0])):
            # plt.plot(rootgrp['crx'][0,i,cutix0:],rootgrp['cry'][0,i,cutix0:],color="C0",linewidth=0.1)
            # fig, axs = plt.subplots(1, 1)
            # print(i)
            # print(rootgrp["bottomiy"][i][j]+1)
            # print(rootgrp["topiy"][i][j]+1)
            # if rootgrp["bottomiy"][i][j]+1 != i-1:
            #     print(i)
            # if rootgrp["topiy"][i][j]+1 != i+1:
            #     print(i)
            # if quantities2d["fna"][i][j] ==0:
            #     print(i,j)
            #     print(rootgrp["rightix"][i][j])
            #     print(rootgrp["leftix"][i][j])
            if rootgrp["leftix"][i][j]+1 != j-1:
                # print("for index", j)   
                # print("neighbor is",rootgrp["leftix"][i][j])
                cutR.append(rootgrp['crx'][0][i][j])
                cutR.append(rootgrp['crx'][2][i][j])
                cutZ.append(rootgrp['cry'][0][i][j])
                cutZ.append(rootgrp['cry'][2][i][j])
                # print(j-1)
                # print(rootgrp["leftix"][i][j])  
                cutiX.append(j)
                cutiY.append(i)
            # if rootgrp["rightix"][i][j]== len(rootgrp['crx'][0][0]):
                # print("-2 exists")
                # print(j)

            if rootgrp["rightix"][i][j]+1 != j+1:
                cutR.append(rootgrp['crx'][1][i][j])
                cutR.append(rootgrp['crx'][3][i][j])
                cutZ.append(rootgrp['cry'][1][i][j])
                cutZ.append(rootgrp['cry'][3][i][j])
                # print(j)
                cutiX.append(j+1)
                cutiY.append(i)
                # print(j-1)
                # print(rootgrp["rightix"][i][j])  
            #     print(rootgrp["rightix"][i][j]+1)      
            # if rootgrp["topix"][i][j] != 
            # print(rootgrp['crx'][1,0,i])
            # print(rootgrp['crx'][0,0,i+1])

            x = [rootgrp['crx'][0,i,j],rootgrp['crx'][2,i,j],rootgrp['crx'][3,i,j],rootgrp['crx'][1,i,j],rootgrp['crx'][0,i,j]]
            y =[rootgrp['cry'][0,i,j],rootgrp['cry'][2,i,j],rootgrp['cry'][3,i,j],rootgrp['cry'][1,i,j],rootgrp['cry'][0,i,j]]
            ix = [j,j+1,j+1,j,j]
            iy =[i,i,i+1,i+1,i]
            
            # plt.plot(x,y)
            color1 = m.to_rgba(plotquantity[i,j])
            # color1 = m2.to_rgba(i)
            if mode=="rectangular":
                axs.fill(ix,iy,color=color1,linewidth=0.01)  
            else:
                axs.fill(x,y,color=color1,linewidth=0.01)  
        
            # axs.fill(x,y,color=color1,linewidth=0.01)  
            # plt.fill_between([x[0],x[0]])
            # plt.plot(x,y,color="black",linewidth=0.1)
            # plt.savefig("firstrun.png",dpi=1000)
            # plt.show()
            # if rootgrp['crx'][1,0,i] != rootgrp['crx'][0,0,i+1] or i==len(R[0])-2:
            #     print(prev)
            #     print(i)
            #     plt.pcolormesh(R[:,prev:i+2],Z[:,prev:i+2],quantities2d["te"][:,prev:i+2])   
            #     prev = i-1
        # if mode=="rectangular":
        #     plt.plot(cutiX,cutiY,color="black",linewidth=0.5)
        # else:
        #     plt.plot(cutR,cutZ,color="black",linewidth=0.5)
    # plt.colorbar()
    divider = make_axes_locatable(plt.gca())
    ax_cb = divider.new_horizontal(size="5%", pad=0.05)#
    label1 = "q"+r'$_{e,||}$'+" (MWm"+r"$^{-2}$"+")"    
    # label1 = "T"+r'$_{e}$'+" (eV)"    
    cb1 = mpl.colorbar.ColorbarBase(ax_cb, cmap=mpl.cm.plasma, orientation='vertical',norm=norm,label=label1)
    plt.xlabel("R (m)")
    plt.ylabel("Z (m)")
    # plt.xlim([1.51,1.87])
    # plt.ylim([-1.61,-1.1])

    # plt.xlabel("ix")
    # plt.ylabel("iy")
    plt.gcf().add_axes(ax_cb)
    fname = quantity+".png"
    plt.savefig(fname,dpi=1000)

def plotUpstream(file,grid="XPT"):
    rootgrp =Dataset(file, "r", format="NETCDF4")
    # print(rootgrp)
    bb = rootgrp['bb']
    dv = np.array(rootgrp['vol'])
    quantities2d = {}
    #grid coordinates:
    quantities2d["r"] = (rootgrp['crx'][0]+rootgrp['crx'][1]+rootgrp['crx'][2]+rootgrp['crx'][3])/4

    quantities2d["z"] = (rootgrp['cry'][0]+rootgrp['cry'][1]+rootgrp['cry'][2]+rootgrp['cry'][3])/4
    quantities2d["z"] = 0.5*quantities2d["z"]
    hx = rootgrp['hx']

    #total magnetic field:
    quantities2d["TotalField"] = np.array(bb[3])
    quantities2d["Bpol"] = np.array(bb[0])
    #length of grid
    s = np.array(hx)*np.abs(np.array(bb[3])/np.array(bb[0]))
    quantities2d["sdiff"] = s

    quantities2d["sdiffpol"] = np.array(hx)
    # Parallel area:
    quantities2d["Area"] = np.array(dv)/s
    #Grid volume
    quantities2d["V"] = dv
    #specific flux ring to focus on
    # print(rootgrp['jsep'][0])
    sep = rootgrp['jsep'][0]
    print(sep)
    ring = sep+5
    quantities2d["ring"] = ring
    #electron density (m^{-3})
    quantities2d["ne"] = rootgrp["ne"]
    #ion density (m^{-3})
    quantities2d["ni"] = rootgrp["na"][1]
    #Conductive electron heat flux (Wm^{-2}):
    fhe_cond = rootgrp['fhe_cond'][0]/np.abs(quantities2d["Area"])
    quantities2d["cond"] = np.abs(fhe_cond)
    #electron temperature (eV):
    te = np.array(rootgrp["te"])
    quantities2d["te"] = te/(1.60*10**(-19))
    print(len(rootgrp["na"]))
    #ion temperature (eV)
    ti = np.array(rootgrp["ti"])
    quantities2d["ti"] = ti/(1.60*10**(-19))
    #artificial impurity radiation (W):
    imprad = np.sum(rootgrp["b2stel_she_bal"],axis=0)+np.sum(rootgrp['eirene_mc_eael_she_bal'],axis=0)
    # imprad = rootgrp["b2stel_she_bal"]
    quantities2d["imprad"] = imprad
    #flow velocity
    vfluid = rootgrp["ua"][1]
    quantities2d["vfluid"] = vfluid

    quantities2d["qpar"] = rootgrp['fhe_cond'][0]+rootgrp['fhe_32'][0]+rootgrp['fhe_52'][0]+rootgrp['fhe_thermj'][0]+rootgrp['fhe_dia'][0]+rootgrp['fhe_ecrb'][0]
    quantities2d["qpar"] =quantities2d["qpar"] +rootgrp['fhe_strange'][0]+rootgrp['fhe_pschused'][0]
    quantities2d["qpar"] = np.abs(quantities2d["qpar"])/quantities2d["Area"]

    dab2 = rootgrp['dab2']
    quantities2d["n0"] = dab2[0]

    midplaneix =  np.argmax(quantities2d["r"][-1])
    Xpointix=0
    if grid=="XPT":
        Xpointix = dataue["ix_cut2"]+1
    else:
        Xpointix = 73
    print(Xpointix)
    Tu = []
    # for i in range(14,21):
    #     heatflux = np.append(np.abs(quantities2d["qpar"][i,midplaneix:111]),np.abs(quantities2d["qpar"][i,133:]))
    #     Spar = np.append(np.abs(quantities2d["sdiff"][i,midplaneix:111]),np.abs(quantities2d["sdiff"][i,133:]))
    #     Spar = cumtrapz(Spar,initial=0)
    #     Tu.append(trapz(heatflux,Spar))

    Rrsep = quantities2d["r"]-0.5*(rootgrp['crx'][2])[8,midplaneix]-0.5*(rootgrp['crx'][3])[8,midplaneix]
    JSEP = 9


    # plt.plot((quantities2d["imprad"])[JSEP:-1,-1])
    # plt.plot((quantities2d["imprad"])[JSEP:-1,-2])
    # plt.savefig("guard.png",dpi=500)
    # plt.show()
    # JSEP = 0
    x = Rrsep[JSEP:-1,midplaneix]
    x = x*1000
    label0 = "q"+r'$_{||,u}$'+" (MWm"+r"$^{-2}$"+")"  
    label1 = "n"+r'$_{e,u}$'+" (m"+r"$^{-3}$"+")"   
    y0 = (quantities2d["qpar"])[JSEP:-1,Xpointix]
    y1 = (quantities2d["ne"])[JSEP:-1,midplaneix]
    y2 = (quantities2d["te"])[JSEP:-1,midplaneix]
    y3 = (quantities2d["qpar"])[JSEP:-1,-2]
    popt, pcov = curve_fit(expFunc, x, y0,p0=[np.amax(y0),0.002])
    
    popt1, pcov = curve_fit(expFunc, x, y1,p0=[np.amax(y1),0.002])
    plt.plot(x,y0,marker="o",linestyle="")
    plt.plot(x,expFunc(x, *popt),label="fit with "+r"$\lambda_{q}$"+"="+str(np.round(popt[1],1))+"mm")
    plt.ylabel(label0)
    plt.xlabel("R-Rsep (mm)")
    plt.legend()

    plt.xlim([-0.07,1.8])
    plt.ylim([0,1.15E9])
    plt.savefig("upstreamqpar.png",dpi=1000)
    plt.show()
    plt.plot(x,y1,marker="o",linestyle="")
    plt.plot(x,expFunc(x, *popt1),label="fit with "+r"$\lambda_{n}$"+"="+str(np.round(popt1[1]/10,0))+"cm")
    plt.ylabel(label1)
    plt.xlabel("R-Rsep (mm)")
    plt.legend()
    plt.xlim([-0.07,1.8])
    plt.ylim([4.74E19,5.02E19])
    plt.savefig("upstreamdens.png",dpi=1000)
    plt.show()
    plt.plot(x,y2,marker="o",linestyle="")
    # plt.plot(x,expFunc(x, *popt1),label="fit with lamda n="+str(np.round(popt1[1],3))+"m")
    plt.ylabel("elec. heat flux at midplane")
    plt.xlabel("R-Rsep (mm)")
    plt.legend()
    plt.savefig("upstreamTemp.png",dpi=1000)
    plt.show()
    plt.plot(x,y3,marker="o",linestyle="")
    plt.ylabel("Tet")
    plt.xlabel("R-Rsep (mm)")
    plt.legend()
    plt.xlim([-0.07,1.8])
    # plt.ylim([0,115])
    plt.savefig("targetTemp.png",dpi=1000)
    plt.show()
    print(pcov)
    print("heat is",trapz(y0[:],x[:])*10)
# plot2d("balFiles//balance.nc","qpar",mode="")
# plot2d("balFiles//balance2MWStandard.nc","te",mode="")
# plot2d("balFiles//balance5MWStandard.nc","qpar",mode="")
# plotUpstream("balFiles//balance2MWStandard.nc",grid="SD")
# plotUpstream("balFiles//balance2MWXPT.nc",grid="XPT")
# plotUpstream("balFiles//balance3MWXPT.nc",grid="XPT")
# plotUpstream("balFiles//balance3MWSN.nc",grid="SN")

