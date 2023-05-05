
from itertools import count
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import trapz,cumtrapz
import re
import matplotlib as mpl
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.optimize import curve_fit
from SharedFunctions import return2d,ImportGridue,plotWALL

plt.rcParams["font.family"] = "serif"
params = {'legend.fontsize': 'small',
         'axes.labelsize': 'small',
         'axes.titlesize':'small',
         'xtick.labelsize':'small',
         'ytick.labelsize':'small',
        #  'figure.figsize': (4,3.2),
         }
plt.rcParams.update(params)





def expFunc(x,A,B):
    return A*np.exp((-x)/B)


# Acess wall loading data from b2plot
LoaddataSmallSep = np.loadtxt("balFiles//wallheatSmallSep",skiprows=1)
LoaddataMedSep = np.loadtxt("balFiles//wallheatMedSep",skiprows=1)
LoaddataLargeSep = np.loadtxt("balFiles//wallheatLargeSep",skiprows=1)
Target2LoadSmallSep = np.sum(LoaddataSmallSep[131:151,2])
# Target2LoadSmallSep = np.sum(LoaddataSmallSep[213:233,2])
Target1LoadSmallSep = np.sum(LoaddataSmallSep[100:120,2])
Target3LoadSmallSep = np.sum(LoaddataSmallSep[167:187,2])
Target2LoadMedSep = np.sum(LoaddataMedSep[131:151,2])
# Target2LoadMedSep = np.sum(LoaddataMedSep[213:233,2])
Target1LoadMedSep = np.sum(LoaddataMedSep[100:120,2])
Target3LoadMedSep = np.sum(LoaddataMedSep[167:187,2])
Target2LoadLargeSep = np.sum(LoaddataLargeSep[131:151,2])
# Target2LoadLargeSep = np.sum(LoaddataLargeSep[213:233,2])
Target1LoadLargeSep = np.sum(LoaddataLargeSep[100:120,2])
Target3LoadLargeSep = np.sum(LoaddataLargeSep[167:187,2])



files = ["balFiles//balanceMastSmallSep.nc","balFiles//balanceMastMedSep.nc","balFiles//balanceMastLargeSep.nc"]
uefiles = ["balFiles//gridueMASTSmallSep","balFiles//gridueMASTMedSep","balFiles//gridueMASTLargeSep"]

fontsize = 8
# plotquantity = np.log10(np.abs(quantities2d0["imprad"]/quantities2d0["V"])+10)

fig, axs = plt.subplots(1, 3)
heat_plate_1 = []
heat_plate_2 = []
heat_plate_3 = []
XptSeparation = []
MaxHeat = []
TargetTemp = []
MaxTemp = []
Rrseps = []
LineWidth = 1
LineColor = "#3E434A"
for a in range(3):
    dataue = ImportGridue(uefiles[a])
    JSEP = dataue["iyseparatrix2"]+1
    JSEP2 = dataue["iyseparatrix3"]+1

    plotWALL("balFiles//inputMAST.dat",axs[a])
    file0 = files[a]
    axs[a].set_aspect('equal')
    axs[a].set_xlim([0.4,0.9])
    axs[a].set_ylim([-1.7,-1])
    if a!=0:
        axs[a].xaxis.set_visible(False)
        axs[a].yaxis.set_visible(False)
    else:
        axs[a].set_xlabel("R (m)")
        axs[a].set_ylabel("Z (m)")
    quantities2d,rootgrp = return2d(file0)
    plotquantity = np.log10(quantities2d["te"])
    norm = mpl.colors.Normalize(vmin=np.amin(plotquantity), vmax=np.amax(plotquantity))
    cmap = cm.plasma
    m = cm.ScalarMappable(norm=norm, cmap=cmap)
    for j in range(0,len(rootgrp['crx'][0][0])):#len(R[0])-1):
        for i in range(0,len(rootgrp['crx'][0])):
            x = [rootgrp['crx'][0,i,j],rootgrp['crx'][2,i,j],rootgrp['crx'][3,i,j],rootgrp['crx'][1,i,j],rootgrp['crx'][0,i,j]]
            y =[rootgrp['cry'][0,i,j],rootgrp['cry'][2,i,j],rootgrp['cry'][3,i,j],rootgrp['cry'][1,i,j],rootgrp['cry'][0,i,j]]
            color1 = m.to_rgba(plotquantity[i,j])
            axs[a].fill(x,y,color=color1,linewidth=0.01)  
    
    # # plot seperatrix 1
    axs[a].plot(dataue['rm'][:dataue["ix_cut3"]+1,dataue["iyseparatrix1"],4],dataue['zm'][:dataue["ix_cut3"]+1,dataue["iyseparatrix1"],4],color=LineColor,linestyle="--",linewidth=LineWidth)
    axs[a].plot(dataue['rm'][dataue["ix_cut3"]+1:dataue["ix_plate3"],dataue["iyseparatrix1"],4],dataue['zm'][dataue["ix_cut3"]+1:dataue["ix_plate3"],dataue["iyseparatrix1"],4],color=LineColor,linestyle="--",linewidth=LineWidth)
    axs[a].plot(dataue['rm'][dataue["ix_plate3"]:dataue["ix_cut4"]+1,dataue["iyseparatrix1"],4],dataue['zm'][dataue["ix_plate3"]:dataue["ix_cut4"]+1,dataue["iyseparatrix1"],4],color=LineColor,linestyle="--",linewidth=LineWidth)
    axs[a].plot(dataue['rm'][dataue["ix_cut4"]+1:,dataue["iyseparatrix1"],4],dataue['zm'][dataue["ix_cut4"]+1:,dataue["iyseparatrix1"],4],color=LineColor,linestyle="--",linewidth=LineWidth)
    # # plot seperatrix 2
    axs[a].plot(dataue['rm'][:dataue["ix_cut3"]+1,dataue["iyseparatrix3"],4],dataue['zm'][:dataue["ix_cut3"]+1,dataue["iyseparatrix3"],4],color=LineColor,linestyle="-.",linewidth=LineWidth)
    axs[a].plot(dataue['rm'][dataue["ix_cut3"]+1:dataue["ix_plate3"],dataue["iyseparatrix3"],4],dataue['zm'][dataue["ix_cut3"]+1:dataue["ix_plate3"],dataue["iyseparatrix3"],4],color=LineColor,linestyle="-.",linewidth=LineWidth)
    axs[a].plot(dataue['rm'][dataue["ix_plate3"]:dataue["ix_cut4"]+1,dataue["iyseparatrix3"],4],dataue['zm'][dataue["ix_plate3"]:dataue["ix_cut4"]+1,dataue["iyseparatrix3"],4],color=LineColor,linestyle="-.",linewidth=LineWidth)
    axs[a].plot(dataue['rm'][dataue["ix_cut4"]+1:,dataue["iyseparatrix3"],4],dataue['zm'][dataue["ix_cut4"]+1:,dataue["iyseparatrix3"],4],color=LineColor,linestyle="-.",linewidth=LineWidth)

    Q1 = quantities2d["qpar"][:,dataue['ix_plate2']-1]*quantities2d["Area"][:,dataue['ix_plate2']-1] 
    Q2 = quantities2d["qpar"][:,dataue['ix_plate3']+1]*quantities2d["Area"][:,dataue['ix_plate3']+1] 
    Q3 = quantities2d["qpar"][:,dataue['ix_plate4']-1]*quantities2d["Area"][:,dataue['ix_plate4']-1] 
    T1 = np.array(quantities2d["qpar"][:,dataue['ix_plate2']-1])
    T2 =  np.array(quantities2d["te"][:,dataue['ix_plate3']+1])
    T3 =  np.array(quantities2d["qpar"][:,dataue['ix_plate4']-1])
    MaxHeat.append(np.max(np.array([ quantities2d["qpar"][:,dataue['ix_plate2']-1],quantities2d["qpar"][:,dataue['ix_plate3']+1],quantities2d["qpar"][:,dataue['ix_plate4']-1]])))
    MaxTemp.append(np.max(np.array([ quantities2d["te"][:,dataue['ix_plate2']-1],quantities2d["te"][:,dataue['ix_plate3']+1],quantities2d["te"][:,dataue['ix_plate4']-1]])))
    TargetTemp.append([T3,T1])
    midplaneix = np.argmax(quantities2d["r"][-1])
    Rrsep = quantities2d["r"]-0.5*(rootgrp['crx'][2])[dataue["iyseparatrix1"],midplaneix]-0.5*(rootgrp['crx'][3])[dataue["iyseparatrix1"],midplaneix]
    Rrseps.append([Rrsep[:,midplaneix],Rrsep[:,midplaneix]])
    heat_plate_1.append(np.sum(Q1))
    heat_plate_2.append(np.sum(Q2))
    heat_plate_3.append(np.sum(Q3))
    XptSeparation.append(Rrsep[JSEP2,midplaneix]-Rrsep[JSEP,midplaneix])
    x = Rrsep[JSEP-2:-1,midplaneix]   
    print("Q1 IS ",np.sum(Q1))
    print("Q2 IS ",np.sum(Q2))
    print("Q3 IS ",np.sum(Q3))

label1 = "log(T"+r'$_{e}$'+")"    
plt.annotate("Target 1",(0.765,-1.32),fontsize=fontsize,color="C0")
plt.annotate("Target 2",(0.715,-1.57),fontsize=fontsize,color="C1")
plt.annotate("Target 3",(0.71,-1.64),fontsize=fontsize,color="C2")
plt.arrow(0.8,-1.34,0.03,-0.05,color="C0")
plt.arrow(0.76,-1.54,0.05,0.05,color="C1")
plt.arrow(0.71,-1.635,-0.025,0.0,color="C2")

cbarax = plt.axes([1.01, 0.2, 0.01, 0.6], facecolor='none')
cb1 = mpl.colorbar.ColorbarBase(cbarax, cmap=mpl.cm.plasma, orientation='vertical',norm=norm,label=label1)

plt.subplots_adjust(wspace=0, hspace=0)
plt.tight_layout(pad=0)
plt.savefig("images//MASTgrids.png",dpi=2000,bbox_inches='tight')
plt.show()


# plot fractional target loads
heat_plate_1 = np.array([Target1LoadSmallSep,Target1LoadMedSep,Target1LoadLargeSep])
heat_plate_2 = np.array([Target2LoadSmallSep,Target2LoadMedSep,Target2LoadLargeSep])
heat_plate_3 = np.array([Target3LoadSmallSep,Target3LoadMedSep,Target3LoadLargeSep])
XptSeparation = np.array(XptSeparation)*1000
totHeat = heat_plate_1+heat_plate_2+heat_plate_3
plt.plot(XptSeparation,heat_plate_1/totHeat,label="Target 1",marker="o",LineWidth=2,linestyle=":")
plt.plot(XptSeparation,heat_plate_2/totHeat,label= "Target 2",marker="^",LineWidth=2,linestyle=":")
plt.plot(XptSeparation,heat_plate_3/totHeat,label="Target 3",marker="s",LineWidth=2,linestyle=":")
plt.xlabel("X point separation (mm)")
plt.ylabel("fraction of outer divertor power deposited")
plt.legend()
plt.savefig("images//MASTPowerSharing.png",dpi=1000,bbox_inches='tight')
plt.show()




import re
import numpy as np
fig, axs = plt.subplots(1, 3)

counter = 0
for fname in ["SmallSep","MedSep","LargeSep"]:
    #get neutral density
    fortFile = "balFiles//fortMAST"+fname+".46"
    dataFort= open(fortFile)
    tline = dataFort.readlines(1)
    linenum = 0
    while True:
        tline = dataFort.readlines(1)
        if  '*eirene data field pdena' in tline[0]:
            break

    Atomdensities = []
    Moldensities = []
    molecules = 0
    while True:
        tline = dataFort.readlines(1)
        co = re.findall("\D+\d+\.\d+\D+\d+",tline[0])

        for i in range(len(co)):
            if molecules:
                Moldensities.append(float(co[i]))
            else:
                Atomdensities.append(float(co[i]))

        if  '*eirene' in tline[0] and molecules==1:
            break
        if  '*eirene data field pdenm' in tline[0]:
            molecules = 1


    densities = np.array(Atomdensities)+np.array(Moldensities)
    densities = densities*10**6

    # get triangle x/y data
    fd = open('balFiles\\fortMAST'+fname+'.33','r')   
    dataTriangImport =  fd.readlines()
    dataTriang = []
    for i in range(len(dataTriangImport)):
        co = re.findall("\D+\d+\.\d+\D+\d+",dataTriangImport[i])
        for j in range(len(co)):
            dataTriang.append(float(co[j]))
    dataTriang = np.array(dataTriang)/100

    datax = dataTriang[0:int(len(dataTriang)/2)]
    datay = dataTriang[int(len(dataTriang)/2):]

    # get triangle indices
    fd = open('balFiles\\fortMAST'+fname+ '.34','r')   
    indices = np.loadtxt(fd,skiprows=1,usecols=(1,2,3))
    axs[counter].set_aspect('equal')
    axs[counter].set_xlim([0.4,0.9])
    axs[counter].set_ylim([-1.7,-1])
    if counter!=0:
        axs[counter].xaxis.set_visible(False)
        axs[counter].yaxis.set_visible(False)
    else:
        axs[counter].set_xlabel("R (m)")
        axs[counter].set_ylabel("Z (m)")
    plotquantity = (densities)
    norm = mpl.colors.Normalize(vmin=0.0, vmax=6.0E17)
    cmap = cm.plasma
    m = cm.ScalarMappable(norm=norm, cmap=cmap)
    for i in range(len(indices)):
        ind0 = int(indices[i][0])-1
        ind1 = int(indices[i][1])-1
        ind2 = int(indices[i][2])-1
        x = [datax[ind0],datax[ind1],datax[ind2]]
        y =[datay[ind0],datay[ind1],datay[ind2]]
        color1 = m.to_rgba(plotquantity[i])
        axs[counter].fill(x,y,color=color1,linewidth=0.01)
    plotWALL("balFiles//inputMAST.dat",axs[counter])
    counter = counter+1

cbarax = plt.axes([1.01, 0.22, 0.01, 0.6], facecolor='none')
cb1 = mpl.colorbar.ColorbarBase(cbarax, cmap=mpl.cm.plasma, orientation='vertical',norm=norm,label="total neutral density (m"+r'$^{-3}$'+")")  

plt.subplots_adjust(wspace=0, hspace=0)
plt.tight_layout(pad=0)
plt.savefig("images//MASTNeutralDens.png",dpi=1000,bbox_inches='tight')
plt.show()
