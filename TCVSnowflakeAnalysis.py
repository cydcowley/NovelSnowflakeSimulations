
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


LineWidth = 1
LineColor = "#3E434A"

def expFunc(x,A,B):
    return A*np.exp((-x)/B)


# Single Null and XPT balance files
# file0 = "balFiles//balanceNoImp.nc"
file0 = "balFiles//balanceSnowflake.nc"
dataue = ImportGridue("balFiles//gridueSnowflake")

#unpack two balance files
quantities2d,rootgrp = return2d(file0)




#define midplane x point indices
midplaneix = int((dataue["ix_cut2"]-dataue["ix_cut1"])*0.75 + dataue["ix_cut1"])
Xpointix=dataue["ix_cut2"]+1
Xpointix0=73
JSEP = dataue["iyseparatrix2"]+1

Rrsep = quantities2d["r"]-0.5*(rootgrp['crx'][2])[dataue["iyseparatrix1"],midplaneix]-0.5*(rootgrp['crx'][3])[dataue["iyseparatrix1"],midplaneix]

x = Rrsep[JSEP:-1,midplaneix]
x = x*1000
y =quantities2d["te"][JSEP:-1,Xpointix0]
plt.plot(y)
plt.show()
popt, pcov = curve_fit(expFunc, x, y,p0=[np.amax(y),0.002])



fig, axs = plt.subplots(1, 1)
axs.plot(x,y,marker="o",color="#EDAD08")
axs.plot(x,expFunc(x, *popt),label="fit with "+r"$\lambda_{q}$"+"="+str(np.round(popt[1],1))+"mm",color="black")
axs.tick_params(axis='y', labelcolor="#EDAD08")
axs.set_ylabel("q"+r'$_{||,u}$'+" (MWm"+r"$^{-2}$"+")")

axs.set_xlabel("R-Rsep (mm)")
axs.legend()
# axs.set_xlim([-0.07,1.8])
# twin1.set_ylim([45,140])
# axs.set_ylim([-0.05,1.15])
plt.savefig("images//SNupstreamqparSnowflake.png",dpi=1000)
plt.show()


fontsize = 8
# plotquantity = np.log10(np.abs(quantities2d0["imprad"]/quantities2d0["V"])+10)

fig, axs = plt.subplots(1, 1)
axs.set_aspect('equal')
cutix0 =Xpointix0

plotquantity = np.log10(np.abs(quantities2d["te"]))
norm = mpl.colors.Normalize(vmin=np.amin(plotquantity), vmax=np.amax(plotquantity))
cmap = cm.plasma
m = cm.ScalarMappable(norm=norm, cmap=cmap)
plotWALL("balFiles//inputSnowflake.dat",axs)
for j in range(0,len(rootgrp['crx'][0][0])):#len(R[0])-1):
    for i in range(0,len(rootgrp['crx'][0])):
        x = [rootgrp['crx'][0,i,j],rootgrp['crx'][2,i,j],rootgrp['crx'][3,i,j],rootgrp['crx'][1,i,j],rootgrp['crx'][0,i,j]]
        y =[rootgrp['cry'][0,i,j],rootgrp['cry'][2,i,j],rootgrp['cry'][3,i,j],rootgrp['cry'][1,i,j],rootgrp['cry'][0,i,j]]
        color1 = m.to_rgba(plotquantity[i,j])
        axs.fill(x,y,color=color1,linewidth=0.01)  

divider = make_axes_locatable(plt.gca())
ax_cb = divider.new_horizontal(size="5%", pad=0.05)#
# label1 = "q"+r'$_{e,||}$'+" (MWm"+r"$^{-2}$"+")"    
label1 = "log(T"+r'$_{e}$'+")"    
# label1 = "log(imp radiation)"    
plt.plot(dataue['rm'][:dataue["ix_cut1"]+1,dataue["iyseparatrix1"],4],dataue['zm'][:dataue["ix_cut1"]+1,dataue["iyseparatrix1"],4],color=LineColor,linestyle="--",linewidth=LineWidth)
plt.plot(dataue['rm'][dataue["ix_cut1"]+1:dataue["ix_cut2"]+1,dataue["iyseparatrix1"],4],dataue['zm'][dataue["ix_cut1"]+1:dataue["ix_cut2"]+1,dataue["iyseparatrix1"],4],color=LineColor,linestyle="--",linewidth=LineWidth)
plt.plot(dataue['rm'][dataue["ix_cut2"]+1:dataue["ix_cut3"]+1,dataue["iyseparatrix1"],4],dataue['zm'][dataue["ix_cut2"]+1:dataue["ix_cut3"]+1,dataue["iyseparatrix1"],4],color=LineColor,linestyle="--",linewidth=LineWidth)
plt.plot(dataue['rm'][dataue["ix_cut3"]+1:dataue["ix_plate3"],dataue["iyseparatrix1"],4],dataue['zm'][dataue["ix_cut3"]+1:dataue["ix_plate3"],dataue["iyseparatrix1"],4],color=LineColor,linestyle="--",linewidth=LineWidth)
plt.plot(dataue['rm'][dataue["ix_plate3"]:dataue["ix_cut4"]+1,dataue["iyseparatrix1"],4],dataue['zm'][dataue["ix_plate3"]:dataue["ix_cut4"]+1,dataue["iyseparatrix1"],4],color=LineColor,linestyle="--",linewidth=LineWidth)
plt.plot(dataue['rm'][dataue["ix_cut4"]+1:,dataue["iyseparatrix1"],4],dataue['zm'][dataue["ix_cut4"]+1:,dataue["iyseparatrix1"],4],color=LineColor,linestyle="--",linewidth=LineWidth)
# # plot seperatrix 2
plt.plot(dataue['rm'][:dataue["ix_cut1"]+1,dataue["iyseparatrix3"],4],dataue['zm'][:dataue["ix_cut1"]+1,dataue["iyseparatrix3"],4],color=LineColor,linestyle="-.",linewidth=LineWidth)
plt.plot(dataue['rm'][dataue["ix_cut1"]+1:dataue["ix_cut2"]+1,dataue["iyseparatrix3"],4],dataue['zm'][dataue["ix_cut1"]+1:dataue["ix_cut2"]+1,dataue["iyseparatrix3"],4],color=LineColor,linestyle="-.",linewidth=LineWidth)
plt.plot(dataue['rm'][dataue["ix_cut2"]+1:dataue["ix_plate3"],dataue["iyseparatrix3"],4],dataue['zm'][dataue["ix_cut2"]+1:dataue["ix_plate3"],dataue["iyseparatrix3"],4],color=LineColor,linestyle="-.",linewidth=LineWidth)
plt.plot(dataue['rm'][dataue["ix_plate3"]:dataue["ix_cut4"]+1,dataue["iyseparatrix3"],4],dataue['zm'][dataue["ix_plate3"]:dataue["ix_cut4"]+1,dataue["iyseparatrix3"],4],color=LineColor,linestyle="-.",linewidth=LineWidth)
plt.plot(dataue['rm'][dataue["ix_cut4"]+1:,dataue["iyseparatrix3"],4],dataue['zm'][dataue["ix_cut4"]+1:,dataue["iyseparatrix3"],4],color=LineColor,linestyle="-.",linewidth=LineWidth)

cb1 = mpl.colorbar.ColorbarBase(ax_cb, cmap=mpl.cm.plasma, orientation='vertical',norm=norm,label=label1)
plt.annotate("PFCs",(1.95,-0.75),fontsize=fontsize,color="black")
plt.annotate("Pump",(1.99,-1.4),fontsize=fontsize,color="black")
plt.arrow(2.16,-0.7,0.16,0.03,color="black")
plt.arrow(2.22,-1.39,0.09,0.0,color="black")

plt.xlabel("R (m)")
plt.ylabel("Z (m)")
# plt.xlim([1.51,1.88])
# plt.ylim([-1.63,-1.1])
plt.gcf().add_axes(ax_cb)
plt.tight_layout()
plt.savefig("images//2dSnowflakeTe.png",dpi=1000,bbox_inches='tight')
plt.show()


fig, axs = plt.subplots(1, 1)
axs.set_aspect('equal')
cutix0 =Xpointix0
plotquantity = quantities2d["vfluid"]/1000#/((quantities2d["te"]+quantities2d["t"]))
norm = mpl.colors.Normalize(vmin=np.amin(plotquantity), vmax=np.amax(plotquantity))
cmap = cm.plasma
m = cm.ScalarMappable(norm=norm, cmap=cmap)
plotWALL("balFiles//inputSnowflake.dat",axs)
for j in range(0,len(rootgrp['crx'][0][0])):#len(R[0])-1):
    for i in range(0,len(rootgrp['crx'][0])):
        x = [rootgrp['crx'][0,i,j],rootgrp['crx'][2,i,j],rootgrp['crx'][3,i,j],rootgrp['crx'][1,i,j],rootgrp['crx'][0,i,j]]
        y =[rootgrp['cry'][0,i,j],rootgrp['cry'][2,i,j],rootgrp['cry'][3,i,j],rootgrp['cry'][1,i,j],rootgrp['cry'][0,i,j]]
        color1 = m.to_rgba(plotquantity[i,j])
        axs.fill(x,y,color=color1,linewidth=0.01)  

divider = make_axes_locatable(plt.gca())
ax_cb = divider.new_horizontal(size="5%", pad=0.05)#
# label1 = "q"+r'$_{e,||}$'+" (MWm"+r"$^{-2}$"+")"    
label1 = "parallel plasma flow velocity (kms"+r'$^{-1}$' +")"    
# label1 = "log(imp radiation)"    
plt.plot(dataue['rm'][:dataue["ix_cut1"]+1,dataue["iyseparatrix1"],4],dataue['zm'][:dataue["ix_cut1"]+1,dataue["iyseparatrix1"],4],color=LineColor,linestyle="--",linewidth=LineWidth)
plt.plot(dataue['rm'][dataue["ix_cut1"]+1:dataue["ix_cut2"]+1,dataue["iyseparatrix1"],4],dataue['zm'][dataue["ix_cut1"]+1:dataue["ix_cut2"]+1,dataue["iyseparatrix1"],4],color=LineColor,linestyle="--",linewidth=LineWidth)
plt.plot(dataue['rm'][dataue["ix_cut2"]+1:dataue["ix_cut3"]+1,dataue["iyseparatrix1"],4],dataue['zm'][dataue["ix_cut2"]+1:dataue["ix_cut3"]+1,dataue["iyseparatrix1"],4],color=LineColor,linestyle="--",linewidth=LineWidth)
plt.plot(dataue['rm'][dataue["ix_cut3"]+1:dataue["ix_plate3"],dataue["iyseparatrix1"],4],dataue['zm'][dataue["ix_cut3"]+1:dataue["ix_plate3"],dataue["iyseparatrix1"],4],color=LineColor,linestyle="--",linewidth=LineWidth)
plt.plot(dataue['rm'][dataue["ix_plate3"]:dataue["ix_cut4"]+1,dataue["iyseparatrix1"],4],dataue['zm'][dataue["ix_plate3"]:dataue["ix_cut4"]+1,dataue["iyseparatrix1"],4],color=LineColor,linestyle="--",linewidth=LineWidth)
plt.plot(dataue['rm'][dataue["ix_cut4"]+1:,dataue["iyseparatrix1"],4],dataue['zm'][dataue["ix_cut4"]+1:,dataue["iyseparatrix1"],4],color=LineColor,linestyle="--",linewidth=LineWidth)
# # plot seperatrix 2
plt.plot(dataue['rm'][:dataue["ix_cut1"]+1,dataue["iyseparatrix3"],4],dataue['zm'][:dataue["ix_cut1"]+1,dataue["iyseparatrix3"],4],color=LineColor,linestyle="-.",linewidth=LineWidth)
plt.plot(dataue['rm'][dataue["ix_cut1"]+1:dataue["ix_cut2"]+1,dataue["iyseparatrix3"],4],dataue['zm'][dataue["ix_cut1"]+1:dataue["ix_cut2"]+1,dataue["iyseparatrix3"],4],color=LineColor,linestyle="-.",linewidth=LineWidth)
plt.plot(dataue['rm'][dataue["ix_cut2"]+1:dataue["ix_plate3"],dataue["iyseparatrix3"],4],dataue['zm'][dataue["ix_cut2"]+1:dataue["ix_plate3"],dataue["iyseparatrix3"],4],color=LineColor,linestyle="-.",linewidth=LineWidth)
plt.plot(dataue['rm'][dataue["ix_plate3"]:dataue["ix_cut4"]+1,dataue["iyseparatrix3"],4],dataue['zm'][dataue["ix_plate3"]:dataue["ix_cut4"]+1,dataue["iyseparatrix3"],4],color=LineColor,linestyle="-.",linewidth=LineWidth)
plt.plot(dataue['rm'][dataue["ix_cut4"]+1:,dataue["iyseparatrix3"],4],dataue['zm'][dataue["ix_cut4"]+1:,dataue["iyseparatrix3"],4],color=LineColor,linestyle="-.",linewidth=LineWidth)

cb1 = mpl.colorbar.ColorbarBase(ax_cb, cmap=mpl.cm.plasma, orientation='vertical',norm=norm,label=label1)
plt.annotate("PFCs",(1.95,-0.75),fontsize=fontsize,color="black")
plt.annotate("Pump",(1.99,-1.4),fontsize=fontsize,color="black")
plt.arrow(2.16,-0.7,0.16,0.03,color="black")
plt.arrow(2.22,-1.39,0.09,0.0,color="black")

plt.xlabel("R (m)")
plt.ylabel("Z (m)")
# plt.xlim([1.51,1.88])
# plt.ylim([-1.63,-1.1])
plt.gcf().add_axes(ax_cb)
plt.tight_layout()
plt.savefig("images//2dSnowflakeFlux.png",dpi=1000,bbox_inches='tight')
plt.show()
