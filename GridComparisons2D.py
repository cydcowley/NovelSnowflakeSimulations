
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

folder0 = "patchData"
dataue = ImportGridue(folder0+"\\gridue")

def expFunc(x,A,B):
    return A*np.exp((-x)/B)


# Single Null and XPT balance files
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

fontsize = 8
fig, axs = plt.subplots(1, 1)
axs.set_aspect('equal')
cutix0 =Xpointix0
plotquantity = quantities2d0["qpar"]

norm = mpl.colors.Normalize(vmin=0, vmax=80)
cmap = cm.plasma
m = cm.ScalarMappable(norm=norm, cmap=cmap)
for j in range(0,len(rootgrp0['crx'][0][0])):#len(R[0])-1):
    for i in range(0,len(rootgrp0['crx'][0])):
        x = [rootgrp0['crx'][0,i,j],rootgrp0['crx'][2,i,j],rootgrp0['crx'][3,i,j],rootgrp0['crx'][1,i,j],rootgrp0['crx'][0,i,j]]
        y =[rootgrp0['cry'][0,i,j],rootgrp0['cry'][2,i,j],rootgrp0['cry'][3,i,j],rootgrp0['cry'][1,i,j],rootgrp0['cry'][0,i,j]]
        color1 = m.to_rgba(plotquantity[i,j])
        axs.fill(x,y,color=color1,linewidth=0.01)  

divider = make_axes_locatable(plt.gca())
ax_cb = divider.new_horizontal(size="5%", pad=0.05)#
label1 = "q"+r'$_{e,||}$'+" (MWm"+r"$^{-2}$"+")"    
# label1 = "T"+r'$_{e}$'+" (eV)"    
cb1 = mpl.colorbar.ColorbarBase(ax_cb, cmap=mpl.cm.plasma, orientation='vertical',norm=norm,label=label1)

plt.annotate("PFCs",(1.6,-1.5),fontsize=fontsize,color="black")
plt.annotate("Pump Chamber",(1.72,-1.32),fontsize=fontsize,color="black")
plt.arrow(1.77,-1.33,0.03,-0.033,color="black")
plt.arrow(1.62,-1.48,0.03,0.033,color="black")
plotWALL("balFiles//input.dat",axs)
plt.xlabel("R (m)")
plt.ylabel("Z (m)")
plt.xlim([1.51,1.9])
plt.ylim([-1.63,-1.1])
plt.gcf().add_axes(ax_cb)
plt.tight_layout()
plt.savefig("images//qpar2dSN.png",dpi=1000,bbox_inches='tight')
plt.show()

fig, axs = plt.subplots(1, 1)
axs.set_aspect('equal')
cutix1 =Xpointix1
plotquantity = quantities2d1["qpar"]
norm = mpl.colors.Normalize(vmin=0, vmax=80)
cmap = cm.plasma

m = cm.ScalarMappable(norm=norm, cmap=cmap)
for j in range(0,len(rootgrp1['crx'][0][0])):#len(R[0])-1):
    for i in range(0,len(rootgrp1['crx'][0])):
        x = [rootgrp1['crx'][0,i,j],rootgrp1['crx'][2,i,j],rootgrp1['crx'][3,i,j],rootgrp1['crx'][1,i,j],rootgrp1['crx'][0,i,j]]
        y =[rootgrp1['cry'][0,i,j],rootgrp1['cry'][2,i,j],rootgrp1['cry'][3,i,j],rootgrp1['cry'][1,i,j],rootgrp1['cry'][0,i,j]]
        color1 = m.to_rgba(plotquantity[i,j])
        axs.fill(x,y,color=color1,linewidth=0.01)  
# # plot seperatrix 1
axs.plot(dataue['rm'][:dataue["ix_cut3"]+1,dataue["iyseparatrix1"],4],dataue['zm'][:dataue["ix_cut3"]+1,dataue["iyseparatrix1"],4],color=LineColor,linestyle="--",linewidth=LineWidth)
axs.plot(dataue['rm'][dataue["ix_cut3"]+1:dataue["ix_plate3"],dataue["iyseparatrix1"],4],dataue['zm'][dataue["ix_cut3"]+1:dataue["ix_plate3"],dataue["iyseparatrix1"],4],color=LineColor,linestyle="--",linewidth=LineWidth)
axs.plot(dataue['rm'][dataue["ix_plate3"]:dataue["ix_cut4"]+1,dataue["iyseparatrix1"],4],dataue['zm'][dataue["ix_plate3"]:dataue["ix_cut4"]+1,dataue["iyseparatrix1"],4],color=LineColor,linestyle="--",linewidth=LineWidth)
axs.plot(dataue['rm'][dataue["ix_cut4"]+1:,dataue["iyseparatrix1"],4],dataue['zm'][dataue["ix_cut4"]+1:,dataue["iyseparatrix1"],4],color=LineColor,linestyle="--",linewidth=LineWidth)
# # plot seperatrix 2
axs.plot(dataue['rm'][:dataue["ix_cut3"]+1,dataue["iyseparatrix3"],4],dataue['zm'][:dataue["ix_cut3"]+1,dataue["iyseparatrix3"],4],color=LineColor,linestyle="-.",linewidth=LineWidth)
axs.plot(dataue['rm'][dataue["ix_cut3"]+1:dataue["ix_plate3"],dataue["iyseparatrix3"],4],dataue['zm'][dataue["ix_cut3"]+1:dataue["ix_plate3"],dataue["iyseparatrix3"],4],color=LineColor,linestyle="-.",linewidth=LineWidth)
axs.plot(dataue['rm'][dataue["ix_plate3"]:dataue["ix_cut4"]+1,dataue["iyseparatrix3"],4],dataue['zm'][dataue["ix_plate3"]:dataue["ix_cut4"]+1,dataue["iyseparatrix3"],4],color=LineColor,linestyle="-.",linewidth=LineWidth)
axs.plot(dataue['rm'][dataue["ix_cut4"]+1:,dataue["iyseparatrix3"],4],dataue['zm'][dataue["ix_cut4"]+1:,dataue["iyseparatrix3"],4],color=LineColor,linestyle="-.",linewidth=LineWidth)

divider = make_axes_locatable(plt.gca())
ax_cb = divider.new_horizontal(size="5%", pad=0.05)#
label1 = "q"+r'$_{e,||}$'+" (MWm"+r"$^{-2}$"+")"    
# label1 = "T"+r'$_{e}$'+" (eV)"    
cb1 = mpl.colorbar.ColorbarBase(ax_cb, cmap=mpl.cm.plasma, orientation='vertical',norm=norm,label=label1)
plotWALL("balFiles//input.dat",axs)
plt.annotate("PFCs",(1.6,-1.5),fontsize=fontsize,color="black")
plt.annotate("Pump Chamber",(1.72,-1.32),fontsize=fontsize,color="black")
plt.annotate("Target 1 ",(1.62,-1.59),fontsize=fontsize,color="black")
plt.annotate("Target 2 ",(1.82,-1.36),fontsize=fontsize,color="black")

plt.arrow(1.77,-1.33,0.03,-0.033,color="black") # pump label
plt.arrow(1.62,-1.48,0.03,0.033,color="black") # PFC label
plt.arrow(1.7,-1.58,0.04,0.023,color="black") # XPT Target 1 label
plt.arrow(1.867,-1.37,-0.014,-0.063,color="black") #XPT Target 2 label

plt.xlabel("R (m)")
plt.ylabel("Z (m)")
plt.xlim([1.51,1.9])
plt.ylim([-1.63,-1.1])
plt.gcf().add_axes(ax_cb)
plt.tight_layout()
plt.savefig("images//qpar2dXPT.png",dpi=1000,bbox_inches='tight')
plt.show()

# fig, axs = plt.subplots(1, 1)
# axs.set_aspect('equal')
# plotquantity = quantities2d["te"]
# norm = mpl.colors.Normalize(vmin=0, vmax=175)
# cmap = cm.plasma
# m = cm.ScalarMappable(norm=norm, cmap=cmap)
# for j in range(cutix0,len(rootgrp['crx'][0][0])):
#     for i in range(0,len(rootgrp['crx'][0])):
#         x = [rootgrp['crx'][0,i,j],rootgrp['crx'][2,i,j],rootgrp['crx'][3,i,j],rootgrp['crx'][1,i,j],rootgrp['crx'][0,i,j]]
#         y =[rootgrp['cry'][0,i,j],rootgrp['cry'][2,i,j],rootgrp['cry'][3,i,j],rootgrp['cry'][1,i,j],rootgrp['cry'][0,i,j]]
#         color1 = m.to_rgba(plotquantity[i,j])
#         axs.fill(x,y,color=color1,linewidth=0.01)  

# divider = make_axes_locatable(plt.gca())
# ax_cb = divider.new_horizontal(size="5%", pad=0.05)# 
# label1 = "T"+r'$_{e}$'+" (eV)"    
# cb1 = mpl.colorbar.ColorbarBase(ax_cb, cmap=mpl.cm.plasma, orientation='vertical',norm=norm,label=label1)
# plt.xlabel("R (m)")
# plt.ylabel("Z (m)")
# plt.xlim([1.51,1.87])
# plt.ylim([-1.61,-1.1])
# plt.gcf().add_axes(ax_cb)
# fname = "te2d.png"
# plt.savefig("images//"+fname,dpi=1000)
# plt.show()



