from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import os
import matplotlib.image as m
import imageio as io
from scipy import interpolate
from PIL import Image
import glob
from scipy.integrate import trapz,cumtrapz
import re
import sys
import matplotlib as mpl
import matplotlib.cm as cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from plot2D import plotWALL, ImportGridue



def plotPatch(dataue,mode="",geometry="XPT"):
    fig, axs = plt.subplots(1, 1)
    LineColor = "#3E434A"
    LineWidth = 1
    xptColor = "black"
    if dataue["iyseparatrix1"] >dataue["iyseparatrix3"]:
        Greater = dataue["iyseparatrix1"]
        Less = dataue["iyseparatrix3"]
        dataue["iyseparatrix1"] = Less
        dataue["iyseparatrix3"] = Greater
    if mode=="rectangular":
        plt.plot(range(len(dataue["rm"])),np.zeros(len(dataue["rm"]))+dataue["iyseparatrix1"]+1,color=LineColor,linestyle="--",linewidth=LineWidth)
        plt.plot(range(len(dataue["rm"])),np.zeros(len(dataue["rm"]))+dataue["iyseparatrix3"]+1,color=LineColor,linestyle="-.",linewidth=LineWidth)
        plt.plot(np.zeros(dataue["iyseparatrix1"]+2)+dataue["ix_cut1"]+1,range(dataue["iyseparatrix1"]+2),color=xptColor,linestyle="-",linewidth=2*LineWidth)
        plt.plot(np.zeros(dataue["iyseparatrix2"]+2)+dataue["ix_cut2"]+1,range(dataue["iyseparatrix2"]+2),color=xptColor,linestyle="-",linewidth=2*LineWidth)
        plt.plot(np.zeros(dataue["iyseparatrix3"]+2)+dataue["ix_cut3"]+1,range(dataue["iyseparatrix3"]+2),color=xptColor,linestyle="-",linewidth=2*LineWidth)
        plt.plot(np.zeros(dataue["iyseparatrix4"]+2)+dataue["ix_cut4"]+1,range(dataue["iyseparatrix4"]+2),color=xptColor,linestyle="-",linewidth=2*LineWidth)
        plt.plot(np.zeros(len(dataue["rm"][0])+1)+dataue["ix_plate3"]+1,range(len(dataue["rm"][0])+1),color=xptColor,linestyle="-",linewidth=2*LineWidth)
        plt.plot(np.zeros(len(dataue["rm"][0])+1)+dataue["ix_plate4"]+2,range(len(dataue["rm"][0])+1),color=xptColor,linestyle="-",linewidth=2*LineWidth)
        plt.plot(np.zeros(len(dataue["rm"][0])+1),range(len(dataue["rm"][0])+1),color=xptColor,linestyle="-",linewidth=2*LineWidth)
        if geometry=="XPT":
            plt.text(0,-1.6,'T'+r'$_{1}$',horizontalalignment='center')
            plt.text(dataue["ix_cut1"]+1,-1.6,'XP'+r'$_{1}$',horizontalalignment='center')
            plt.text(dataue["ix_cut1"]+1,-3,'Left',horizontalalignment='center')
            plt.text(dataue["ix_cut2"]+1,-1.6,'XP'+r'$_{1}$',horizontalalignment='center')
            plt.text(dataue["ix_cut2"]+1,-3,'Right',horizontalalignment='center')
            plt.text(dataue["ix_cut3"]-2,-1.6,'XP'+r'$_{2}$',horizontalalignment='center')
            plt.text(dataue["ix_cut3"]-2,-3,'Left',horizontalalignment='center')
            plt.text(dataue["ix_cut4"]+3,-1.6,'XP'+r'$_{2}$',horizontalalignment='center')
            plt.text(dataue["ix_cut4"]+3,-3,'Right',horizontalalignment='center')
            plt.text(dataue["ix_plate3"]+1,-1.6,'T'r'$_{2}$'+"/T"+r'$_{3}$',horizontalalignment='center')
            plt.text(dataue["ix_plate4"]+4,-1.6,'T'+r'$_{4}$',horizontalalignment='center')
        else:
            plt.text(0,-1.6,'T'+r'$_{1}$',horizontalalignment='center')
            plt.text(dataue["ix_cut1"]+1,-1.6,'XP'+r'$_{1}$',horizontalalignment='center')
            plt.text(dataue["ix_cut1"]+1,-3,'Left',horizontalalignment='center')
            plt.text(dataue["ix_cut2"]+1,-1.6,'XP'+r'$_{2}$',horizontalalignment='center')
            plt.text(dataue["ix_cut2"]+1,-3,'Left',horizontalalignment='center')
            plt.text(dataue["ix_cut3"]+1,-1.6,'XP'+r'$_{2}$',horizontalalignment='center')
            plt.text(dataue["ix_cut3"]+1,-3,'Right',horizontalalignment='center')
            plt.text(dataue["ix_cut4"]+1,-1.6,'XP'+r'$_{1}$',horizontalalignment='center')
            plt.text(dataue["ix_cut4"]+1,-3,'Right',horizontalalignment='center')
            plt.text(dataue["ix_plate3"]+1,-1.6,'T'r'$_{2}$'+"/T"+r'$_{3}$',horizontalalignment='center')
            plt.text(dataue["ix_plate4"]+1,-1.6,'T'+r'$_{4}$',horizontalalignment='center')
    
            
    elif geometry=="XPT":
        axs.set_aspect('equal')       
        # plot seperatrix 1
        plt.plot(dataue['rm'][:dataue["ix_cut3"]+1,dataue["iyseparatrix1"],4],dataue['zm'][:dataue["ix_cut3"]+1,dataue["iyseparatrix1"],4],color=LineColor,linestyle="--",linewidth=LineWidth)
        plt.plot(dataue['rm'][dataue["ix_cut3"]+1:dataue["ix_plate3"],dataue["iyseparatrix1"],4],dataue['zm'][dataue["ix_cut3"]+1:dataue["ix_plate3"],dataue["iyseparatrix1"],4],color=LineColor,linestyle="--",linewidth=LineWidth)
        plt.plot(dataue['rm'][dataue["ix_plate3"]:dataue["ix_cut4"]+1,dataue["iyseparatrix1"],4],dataue['zm'][dataue["ix_plate3"]:dataue["ix_cut4"]+1,dataue["iyseparatrix1"],4],color=LineColor,linestyle="--",linewidth=LineWidth)
        plt.plot(dataue['rm'][dataue["ix_cut4"]+1:,dataue["iyseparatrix1"],4],dataue['zm'][dataue["ix_cut4"]+1:,dataue["iyseparatrix1"],4],color=LineColor,linestyle="--",linewidth=LineWidth)
        # # plot seperatrix 2
        plt.plot(dataue['rm'][:dataue["ix_cut3"]+1,dataue["iyseparatrix3"],4],dataue['zm'][:dataue["ix_cut3"]+1,dataue["iyseparatrix3"],4],color=LineColor,linestyle="-.",linewidth=LineWidth)
        plt.plot(dataue['rm'][dataue["ix_cut3"]+1:dataue["ix_plate3"],dataue["iyseparatrix3"],4],dataue['zm'][dataue["ix_cut3"]+1:dataue["ix_plate3"],dataue["iyseparatrix3"],4],color=LineColor,linestyle="-.",linewidth=LineWidth)
        plt.plot(dataue['rm'][dataue["ix_plate3"]:dataue["ix_cut4"]+1,dataue["iyseparatrix3"],4],dataue['zm'][dataue["ix_plate3"]:dataue["ix_cut4"]+1,dataue["iyseparatrix3"],4],color=LineColor,linestyle="-.",linewidth=LineWidth)
        plt.plot(dataue['rm'][dataue["ix_cut4"]+1:,dataue["iyseparatrix3"],4],dataue['zm'][dataue["ix_cut4"]+1:,dataue["iyseparatrix3"],4],color=LineColor,linestyle="-.",linewidth=LineWidth)
        # #plot x point 1
        plt.plot(dataue['rm'][0,:,4],dataue['zm'][0,:,4],color=xptColor,linestyle="-",linewidth=2*LineWidth)
        plt.plot(dataue['rm'][-1,:,4],dataue['zm'][-1,:,4],color=xptColor,linestyle="-",linewidth=2*LineWidth)
        plt.plot(dataue['rm'][dataue["ix_plate3"],:,4],dataue['zm'][dataue["ix_plate3"],:,4],color=xptColor,linestyle="-",linewidth=2*LineWidth)
        plt.plot(dataue['rm'][dataue["ix_plate2"],:,4],dataue['zm'][dataue["ix_plate2"],:,4],color=xptColor,linestyle="-",linewidth=2*LineWidth)
        plt.plot(dataue['rm'][dataue["ix_cut1"],:dataue["iyseparatrix1"]+1,4],dataue['zm'][dataue["ix_cut1"],:dataue["iyseparatrix1"]+1,4],color=xptColor,linestyle="-",linewidth=2*LineWidth)
        plt.plot(dataue['rm'][dataue["ix_cut2"],:dataue["iyseparatrix2"]+1,4],dataue['zm'][dataue["ix_cut2"],:dataue["iyseparatrix2"]+1,4],color=xptColor,linestyle="-",linewidth=2*LineWidth)
        plt.plot(dataue['rm'][dataue["ix_cut3"],:dataue["iyseparatrix3"]+1,4],dataue['zm'][dataue["ix_cut3"],:dataue["iyseparatrix3"]+1,4],color=xptColor,linestyle="-",linewidth=2*LineWidth)
        plt.plot(dataue['rm'][dataue["ix_cut4"],:dataue["iyseparatrix4"]+1,4],dataue['zm'][dataue["ix_cut4"],:dataue["iyseparatrix4"]+1,4],color=xptColor,linestyle="-",linewidth=2*LineWidth)
    
    elif geometry == "DDN":

        axs.set_aspect('equal')
        # plot seperatrix 1
        plt.plot(dataue['rm'][:dataue["ix_cut2"]+1,dataue["iyseparatrix1"],4],dataue['zm'][:dataue["ix_cut2"]+1,dataue["iyseparatrix1"],4],color=LineColor,linestyle="--",linewidth=LineWidth)
        plt.plot(dataue['rm'][dataue["ix_cut2"]+1:dataue["ix_plate3"],dataue["iyseparatrix1"],4],dataue['zm'][dataue["ix_cut2"]+1:dataue["ix_plate3"],dataue["iyseparatrix1"],4],color=LineColor,linestyle="--",linewidth=LineWidth)
        plt.plot(dataue['rm'][dataue["ix_plate3"]:dataue["ix_cut3"]+1,dataue["iyseparatrix1"],4],dataue['zm'][dataue["ix_plate3"]:dataue["ix_cut3"]+1,dataue["iyseparatrix1"],4],color=LineColor,linestyle="--",linewidth=LineWidth)
        plt.plot(dataue['rm'][dataue["ix_cut3"]+1:,dataue["iyseparatrix1"],3],dataue['zm'][dataue["ix_cut3"]+1:,dataue["iyseparatrix1"],3],color=LineColor,linestyle="--",linewidth=LineWidth)
        
         # # plot seperatrix 2
        plt.plot(dataue['rm'][:dataue["ix_plate3"],dataue["iyseparatrix3"],4],dataue['zm'][:dataue["ix_plate3"],dataue["iyseparatrix3"],4],color=LineColor,linestyle="-.",linewidth=LineWidth)
        plt.plot(dataue['rm'][dataue["ix_plate3"]+1:,dataue["iyseparatrix3"],4],dataue['zm'][dataue["ix_plate3"]+1:,dataue["iyseparatrix3"],4],color=LineColor,linestyle="-.",linewidth=LineWidth)
        # #plot x point 1
        plt.plot(dataue['rm'][dataue["ix_cut1"],:dataue["iyseparatrix1"]+1,4],dataue['zm'][dataue["ix_cut1"],:dataue["iyseparatrix1"]+1,4],color=xptColor,linestyle="-",linewidth=2*LineWidth)
        plt.plot(dataue['rm'][dataue["ix_cut2"],:dataue["iyseparatrix2"]+1,4],dataue['zm'][dataue["ix_cut2"],:dataue["iyseparatrix2"]+1,4],color=xptColor,linestyle="-",linewidth=2*LineWidth)
        plt.plot(dataue['rm'][dataue["ix_cut3"],:dataue["iyseparatrix3"]+1,4],dataue['zm'][dataue["ix_cut3"],:dataue["iyseparatrix3"]+1,4],color=xptColor,linestyle="-",linewidth=2*LineWidth)
        plt.plot(dataue['rm'][dataue["ix_cut4"],:dataue["iyseparatrix4"]+1,4],dataue['zm'][dataue["ix_cut4"],:dataue["iyseparatrix4"]+1,4],color=xptColor,linestyle="-",linewidth=2*LineWidth)
        plt.plot(dataue['rm'][0,:,4],dataue['zm'][0,:,4],color=xptColor,linestyle="-",linewidth=2*LineWidth)
        plt.plot(dataue['rm'][-1,:,4],dataue['zm'][-1,:,4],color=xptColor,linestyle="-",linewidth=2*LineWidth)
        plt.plot(dataue['rm'][dataue["ix_plate3"],:,4],dataue['zm'][dataue["ix_plate3"],:,4],color=xptColor,linestyle="-",linewidth=2*LineWidth)
        plt.plot(dataue['rm'][dataue["ix_plate2"],:,4],dataue['zm'][dataue["ix_plate2"],:,4],color=xptColor,linestyle="-",linewidth=2*LineWidth)
        # plt.plot(dataue['rm'][-1,:dataue["iyseparatrix1"]+1,4],dataue['zm'][dataue["ix_cut1"],:dataue["iyseparatrix1"]+1,4],color=xptColor,linestyle="-",linewidth=2*LineWidth)
    else:
        axs.set_aspect('equal')       
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
    
        plt.plot(dataue['rm'][0,:,4],dataue['zm'][0,:,4],color=xptColor,linestyle="-",linewidth=2*LineWidth)
        plt.plot(dataue['rm'][-1,:,4],dataue['zm'][-1,:,4],color=xptColor,linestyle="-",linewidth=2*LineWidth)
        plt.plot(dataue['rm'][dataue["ix_plate3"],:,4],dataue['zm'][dataue["ix_plate3"],:,4],color=xptColor,linestyle="-",linewidth=2*LineWidth)
        plt.plot(dataue['rm'][dataue["ix_plate2"],:,4],dataue['zm'][dataue["ix_plate2"],:,4],color=xptColor,linestyle="-",linewidth=2*LineWidth)
    for j in range(0,len(dataue["rm"])):
        for i in range(0,len(dataue["rm"][0])):
            x, y = 0,0
            if mode=="rectangular":
                x = [j,j,j+1,j+1,j]
                y =[i,i+1,i+1,i,i]
            else:
                x = [dataue["rm"][j,i,1],dataue["rm"][j,i,3],dataue["rm"][j,i,4],dataue["rm"][j,i,2],dataue["rm"][j,i,1]]
                y = [dataue["zm"][j,i,1],dataue["zm"][j,i,3],dataue["zm"][j,i,4],dataue["zm"][j,i,2],dataue["zm"][j,i,1]]
            
            if j<=dataue["ix_cut1"]:
                if i<=dataue["iyseparatrix1"]:
                    axs.fill(x,y,color="#F78BCC",linewidth=0.0) 
                elif i<=dataue["iyseparatrix3"]:
                    axs.fill(x,y,color="#CD70D4",linewidth=0.0) 
                else:
                    axs.fill(x,y,color="#A92BAD",linewidth=0.0)   
            elif j<=dataue["ix_cut2"]:
                            
                if i<=dataue["iyseparatrix1"]:
                    axs.fill(x,y,color="#19ABB5",linewidth=0.0) 
                elif i<=dataue["iyseparatrix3"]:
                    axs.fill(x,y,color="#1D6996",linewidth=0.0) 
                else:
                    axs.fill(x,y,color="#20527A",linewidth=0.0) 
           
            elif j<=dataue["ix_cut3"] and (geometry=="XPT"or  geometry=="Snowflake"):
                if i<=dataue["iyseparatrix1"]:
                    axs.fill(x,y,color="#73AF48",linewidth=0.0) 
                elif i<=dataue["iyseparatrix3"]:
                    axs.fill(x,y,color="#3F6E3D",linewidth=0.0) 
                else:
                    axs.fill(x,y,color="#21592C",linewidth=0.0)
            elif j<=dataue["ix_plate3"] and geometry=="DDN":
                if i<=dataue["iyseparatrix1"]:
                    axs.fill(x,y,color="#73AF48",linewidth=0.0) 
                elif i<=dataue["iyseparatrix3"]:
                    axs.fill(x,y,color="#3F6E3D",linewidth=0.0) 
                else:
                    axs.fill(x,y,color="#21592C",linewidth=0.0)
                       
            elif j<=dataue["ix_plate3"] and (geometry=="XPT"or  geometry=="Snowflake"):
                if i<=dataue["iyseparatrix1"]:
                    axs.fill(x,y,color="#F5E084",linewidth=0.0) 
                elif i<=dataue["iyseparatrix3"]:
                    axs.fill(x,y,color="#EDAD08",linewidth=0.0) 
                else:
                    axs.fill(x,y,color="#CC7D1D",linewidth=0.0)
            
            elif j<=dataue["ix_cut3"] and geometry=="DDN":
                if i<=dataue["iyseparatrix1"]:
                    axs.fill(x,y,color="#F5E084",linewidth=0.0) 
                elif i<=dataue["iyseparatrix3"]:
                    axs.fill(x,y,color="#EDAD08",linewidth=0.0) 
                else:
                    axs.fill(x,y,color="#CC7D1D",linewidth=0.0)
            
            elif j<=dataue["ix_cut4"]:
                if i<=dataue["iyseparatrix1"]:
                    axs.fill(x,y,color="#CC503E",linewidth=0.0) 
                elif i<=dataue["iyseparatrix3"]:
                    axs.fill(x,y,color="#B01010",linewidth=0.0) 
                else:
                    axs.fill(x,y,color="#991717",linewidth=0.0)
            else:
                if i<=dataue["iyseparatrix1"]:
                    axs.fill(x,y,color="#87377B",linewidth=0.0) 
                elif i<=dataue["iyseparatrix3"]:
                    axs.fill(x,y,color="#6F4070",linewidth=0.0) 
                else:
                    axs.fill(x,y,color="#3E0E45",linewidth=0.0)   
    fontsize = 9
    if mode=="rectangular":
        plt.xlabel("ix")
        plt.ylabel("iy")
        if geometry=="XPT":
            plt.ylim([-4,35])
            plt.xlim([-4,189])
        else:
            plt.ylim([-4,33])
            plt.xlim([-2,102])     
    else:
        axs.set_xticks([])
        axs.set_yticks([])
        # plt.xlabel("R (m)")
        # plt.ylabel("Z (m)")
        if geometry=="XPT":
            plt.annotate("XP"+r"$_{1}$",(1.48,-1.25),fontsize=fontsize,color="black")
            plt.annotate("XP"+r"$_{2}$",(1.63,-1.51),fontsize=fontsize,color="black")
            plt.annotate("T"+r"$_{3}$",(1.82,-1.61),fontsize=fontsize,color="black")
            plt.annotate("T"+r"$_{2}$",(1.82,-1.4),fontsize=fontsize,color="black")
            plt.annotate("T"+r"$_{1}$",(1.34,-1.12),fontsize=fontsize,color="black")
            plt.annotate("T"+r"$_{4}$",(1.7,-1.59),fontsize=fontsize,color="black")
            plt.ylim([-1.62,-0.48])
            plt.xlim([1.31,1.88])

        elif geometry =="DDN":
            plt.annotate("XP"+r"$_{1}$",(1.5,-0.7),fontsize=fontsize,color="black")
            plt.annotate("XP"+r"$_{2}$",(1.5,0.54),fontsize=fontsize,color="black")
            plt.annotate("T"+r"$_{3}$",(1.49,1.39),fontsize=fontsize,color="black")
            plt.annotate("T"+r"$_{2}$",(1.07,1.39),fontsize=fontsize,color="black")
            plt.annotate("T"+r"$_{1}$",(1.02,-1.42),fontsize=fontsize,color="black")
            plt.annotate("T"+r"$_{4}$",(1.63,-1.43),fontsize=fontsize,color="black")
            plt.ylim([-1.5,1.5])
            plt.xlim([0.98,2.4])
    plt.tight_layout()
    fname = "patch"+mode+geometry+".png"
    plt.savefig("images/"+fname,dpi=1000,bbox_inches='tight')
# folder0 = "D:\\my stuff\\PhD\\collaboratory\\XPT_Building\\SPARCSNXPT_0-5mmSeparation\\FineGrid"
folder0 = "D:\\my stuff\\PhD\\collaboratory\\Paper1\\patchData"

# dataue = ImportGridue(folder0+"\\gridue")
# plotPatch(dataue,mode="normal",geometry="XPT")
# dataue = ImportGridue(folder0+"\\gridueSnowflake")
# plotPatch(dataue,mode="normal",geometry="Snowflake")
dataue = ImportGridue(folder0+"\\gridueDDN")
plotPatch(dataue,mode="normal",geometry="DDN")
# dataue = ImportGridue(folder0+"\\gridueDDN")
# plotPatch(dataue,mode="rectangular",geometry="DDN")