from typing import List
# from plot2D import plotWALL
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
from SharedFunctions import return2d
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
plt.rcParams["font.family"] = "serif"
params = {'legend.fontsize': 'medium',
         'axes.labelsize': 'medium',
         'axes.titlesize':'medium',
         'xtick.labelsize':'medium',
         'ytick.labelsize':'medium',
         'figure.figsize': (8.3,5),
         }
plt.rcParams.update(params)

fig, axs = plt.subplots(1, 4)
for i in range(4):
    axs[i].set_aspect('equal')

colors = ["#0AC48D","#00CFDE","#00A838","00A838"]
colors = ["#004C4C","#00CDC8","#006666","#009E9A"]
colors = ["C0","C0","C0","C0"]
targetColor = "C2"
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


def orderData(data):
    datanew = data[0:2]
    data = np.delete(data,range(0,2),axis=0)
    for i in range(len(data)-1):

        for j in range(0,len(data),2):
            if j>=len(data):
                break
            if datanew[0][0] ==data[j+1][0] and datanew[0][1] ==data[j+1][1]:
                datanew= np.append(data[j:j+2],datanew,axis=0)
                data = np.delete(data,range(j,j+2),axis=0)
            
            elif datanew[-1][0] ==data[j][0] and datanew[-1][1] ==data[j][1]:
                datanew= np.append(datanew,data[j:j+2],axis=0)
                data = np.delete(data,range(j,j+2),axis=0)

    return datanew


def plotContours():
    data1 = np.loadtxt(Dir+"\\contour1.ogr",delimiter=" ",skiprows=1)
    data1 = np.array(data1)
    datanew = orderData(data1)
    data1 = datanew.copy()/1000
    data2 = np.loadtxt(Dir+"\\contour2.ogr",delimiter=" ",skiprows=1)
    data2 = np.array(data2)
    datanew = orderData(data2)
    data2 = datanew.copy()/1000
    data3 = np.loadtxt(Dir+"\\contour3.ogr",delimiter=" ",skiprows=1)
    data3 = np.array(data3)
    datanew = orderData(data3)
    data3 = datanew.copy()/1000
    data4 = np.loadtxt(Dir+"\\contour4.ogr",delimiter=" ",skiprows=1)
    data4 = np.array(data4)
    datanew = orderData(data4)
    data4 = datanew.copy()/1000

    for i in range(4):
        axs[i].plot([data2[0,0],data3[-1,0]],[data2[0,1],data3[-1,1]],color=targetColor,linewidth=1)
        axs[i].plot([data1[0,0],data2[-1,0]],[data1[0,1],data2[-1,1]],color=targetColor,linewidth=1)
        axs[i].plot([data3[0,0],data4[-1,0]],[data3[0,1],data4[-1,1]],color=targetColor,linewidth=1)
        axs[i].plot([data4[0,0],data1[-1,0]],[data4[0,1],data1[-1,1]],color=targetColor,linewidth=1)
        axs[i].plot(data1[:,0],data1[:,1],color=colors[0])
        axs[i].plot(data2[:,0],data2[:,1],color=colors[1])
        axs[i].plot(data3[:,0],data3[:,1],color=colors[2])
        axs[i].plot(data4[:,0],data4[:,1],color=colors[3])

    # for i in range(2,4):
    #     axs[i].plot([data2[0,0],data3[-1,0]],[data2[0,1],data3[-1,1]],color=targetColor,linewidth=1)
    #     axs[i].plot([data1[0,0],data2[-1,0]],[data1[0,1],data2[-1,1]],color=targetColor,linewidth=1)
    #     axs[i].plot([data3[0,0],data4[-1,0]],[data3[0,1],data4[-1,1]],color=targetColor,linewidth=1)
    #     axs[i].plot([data4[0,0],data1[-1,0]],[data4[0,1],data1[-1,1]],color=targetColor,linewidth=1)
    #     axs[i].plot(data1[:,0],data1[:,1],color="black",linewidth=0.15)
    #     axs[i].plot(data2[:,0],data2[:,1],color="black",linewidth=0.15)
    #     axs[i].plot(data3[:,0],data3[:,1],color="black",linewidth=0.15)
    #     axs[i].plot(data4[:,0],data4[:,1],color="black",linewidth=0.15)


# plotWALL("D:\\my stuff\\PhD\\collaboratory\\XPT_Building\\templates\\input.dat.template",axs)



folder0 = "D:\\my stuff\\PhD\\collaboratory\\XPT_Building\\SPARCSNXPT_0-5mmSeparation\\FineGrid"
data = ImportGridue(folder0+"\\gridue")

for i in range(len(data["rm"])):
    for j in range(len(data["rm"][0])):
        r = [data["rm"][i][j][1],data["rm"][i][j][3],data["rm"][i][j][4],data["rm"][i][j][2]]
        z = [data["zm"][i][j][1],data["zm"][i][j][3],data["zm"][i][j][4],data["zm"][i][j][2]]
        axs[1].plot(r,z,linewidth=0.15,color="black")
        # plt.plot(data["rm"][i][j][1],data["zm"][i][j][1],linestyle="",marker="o")


file0 = "balFiles//balanceXPTTest.nc"

#unpack two balance files
quantities2d0,rootgrp0 = return2d(file0)


#define midplane x point indices


plotquantity = np.log10(quantities2d0["te"])
norm = mpl.colors.Normalize(vmin=np.amin(plotquantity), vmax=np.amax(plotquantity))
cmap = cm.plasma
m = cm.ScalarMappable(norm=norm, cmap=cmap)
for j in range(0,len(rootgrp0['crx'][0][0])):#len(R[0])-1):
    for i in range(0,len(rootgrp0['crx'][0])):
        x = [rootgrp0['crx'][0,i,j],rootgrp0['crx'][2,i,j],rootgrp0['crx'][3,i,j],rootgrp0['crx'][1,i,j],rootgrp0['crx'][0,i,j]]
        y =[rootgrp0['cry'][0,i,j],rootgrp0['cry'][2,i,j],rootgrp0['cry'][3,i,j],rootgrp0['cry'][1,i,j],rootgrp0['cry'][0,i,j]]
        color1 = m.to_rgba(plotquantity[i,j])
        axs[3].fill(x,y,color=color1,linewidth=0.01)  
label1 = "log (T"+r'$_{e}$'+")"  
# divider = make_axes_locatable(axs[3])
# ax_cb = divider.append_axes("right", size="5%", pad=-0.1) 
cbarax = plt.axes([0.92, 0.1, 0.01, 0.2], facecolor='none')
cb1 = mpl.colorbar.ColorbarBase(cbarax, cmap=mpl.cm.plasma, orientation='vertical',norm=norm,label=label1)

fd = open('balFiles\\fort.33','r')    
d = np.loadtxt(fd,skiprows=1)
data = d.flatten()
datax = data[0:int(len(data)/2)]/100
datay = data[int(len(data)/2):]/100
print(len(datay))
print(len(datax))
fd.close()

fd = open('balFiles\\fort.34','r')   
indices = np.loadtxt(fd,skiprows=1,usecols=(1,2,3))
# data = np.array([d["col1"],d["col2"],d["col3"],d["col4"]])
# print(d["col2"])
for i in range(len(indices)):
    ind0 = int(indices[i][0])-1
    ind1 = int(indices[i][1])-1
    ind2 = int(indices[i][2])-1
    # print(ind0,ind1,ind2)
    axs[2].plot([datax[ind0],datax[ind1],datax[ind2]],[datay[ind0],datay[ind1],datay[ind2]],color="black",linewidth=0.15)

Dir = "D:\\my stuff\\PhD\\collaboratory\\XPT_Building\\SPARC_LSNX3\\"
fileIn = Dir+"\\contour1.ogr"
colornum = 1
plotContours()


axs[0].set_xlabel("R (m)")
axs[0].set_ylabel("Z (m)")
fontsize = 7
axs[0].annotate("contour 1",(2,0.96),fontsize=fontsize,color=colors[0])
axs[0].annotate("contour 2",(1.9,-1.52),fontsize=fontsize,color=colors[1])
axs[0].annotate("contour 3",(1.45,-1.67),fontsize=fontsize,color=colors[2])
axs[0].annotate("contour 4",(1.22,-1.38),fontsize=fontsize,color=colors[3])
axs[0].annotate("targets",(1.4,-1.52),fontsize=fontsize,color=targetColor)
axs[0].annotate("a)",(2.1,1.15),fontsize=12,color="black")
axs[1].annotate("b)",(2.1,1.15),fontsize=12,color="black")
axs[2].annotate("c)",(2.1,1.15),fontsize=12,color="black")
axs[3].annotate("d)",(2.1,1.15),fontsize=12,color="black")
for i in range(1,4):
    axs[i].set_yticklabels([])
    axs[i].set_xticklabels([])
# chartBox = axs[4].get_position()
# axs[4].set_position([-chartBox.x0, chartBox.y0,
#                  chartBox.width*0.1,
#                  chartBox.height * 0.1])
  
plt.subplots_adjust(wspace=0, hspace=0)
plt.tight_layout(pad=0)
plt.savefig("images/Process.png",dpi=800, bbox_inches='tight')
plt.show()

