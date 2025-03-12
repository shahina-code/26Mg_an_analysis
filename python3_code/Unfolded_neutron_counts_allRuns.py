import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from scipy.integrate import simps
from numpy import trapz

from scipy.optimize import curve_fit
from scipy.stats import rayleigh



from scipy.optimize import curve_fit
from scipy import integrate
import matplotlib.ticker as mtick
from matplotlib.ticker import FormatStrFormatter
import matplotlib as mpl
from scipy import integrate


def gaussian1(x,  norm1,mean1,sigma1):
    return (norm1/(sigma1*np.sqrt(2*np.pi)))*np.exp(-np.power(x - mean1, 2)/(2*np.power(sigma1, 2)))
    
def gaussian2(x,  norm2,mean2,sigma2):
    return (norm2/(sigma2*np.sqrt(2*np.pi)))*np.exp(-np.power(x - mean2, 2)/(2*np.power(sigma2, 2)))

def fitFunc(x, norm1,mean1,sigma1,norm2,mean2,sigma2):
    return  gaussian1(x, norm1,mean1,sigma1) +  gaussian2(x,  norm2,mean2,sigma2)
    
    
total_counts_25Mg_n1=[]
total_counts_25Mg_n0=[]
total_counts_13C = []
total_counts_11B = []
total_counts_17O_n0=[]
total_counts_17O_n1 =[]


data_run = np.loadtxt("../25Mg_a_n_columns.txt")
runNumber = data_run[:,0][0:108]
charge = data_run[:,1][0:108]
beamEnergy = data_run[:,2][0:108]


for i in range(len(runNumber)):
    data_unfolded = np.loadtxt("../Unfolded_output_files/%d.txt" %(runNumber[i]))
    bins_unfolded = data_unfolded[:,0]
    counts_unfolded = data_unfolded[:,1]
    
    bins_unfolded = np.array(bins_unfolded)
    counts_unfolded = np.array(counts_unfolded)
    
    
    neutron_data = pd.read_csv("../Unfolded_count_files/%d.txt" %(runNumber[i]),header=None,sep=" ")
    df = pd.DataFrame(neutron_data)
    df.columns = ["run","beam-energy","reaction","blah","neutron_energy"]
    neutron_energy = df.values[:,4]
    reactions = df.values[:,2]
    
    
    
    
    #Masking on the region of the 17O(a,n1) and 11B(a,n0) peak to fit it with a double gaussian
    E_min_11B = neutron_energy[8] -0.25
    E_max_11B = neutron_energy[12] + 0.4

    mask1_11B = bins_unfolded > E_min_11B
    bins_unfolded_mask1_11B = bins_unfolded[mask1_11B]
    counts_unfolded_mask1_11B = counts_unfolded[mask1_11B]
    
    mask2_11B = bins_unfolded_mask1_11B < E_max_11B
    bins_unfolded_masked_11B = bins_unfolded_mask1_11B[mask2_11B]
    counts_unfolded_masked_11B = counts_unfolded_mask1_11B[mask2_11B]
    
    mask3_11B = counts_unfolded_masked_11B >0
    bins_unfolded_masked_pos_11B = bins_unfolded_masked_11B[mask3_11B]
    counts_unfolded_masked_pos_11B= counts_unfolded_masked_11B[mask3_11B]
    total_counts_11B.append(np.sum(counts_unfolded_masked_pos_11B))
    
    # Fitting Parameters for fitting double gaussian to the  17O(a,n1) and 11B(a,n0)  peak
    
    norm1_min_11B=0.01
    norm1_max_11B=(np.max(counts_unfolded_masked_pos_11B)+1)*10
    mean1_min_11B= neutron_energy[8] -0.4
    mean1_max_11B= neutron_energy[8] + 0.4
    sigma1_min_11B=0.01
    sigma1_max_11B=0.3
    
    norm2_min_11B=0.01
    norm2_max_11B=np.max(counts_unfolded_masked_pos_11B)*10
    mean2_min_11B= neutron_energy[12] -0.4
    mean2_max_11B= neutron_energy[12] + 0.4
    sigma2_min_11B=0.01
    sigma2_max_11B=0.3
    
    """
    popt_11B, pcov_11B = curve_fit(fitFunc,bins_unfolded_masked_pos_11B,counts_unfolded_masked_pos_11B,method="trf",max_nfev=10000,bounds=([norm1_min_11B,mean1_min_11B,sigma1_min_11B,norm2_min_11B,mean2_min_11B,sigma2_min_11B],[norm1_max_11B,mean1_max_11B,sigma1_max_11B,norm2_max_11B,mean2_max_11B,sigma2_max_11B]))
    
    x_11B = np.linspace(E_min_11B, E_max_11B, 2000)
    y_11B = fitFunc(x_11B, *popt_11B)
    perr_11B =  np.sqrt(np.diag(pcov_11B))
   
   
    Counts_17O_n1 = popt_11B[0]
    total_counts_17O_n1.append(Counts_17O_n1)
    Counts_11B = popt_11B[3]
    """
    
    
    
    
    # Emin and Emax for Integrating the 25Mg(a,n1) peak
    E_min_25Mg_n1 = neutron_energy[2]-0.3 #2.43
    
    E_max_25Mg_n1 = neutron_energy[2]+0.4 #3.46


    mask1_25Mg_n1 = bins_unfolded > E_min_25Mg_n1
    bins_unfolded_mask1_25Mg_n1 = bins_unfolded[mask1_25Mg_n1]
    counts_unfolded_mask1_25Mg_n1 = counts_unfolded[mask1_25Mg_n1]
    
    mask2_25Mg_n1 = bins_unfolded_mask1_25Mg_n1 < E_max_25Mg_n1
    bins_unfolded_masked_25Mg_n1 = bins_unfolded_mask1_25Mg_n1[mask2_25Mg_n1]
    counts_unfolded_masked_25Mg_n1 = counts_unfolded_mask1_25Mg_n1[mask2_25Mg_n1]
    
    mask3_25Mg_n1 = counts_unfolded_masked_25Mg_n1 >0
    bins_unfolded_masked_pos_25Mg_n1 = bins_unfolded_masked_25Mg_n1[mask3_25Mg_n1]
    counts_unfolded_masked_pos_25Mg_n1= counts_unfolded_masked_25Mg_n1[mask3_25Mg_n1]
    total_counts_25Mg_n1.append(np.sum(counts_unfolded_masked_pos_25Mg_n1))
    
    
    # Fitting Parameters for fitting double gaussian to the 17O(a,n0) and 25Mg(a,n1) peak
    
    norm1_min_25Mg_n1=0.01
    norm1_max_25Mg_n1=np.max(counts_unfolded_masked_pos_25Mg_n1)*10
    mean1_min_25Mg_n1= neutron_energy[7] -0.4
    mean1_max_25Mg_n1= neutron_energy[7] + 0.4
    sigma1_min_25Mg_n1=0.01
    sigma1_max_25Mg_n1=0.3
    
    norm2_min_25Mg_n1=0.01
    norm2_max_25Mg_n1=np.max(counts_unfolded_masked_pos_25Mg_n1)*10
    mean2_min_25Mg_n1= neutron_energy[2] -0.4
    mean2_max_25Mg_n1= neutron_energy[2] + 0.4
    sigma2_min_25Mg_n1=0.01
    sigma2_max_25Mg_n1=0.3
    
    """
    popt_25Mg_n1, pcov_25Mg_n1 = curve_fit(fitFunc,bins_unfolded_masked_pos_25Mg_n1,counts_unfolded_masked_pos_25Mg_n1,method="trf",max_nfev=10000,bounds=([norm1_min_25Mg_n1,mean1_min_25Mg_n1,sigma1_min_25Mg_n1,norm2_min_25Mg_n1,mean2_min_25Mg_n1,sigma2_min_25Mg_n1],[norm1_max_25Mg_n1,mean1_max_25Mg_n1,sigma1_max_25Mg_n1,norm2_max_25Mg_n1,mean2_max_25Mg_n1,sigma2_max_25Mg_n1]))
    
    x_25Mg_n1 = np.linspace(E_min_25Mg_n1, E_max_25Mg_n1, 2000)
    y_25Mg_n1 = fitFunc(x_25Mg_n1, *popt_25Mg_n1)
    perr_25Mg_n1 =  np.sqrt(np.diag(pcov_25Mg_n1))
    
    Counts_17O_n0 = popt_25Mg_n1[0]
    total_counts_17O_n0.append(Counts_17O_n0)
    Counts_25Mg_n1 = popt_25Mg_n1[3]
    """
    
    
    
    #Masking on the region of the 13C(a,n) and 25Mg(a,n0) peak to fit it with a double gaussian
    E_min_13C = neutron_energy[0]-0.6 #2.43
    
    E_max_13C = neutron_energy[0]+0.3 #3.46


    mask1_13C = bins_unfolded > E_min_13C
    bins_unfolded_mask1_13C = bins_unfolded[mask1_13C]
    counts_unfolded_mask1_13C= counts_unfolded[mask1_13C]
    
    mask2_13C = bins_unfolded_mask1_13C < E_max_13C
    bins_unfolded_masked_13C = bins_unfolded_mask1_13C[mask2_13C]
    counts_unfolded_masked_13C = counts_unfolded_mask1_13C[mask2_13C]
    
    mask3_13C = counts_unfolded_masked_13C >0
    bins_unfolded_masked_pos_13C = bins_unfolded_masked_13C[mask3_13C]
    counts_unfolded_masked_pos_13C= counts_unfolded_masked_13C[mask3_13C]
    total_counts_13C.append(np.sum(counts_unfolded_masked_pos_13C))
    
    
    E_min_25Mg_n0 = neutron_energy[1]-0.3 #2.43
    
    E_max_25Mg_n0 = neutron_energy[1]+0.3 #3.46


    mask1_25Mg_n0 = bins_unfolded > E_min_25Mg_n0
    bins_unfolded_mask1_25Mg_n0 = bins_unfolded[mask1_25Mg_n0]
    counts_unfolded_mask1_25Mg_n0 = counts_unfolded[mask1_25Mg_n0]
    
    mask2_25Mg_n0 = bins_unfolded_mask1_25Mg_n0 < E_max_25Mg_n0
    bins_unfolded_masked_25Mg_n0 = bins_unfolded_mask1_25Mg_n0[mask2_25Mg_n0]
    counts_unfolded_masked_25Mg_n0 = counts_unfolded_mask1_25Mg_n0[mask2_25Mg_n0]
    
    mask3_25Mg_n0 = counts_unfolded_masked_25Mg_n0 >0
    bins_unfolded_masked_pos_25Mg_n0 = bins_unfolded_masked_25Mg_n0[mask3_25Mg_n0]
    counts_unfolded_masked_pos_25Mg_n0 = counts_unfolded_masked_25Mg_n0[mask3_25Mg_n0]
    total_counts_25Mg_n0.append(np.sum(counts_unfolded_masked_pos_25Mg_n0))

    
    # Fitting Parameters for fitting double gaussian to the 13C(a,n) and 25Mg(a,n0) peak
    """
    norm1_min=0.01
    norm1_max=np.max(counts_unfolded_masked_pos)*10
    mean1_min= neutron_energy[0] -0.4
    mean1_max= neutron_energy[0] + 0.4
    sigma1_min=0.01
    sigma1_max=0.3
    
    norm2_min=0.01
    norm2_max=np.max(counts_unfolded_masked_pos)*10
    mean2_min= neutron_energy[1] -0.4
    mean2_max= neutron_energy[1] + 0.4
    sigma2_min=0.01
    sigma2_max=0.3
    
    
    
    
    popt, pcov = curve_fit(fitFunc,bins_unfolded_masked_pos,counts_unfolded_masked_pos,method="trf",max_nfev=10000,bounds=([norm1_min,mean1_min,sigma1_min,norm2_min,mean2_min,sigma2_min],[norm1_max,mean1_max,sigma1_max,norm2_max,mean2_max,sigma2_max]))
    
    x = np.linspace(E_min_fit, E_max_fit, 2000)
    y = fitFunc(x, *popt)
    perr =  np.sqrt(np.diag(pcov))
    
    Counts_13C = popt[0]
    total_counts_13C.append(Counts_13C)
    Counts_25Mg_n0 = popt[3]
    """
    
    
   
    
    fig,ax=plt.subplots(1,1,figsize=(10,10))
    fig.set_dpi(150)
    fig.set_size_inches(10,6)
    #Axis ticks
    ax.tick_params(which='major', direction='in', width=2, length=10, color='k', pad=15)
    ax.tick_params(which='minor', direction='in', width=1.5, length=5, color='k', pad=15)
    ax.xaxis.set_minor_locator(mtick.AutoMinorLocator(5))
    ax.yaxis.set_minor_locator(mtick.AutoMinorLocator(5))
    ax.tick_params(labelsize=15)
    ax.patch.set_edgecolor('k')
    ax.patch.set_linewidth(2)
    
    plt.step(bins_unfolded, counts_unfolded)
    #plt.plot(x,y,label="Fit",linewidth=3)
    #plt.plot(x_25Mg_n1,y_25Mg_n1,label="Fit",linewidth=3)
    #plt.plot(x_11B,y_11B,label="Fit",linewidth=3)
    plt.xlim(0,8)
    
    plt.title(" Target: 25Mg, Run_%d , E_beam: %0.2f keV" %(runNumber[i],beamEnergy[i]),fontsize=15)
    
    """
    if (runNumber[i] <1711):
        
        plt.title(" Target: 25Mg, Run_%d , E_beam: %0.2f keV" %(runNumber[i],beamEnergy[i]),fontsize=15)
    else:
        plt.title(" Target: 24Mg, Run_%d , E_beam: %0.2f keV" %(runNumber[i],beamEnergy[i]),fontsize=15)
   """

    #plt.axvline(x=neutron_energy[0],color="blue",linewidth=2,linestyle="--",label="13C(a,n0) En: %0.3f MeV" %(neutron_energy[0]) )
    
    #plt.axvline(x=neutron_energy[4],color="blue",label="26Mg(a,n0) En: %0.3f MeV" %(neutron_energy[4]) )
    #plt.axvline(x=neutron_energy[5],color="green",label="26Mg(a,n1) En: %0.3f MeV" %(neutron_energy[5]) )
    #plt.axvline(x=neutron_energy[7],color="black",linewidth=2,linestyle="--",label="17O(a,n0) En: %0.3f MeV"%(neutron_energy[7]) )
    #plt.axvline(x=neutron_energy[8],color="brown",linewidth=2,linestyle="--",label="17O(a,n1) En: %0.3f MeV"%(neutron_energy[8]) )
    #plt.axvline(x=neutron_energy[9],color="yellow",linewidth=2,linestyle="--",label="18O(a,n0) En: %0.3f MeV"%(neutron_energy[9]) )
    #plt.axvline(x=neutron_energy[11],color="red",linewidth=2,linestyle="--",label="10B(a,n0) En: %0.3f MeV"%(neutron_energy[11]) )
    #plt.axvline(x=neutron_energy[12],color="deeppink",linewidth=2,linestyle="--",label="11B(a,n0) En: %0.3f MeV"%(neutron_energy[12]) )
    
    #plt.axvline(x=neutron_energy[0],color="red",linewidth=2,linestyle="--",label="13C(a,n0) En: %0.3f MeV" %(neutron_energy[0]) )
    #plt.axvline(x=neutron_energy[1],color="purple",linewidth=2,linestyle="--",label="25Mg(a,n0) En: %0.3f MeV" %(neutron_energy[1]) )
    ##plt.axvline(x=neutron_energy[2],color="green",linewidth=2,linestyle="--",label="25Mg(a,n1) En: %0.3f MeV" %(neutron_energy[2]) )
    #plt.axvline(x=neutron_energy[4],color="black",label="17O(a,n0) En: %0.3f MeV"%(neutron_energy[4]) )
    #plt.axvline(x=neutron_energy[3],color="red",label="25Mg(a,n2) En: %0.3f MeV" %(neutron_energy[3]) )
    #plt.axvline(x=neutron_energy[8]+0.5,color="black",label="10B(a,n0) En: %0.3f MeV"%(neutron_energy[8]+0.5) )
    #plt.axvline(x=neutron_energy[9]+0.5,color="blue",label="11B(a,n0) En: %0.3f MeV"%(neutron_energy[9]+0.5) )
    #plt.axvline(x=neutron_energy[10],color="blue",label="7Li(p,n0) En: %0.3f MeV"%(neutron_energy[10]) )
    #plt.axvline(x=neutron_energy[11],color="blue",label="7Li(p,n1) En: %0.3f MeV"%(neutron_energy[11]) )
    
    

    plt.xlim(0,np.max(bins_unfolded))
    plt.ylim(0,1.2*np.max(counts_unfolded))

    
    plt.ylabel("Counts/50 keVee",fontsize=15)
    plt.xlabel("$E_{n}$ (MeV)",fontsize=15)
    
    #plt.show()

    
    
    #ax.fill_between(neutron_energy_masked_13C_1,counts_unfolded_masked_13C_1,step="pre",facecolor="blue",alpha=0.5,label="13C(a,n0) En: %0.3f MeV" %(neutron_energy[0]))
    #ax.fill_between(neutron_energy_25Mg_n1_masked_1, counts_unfolded_25Mg_n1_masked_1 , step="pre",facecolor='green', alpha=.5,label="25Mg(a,n1) En: %0.3f MeV" %(neutron_energy[2]))
    #ax.fill_between(neutron_energy_masked_17O_1,counts_unfolded_masked_17O_1,step="pre",facecolor="orange",alpha=0.5,label="17O(a,n0) En: %0.3f MeV"%(neutron_energy[7]))
    
    #plt.legend()
    plt.savefig("Unfolded_spectrums_OU_1/Unfolded_spectrum_run_%d.png" %(runNumber[i]),dpi=300)
    plt.close()





# Writing arrays to a file

"""
output_file = "counts_unfolding.txt"

with open(output_file, 'w') as file:
    for i in range(len(total_counts_25Mg_n1)):
        file.write("%d \t  %d \t %f \t %d \t %d \t %d \t %d \n" %(runNumber[i],charge[i],beamEnergy[i],total_counts_25Mg_n1[i],total_counts_25Mg_n0[i],total_counts_13C[i],total_counts_11B[i]))
file.close()
"""

