
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
from __future__ import print_function
import seaborn as sns
from matplotlib.ticker import NullFormatter, MaxNLocator
import matplotlib.ticker as ticker
import plotly.graph_objects as go
import scipy as sp
from itertools import chain
import matplotlib as mpl
from matplotlib.lines import Line2D
from scipy import spatial
from scipy.spatial import ConvexHull
from scipy.optimize import curve_fit
from matplotlib import path
from scipy.stats import probplot,shapiro, sem
from scipy.interpolate import make_interp_spline
from mpl_toolkits.mplot3d.axes3d import Axes3D
from matplotlib.ticker import (AutoMinorLocator, MultipleLocator)
from matplotlib import cm
from numpy import linspace
import pylab
import os
from scipy.ndimage import gaussian_filter


# # initialization steps

# In[2]:


#protamine details
#need to remove frames and simulations that are not necessary
#only show equilibrated part
#this is for salmon
salmon_rg = pd.read_csv("~/protamine_simulation/jupyter_nb/salmon_rg.csv").iloc[:,:3]
salmon_ree = pd.read_csv("~/protamine_simulation/jupyter_nb/salmon_ree.csv")
salmon_rg_ree=pd.concat([salmon_rg,salmon_ree['Ree']],axis=1)
salmon_rg_ree=salmon_rg_ree[salmon_rg_ree.sim!=2]
salmon_rg_ree=salmon_rg_ree[salmon_rg_ree.sim!=11]
#this step is because some salmon simulations have more than 200ns
temp_df=salmon_rg_ree.copy()
for sim in temp_df.sim.unique():
    n=temp_df[temp_df.sim==sim].frames.max()-100000
    if n>0:
        temp_df.drop(temp_df[temp_df.sim==sim].tail(n).index,inplace=True)
salmon_rg_ree=temp_df.copy()
del temp_df
#this step is taking last 40000 frames i.e. equilibrated portion
blank_df = pd.DataFrame()
for sim in salmon_rg_ree.sim.unique():
    temp_df = salmon_rg_ree[salmon_rg_ree.sim==sim].iloc[-40000:,:]
    blank_df = blank_df.append(temp_df)
salmon_rg_ree = blank_df.copy()
del blank_df
del temp_df
del salmon_rg, salmon_ree

#need to remove frames and simulations that are not necessary
#only show equilibrated part
#this is for p1
p1_rg = pd.read_csv("~/protamine_simulation/jupyter_nb/p1_rg.csv")
p1_ree = pd.read_csv("~/protamine_simulation/jupyter_nb/p1_ree.csv")
p1_rg_ree=pd.concat([p1_rg,p1_ree['Ree']],axis=1)
p1_rg_ree=p1_rg_ree[p1_rg_ree.sim!=15]
p1_rg_ree=p1_rg_ree[p1_rg_ree.sim!=13]
blank_df = pd.DataFrame()
for sim in p1_rg_ree.sim.unique():
    temp_df = p1_rg_ree[p1_rg_ree.sim==sim].iloc[-40000:,:]
    blank_df = blank_df.append(temp_df)
p1_rg_ree = blank_df.copy()
del blank_df
del temp_df
del p1_rg, p1_ree

#need to remove frames and simulations that are not necessary
#only show equilibrated part
#this is for bull
bull_rg = pd.read_csv("~/protamine_simulation/jupyter_nb/bull_Rg_master_out.csv")
bull_ree = pd.read_csv("~/protamine_simulation/jupyter_nb/bull_Ree_master_out.csv")
bull_rg_ree=pd.concat([bull_rg,bull_ree['Ree']],axis=1)
bull_rg_ree=bull_rg_ree[bull_rg_ree.sim!=14]
bull_rg_ree=bull_rg_ree[bull_rg_ree.sim!=15]
blank_df = pd.DataFrame()
for sim in bull_rg_ree.sim.unique():
    temp_df = bull_rg_ree[bull_rg_ree.sim==sim].iloc[-40000:,:]
    blank_df = blank_df.append(temp_df)
bull_rg_ree = blank_df.copy()

del blank_df
del temp_df
del bull_rg, bull_ree


# In[3]:


#human P1 protamine ESFF1 13 simulations
#drop the same trajectories as ff14sb 13 simulations
#use equilibrated region same region as ff14sb

p1_rg_esff1= pd.read_csv('/home/hshadman/protamine_simulation/ESFF1/protamine_human/P1/explicit_solvent/more_nosalt_runs/explicit_continued_from_implicit/Rg_Ree/p1_Rg_master_out.csv')
p1_ree_esff1= pd.read_csv('/home/hshadman/protamine_simulation/ESFF1/protamine_human/P1/explicit_solvent/more_nosalt_runs/explicit_continued_from_implicit/Rg_Ree/p1_Ree_master_out.csv')
p1_rg_ree_esff1=pd.concat([p1_rg_esff1,p1_ree_esff1['Ree']],axis=1)
p1_rg_ree_esff1=p1_rg_ree_esff1[p1_rg_ree_esff1.sim!=15]
p1_rg_ree_esff1=p1_rg_ree_esff1[p1_rg_ree_esff1.sim!=13]
blank_df = pd.DataFrame()
for sim in p1_rg_ree_esff1.sim.unique():
    temp_df = p1_rg_ree_esff1[p1_rg_ree_esff1.sim==sim].iloc[-40000:,:]
    blank_df = blank_df.append(temp_df)
p1_rg_ree_esff1 = blank_df.copy()
del blank_df
del temp_df
del p1_rg_esff1, p1_ree_esff1


# In[4]:


#the long simulations 
salmon_rg_oldff=pd.read_csv('~/protamine_simulation/jupyter_nb/salmon_Rg_eighth_oldff.csv')                                                                    
salmon_ree_oldff=pd.read_csv('~/protamine_simulation/jupyter_nb/salmon_Ree_eighth_oldff.csv')                                                                  
salmon_rg_ree_oldff=pd.concat([salmon_rg_oldff,salmon_ree_oldff['Ree']],axis=1)                                                  
salmon_rg_ree_oldff['ratio']=salmon_rg_ree_oldff.Ree.values**2/salmon_rg_ree_oldff.Rg.values**2                                  

p1_rg_oldff=pd.read_csv('~/protamine_simulation/jupyter_nb/p1_Rg_eighth_oldff.csv')                                                                    
p1_ree_oldff=pd.read_csv('~/protamine_simulation/jupyter_nb/p1_Ree_eighth_oldff.csv')                                                                  
p1_rg_ree_oldff=pd.concat([p1_rg_oldff,p1_ree_oldff['Ree']],axis=1)                                                  
p1_rg_ree_oldff['ratio']=p1_rg_ree_oldff.Ree.values**2/p1_rg_ree_oldff.Rg.values**2                                  

bull_rg_oldff=pd.read_csv('~/protamine_simulation/jupyter_nb/bull_Rg_fourth_oldff.csv')                                                                    
bull_ree_oldff=pd.read_csv('~/protamine_simulation/jupyter_nb/bull_Ree_fourth_oldff.csv')                                                                  
bull_rg_ree_oldff=pd.concat([bull_rg_oldff,bull_ree_oldff['Ree']],axis=1)     
bull_rg_ree_oldff['ratio']=bull_rg_ree_oldff.Ree.values**2/bull_rg_ree_oldff.Rg.values**2                                  

#use below two dataframes for plotting                                                                               

del salmon_rg_oldff, salmon_ree_oldff,p1_rg_oldff, p1_ree_oldff, bull_rg_oldff, bull_ree_oldff
  


# In[5]:


#protamine details (only for salmon dataframe i added Rg/Rg_mean)
salmon_rg_ree_ratheatmap=salmon_rg_ree.copy()
salmon_rg_ree_ratheatmap['ratio']=salmon_rg_ree_ratheatmap.Ree.values**2/salmon_rg_ree_ratheatmap.Rg.values**2
del salmon_rg_ree

test=salmon_rg_ree_ratheatmap.copy()
temp_rg_mean=[]
for sim in test.sim.unique():
    temp_rg_mean.append(list(np.repeat(test[test.sim==sim]['Rg'].mean(),
                                 test[test.sim==sim]['Rg'].shape[0])))
test['Rg_mean']=list(chain.from_iterable(temp_rg_mean))
test['Rg/Rg_mean']=test['Rg'].values/test['Rg_mean']    
salmon_rg_ree_ratheatmap=test.copy()
del test, temp_rg_mean


# In[6]:


#protamine details
p1_rg_ree_ratheatmap=p1_rg_ree.copy()
p1_rg_ree_ratheatmap['ratio']=p1_rg_ree_ratheatmap.Ree.values**2/p1_rg_ree_ratheatmap.Rg.values**2
del p1_rg_ree
p1_rg_ree_ratheatmap


# In[7]:


#protamine details
p1_rg_ree_esff1_ratheatmap=p1_rg_ree_esff1.copy()
p1_rg_ree_esff1_ratheatmap['ratio']=p1_rg_ree_esff1_ratheatmap.Ree.values**2/p1_rg_ree_esff1_ratheatmap.Rg.values**2
del p1_rg_ree_esff1
p1_rg_ree_esff1_ratheatmap


# In[8]:


#protamine details
bull_rg_ree_ratheatmap=bull_rg_ree.copy()
bull_rg_ree_ratheatmap['ratio']=bull_rg_ree_ratheatmap.Ree.values**2/bull_rg_ree_ratheatmap.Rg.values**2
del bull_rg_ree
bull_rg_ree_ratheatmap


# In[9]:


#salmon tensors received and converted to moments at below directory. only equilibrated region for moments for time
salmon_moments=pd.read_csv('/home/hshadman/protamine_simulation/jupyter_nb/salmon_moments_equilibrated.csv')
salmon_moments.insert(0,salmon_rg_ree_ratheatmap.columns[0],
                      salmon_rg_ree_ratheatmap[salmon_rg_ree_ratheatmap.columns[0]].values)
salmon_moments.insert(1,salmon_rg_ree_ratheatmap.columns[1],
                      salmon_rg_ree_ratheatmap[salmon_rg_ree_ratheatmap.columns[1]].values)
salmon_moments.insert(2,salmon_rg_ree_ratheatmap.columns[2],
                      salmon_rg_ree_ratheatmap[salmon_rg_ree_ratheatmap.columns[2]].values)
salmon_R1=salmon_moments.R1.values
salmon_R2=salmon_moments.R2.values
salmon_R3=salmon_moments.R3.values
salmon_Rg=(salmon_R1+salmon_R2+salmon_R3)**0.5
salmon_moments['asphericity']=(salmon_R1) - ((0.5) * ((salmon_R2) + (salmon_R3)))
salmon_moments['acylindricity']=(salmon_R2) - (salmon_R3)
salmon_moments['RSA'] = ((salmon_moments.asphericity.values**2) + (0.75*(salmon_moments.acylindricity.values**2)))/(salmon_Rg**4)
del salmon_R1, salmon_R2, salmon_R3, salmon_Rg


#p1 tensors received and converted to moments at below directory. only equilibrated region for moments for time
p1_moments=pd.read_csv('/home/hshadman/protamine_simulation/jupyter_nb/p1_moments_equilibrated.csv')
p1_moments.insert(0,p1_rg_ree_ratheatmap.columns[0],
                      p1_rg_ree_ratheatmap[p1_rg_ree_ratheatmap.columns[0]].values)
p1_moments.insert(1,p1_rg_ree_ratheatmap.columns[1],
                      p1_rg_ree_ratheatmap[p1_rg_ree_ratheatmap.columns[1]].values)
p1_moments.insert(2,p1_rg_ree_ratheatmap.columns[2],
                      p1_rg_ree_ratheatmap[p1_rg_ree_ratheatmap.columns[2]].values)

p1_R1=p1_moments.R1.values
p1_R2=p1_moments.R2.values
p1_R3=p1_moments.R3.values
p1_Rg=(p1_R1+p1_R2+p1_R3)**0.5
p1_moments['asphericity']=(p1_R1) - ((0.5) * ((p1_R2) + (p1_R3)))
p1_moments['acylindricity']=(p1_R2) - (p1_R3)
p1_moments['RSA'] = ((p1_moments.asphericity.values**2) + (0.75*(p1_moments.acylindricity.values**2)))/(p1_Rg**4)
del p1_R1, p1_R2, p1_R3, p1_Rg


#bull tensors received and converted to moments at below directory. only equilibrated region for moments for time
bull_moments=pd.read_csv('/home/hshadman/protamine_simulation/jupyter_nb/bull_moments_equilibrated.csv')
bull_moments.insert(0,bull_rg_ree_ratheatmap.columns[0],
                      bull_rg_ree_ratheatmap[bull_rg_ree_ratheatmap.columns[0]].values)
bull_moments.insert(1,bull_rg_ree_ratheatmap.columns[1],
                      bull_rg_ree_ratheatmap[bull_rg_ree_ratheatmap.columns[1]].values)
bull_moments.insert(2,bull_rg_ree_ratheatmap.columns[2],
                      bull_rg_ree_ratheatmap[bull_rg_ree_ratheatmap.columns[2]].values)
bull_R1=bull_moments.R1.values
bull_R2=bull_moments.R2.values
bull_R3=bull_moments.R3.values
bull_Rg=(bull_R1+bull_R2+bull_R3)**0.5
bull_moments['asphericity']=(bull_R1) - ((0.5) * ((bull_R2) + (bull_R3)))
bull_moments['acylindricity']=(bull_R2) - (bull_R3)
bull_moments['RSA'] = ((bull_moments.asphericity.values**2) + (0.75*(bull_moments.acylindricity.values**2)))/(bull_Rg**4)
del bull_R1, bull_R2, bull_R3, bull_Rg


# In[10]:


salmon_moments


# In[11]:


bull_moments


# In[12]:


p1_moments


# In[13]:


bull_moments


# In[14]:


#BE CAREFUL with working directory for this first part of code
#-------
#only SAW_equil_chain_rg_ree was changed to SAW_SAW_equil_chain_rg_ree
# BE CAREFUL when changing the RW to SAW and RW
cwd= os.getcwd()
os.chdir('/home/hshadman/polym_sep/SAW/varying_chain_lengths/')
chain_lengths=[]
for file in os.listdir():
    if file.split('_')[0]=='chain':
        chain_lengths.append(int(file.split('_')[2]))
os.chdir(cwd)
#--------
folders=[]
chain_lengths=sorted(chain_lengths)
for i in chain_lengths:
    folders.append('/home/hshadman/polym_sep/SAW/varying_chain_lengths/chain_length_'+str(i)+'/varying_epas_chain_length_'+str(i)+'_master_out.csv')
j=0
for file in folders:
    test=pd.read_csv(file)
    test['Rend2']=test.Rendx+test.Rendy+test.Rendz
    test['Rg2']=test.Rgx+test.Rgy+test.Rgz
    test['ratio']=test.Rend2.values/test.Rg2.values
    test['asphericity']=test.Rgx.values-(0.5*(test.Rgy.values+test.Rgz.values))
    test['acylindricity']=test.Rgy.values-test.Rgz.values
    test['RSA']=((test.asphericity.values**2+(0.75*test.acylindricity.values**2))/(test.Rg2.values)**2)**0.5

    epas_considered = test.epas.unique()
    blank_df = pd.DataFrame()
    for epas in test.epas.unique():
        if epas in epas_considered:
            frames_number=len(test[test.epas==epas].index)
            equil_frames=int(0.90*frames_number)
            temp_df = test[test.epas==epas].iloc[-equil_frames:,:]
            blank_df = blank_df.append(temp_df)
    equil_test=blank_df.copy()
    equil_test=equil_test.drop(['frames','econf'],axis=1)
    if j==0:
        SAW_equil_chain_rg_ree=equil_test.copy()
        for epas_val in equil_test.epas.unique():
            if epas_val==equil_test.epas.unique()[0]:
                running_df=pd.DataFrame(data=equil_test[equil_test.epas==epas_val].mean(axis=0)).T
            elif epas_val!=equil_test.epas.unique()[0]:
                running_df=running_df.append(pd.DataFrame(data=equil_test[equil_test.epas==epas_val].mean(axis=0)).T)            
        summary_df_SAW=running_df.copy()
    else:
        SAW_equil_chain_rg_ree=SAW_equil_chain_rg_ree.append(equil_test)
        for epas_val in equil_test.epas.unique():
            if epas_val==equil_test.epas.unique()[0]:
                running_df=pd.DataFrame(data=equil_test[equil_test.epas==epas_val].mean(axis=0)).T
            elif epas_val!=equil_test.epas.unique()[0]:
                running_df=running_df.append(pd.DataFrame(data=equil_test[equil_test.epas==epas_val].mean(axis=0)).T)            
        summary_df_SAW=pd.concat([summary_df_SAW,running_df],axis=0)    
    del blank_df, temp_df, test, equil_test    
    j+=1
del running_df, equil_frames
summary_df_SAW=summary_df_SAW.round({'epas':2})
SAW_equil_chain_rg_ree=SAW_equil_chain_rg_ree.round({'epas':2})


# In[15]:


#BE CAREFUL with working directory for this first part of code
#-------
#only RW_equil_chain_rg_ree was changed to RW_RW_equil_chain_rg_ree
# BE CAREFUL when changing the RW to SAW and RW
cwd= os.getcwd()
os.chdir('/home/hshadman/polym_sep/random_walk/varying_chain_lengths/')
chain_lengths=[]
for file in os.listdir():
    if file.split('_')[0]=='chain':
        chain_lengths.append(int(file.split('_')[2]))
os.chdir(cwd)
#--------
folders=[]
chain_lengths=sorted(chain_lengths)
for i in chain_lengths:
    folders.append('/home/hshadman/polym_sep/random_walk/varying_chain_lengths/chain_length_'+str(i)+'/varying_epas_chain_length_'+str(i)+'_master_out.csv')
j=0
for file in folders:
    test=pd.read_csv(file)
    test['Rend2']=test.Rendx+test.Rendy+test.Rendz
    test['Rg2']=test.Rgx+test.Rgy+test.Rgz
    test['ratio']=test.Rend2.values/test.Rg2.values
    test['asphericity']=test.Rgx.values-(0.5*(test.Rgy.values+test.Rgz.values))
    test['acylindricity']=test.Rgy.values-test.Rgz.values
    test['RSA']=((test.asphericity.values**2+(0.75*test.acylindricity.values**2))/(test.Rg2.values)**2)**0.5

    epas_considered = test.epas.unique()
    blank_df = pd.DataFrame()
    for epas in test.epas.unique():
        if epas in epas_considered:
            frames_number=len(test[test.epas==epas].index)
            equil_frames=int(0.90*frames_number)
            temp_df = test[test.epas==epas].iloc[-equil_frames:,:]
            blank_df = blank_df.append(temp_df)
    equil_test=blank_df.copy()
    equil_test=equil_test.drop(['frames','econf'],axis=1)
    if j==0:
        RW_equil_chain_rg_ree=equil_test.copy()
        for epas_val in equil_test.epas.unique():
            if epas_val==equil_test.epas.unique()[0]:
                running_df=pd.DataFrame(data=equil_test[equil_test.epas==epas_val].mean(axis=0)).T
            elif epas_val!=equil_test.epas.unique()[0]:
                running_df=running_df.append(pd.DataFrame(data=equil_test[equil_test.epas==epas_val].mean(axis=0)).T)            
        summary_df_RW=running_df.copy()
    else:
        RW_equil_chain_rg_ree=RW_equil_chain_rg_ree.append(equil_test)
        for epas_val in equil_test.epas.unique():
            if epas_val==equil_test.epas.unique()[0]:
                running_df=pd.DataFrame(data=equil_test[equil_test.epas==epas_val].mean(axis=0)).T
            elif epas_val!=equil_test.epas.unique()[0]:
                running_df=running_df.append(pd.DataFrame(data=equil_test[equil_test.epas==epas_val].mean(axis=0)).T)            
        summary_df_RW=pd.concat([summary_df_RW,running_df],axis=0)    
    del blank_df, temp_df, test, equil_test    
    j+=1
del running_df, equil_frames
summary_df_RW=summary_df_RW.round({'epas':2})
RW_equil_chain_rg_ree=RW_equil_chain_rg_ree.round({'epas':2})


# In[16]:


#PEI code
pei_ratheatmap=pd.read_csv("/home/hshadman/caleb_cyu1/jupyter_nb/pei_Rg_Ree_master_out.csv")
pei_ratheatmap['ratio']=pei_ratheatmap.Ree.values**2/pei_ratheatmap.Rg.values**2

test=pei_ratheatmap.copy()

temp_rg_mean=[]
for sim in test.proton.unique():
    temp_rg_mean.append(list(np.repeat(test[test.proton==sim]['Rg'].mean(),
                                 test[test.proton==sim]['Rg'].shape[0])))
    
test['Rg_mean']=list(chain.from_iterable(temp_rg_mean))
test['Rg/Rg_mean']=test['Rg'].values/test['Rg_mean']

pei_ratheatmap=test.copy()
del test, temp_rg_mean


# In[17]:


pei_ratheatmap


# # formula from 2012 macromolecules elder jayaraman paper
# $$
# asphericity=(R_1)-\frac{1}{2}(R_2+R_3) \\  
# acylindricity=(R_2-R_3) \\
# RSA=\frac{asphericity^2 + 0.75acylindricity^2}{R_g^4}
# $$

# In[18]:


#double check RSA definitions!!! (have been double checked)
pei_moments=pd.read_csv("/home/hshadman/caleb_cyu1/jupyter_nb/pei_moments.csv")
pei_tensor=pd.read_csv("/home/hshadman/caleb_cyu1/pei_tensor/pei_tensor_master_out.csv")
pei_moments.insert(0,pei_tensor.columns[0],pei_tensor[pei_tensor.columns[0]])
pei_moments.insert(1,pei_tensor.columns[1],pei_tensor[pei_tensor.columns[1]])
pei_moments.insert(2,pei_tensor.columns[2],pei_tensor[pei_tensor.columns[2]])
pei_R1=pei_moments.R1.values
pei_R2=pei_moments.R2.values
pei_R3=pei_moments.R3.values
pei_Rg=pei_moments.Rg.values
pei_moments['asphericity']=(pei_R1) - ((0.5) * ((pei_R2) + (pei_R3)))
pei_moments['acylindricity']=(pei_R2) - (pei_R3)
pei_moments['RSA'] = ((pei_moments.asphericity.values**2) + (0.75*(pei_moments.acylindricity.values**2)))/(pei_Rg**4)
del pei_tensor, pei_R1, pei_R2, pei_R3, pei_Rg
pei_moments


# In[19]:


(pei_moments.R1+pei_moments.R2+pei_moments.R3)**0.5


# In[20]:


GW_ind_runs_chainlen25 = pd.read_csv('/home/hshadman/polym_sep/gaussian_chain/code/gaussian_single_chain_length_25_6_runs.csv')
GW_ind_runs_chainlen100 = pd.read_csv('/home/hshadman/polym_sep/gaussian_chain/code/gaussian_single_chain_length_100_6_runs.csv')
GW_ind_runs_chainlen100_2 = pd.read_csv('/home/hshadman/polym_sep/gaussian_chain/code/2nd_time_gaussian_single_chain_length_100_9_runs.csv')
GW_ind_runs_chainlen100_3 = pd.read_csv('/home/hshadman/polym_sep/gaussian_chain/code/3rd_time_gaussian_single_chain_length_100_6_runs.csv')
GW_ind_runs_chainlen150 = pd.read_csv('/home/hshadman/polym_sep/gaussian_chain/code/gaussian_single_chain_length_150_6_runs.csv')

#add rg/rg_mean column for each dataframe saved

#first the GW_ind_runs_chainlen25 (2 instances of dataframe in each chunk)
test=GW_ind_runs_chainlen25.copy()
test['Rg']=test.Rg2.values**0.5
temp_rg_mean=[]
for run_num in test.run_number.unique():
    temp_rg_mean.append(list(np.repeat(test[test.run_number==run_num]['Rg'].mean(),
                                 test[test.run_number==run_num]['Rg'].shape[0])))

test['Rg_mean']=list(chain.from_iterable(temp_rg_mean))
test['Rg/Rg_mean']=test['Rg'].values/test['Rg_mean']

GW_ind_runs_chainlen25=test.copy()
del test, temp_rg_mean

#second the GW_ind_runs_chainlen100 (2 instances of dataframe in each chunk)
test=GW_ind_runs_chainlen100.copy()
test['Rg']=test.Rg2.values**0.5
temp_rg_mean=[]
for run_num in test.run_number.unique():
    temp_rg_mean.append(list(np.repeat(test[test.run_number==run_num]['Rg'].mean(),
                                 test[test.run_number==run_num]['Rg'].shape[0])))

test['Rg_mean']=list(chain.from_iterable(temp_rg_mean))
test['Rg/Rg_mean']=test['Rg'].values/test['Rg_mean']

GW_ind_runs_chainlen100=test.copy()
del test, temp_rg_mean

#second second the GW_ind_runs_chainlen100 (2 instances of dataframe in each chunk)
test=GW_ind_runs_chainlen100_2.copy()
test['Rg']=test.Rg2.values**0.5
temp_rg_mean=[]
for run_num in test.run_number.unique():
    temp_rg_mean.append(list(np.repeat(test[test.run_number==run_num]['Rg'].mean(),
                                 test[test.run_number==run_num]['Rg'].shape[0])))

test['Rg_mean']=list(chain.from_iterable(temp_rg_mean))
test['Rg/Rg_mean']=test['Rg'].values/test['Rg_mean']

GW_ind_runs_chainlen100_2=test.copy()

#second second second the GW_ind_runs_chainlen100 (2 instances of dataframe in each chunk)
test=GW_ind_runs_chainlen100_3.copy()
test['Rg']=test.Rg2.values**0.5
temp_rg_mean=[]
for run_num in test.run_number.unique():
    temp_rg_mean.append(list(np.repeat(test[test.run_number==run_num]['Rg'].mean(),
                                 test[test.run_number==run_num]['Rg'].shape[0])))

test['Rg_mean']=list(chain.from_iterable(temp_rg_mean))
test['Rg/Rg_mean']=test['Rg'].values/test['Rg_mean']

GW_ind_runs_chainlen100_3=test.copy()



#update run number in second  and 3rd dataframes of chain length 100
GW_ind_runs_chainlen100_2['run_number']= GW_ind_runs_chainlen100_2['run_number'] + GW_ind_runs_chainlen100['run_number'].max()
GW_ind_runs_chainlen100_3['run_number']= GW_ind_runs_chainlen100_3['run_number'] + GW_ind_runs_chainlen100_2['run_number'].max()
#combine the three dataframes
GW_ind_runs_chainlen100 = pd.concat([GW_ind_runs_chainlen100,
                                     GW_ind_runs_chainlen100_2,
                                    GW_ind_runs_chainlen100_3],axis=0)
del test, temp_rg_mean,GW_ind_runs_chainlen100_2,GW_ind_runs_chainlen100_3



#third the GW_ind_runs_chainlen150 (2 instances of dataframe in each chunk)
test=GW_ind_runs_chainlen150.copy()
test['Rg']=test.Rg2.values**0.5
temp_rg_mean=[]
for run_num in test.run_number.unique():
    temp_rg_mean.append(list(np.repeat(test[test.run_number==run_num]['Rg'].mean(),
                                 test[test.run_number==run_num]['Rg'].shape[0])))

test['Rg_mean']=list(chain.from_iterable(temp_rg_mean))
test['Rg/Rg_mean']=test['Rg'].values/test['Rg_mean']

GW_ind_runs_chainlen150=test.copy()
del test, temp_rg_mean


# In[21]:


GW_equil_chain_rg_ree = pd.read_csv('/home/hshadman/polym_sep/gaussian_chain/code/gaussian_chain_master_out.csv')[['chain_length','Rg2','Rend2','ratio']]
GW_moments = pd.read_csv('/home/hshadman/polym_sep/gaussian_chain/code/gaussian_chain_master_out.csv')[['chain_length','Rg2','R1','R2','R3','asphericity','acylindricity','RSA']]
#GW_moments already has asphericity, acylindricity and RSA from csv file but we are replacing that and recalcualting

GW_R1=GW_moments.R1.values
GW_R2=GW_moments.R2.values
GW_R3=GW_moments.R3.values
GW_Rg=GW_moments.Rg2.values**0.5
GW_moments['asphericity']=(GW_R1) - ((0.5) * ((GW_R2) + (GW_R3)))
GW_moments['acylindricity']=(GW_R2) - (GW_R3)
GW_moments['RSA'] = ((GW_moments.asphericity.values**2) + (0.75*(GW_moments.acylindricity.values**2)))/(GW_Rg**4)
del GW_R1, GW_R2, GW_R3, GW_Rg


# In[22]:


GW_moments


# In[23]:


#adding a column of Rg/Rg_mean for GW ONLY
#careful about changing anything here
#only GW don't do this for RW and SAW having epas and chainlength complicates things
test=GW_equil_chain_rg_ree.copy()
test['Rg']=test.Rg2.values**0.5
temp_rg_mean=[]
for chain_length in test.chain_length.unique():
    temp_rg_mean.append(list(np.repeat(test[test.chain_length==chain_length]['Rg'].mean(),
                                 test[test.chain_length==chain_length]['Rg'].shape[0])))
    
test['Rg_mean']=list(chain.from_iterable(temp_rg_mean))
test['Rg/Rg_mean']=test['Rg'].values/test['Rg_mean']

GW_equil_chain_rg_ree=test.copy()
del test, temp_rg_mean


# In[24]:


#ab40_ff14SB (disordered)
#CAREFUL if using frame numbers -- they repeat
ab40_ff14sb_rg_noimage = pd.read_csv("/home/hshadman/trajectories_folder_Li_zhengxin/ab40-ff14SB/analysis/Rg_ree_data/ab40_ff14SB_Rg_noimage_master_out.csv")
ab40_ff14sb_ree_noimage = pd.read_csv("/home/hshadman/trajectories_folder_Li_zhengxin/ab40-ff14SB/analysis/Rg_ree_data/ab40_ff14SB_Ree_noimage_master_out.csv")
ab40_ff14sb_rg_ree_ratheatmap_noimage = pd.concat([ab40_ff14sb_rg_noimage.frames,
                                           ab40_ff14sb_rg_noimage.sim,
                                          ab40_ff14sb_rg_noimage.Rg,
                                          ab40_ff14sb_ree_noimage.Ree],axis=1)
ab40_ff14sb_rg_ree_ratheatmap_noimage['ratio']=ab40_ff14sb_rg_ree_ratheatmap_noimage.Ree**2/ab40_ff14sb_rg_ree_ratheatmap_noimage.Rg**2

test=ab40_ff14sb_rg_ree_ratheatmap_noimage.copy()

temp_rg_mean=[]
for sim in test.sim.unique():
    temp_rg_mean.append(list(np.repeat(test[test.sim==sim]['Rg'].mean(),
                                 test[test.sim==sim]['Rg'].shape[0])))
    
test['Rg_mean']=list(chain.from_iterable(temp_rg_mean))
test['Rg/Rg_mean']=test['Rg'].values/test['Rg_mean']

ab40_ff14sb_rg_ree_ratheatmap_noimage=test.copy()
del test, temp_rg_mean

del ab40_ff14sb_rg_noimage, ab40_ff14sb_ree_noimage


# In[25]:


ab40_ff14sb_rg_ree_ratheatmap_noimage


# In[26]:


ab40_ff14sb_moments=pd.read_csv("/home/hshadman/trajectories_folder_Li_zhengxin/jupyter_nb/ab40_ff14sb_rg_noimage_moments.csv")

ab40_ff14sb_R1=ab40_ff14sb_moments.R1.values
ab40_ff14sb_R2=ab40_ff14sb_moments.R2.values
ab40_ff14sb_R3=ab40_ff14sb_moments.R3.values
ab40_ff14sb_Rg=(ab40_ff14sb_R1+ab40_ff14sb_R2+ab40_ff14sb_R3)**0.5
ab40_ff14sb_moments['asphericity']=(ab40_ff14sb_R1) - ((0.5) * ((ab40_ff14sb_R2) + (ab40_ff14sb_R3)))
ab40_ff14sb_moments['acylindricity']=(ab40_ff14sb_R2) - (ab40_ff14sb_R3)
ab40_ff14sb_moments['RSA'] = ((ab40_ff14sb_moments.asphericity.values**2) + (0.75*(ab40_ff14sb_moments.acylindricity.values**2)))/(ab40_ff14sb_Rg**4)
del ab40_ff14sb_R1, ab40_ff14sb_R2, ab40_ff14sb_R3, ab40_ff14sb_Rg
ab40_ff14sb_moments


# In[27]:


#ab40_ff14SB (disordered)
#CAREFUO if using frame numbers -- they repeat
ab40_ff14sb_rg = pd.read_csv("/home/hshadman/trajectories_folder_Li_zhengxin/ab40-ff14SB/analysis/Rg_ree_data/ab40_ff14SB_Rg_master_out.csv")
ab40_ff14sb_ree = pd.read_csv("/home/hshadman/trajectories_folder_Li_zhengxin/ab40-ff14SB/analysis/Rg_ree_data/ab40_ff14SB_Ree_master_out.csv")
ab40_ff14sb_rg_ree_ratheatmap = pd.concat([ab40_ff14sb_rg.frames,
                                           ab40_ff14sb_rg.sim,
                                          ab40_ff14sb_rg.Rg,
                                          ab40_ff14sb_ree.Ree],axis=1)
ab40_ff14sb_rg_ree_ratheatmap['ratio']=ab40_ff14sb_rg_ree_ratheatmap.Ree**2/ab40_ff14sb_rg_ree_ratheatmap.Rg**2

test=ab40_ff14sb_rg_ree_ratheatmap.copy()

temp_rg_mean=[]
for sim in test.sim.unique():
    temp_rg_mean.append(list(np.repeat(test[test.sim==sim]['Rg'].mean(),
                                 test[test.sim==sim]['Rg'].shape[0])))
    
test['Rg_mean']=list(chain.from_iterable(temp_rg_mean))
test['Rg/Rg_mean']=test['Rg'].values/test['Rg_mean']

ab40_ff14sb_rg_ree_ratheatmap=test.copy()
del test, temp_rg_mean

del ab40_ff14sb_rg, ab40_ff14sb_ree


# In[28]:


ab40_ff14sb_rg_ree_ratheatmap


# In[29]:


#tauF4_esff1 (disordered)
#CAREFUO if using frame numbers -- they repeat
tauF4_esff1_rg = pd.read_csv("/home/hshadman/trajectories_folder_Li_zhengxin/tauF4-ESFF1/analysis/Rg_ree_data/tauF4_ESFF1_Rg_master_out.csv")
tauF4_esff1_ree = pd.read_csv("/home/hshadman/trajectories_folder_Li_zhengxin/tauF4-ESFF1/analysis/Rg_ree_data/tauF4_ESFF1_Ree_master_out.csv")
tauF4_esff1_rg_ree_ratheatmap = pd.concat([tauF4_esff1_rg.frames,
                                           tauF4_esff1_rg.sim,
                                          tauF4_esff1_rg.Rg,
                                          tauF4_esff1_ree.Ree],axis=1)
tauF4_esff1_rg_ree_ratheatmap['ratio']=tauF4_esff1_rg_ree_ratheatmap.Ree**2/tauF4_esff1_rg_ree_ratheatmap.Rg**2

test=tauF4_esff1_rg_ree_ratheatmap.copy()

temp_rg_mean=[]
for sim in test.sim.unique():
    temp_rg_mean.append(list(np.repeat(test[test.sim==sim]['Rg'].mean(),
                                 test[test.sim==sim]['Rg'].shape[0])))
    
test['Rg_mean']=list(chain.from_iterable(temp_rg_mean))
test['Rg/Rg_mean']=test['Rg'].values/test['Rg_mean']

tauF4_esff1_rg_ree_ratheatmap=test.copy()
del test, temp_rg_mean

del tauF4_esff1_rg, tauF4_esff1_ree


# In[30]:


tauF4_esff1_rg_ree_ratheatmap


# In[31]:


tauF4_esff1_moments=pd.read_csv("/home/hshadman/trajectories_folder_Li_zhengxin/jupyter_nb/tauF4_esff1_rg_moments.csv")

tauF4_esff1_R1=tauF4_esff1_moments.R1.values
tauF4_esff1_R2=tauF4_esff1_moments.R2.values
tauF4_esff1_R3=tauF4_esff1_moments.R3.values
tauF4_esff1_Rg=(tauF4_esff1_R1+tauF4_esff1_R2+tauF4_esff1_R3)**0.5
tauF4_esff1_moments['asphericity']=(tauF4_esff1_R1) - ((0.5) * ((tauF4_esff1_R2) + (tauF4_esff1_R3)))
tauF4_esff1_moments['acylindricity']=(tauF4_esff1_R2) - (tauF4_esff1_R3)
tauF4_esff1_moments['RSA'] = ((tauF4_esff1_moments.asphericity.values**2) + (0.75*(tauF4_esff1_moments.acylindricity.values**2)))/(tauF4_esff1_Rg**4)
del tauF4_esff1_R1, tauF4_esff1_R2, tauF4_esff1_R3, tauF4_esff1_Rg
tauF4_esff1_moments


# In[32]:


#FKBP12_ESFF1 (ordered)

FKBP12_ESFF1_rg = pd.read_csv("/home/hshadman/trajectories_folder_Li_zhengxin/FKBP12/FKBP12/ESFFx/analysis/Rg_ree_data/FKBP12_ESFF1_Rg_master_out.csv")
FKBP12_ESFF1_ree = pd.read_csv("/home/hshadman/trajectories_folder_Li_zhengxin/FKBP12/FKBP12/ESFFx/analysis/Rg_ree_data/FKBP12_ESFF1_Ree_master_out.csv")
FKBP12_ESFF1_rg_ree_ratheatmap = pd.concat([FKBP12_ESFF1_rg.frames,
                                           FKBP12_ESFF1_rg.sim,
                                          FKBP12_ESFF1_rg.Rg,
                                          FKBP12_ESFF1_ree.Ree],axis=1)
FKBP12_ESFF1_rg_ree_ratheatmap['ratio']=FKBP12_ESFF1_rg_ree_ratheatmap.Ree**2/FKBP12_ESFF1_rg_ree_ratheatmap.Rg**2
test=FKBP12_ESFF1_rg_ree_ratheatmap.copy()

temp_rg_mean=[]
for sim in test.sim.unique():
    temp_rg_mean.append(list(np.repeat(test[test.sim==sim]['Rg'].mean(),
                                 test[test.sim==sim]['Rg'].shape[0])))
    
test['Rg_mean']=list(chain.from_iterable(temp_rg_mean))
test['Rg/Rg_mean']=test['Rg'].values/test['Rg_mean']

FKBP12_ESFF1_rg_ree_ratheatmap=test.copy()
del test, temp_rg_mean

del FKBP12_ESFF1_rg, FKBP12_ESFF1_ree


# In[33]:


FKBP12_ESFF1_rg_ree_ratheatmap


# In[34]:


FKBP12_ESFF1_moments=pd.read_csv("/home/hshadman/trajectories_folder_Li_zhengxin/jupyter_nb/FKBP12_ESFF1_rg_moments.csv")

FKBP12_ESFF1_R1=FKBP12_ESFF1_moments.R1.values
FKBP12_ESFF1_R2=FKBP12_ESFF1_moments.R2.values
FKBP12_ESFF1_R3=FKBP12_ESFF1_moments.R3.values
FKBP12_ESFF1_Rg=(FKBP12_ESFF1_R1+FKBP12_ESFF1_R2+FKBP12_ESFF1_R3)**0.5
FKBP12_ESFF1_moments['asphericity']=(FKBP12_ESFF1_R1) - ((0.5) * ((FKBP12_ESFF1_R2) + (FKBP12_ESFF1_R3)))
FKBP12_ESFF1_moments['acylindricity']=(FKBP12_ESFF1_R2) - (FKBP12_ESFF1_R3)
FKBP12_ESFF1_moments['RSA'] = ((FKBP12_ESFF1_moments.asphericity.values**2) + (0.75*(FKBP12_ESFF1_moments.acylindricity.values**2)))/(FKBP12_ESFF1_Rg**4)
del FKBP12_ESFF1_R1, FKBP12_ESFF1_R2, FKBP12_ESFF1_R3, FKBP12_ESFF1_Rg
FKBP12_ESFF1_moments


# In[35]:


#ubiquitin_ESFF1 (ordered)

ubiquitin_ESFF1_rg = pd.read_csv("/home/hshadman/trajectories_folder_Li_zhengxin/ubiquitin-ESFF1/combined_traj1/Rg_ree_data/ubiquitin_ESFF1_Rg_master_out.csv")
ubiquitin_ESFF1_ree = pd.read_csv("/home/hshadman/trajectories_folder_Li_zhengxin/ubiquitin-ESFF1/combined_traj1/Rg_ree_data/ubiquitin_ESFF1_Ree_master_out.csv")
ubiquitin_ESFF1_rg_ree_ratheatmap = pd.concat([ubiquitin_ESFF1_rg.frames,
                                           ubiquitin_ESFF1_rg.sim,
                                          ubiquitin_ESFF1_rg.Rg,
                                          ubiquitin_ESFF1_ree.Ree],axis=1)
ubiquitin_ESFF1_rg_ree_ratheatmap['ratio']=ubiquitin_ESFF1_rg_ree_ratheatmap.Ree**2/ubiquitin_ESFF1_rg_ree_ratheatmap.Rg**2

test=ubiquitin_ESFF1_rg_ree_ratheatmap.copy()

temp_rg_mean=[]
for sim in test.sim.unique():
    temp_rg_mean.append(list(np.repeat(test[test.sim==sim]['Rg'].mean(),
                                 test[test.sim==sim]['Rg'].shape[0])))
    
test['Rg_mean']=list(chain.from_iterable(temp_rg_mean))
test['Rg/Rg_mean']=test['Rg'].values/test['Rg_mean']

ubiquitin_ESFF1_rg_ree_ratheatmap=test.copy()
del test, temp_rg_mean

del ubiquitin_ESFF1_rg, ubiquitin_ESFF1_ree


# In[36]:


ubiquitin_ESFF1_rg_ree_ratheatmap


# In[37]:


ubiquitin_ESFF1_moments= pd.read_csv("/home/hshadman/trajectories_folder_Li_zhengxin/jupyter_nb/ubiquitin_ESFF1_rg_moments.csv")

ubiquitin_ESFF1_R1=ubiquitin_ESFF1_moments.R1.values
ubiquitin_ESFF1_R2=ubiquitin_ESFF1_moments.R2.values
ubiquitin_ESFF1_R3=ubiquitin_ESFF1_moments.R3.values
ubiquitin_ESFF1_Rg=(ubiquitin_ESFF1_R1+ubiquitin_ESFF1_R2+ubiquitin_ESFF1_R3)**0.5
ubiquitin_ESFF1_moments['asphericity']=(ubiquitin_ESFF1_R1) - ((0.5) * ((ubiquitin_ESFF1_R2) + (ubiquitin_ESFF1_R3)))
ubiquitin_ESFF1_moments['acylindricity']=(ubiquitin_ESFF1_R2) - (ubiquitin_ESFF1_R3)
ubiquitin_ESFF1_moments['RSA'] = ((ubiquitin_ESFF1_moments.asphericity.values**2) + (0.75*(ubiquitin_ESFF1_moments.acylindricity.values**2)))/(ubiquitin_ESFF1_Rg**4)
del ubiquitin_ESFF1_R1, ubiquitin_ESFF1_R2, ubiquitin_ESFF1_R3, ubiquitin_ESFF1_Rg
ubiquitin_ESFF1_moments


# In[38]:


#lush_ESFF1 (ordered)
lush_ESFF1_rg = pd.read_csv("/home/hshadman/trajectories_folder_Li_zhengxin/traj_LUSH/combined_traj1/Rg_ree_data/LUSH_ESFF1_Rg_master_out.csv")
lush_ESFF1_ree = pd.read_csv("/home/hshadman/trajectories_folder_Li_zhengxin/traj_LUSH/combined_traj1/Rg_ree_data/LUSH_ESFF1_Ree_master_out.csv")
lush_ESFF1_rg_ree_ratheatmap = pd.concat([lush_ESFF1_rg.frames,
                                           lush_ESFF1_rg.sim,
                                          lush_ESFF1_rg.Rg,
                                          lush_ESFF1_ree.Ree],axis=1)
lush_ESFF1_rg_ree_ratheatmap['ratio']=lush_ESFF1_rg_ree_ratheatmap.Ree**2/lush_ESFF1_rg_ree_ratheatmap.Rg**2

test=lush_ESFF1_rg_ree_ratheatmap.copy()

temp_rg_mean=[]
for sim in test.sim.unique():
    temp_rg_mean.append(list(np.repeat(test[test.sim==sim]['Rg'].mean(),
                                 test[test.sim==sim]['Rg'].shape[0])))
    
test['Rg_mean']=list(chain.from_iterable(temp_rg_mean))
test['Rg/Rg_mean']=test['Rg'].values/test['Rg_mean']

lush_ESFF1_rg_ree_ratheatmap=test.copy()
del test, temp_rg_mean

del lush_ESFF1_rg, lush_ESFF1_ree


# In[39]:


lush_ESFF1_rg_ree_ratheatmap


# In[40]:


lush_ESFF1_moments=pd.read_csv("/home/hshadman/trajectories_folder_Li_zhengxin/jupyter_nb/lush_ESFF1_rg_moments.csv")

lush_ESFF1_R1=lush_ESFF1_moments.R1.values
lush_ESFF1_R2=lush_ESFF1_moments.R2.values
lush_ESFF1_R3=lush_ESFF1_moments.R3.values
lush_ESFF1_Rg=(lush_ESFF1_R1+lush_ESFF1_R2+lush_ESFF1_R3)**0.5
lush_ESFF1_moments['asphericity']=(lush_ESFF1_R1) - ((0.5) * ((lush_ESFF1_R2) + (lush_ESFF1_R3)))
lush_ESFF1_moments['acylindricity']=(lush_ESFF1_R2) - (lush_ESFF1_R3)
lush_ESFF1_moments['RSA'] = ((lush_ESFF1_moments.asphericity.values**2) + (0.75*(lush_ESFF1_moments.acylindricity.values**2)))/(lush_ESFF1_Rg**4)
del lush_ESFF1_R1, lush_ESFF1_R2, lush_ESFF1_R3, lush_ESFF1_Rg
lush_ESFF1_moments


# # initialization complete




def fA_using_cdist(protein_var,GW_var,upto_ind_run,every_ith_snap,radius_):

    po_x=GW_var[GW_var.run_number<=upto_ind_run]['Rg/Rg_mean'].values
    po_y=GW_var[GW_var.run_number<=upto_ind_run]['ratio'].values

    po_x=(po_x-np.mean(po_x))/np.std(po_x)    
    po_y=(po_y-np.mean(po_y))/np.std(po_y)
    
    GW_points=np.c_[po_x, po_y]

    plt.scatter(po_x,po_y,color='black')

    pro_x=protein_var['Rg/Rg_mean'].values[::every_ith_snap]
    pro_y=protein_var['ratio'].values[::every_ith_snap]
    
    pro_x=(pro_x-np.mean(pro_x))/np.std(pro_x)
    pro_y=(pro_y-np.mean(pro_y))/np.std(pro_y)

    protein_points=np.c_[pro_x, pro_y]
    plt.scatter(po_x,po_y,color='black')
    plt.scatter(pro_x,pro_y,color='magenta',alpha=0.05)

    tree_GW=spatial.cKDTree(GW_points)
    tree_protein=spatial.cKDTree(protein_points)


    GW_points_within_range=np.array([])
    not_in_range=[]
    j = 0
    for point in protein_points:
        checking_GW_points=tree_GW.query_ball_point(point,radius_)
        if not checking_GW_points:
            not_in_range.append(point)
        else:
            GW_points_within_range=np.append(GW_points_within_range,checking_GW_points)
            GW_points_within_range=np.unique(GW_points_within_range)
        del checking_GW_points
        j+=1
        if j%100000==0:
            print(f'{j} protein snapshots completed')

    print(f'{((len(not_in_range)/protein_points.shape[0])*100)}% of protein snapshots not in range of GW')
    fA_value = GW_points_within_range.shape[0]/GW_points.shape[0]
    
    return fA_value
    
fA_using_cdist(GW_in,GW_var,upto_ind_run,every_ith_snap,radius_)

    
