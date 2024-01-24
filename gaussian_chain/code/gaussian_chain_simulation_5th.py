#this is a copy of the gaussian_chain jupyter notebook (partial copy)
#i copied it over because i want to submit as a job one gaussian chain simulation

import numpy as np
import math
import pandas as pd
import random as rd
from numpy import linalg as la
import matplotlib.pyplot as plt
def Ree2(x,y,z):
    return ((x[0]-x[len(x)-1])**2+(y[0]-y[len(y)-1])**2+(z[0]-z[len(z)-1])**2)
def Rgx2(x,chain_length): 
    Rgx = np.sum((x - np.mean(x))**2)/chain_length
    return Rgx
def Rgy2(y,chain_length): 
    Rgy = np.sum((y - np.mean(y))**2)/chain_length
    return Rgy
def Rgz2(z,chain_length): 
    Rgz = np.sum((z - np.mean(z))**2)/chain_length
    return Rgz
def calculate_tensors(x,y,z,chain_length):
    xx = np.sum((x-np.mean(x))*(x-np.mean(x)))/chain_length
    xy = np.sum((x-np.mean(x))*(y-np.mean(y)))/chain_length
    xz = np.sum((x-np.mean(x))*(z-np.mean(z)))/chain_length
    yx = np.sum((y-np.mean(y))*(x-np.mean(x)))/chain_length
    yy = np.sum((y-np.mean(y))*(y-np.mean(y)))/chain_length 
    yz = np.sum((y-np.mean(y))*(z-np.mean(z)))/chain_length 
    zx = np.sum((z-np.mean(z))*(x-np.mean(x)))/chain_length 
    zy = np.sum((z-np.mean(z))*(y-np.mean(y)))/chain_length
    zz = np.sum((z-np.mean(z))*(z-np.mean(z)))/chain_length
    
    Sij = np.array([[xx,xy,xz],
                   [yx,yy,yz],
                   [zx,zy,zz]])
    return Sij
def diagonalize_tensor(rg_rensors):
    tensor_df = pd.DataFrame()
    for i in rg_tensors.index:# written by me
        temp_mat = np.zeros((3,3))
        temp_mat[0,:] = rg_tensors.loc[i,['XX','XY','XZ']].values
        temp_mat[1,:] = rg_tensors.loc[i,['XY','YY','YZ']].values
        temp_mat[2,:] = rg_tensors.loc[i,['XZ','YZ','ZZ']].values
        temp_mat = pd.DataFrame(np.sort(la.eig(temp_mat)[0].real)[::-1]).T
        temp_mat.columns = ['R1','R2','R3']
        tensor_df = tensor_df.append(temp_mat)
        if i%100000==0:
            print('Step: ', i,' complete')
    tensor_df.reset_index(drop=True)#.to_csv('gaussian_chain_moments.csv',index=False)
    return tensor_df.reset_index(drop=True)

def total_fixed_length(x,y,z):
    j=0
    seg_length=[]
    for chain_ind in np.arange(0,x.shape[0]):
        if x.shape[0]==y.shape[0] and y.shape[0]==z.shape[0]:
            bead_coor = np.array([x[chain_ind],y[chain_ind],z[chain_ind]])
            if j==0:
                no_change=1
            elif j>0 and j<=x.shape[0]:
                seg_length.append(((x[chain_ind]-x[chain_ind-1])**2+(y[chain_ind]-y[chain_ind-1])**2+(z[chain_ind]-z[chain_ind-1])**2)**0.5)

        else:
            print('error')
        j+=1
    return sum(seg_length)
    
#gaussian chain code accounting for interval and nosnaps

def gaussian_chain_3d(chain_length,nosnaps,interval,mu,sigma):
    global Rend2,Rg2, shape_ratio, rg_tensors, x, y, z
    chain_length=chain_length
    #no of steps
    x = np.zeros(chain_length)
    y = np.zeros(chain_length)
    z = np.zeros(chain_length)
    nosnaps = nosnaps
    interval=interval
    step_size=1
    snapshot=0
    mu = mu
    sigma = sigma #kuhn length
    delta_x=[]
    delta_x2 = []
    Rend2 = []
    rg_tensors=pd.DataFrame(columns=['XX','XY','XZ','YX','YY','YZ','ZX','ZY','ZZ'])
    Rg2 = []
    shape_ratio=[]
    mean_sq_disp_i_master = []
    while snapshot<(nosnaps*interval):
        mean_sq_disp_i = []
        for i in range(1,chain_length,1):
            x[i] = x[i-1] + rd.gauss(mu,sigma)
            y[i] = y[i-1] + rd.gauss(mu,sigma)
            z[i] = z[i-1] + rd.gauss(mu,sigma)        
            mean_sq_disp_i.append((x[i]-x[0])**2)
        if snapshot in np.arange(0,(nosnaps*interval),interval):
            Rend2.append(Ree2(x,y,z))
            Rg2.append(Rgx2(x,chain_length)+Rgy2(y,chain_length)+Rgz2(z,chain_length))
            shape_ratio.append(Ree2(x,y,z)/(Rgx2(x,chain_length)+Rgy2(y,chain_length)+Rgz2(z,chain_length)))
            rg_tensors.loc[len(rg_tensors)]=calculate_tensors(x,y,z,chain_length).flatten()
        snapshot = snapshot + 1   
    print(f'chain_length={chain_length}, Nosnaps={nosnaps}, interval= {interval}')   
    master_out=pd.DataFrame(np.array([Rg2,
                       Rend2,
                       shape_ratio]).T,columns=['Rg2','Rend2','ratio'])
    master_out=pd.concat([master_out,diagonalize_tensor(rg_tensors)],axis=1)#be careful here if messing w/ master_out
    master_out.insert(0,'chain_length',np.repeat(chain_length,nosnaps))
    return master_out
#this funciton is same as gaussian_chain_3d function but intended for multiple independent runs of single chain length
def gaussian_chain_length_single_chain_length_multiple_runs(chain_length_choice,no_of_runs,nosnaps,interval,mu,sigma):
    global final_master
    for k in range(1,no_of_runs+1):
        chain_lengths=[chain_length_choice]
        nosnaps=nosnaps
        interval=interval
        mu=mu
        sigma=sigma
        j=0
        for chain in chain_lengths:
            if j==0:
                master=gaussian_chain_3d(chain,nosnaps,interval,mu,sigma)
            else:
                running=gaussian_chain_3d(chain,nosnaps,interval,mu,sigma)
                master=pd.concat([master,running],axis=0,ignore_index=True)
            j+=1
        if k==1:   
            final_master=master.copy()
            final_master.insert(0,'run_number',np.repeat(k,final_master.shape[0]))
            print('1 run complete')
        elif k>1:
            master.insert(0,'run_number',np.repeat(k,master.shape[0]))
            final_master=pd.concat([final_master,master],axis=0,ignore_index=True)
            print(f'{k} runs complete')            
    return final_master

#run gaussian chain simulation
#specify no of runs, chain length etc. as demonstrated below 
chain_length_choice=3000
no_of_runs=64
no_snapshots=30000
gaussian_chain_length_single_chain_length_multiple_runs(chain_length_choice,no_of_runs,no_snapshots,10,0,1)

final_master.to_csv(f"reviewer_gaussian_chain_length_{chain_length_choice}_{no_snapshots}_{no_of_runs}runs.csv",
                    index=False)

