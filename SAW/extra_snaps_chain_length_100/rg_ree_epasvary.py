#!/user/bin/python
import os
import string
import csv
import math
import numpy as np
import pandas as pd
chain_length=int(os.getcwd().split('_')[len(os.getcwd().split('_'))-1])

epas=[epas_val.split('_')[len(epas_val.split('_'))-1] for epas_val in list(filter(lambda x:x[0:4]=='epas',os.listdir()))]


folders=[]
for i in range(len(epas)):
    folders.append('./epas_'+epas[i]+'/')

j=0
for running_epas in folders:
    ind_run_list=os.listdir(running_epas)
    for ind_run in ind_run_list:
        run_number = int(ind_run.split('_')[0][0])
        epas_val = float(running_epas.split('_')[1][0:-1])
        sim = running_epas+ind_run+'/Monitor.dat'
        i = 0
        input_file=pd.read_csv(sim,delim_whitespace=True,skiprows=[0], header=None)
        input_file.columns=['frames','econf','Rgx', 'Rgy','Rgz','Rendx','Rendy','Rendz']
        length=int(input_file.frames.max())
        output=input_file.loc[:,'frames'].to_frame().iloc[:length,:]
        output['econf']=input_file.iloc[int(i*length):int((i+1)*length),:].econf.values
        output['Rgx']=input_file.iloc[int(i*length):int((i+1)*length),:].Rgx.values
        output['Rgy']=input_file.iloc[int(i*length):int((i+1)*length),:].Rgy.values
        output['Rgz']=input_file.iloc[int(i*length):int((i+1)*length),:].Rgz.values
        output['Rendx']=input_file.iloc[int(i*length):int((i+1)*length),:].Rendx.values
        output['Rendy']=input_file.iloc[int(i*length):int((i+1)*length),:].Rendy.values
        output['Rendz']=input_file.iloc[int(i*length):int((i+1)*length),:].Rendz.values
        i=i+1
        output.insert(1,'chain_length',np.repeat(chain_length,len(output.index)))
        output.insert(1,'epas',np.repeat(epas_val,len(output.index)))
        output.insert(1,'run_number',np.repeat(run_number,len(output.index)))
        if j==0:
            master_output=output.copy()
        else:
            master_output=master_output.append(output)
        j=j+1
master_output.to_csv('extra_snaps_'+str(chain_length)+'_master_out.csv',index=False)
