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
    folders.append('epas_'+epas[i]+'/Monitor.dat')
j=0

for sim in folders:    
    inputfile1 = open(sim, "r+") # open the file
    file1_lines = inputfile1.readlines()
    skiprow_index=[]
    bond_list=[]
    i=0
    for line in file1_lines:
        if line[0]!='#': 
            frame=int(line[0:12])
        i=i+1
    i=0
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
    output.insert(1,'epas',np.repeat(epas[j],len(output.index)))
    if j==0:
        master_output=output.copy()
    else:
        master_output=master_output.append(output)
    j=j+1
master_output.to_csv('varying_epas_chain_length_'+str(chain_length)+'_master_out.csv',index=False)

