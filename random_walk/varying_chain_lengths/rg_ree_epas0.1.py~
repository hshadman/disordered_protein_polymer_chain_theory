#!/user/bin/python

import os
import string
import csv
import math
import numpy as np
import pandas as pd

chain_lengths=['25','50','75','100','125','150']
folders=[]
for i in range(len(chain_lengths)):
    folders.append('../chain_length_'+chain_lengths[i]+'/epas_0.1/Monitor.dat')
j=0

for sim in folders:
    
    inputfile1 = open(sim, "r+") # open the file
    file1_lines = inputfile1.readlines()

    skiprow_index=[]
    bond_list=[]
    i=0
    for line in file1_lines:
        if line[0]=='#':
            a=2
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
    output.insert(1,'chain_length',np.repeat(chain_lengths[j],len(output.index)))
    if j==0:
        master_output=output.copy()
    else:
        master_output=master_output.append(output)
    j=j+1
master_output.to_csv('varying_chain_lengths_master_out_epas=0.1.csv',index=False)




