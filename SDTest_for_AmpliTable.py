#!/usr/bin/python3 
from  scipy import stats
import os
import sys
import re
import numpy as np

z_cutoff=int(sys.argv[3])

#process sample you want to test with base line
fi_clinical=open(sys.argv[1])
line=fi_clinical.readline()
line=fi_clinical.readline()
ID2Soft2reads_cl={}
ID2total_reads_cl={}
while(line):
    line=line.rstrip()
    arr=line.split('\t')
    if not arr[1] in ID2Soft2reads_cl:
        ID2Soft2reads_cl[arr[1]]={}
    if not arr[1] in ID2total_reads_cl:
        ID2total_reads_cl[arr[1]]=arr[4]
    if not arr[2] in ID2Soft2reads_cl[arr[1]]:
        ID2Soft2reads_cl[arr[1]][arr[2]]=arr[3]
    
    line=fi_clinical.readline()

ID2So2per={}
for key in ID2Soft2reads_cl:
    for i in ID2Soft2reads_cl[key]:
        per=float(ID2Soft2reads_cl[key][i])/float(ID2total_reads_cl[key])
        #print(key+'\t'+i+'\t'+ID2Soft2reads_cl[key][i]+'\t'+str(per))
        #print(ID2total_reads_cl[key])
        ID2So2per.setdefault(key,{}).setdefault(i,per)


#open the baseline amplicons 
fi_base=open(sys.argv[2])
Amp2So2Sam2read={}
Amp2So2Sam2total={}
line_b=fi_base.readline()
line_b=fi_base.readline()
while(line_b):
    line_b=line_b.rstrip()
    arr_b=line_b.split('\t')
    if not arr_b[1] in Amp2So2Sam2read:
        Amp2So2Sam2read.setdefault(arr_b[1],{})
        Amp2So2Sam2total.setdefault(arr_b[1],{})
    if not arr_b[2] in Amp2So2Sam2read[arr_b[1]]:
        Amp2So2Sam2read[arr_b[1]].setdefault(arr_b[2],{})
        Amp2So2Sam2total[arr_b[1]].setdefault(arr_b[2],{})
    if not arr_b[0] in Amp2So2Sam2read[arr_b[1]][arr_b[2]]:
        Amp2So2Sam2read[arr_b[1]][arr_b[2]].setdefault(arr_b[0],float(arr_b[3])/float(arr_b[4]))
        Amp2So2Sam2total[arr_b[1]][arr_b[2]].setdefault(arr_b[0],float(arr_b[3]))
    line_b=fi_base.readline()


Amp2So2Mean={}
Amp2So2SD={}
for amp in Amp2So2Sam2read:
    Amp2So2Mean.setdefault(amp,{})
    Amp2So2SD.setdefault(amp,{})
    for so in Amp2So2Sam2read[amp]:
        li=[]
        for sam in Amp2So2Sam2read[amp][so]:
            li.append(Amp2So2Sam2read[amp][so][sam])
        Mean=np.mean(li)
        #print(Mean)
        SD=np.std(li)
        Amp2So2Mean[amp].setdefault(so,Mean)
        Amp2So2SD[amp].setdefault(so,SD)

Amp2So2Zscr={}
for idd in ID2So2per:
    Amp2So2Zscr.setdefault(idd,{})
    if idd in Amp2So2SD:
        for so in ID2So2per[idd]:
            if(Amp2So2SD[idd][so]==0):
                Amp2So2Zscr[idd][so]=int(ID2Soft2reads_cl[idd][so])
                #print(idd+'\t'+so)
            else:
                Z_scr=float(ID2So2per[idd][so]-Amp2So2Mean[idd][so])/float(Amp2So2SD[idd][so])
                Amp2So2Zscr[idd][so]=Z_scr
                
for key in Amp2So2Zscr:
    for key2 in Amp2So2Zscr[key]:
        if(Amp2So2Zscr[key][key2]>z_cutoff):
            print(key+'\t'+key2+'\t'+str(Amp2So2Zscr[key][key2]))







