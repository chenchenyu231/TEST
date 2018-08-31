#!/usr/bin/python3

from  scipy import stats
import os
import sys
import re
import numpy as np


cut=int(sys.argv[3])

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
#for key2 in ID2So2per:
#    for k in ID2So2per[key2]:
#        print(ID2So2per[key2][k])



fi_base=open(sys.argv[2])
Amp2So2Sam2read={}
#Amp2So2Sam2total={}
line_b=fi_base.readline()
line_b=fi_base.readline()
while(line_b):
    #print(line_b,end='')
    line_b=line_b.rstrip()
    arr_b=line_b.split('\t')
    if not arr_b[1] in Amp2So2Sam2read:
        Amp2So2Sam2read.setdefault(arr_b[1],{})
        #Amp2So2Sam2total.setdefault(arr_b[1],{})
    if not arr_b[2] in Amp2So2Sam2read[arr_b[1]]:
        Amp2So2Sam2read[arr_b[1]].setdefault(arr_b[2],{})
        #Amp2So2Sam2total[arr_b[1]].setdefault(arr_b[2],{})
    if not arr_b[0] in Amp2So2Sam2read[arr_b[1]][arr_b[2]]:
        Amp2So2Sam2read[arr_b[1]][arr_b[2]].setdefault(arr_b[0],float(arr_b[3])/float(arr_b[4]))
                                         
    line_b=fi_base.readline()
#print(Amp2So2Sam2read['GENEID_POLE_Pool_3_ID_AMPL7154392346']['26-30S']['A_424_56_ampli'])


Amp2So2Median={}
Amp2So2Mad={}
for amp in Amp2So2Sam2read:
    Amp2So2Median.setdefault(amp,{})
    Amp2So2Mad.setdefault(amp,{})
    for so in Amp2So2Sam2read[amp]:
        #print(Amp2So2Sam2read[amp][so])
        li=[]
        for sam in Amp2So2Sam2read[amp][so]:
            li.append(Amp2So2Sam2read[amp][so][sam])
        #print(li)
        #print(len(li))
        Med=np.median(li)
        #print(Med)
        li2=[]
        for m in li:
            li2.append(abs(m-Med))
        Mad=np.median(li2)
        #print(li2)
        Amp2So2Median[amp].setdefault(so,Med)
        Amp2So2Mad[amp].setdefault(so,Mad)


Amp2So2MADScr={}
for idd in ID2So2per:
    Amp2So2MADScr.setdefault(idd,{})
    if idd in Amp2So2Median:
        for so in ID2So2per[idd]:
            if(Amp2So2Mad[idd][so]==0):
                Amp2So2MADScr[idd][so]=int(ID2Soft2reads_cl[idd][so])
            else:
                MADscr=float(ID2So2per[idd][so]-Amp2So2Median[idd][so])/float(Amp2So2Mad[idd][so])
                Amp2So2MADScr[idd][so]=MADscr


for key in Amp2So2MADScr:
    for key2 in Amp2So2MADScr[key]:
        if(Amp2So2MADScr[key][key2]>cut):
            print(key+'\t'+key2+'\t'+str(Amp2So2MADScr[key][key2]))

#print(key)

for k in Amp2So2MADScr['GENEID_RB1_Pool_3_ID_AMPL7154413607']:
    print(k+'\t'+str(Amp2So2MADScr['GENEID_RB1_Pool_3_ID_AMPL7154413607'][k]))


#print(Amp2So2MADScr['GENEID_RB1_Pool_3_ID_AMPL7154413607']['26-30S'])


