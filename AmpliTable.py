#!/usr/bin/python3
import re
import sys
from operator import itemgetter
import os

fo=open('AmpliTable_Five.txt','w')
fo.write('Sample\tReadID\tSoft_len\tReads\tdepth\n')

for ind_fi in range(1,len(sys.argv)):
    #print(ind_fi)
    bam_file_prefix = os.path.basename(sys.argv[ind_fi])
    fi_ampli=open(sys.argv[ind_fi])
    print(sys.argv[ind_fi])
    Id2InfoNum={}
    Id2dep={}
    line=fi_ampli.readline()
    while(line):
        arr=line.rstrip().split('\t')
    
        if not arr[2] in Id2dep:
            Id2dep.setdefault(arr[2],1)
        else:
            Id2dep[arr[2]]+=1
        ID=arr[2]
        Info=arr[1]
        find=re.findall(r'(\d+)S+',Info)
        #print(find)
        if not ID in Id2InfoNum:
            Id2InfoNum.setdefault(ID,{})
                                                
        if(len(find)>0):
            for num in find:
                if not num in Id2InfoNum[ID]:
                    Id2InfoNum[ID].setdefault(num,1)
                else:
                    Id2InfoNum[ID][num]+=1
                                                                                                                                        
        line=fi_ampli.readline()
    print('fin read')
    for key2 in Id2InfoNum:
        Id2InfoNum[key2] = {int(k):int(v) for k,v in Id2InfoNum[key2].items()}
        for num in range(1,500):
            if num not in Id2InfoNum[key2]:
                Id2InfoNum[key2].setdefault(num,0)
        #print(key2,Id2InfoNum[key2])
        count_min=0
        count_large=0
        count=0
        for key,value in sorted(Id2InfoNum[key2].items()):
            if(key<=10):
                count_min+=value
            if(key==10):
                fo.write(bam_file_prefix+'\t'+key2+'\t<10S\t'+str(count_min)+'\t'+str(Id2dep[key2])+'\n')
            if(key>10 and key<=80):
                count+=value
                if(key%5==0):
                    fo.write(bam_file_prefix+'\t'+key2+'\t'+str(key-4)+'-'+str(key)+'S'+'\t'+str(count)+'\t'+str(Id2dep[key2])+'\n')
                    count=0
            if(key>80):
                count_large+=value
            if(key==499):
                fo.write(bam_file_prefix+'\t'+key2+'\t>80S\t'+str(count_large)+'\t'+str(Id2dep[key2])+'\n')

    ind_fi+=1
