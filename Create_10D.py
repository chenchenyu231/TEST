import re
import subprocess as sb
import sys
import re
import pandas as pd
import os
import pybedtools as bt
import pysam


fi=open('/home/chenyu/Softclip/B_589_66_new.sam')
li=fi.readline()
dic={}
while(li):
    arr=li.rstrip().split('\t')
    col_fea=arr[5]
    #print(col_fea,end='')
    find=re.findall(r'(\d+)I+',col_fea)
    find = [int(i) for i in find]
    find.sort()
    #print(find,end='')
    if(len(find)>0):
        if(find[-1]>=15):
            if not arr[3] in dic:
                dic.setdefault(arr[3], {}).setdefault("loci", arr[2])
                print(arr[0]+'\t'+arr[2]+'\t'+arr[3]+'\t'+arr[5])
                                                                                                        
    li=fi.readline()
    #print(dic.keys)

