#!/usr/bin/python3
from operator import itemgetter  
import sys
import re

Aml2Num={}
Aml2So={}
f=sys.argv[1]
fi=open(f)
line=fi.readline()
while(line):
    arr=line.rstrip().split('\t')
    ID=arr[2]
    Soft=arr[1]
    check=re.sub(r'\d','',Soft)
    string=re.sub(r'\D','_',Soft)
    num=string.split('_')

    if(check[0]=='S' or  check[-1]=='S'):
        if not ID in Aml2Num:
            Aml2Num[ID]=0
        else:
            Aml2Num[ID]+=1
        if(check[0]=='S' and  check[-1]=='S'):
            if(int(num[0]) >10 and int(num[0]) <80):
                if not ID in Aml2So:
                    Aml2So[ID]=0
                else:
                    Aml2So[ID]+=1
            if(int(num[-2]) >10 and int(num[-2]) <80):
                if not ID in Aml2So:
                    Aml2So[ID]=0
                else:
                    Aml2So[ID]+=1
        elif(check[0]=='S'):
            if(int(num[0]) >10 and int(num[0]) <80):
                if not ID in Aml2So:
                    Aml2So[ID]=0
                else:
                    Aml2So[ID]+=1
        else:
            if(int(num[-2]) >10 and int(num[-2]) <80):
                if not ID in Aml2So:
                    Aml2So[ID]=0
                else:
                    Aml2So[ID]+=1







    line=fi.readline()


for key, value in sorted(Aml2Num.items(), key = itemgetter(1), reverse = True):
    print('%s\t%s\t%.5f' %(key,value,round(float(Aml2So[key])/float(value),5)))
            
