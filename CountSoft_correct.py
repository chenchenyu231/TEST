#!/usr/bin/python
from __future__ import print_function
import sys
import re
import subprocess as sb


com='samtools view -h '+ sys.argv[1]
proc=sb.Popen(args=com,shell=True,stdout=sb.PIPE,stderr=sb.PIPE)
(stdout_get,stderr_get)=proc.communicate()

arr=stdout_get.split('\n')
total=0
dic={}
for i in arr:
    if not re.match('^@',i):
        col=i.split('\t')
        if(len(col)>5):
            total+=1
            col_clip=col[5]
            check=re.sub(r'\d','',col_clip)
            string=re.sub(r'\D','_',col_clip)
            num=string.split('_')
            if(check[0]=='S' or check[-1]=='S'):
                if(num[0] in dic):
                    dic[num[0]]+=1
                else:
                    dic[num[0]]=1
                if(num[-2] in dic):
                    dic[num[-2]]+=1
                else:
                    dic[num[-2]]=1


fi_na=sys.argv[1].split('/')
dic = {int(k):int(v) for k,v in dic.items()}

print(fi_na[-1])
#for key, value in sorted(dic.iteritems(),reverse=True, key=lambda (k,v): (v,k)):
#    print("%s: %s\t %.3f" % (key, value,round((float(value)*100)/total/2,3)))

for key in sorted(dic.iterkeys()):
    print( "%s: %s\t %.5f" % (key, dic[key],round((float(dic[key])*100)/total/2,6)))


