#!/usr/bin/python
from __future__ import print_function
import sys
import re
import subprocess as sb

for k in range(1,len(sys.argv)):
    com='samtools view -h '+ sys.argv[k]
    proc=sb.Popen(args=com,shell=True,stdout=sb.PIPE,stderr=sb.PIPE)
    (stdout_get,stderr_get)=proc.communicate()

    arr=stdout_get.split('\n')
    #fow_num={}
    #rev_num={}
    dic={}
    #total_reads=0
    #both_Soreads=0
    #fow_Sofreads=0
    #rev_Sofreads=0
    all_so=0
    total_80to10=0
    for i in arr:
        if not re.match('^@',i):
            col=i.split('\t')
            if(len(col)>5):
                #total_reads +=1
                col_clip=col[5]
                #print(col_clip)
                check=re.sub(r'\d','',col_clip)
                #print(check)
                string=re.sub(r'\D','_',col_clip)
                num=string.split('_')
                if(check[0]=='S' or check[-1]=='S'):
                    all_so+=1
                    if(int(num[0])<=80 and int(num[0])>=10):
                        if num[0] in dic:
                            dic[num[0]] +=1
                        else:
                            dic[num[0]] =1

                        #for m in range(0,5):
                            #print(col[m],'\t')
                        #print('\n')
                        total_80to10 +=1
                    elif(int(num[-2])<=80 and int(num[-2])>=10):
                        #for m in range(0,5):
                            #print(col[m],'\t')
                        #print('\n')
                        if num[-2] in dic:
                            dic[num[-2]] +=1
                        else:
                            dic[num[-2]] =1

                        total_80to10+=1
                    #if num[0] in fow_num:
                        #fow_num[num[0]] +=1
                    #else:
                        #fow_num[num[0]] =1
                    #if num[-2] in rev_num:
                        #rev_num[num[-2]] +=1
                    #else:
                        #rev_num[num[-2]] =1
                    #both_Soreads +=1

                '''elif(check[0]=='S'):
                    if num[0] in fow_num:
                        fow_num[num[0]] +=1
                    else:
                        fow_num[num[0]] =1
                    fow_Sofreads +=1
                elif(check[-1]=='S'):
                    if num[-2] in rev_num:
                        rev_num[num[-2]] +=1
                    else:
                        rev_num[num[-2]] =1
                    rev_Sofreads +=1
                #print(fow_num)
                #print(rev_num)'''
    fi_na=sys.argv[k].split('/')
    #print(total_reads)
    #no_clip=total_reads-both_Soreads-fow_Sofreads-rev_Sofreads
    #print(no_clip)
    #clippp=both_Soreads+fow_Sofreads+rev_Sofreads
    #print(fi_na[-1]+'\tBoth_Socl_reads:'+str(both_Soreads)+'\tforw_Socl_reads:'+str(fow_Sofreads)+'\trev_Socl_reads:'+str(rev_Sofreads)+'\tnoSoCl:'+str(no_clip))
    #per_noClip=no_clip/total_reads
    #print('no SoftClip: %.3f' %(round(float(no_clip)*100/clippp,3)))
    #print(clippp)
    #print('\n')
    dic = {int(k):int(v) for k,v in dic.items()}
    #rev_num = {int(k):int(v) for k,v in rev_num.items()}
    print(fi_na[-1])
    print(total_80to10)
    #for k in dic:
    #    print(str(k)+'\t'+str(dic[k]))
    for key, value in sorted(dic.iteritems(),reverse=True, key=lambda (k,v): (v,k)):
        print("%s: %s\t %.3f" % (key, value,round((float(value)*100/all_so),3)))
    #print('\n')
    #long_clin sorted(fow_num.iterkeys(),reverse=True)[:5]:
    #p=0
    #                                                       
    #    print "%s: %s\t %.3f" % (key, fow_num[key],round((float(fow_num[key])*100/clippp),3))
        #if(key>100):
            #long_clip=long_clip+fow_num[key]
    #print('fow >100:'+str(round(long_clip/clippp,2)))
    #print(long_clip)
    #print('='*50)
    #for key, value in sorted(rev_num.iteritems(),reverse=True, key=lambda (k,v): (v,k))[:5]:
    #    print "%s: %s\t %.3f" % (key, value,(round(float(value)*100/clippp,3)))
    #print('\n')
    #long_clip_2=0
    #for key in sorted(rev_num.iterkeys(),reverse=True)[:5]:
    #    print "%s: %s\t %.3f" % (key, rev_num[key],round((float(rev_num[key])*100/clippp),3))
        #if(key>100):
            #long_clip_2=long_clip_2+rev_num[key]

    #print('rev >100:'+str(round(long_clip_2/clippp,2)))
    #print(long_clip_2)
    #print('fow\n')
    #for key, value in sorted(fow_num.iteritems(),reverse=True, key=lambda (k,v): (v,k))[:80]:
    #    print "%s: %s\t %.3f" % (key, value,round((float(value)*100/clippp),3))
    #print('rev')
    #for key, value in sorted(rev_num.iteritems(),reverse=True, key=lambda (k,v): (v,k))[:80]:
    #    print "%s: %s\t %.3f" % (key, value,(round(float(value)*100/clippp,3)))
#'''

