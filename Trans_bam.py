#!/usr/bin/python
import subprocess as sb
import sys
import re
import os
import argparse

if __name__ == "__main__":
        parser = argparse.ArgumentParser()
        parser.add_argument("file", nargs='+',help="Bam files to trans smaller")       
        args=parser.parse_args()


for i in args.file:
    name=i.split("/")
    coma='samtools view -h '+i
    fo_name=name[-1]+"_thin.bam"
    fo=open(fo_name,"w")
    #print(coma)
    proc = sb.Popen(args=coma, shell=True, stdout=sb.PIPE, stderr=sb.PIPE)
    (stdout_get, stderr_get) = proc.communicate()
    arr=str(stdout_get).split("\n")
    #del arr[-1]
    for k in arr:
        if(re.match('^@',k)):
            fo.write(k+"\n")
        else:
            col=k.split("\t")
            if(len(col)<12):
                fo.write(k)
            else:
                for m in range(0,11):
                    fo.write(col[m]+"\t")
                if(k!=arr[-1]):
                    fo.write("\n")
    out_bam=name[-1]+"_trans.bam"

    commb='samtools view -b '+fo_name+" > "+out_bam
    commc='samtools index -b '+out_bam
    os.system(commb)
    os.system(commc)
    commd='rm '+fo_name
    #os.system(commd)


