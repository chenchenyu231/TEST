#!/usr/bin/python3

import subprocess as sb
import sys
import re
import pandas as pd
import os
import pybedtools as bt
import pysam


fi_bed=sys.argv[1]
fi_bam=sys.argv[2]
pre_bam=re.search(r'(\S+)\.bam',os.path.basename(fi_bam)).groups()[0]


bam_bt=bt.BedTool(fi_bam)
bed_bt=bt.BedTool(fi_bed)
bam_intersect=bam_bt.intersect(bed_bt,bed=True,wb=True)

Id2Cig={}
for k_line in bam_bt:
    k_line=str(k_line)
    col=k_line.rstrip().split('\t')
    #print(col[0]+'\t'+col[5])
    Id2Cig.setdefault(col[0],col[5])
    #bam_bt.head()
    #bed_bt.head()
    #bam_intersect.head()

read_id_with_cov = {}
target_id_dict = {}
count=0
for in_line in bam_intersect:
    in_line = str(in_line)
    fields = in_line.rstrip("\n").split("\t")
    read_chr   = fields[0]
    read_start = int(fields[1])
    read_end   = int(fields[2])
    read_id    = fields[3]

    target_chr   = fields[12]
    target_start = int(fields[13])
    target_end   = int(fields[14]) 
    target_id    = fields[15]
    target_str   = target_chr+":"+str(target_start)+"-"+str(target_end)+";"+str(target_id)

    covearge = float(read_end-read_start)/float(target_end-target_start)
    if not read_id in read_id_with_cov:
        read_id_with_cov.setdefault(read_id, {}).setdefault("targetID", target_str)
        read_id_with_cov.setdefault(read_id, {}).setdefault("Coverage", covearge)
        target_id_dict.setdefault(target_str, None)

    else:
        if float(covearge) > float(read_id_with_cov[read_id]["Coverage"]):
            read_id_with_cov[read_id]["Coverage"] = float(covearge)
            read_id_with_cov[read_id]["targetID"] = target_str
            target_id_dict.setdefault(target_str, None)
    count+=1

for key in read_id_with_cov.keys():
    ID=read_id_with_cov[key]['targetID'].split(';')
    print(key+'\t'+Id2Cig[key]+'\t'+ID[1])

