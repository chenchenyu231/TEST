
# coding: utf-8

import pandas as pd
import re 
import pysam
import numpy as np
from scipy import stats
import subprocess
from tempfile import SpooledTemporaryFile as tempfile
from timeit import default_timer as timer
from collections import Counter
from itertools import combinations
import os
import getopt
import sys
import pybedtools as bt
import timeit

bed_file = "/home/chenyu/Igv_pra/Bed_file/4Pool_ONCOv2_Designed_plus_PGx_20170704.bed"
in_file = sys.argv[1]
bam_file_prefix = re.search(r'(\S+)\.bam',os.path.basename(in_file)).groups()[0]

def get_bedpe(target_id_in):
    
    (target_chr, target_start, target_end ) = re.search(r'^(chr\S+?):(\d+?)-(\d+);',target_id_in).groups()
    five_primer_start = int(target_start)-26
    five_primer_end = int(target_start)-1
    three_primer_start = int(target_end)
    three_primer_end = int(target_end)+25
    out =  "\t".join([target_chr,str(five_primer_start),str(five_primer_end),target_chr,str(three_primer_start),str(three_primer_end)])+"\n"
    return out


bam_bt = bt.BedTool(in_file)
target_bt = bt.BedTool(bed_file)
bam_target_intersect = bam_bt.intersect(target_bt,bed=True,wb=True)


timer_start = timeit.default_timer()
read_id_with_cov = {}
target_id_dict = {}
count = 0
for in_line in bam_target_intersect:
    in_line = str(in_line)
    count += 1
    if (count % 5000) == 0:
        timer_end = timeit.default_timer()
  
    
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
    #print 
    if not read_id in read_id_with_cov:
        read_id_with_cov.setdefault(read_id, {}).setdefault("targetID", target_str)
        read_id_with_cov.setdefault(read_id, {}).setdefault("Coverage", covearge)
        target_id_dict.setdefault(target_str, None)
    else:
        if float(covearge) > float(read_id_with_cov[read_id]["Coverage"]):
            read_id_with_cov[read_id]["Coverage"] = float(covearge)
            read_id_with_cov[read_id]["targetID"] = target_str
            target_id_dict.setdefault(target_str, None)
            
timer_end = timeit.default_timer()


timer_start = timeit.default_timer()
samfile = pysam.AlignmentFile(in_file, "rb")
outfile = pysam.AlignmentFile(bam_file_prefix+"_add_ampID.sam", "w", template=samfile)



for read in samfile:
   
    read_name = read.query_name
    if read_name in read_id_with_cov:
        
        read.tags += [("ZL:Z:",read_id_with_cov[read_name]['targetID'])]

    else:
        read.tags += [("ZL:Z:","NA")]
	
	outfile.write(read)

samfile.close()
outfile.close()

print timer_end - timer_start      



