#!/usr/bin/python
from __future__ import print_function
import subprocess as sb
import sys
import re
import os
import argparse




com='samtools view -h '+sys.argv[1]
#print(com)
proc = sb.Popen(args=com, shell=True, stdout=sb.PIPE, stderr=sb.PIPE)
(stdout_get, stderr_get) = proc.communicate()
arr=str(stdout_get).split("\n")
for i in arr:
    if not re.match('^@',i):
        col=i.split('\t')
        if(len(col)>5):
            for k in range(0,12):
                print(col[k],end='\t')
            print('\n',end='')


