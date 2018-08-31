#!/usr/bin/python3 -w
import sys
import argparse
import re
import os
if __name__ == "__main__":
        parser = argparse.ArgumentParser()
        parser.add_argument("file", nargs='+',help="Files to get snapshot")
        parser.add_argument("-o_dir","--directory",dest="dir",default=".",help="Directory to store snapshot")
        parser.add_argument("-loci","--the chr locus",dest="loci",help="chromosome loci taking shot")
        parser.add_argument("-v","--input variant file",dest="vcf",help="input the variant ID file")
        parser.add_argument("-d","--num of downsample reads ",dest="downsample",help="num of downsample reads")
        args=parser.parse_args()

#reads num
if(args.downsample):
    reads=args.downsample
else:
    reads=3000
#give path to prefs.properties
fi_pro=open("/home/chenyu/igv/prefs.properties")
fo_new=open("/home/chenyu/igv/tmp",'w')
line_pro=fi_pro.readline()
while(line_pro):
    if(re.match('^SAM.MAX_LEVELS',line_pro)):
        fo_new.write("SAM.MAX_LEVELS=")
        fo_new.write(str(reads))
        fo_new.write("\n")
    elif(re.match('^SAM.GROUP_OPTION',line_pro)):
        fo_new.write('SAM.GROUP_OPTION=STRAND\n')
    elif(re.match('^SAM.COLOR_BY',line_pro)):
        fo_new.write('SAM.COLOR_BY=READ_STRAND\n')
    elif(re.match('^SAM.SORT_OPTION',line_pro)):
        fo_new.write('SAM.SORT_OPTION=STRAND\n')
    else:
        fo_new.write(line_pro)
    line_pro=fi_pro.readline()
fi_pro.close()
fo_new.close()
Comma='mv /home/chenyu/igv/prefs.properties /home/chenyu/igv/prefs.properties.1'
Commc='rm /home/chenyu/igv/prefs.properties.1'
Commb='mv /home/chenyu/igv/tmp /home/chenyu/igv/prefs.properties'

os.system(Comma)

os.system(Commc)
os.system(Commb)


def Snapshot(fi,chor,chr_loci,mut='',di='.'):
    fo=open('igv_batchtry.txt','w')
    fo.write("new\n")
    fo.write("genome  http://s3.amazonaws.com/igv.broadinstitute.org/genomes/hg19.genome\n")
    fo.write("maxPanelHeight 8000\n")
    for i in fi:
        fo.write("load ")
        fo.write(i)
        fo.write("\n")
    fo.write("snapshotDirectory ")
    fo.write(di+"\n")
    chr_range=[200,40,20,10]
    for i in chr_range:
        fo.write("goto chr"+chor)
        fo.write(":")
        ran_min=int(chr_loci)-int(i)/2
        fo.write(str(int(ran_min)))
        fo.write("-")
        ran_max=int(chr_loci)+int(i)/2
        #print("ran_max"+str(int(ran_max)))
        fo.write(str(int(ran_max))+"\n")
        fo.write("sort strand\n")
        fo.write("squish\n")
        fo.write("snapshot ")
        nam=args.file[0].split("/")
        nam=nam[-1].split(".")
        nam=nam[0]+str(chor)+":"+str(chr_loci)+"_"+mut
        fo.write(nam)
        fo.write("_r")
        fo.write(str(i))
        fo.write(".png\n")
    fo.write("\nexit")

command='xvfb-run --server-args="-screen 0 1920x1080x24" java -Xmx20G -jar  ~/tools/IGV_2.4.13/igv.jar -b igv_batchtry.txt'

if(args.loci):
    if(re.match('^chr',args.loci)):
        if(args.loci[5]==":"):
            chor=args.loci[3:5]
            chr_loci=args.loci[6:]
        else:
            chor=args.loci[3:4]
            chr_loci=args.loci[5:]
    else:
        if(args.loci[2]==":"):
            chor=args.loci[0:2]
            chr_loci=args.loci[3:]
        else:
            chor=args.loci[0:1]
            chr_loci=args.loci[2:]
    if(args.dir):
        Snapshot(args.file,chor,chr_loci,'',args.dir)
    else:
        Snapshot(args.file,chor,chr_loci)
    os.system(command)

var2run={}
if(args.vcf):
    fi_vaf=open(args.vcf)
    line_vaf=fi_vaf.readline()
    arr_head=line_vaf.split("\t")
    for i in range(0,len(arr_head)):
        arr_head[i]=arr_head[i].lower()
        if(arr_head[i]=='variant_id'):
            col_var=i
        if(arr_head[i]=='run_name'):
            col_run=i
    line_vaf=fi_vaf.readline()
    while(line_vaf):
        arr=line_vaf.split("\t")
        var=arr[col_var]
        run=arr[col_run]
        var2run[var]=run
        line_vaf=fi_vaf.readline()
    
    for i in var2run.keys():
        tmp=i.split("_")
        chor=tmp[1]
        chor=chor[3:]
        print("chor:"+chor)
        chr_loci=tmp[2]
        print("chr_loci:"+chr_loci)
        mut=tmp[3].replace('>','=')
        if(args.dir):
            Snapshot(args.file,chor,chr_loci,mut,args.dir)
        else:
            Snapshot(args.file,chor,chr_loci,mut)
        os.system(command)








