#!/usr/bin/python3 
import sys
import argparse
import re
import os

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("file", nargs='+',help="Files to get snapshot")
    parser.add_argument("-dir","--directory",dest="dir",default=".",help="Directory to store snapshot") 
    parser.add_argument("-l","--the chr locus",dest="loci",help="chromosome loci taking shot")
    parser.add_argument("-v","--input variant file",dest="vcf",help="input the variant ID file")
    parser.add_argument("-b","--input bed file",dest="bed",help="input the bed file")
    parser.add_argument("-d","--num of downsample reads ",dest="downsample",help="num of downsample reads")    
    args=parser.parse_args()

#reads num
if(args.downsample):
    reads=args.downsample
else:
    reads=3000
#give path to prefs.properties
#remember to go to igv preference check 'show center line' or go to prefs.properties write :'SAM.SHOW_CENTER_LINE=true' 
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


#determine the maxPanelHeight(depends on the number of amplicon reads in the region

fo=open('igv_batchtry.txt','w')
fo.write("new\n")
fo.write("genome  http://s3.amazonaws.com/igv.broadinstitute.org/genomes/hg19.genome\n")
fo.write("maxPanelHeight ")
height=int(reads)*2
fo.write(str(height)+"\n")
if(args.bed):
    fo.write("load ")
    fo.write(args.bed+'\n')
for i in args.file:
    fo.write("load ")
    fo.write(i)
    fo.write("\n")
fo.write("snapshotDirectory ")
if(args.dir):
    fo.write(args.dir+"\n")
else:
    fo.write(".\n")
    


def Snapshot(chor,chr_loci,mut='',flag=0,ran=0):
    if(ran==0):
        chr_range=[200,40,20,10]
    elif(ran==10):
        chr_range=[200,10]
    elif(ran==20):
        chr_range=[200,20]
    elif(ran==40):
        chr_range=[200,40]
    for i in chr_range:
        fo.write("goto chr"+chor)
        fo.write(":")
        ran_min=int(chr_loci)-int(i)/2
        fo.write(str(int(ran_min)))
        fo.write("-")
        ran_max=int(chr_loci)+int(i)/2
        #print("ran_max"+str(int(ran_max)))
        fo.write(str(int(ran_max))+"\n")
        if(flag==1):
            fo.write("sort strand\n")
            fo.write("squish\n")
            flag=0
        fo.write("snapshot ")
        nam=args.file[0].split("/")
        nam=nam[-1].split(".")
        nam=nam[0]+str(chor)+":"+str(chr_loci)+"_"+mut
        fo.write(nam)
        fo.write("_r")
        fo.write(str(i))
        fo.write(".png\n")




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
    Snapshot(chor,chr_loci,'',1)



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
        #print("chor:"+chor)
        chr_loci=tmp[2]
        #print("chr_loci:"+chr_loci)
        mut=tmp[3].replace('>','=')
        if(len(mut)>=10 and len(mut)<20):
            Snapshot(chor,chr_loci,mut,1,20)
        elif(len(mut)>=20):
            Snapshot(chor,chr_loci,mut,1,40)
        else:
            Snapshot(chor,chr_loci,mut,1,10)


        
        

fo.write("\nexit")
fo.close()

os.system(command)





