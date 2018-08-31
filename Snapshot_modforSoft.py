#!/usr/bin/python3 
import sys
import argparse
import re
import os
import numpy as np
import pandas as pd

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-b",'--bed',dest='bed',help="INPUT MERGED BED FILE")
    parser.add_argument("file", nargs='+',help="Files to get snapshot(bam)")
    parser.add_argument("-dir","--directory",dest="dir",default=".",help="Directory to store snapshot")
    parser.add_argument("-l","--the chr locus",dest="loci",help="chromosome loci taking shot")
    parser.add_argument("-S","--for softclip",dest="soft",help="for softclip check")
    parser.add_argument("-v","--input variant list",dest="vli",help="input the variant list")
    parser.add_argument("-d","--num of downsample reads ",dest="downsample",help="num of downsample reads(defalut 3000)")
    args=parser.parse_args()
if not args.bed:
    sys.exit('ERROR: Please input MERGED BED FILE')
else:
    if not re.match('(.+)bed$',args.bed):
        sys.exit('ERROR: Please input MERGED BED FILE')
if not args.loci and not args.vli and not args.soft:
    sys.exit('ERROR: Please give a chr locus or input variant list file or softclip check')
#reads num(defalut 3000)
if(args.downsample):
    reads=args.downsample
else:
    reads=3000
#Parameters

#igv_session.xml dir
sess_dir='./igv_session.xml'
#Defalut panel Height
df_PanHeight=int(reads)+1000
#Defalut multi-amplicon panel height
df_MulPanHeight=int(reads)*2+2000
#Default multi downsample 
df_MulDownsample=int(reads)*2

#recreate the origin session file
fi_ses=open(sess_dir)
fo_new=open('new_sesion.xml','w')
line_ses=fi_ses.readline()
count=1
tab=' '*4
name=args.bed
spi_nam=name.split('/')
while(line_ses):
    if(count==3):
        fo_new.write(tab+'<Resources>\n')
        fo_new.write(tab*2+'<Resource path="'+name+'"/>\n')
        fo_new.write(tab+'</Resources>\n')
    if(count==7):
        fo_new.write(tab*2+'<Track altColor="0,0,178" autoScale="false" clazz="org.broad.igv.track.FeatureTrack" color="0,0,178" colorScale="ContinuousColorScale;0.0;211.0;255,255,255;0,0,178" displayMode="EXPAND" featureVisibilityWindow="-1" fontSize="10" id="'+name+'" name="'+spi_nam[-1]+'" renderer="BASIC_FEATURE" sortable="false" visible="true" windowFunction="count">\n')
        fo_new.write(tab*3+'<DataRange baseline="0.0" drawBaseline="true" flipAxis="false" maximum="211.0" minimum="0.0" type="LINEAR"/>\n')
        fo_new.write(tab*2+'</Track>\n')
    fo_new.write(line_ses)
    line_ses=fi_ses.readline()
    count+=1
fi_ses.close()
fo_new.close()



#determine the Downsample num(depends on the number of amplicon reads in the region
#trans bed file to array


fi_bed=open(args.bed)
line_bed=fi_bed.readline()
col=line_bed.split('\t')
num=col[3].split('&')
cho=col[0]
arr_bed=[[int(cho[3:]),int(col[1]),int(col[2]),len(num)]]
line_bed=fi_bed.readline()

while(line_bed):
    col=line_bed.split('\t')
    num=col[3].split('&')
    cho=col[0]
    arr_bed=np.append(arr_bed,[[cho[3:],int(col[1]),int(col[2]),len(num)]],axis=0)
    line_bed=fi_bed.readline()
    fr=pd.DataFrame(arr_bed)
    df=fr.rename(columns={0: "chr", 1: "str",2:'end',3:'num'})
    df[['str','end','num']]=df[['str','end','num']].astype(dtype='int32')
print('Finish Bed Process')

#writing batch file
fo=open('igv_batchtry.txt','w')
fo.write("new\n")
fo.write("genome  http://s3.amazonaws.com/igv.broadinstitute.org/genomes/hg19.genome\n")
fo.write("maxPanelHeight ")
fo.write(str(df_PanHeight)+"\n")
fo.write("load ")
fo.write("new_sesion.xml\n")
#set preference here 
fo.write("preference SAM.SHOW_CENTER_LINE=true\n")
fo.write("preference SAM.GROUP_OPTION=STRAND\n")
fo.write("preference SORT_OPTION=BASE\n")
fo.write("preference SAM.COLOR_BY=READ_STRAND\n")
fo.write('load ')
fo.write('/home/chenyu/Igv_pra/Bed_file/4Pool_ONCOv2_Designed_plus_PGx_20170704.bed\n')
#load bed
#if(args.bed):
#    fo.write("load ")
#    fo.write(args.bed+'\n')
#load bam
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
    #ran control bps of photo
    print('in Snapshot')
    if(ran==0):
        chr_range=[200,10]
    elif(ran==10):
        chr_range=[200,10]
    elif(ran==20):
        chr_range=[200,20]
    elif(ran==40):
        chr_range=[200,40]
    elif(ran==100):
        chr_range=[200,300]
    #go to loci 
    #range change(for softclip)
    for i in chr_range:
        fo.write("goto chr"+chor)
        fo.write(":")
        ran_min=int(chr_loci)-50
        fo.write(str(int(ran_min)))
        fo.write("-")
        ran_max=int(chr_loci)+int(i)
        fo.write(str(int(ran_max))+"\n")
        #sort,squish when change var 
        if(flag==1):
            fo.write("sort base\n")
            fo.write("squish\n")
            flag=0
        fo.write("snapshot ")
        #snapshot file name
        nam=args.file[0].split("/")
        nam=nam[-1].split(".")
        nam=nam[0]+"_"+str(chor)+":"+str(chr_loci)+"_"+mut
        fo.write(nam)
        fo.write("_r")
        fo.write(str(i))
        fo.write(".png\n")

command='xvfb-run --server-args="-screen 0 1920x1080x24" java -Xmx20G -jar  ~/tools/IGV_2.4.13/igv.jar -b igv_batchtry.txt'
#devide multi amplicon and mono 
ampli_multi=np.array([])
ampli_one=np.array([])

#for softclip
if(args.soft):
    print('Soft mode')
    fi_so=open(args.soft)
    line_so=fi_so.readline()
    while(line_so):
        arr_so=line_so.rstrip().split('\t')
        col_cho=arr_so[1]
        #print(arr_so)
        chor=col_cho[3:]
        chr_loci=arr_so[2]
        locii=df['chr']==chor
        #print(chor+'\t'+chr_loci)
        strr=df['str']<=int(chr_loci)
        endd=df['end']>=int(chr_loci)
        df2=df[(locii & strr & endd)]
        if(df2.empty==True):
            print('ERROR:chromosome loci is not on the design panel')
            print(chor+'\t'+chr_loci)
        amp=int(df2['num'])
        #print(amp)
        if(amp==1):
            if(len(ampli_one)==0):
                ampli_one=np.array([[chor,chr_loci,'',1,100]])
            else:
                ampli_one=np.append(ampli_one,[[chor,chr_loci,'',1,100]],axis=0)
        else:
            if(len(ampli_multi)==0):
                ampli_multi=np.array([[chor,chr_loci,'',1,100]])
            else:
                ampli_multi=np.append(ampli_multi,[[chor,chr_loci,'',1,100]],axis=0)
        line_so=fi_so.readline()



#detect which format of loci you input
if(args.loci):
    print('loci mode')
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
    #set and find the amplicon in list
    locii=df['chr']==chor
    strr=df['str']<=int(chr_loci)
    endd=df['end']>=int(chr_loci)
    df2=df[(locii & strr & endd)]
    if(df2.empty==True):
        print('ERROR:chromosome loci is not on the design panel')
    amp=int(df2['num'])
    if(amp==1):
        if(len(ampli_one)==0):
            ampli_one=np.array([[chor,chr_loci,'',1,0]])
        else:
            ampli_one=np.append(ampli_one,[[chor,chr_loci,'',1,0]],axis=0)         
    else:
        if(len(ampli_multi)==0):
            ampli_multi=np.array([[chor,chr_loci,'',1,0]])
        else:
            ampli_multi=np.append(ampli_multi,[[chor,chr_loci,'',1,0]],axis=0)
            
var2run={}
#input Variant list  file
if(args.vli):
    print('vli mode')
    fi_vaf=open(args.vli)
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
        chr_loci=tmp[2]
        mut=tmp[3].replace('>','=')
        locii=df['chr']==chor
        strr=df['str']<=int(chr_loci)
        endd=df['end']>=int(chr_loci)
        df2=df[(locii & strr & endd)]
        if(df2.empty==True):
            print('ERROR:chrosome loci is not on the design panel')
        if(len(df2)>1):
            sys.exit('ERROR:the Bed file you input is Should BE A MERGED BED FILE')
        amp=int(df2['num'])
        #mut detect mutation len 
        if(len(mut)>=10 and len(mut)<20):
            if(amp==1):
                if(len(ampli_one)==0):
                    ampli_one=np.array([[chor,chr_loci,mut,1,20]])
                else:
                    ampli_one=np.append(ampli_one,[[chor,chr_loci,mut,1,20]],axis=0)
            else:
                if(len(ampli_multi)==0):
                    ampli_multi=np.array([[chor,chr_loci,mut,1,20]])
                else:
                    ampli_multi=np.append(ampli_multi,[[chor,chr_loci,mut,1,20]],axis=0)
        elif(len(mut)>=20):
            if(amp==1):
                if(len(ampli_one)==0):
                    ampli_one=np.array([[chor,chr_loci,mut,1,40]])
                else:
                    ampli_one=np.append(ampli_one,[[chor,chr_loci,mut,1,40]],axis=0)
            else:
                if(len(ampli_multi)==0):
                    ampli_multi=np.array([[chor,chr_loci,mut,1,40]])
                else:
                    ampli_multi=np.append(ampli_multi,[[chor,chr_loci,mut,1,40]],axis=0)
        else:
            if(amp==1):
                if(len(ampli_one)==0):
                    ampli_one=np.array([[chor,chr_loci,mut,1,10]])
                else:
                    ampli_one=np.append(ampli_one,[[chor,chr_loci,mut,1,10]],axis=0)
            else:
                if(len(ampli_multi)==0):
                    ampli_multi=np.array([[chor,chr_loci,mut,1,10]])
                else:
                    ampli_multi=np.append(ampli_multi,[[chor,chr_loci,mut,1,10]],axis=0)

df_one=pd.DataFrame(ampli_one)
df_multi=pd.DataFrame(ampli_multi)
#set downsample
if not len(ampli_one)==0:
    fo.write('preference SAM.MAX_LEVELS ')
    fo.write(str(reads)+'\n')
    for i in range(0,len(ampli_one)):
        chor=df_one.loc[i,0]
        chr_loci=df_one.loc[i,1]
        mut=df_one.loc[i,2]
        fla=df_one.loc[i,3]
        num=df_one.loc[i,4]
        Snapshot(chor,int(chr_loci),mut,int(fla),int(num))

if not len(ampli_multi)==0:
#multi 
    fo.write('maxPanelHeight ')

    fo.write(str(df_MulPanHeight)+"\n")

    fo.write('preference SAM.MAX_LEVELS ')
    fo.write(str(df_MulDownsample)+'\n')
    fo.write('sort base\n')

    for i in range(0,len(ampli_multi)):
        chor=df_multi.loc[i,0]
        chr_loci=df_multi.loc[i,1]
        mut=df_multi.loc[i,2]
        fla=df_multi.loc[i,3]
        num=df_multi.loc[i,4]
        Snapshot(chor,int(chr_loci),mut,int(fla),int(num))
                        

fo.write("exit")
fo.close()


os.system(command)

os.system('rm ./new_sesion.xml')
