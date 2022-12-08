#!/usr/bin/python3
import cgi
import cgitb; cgitb.enable()  # for troubleshooting
import os
import gzip
import subprocess
import textwrap

rootdir=os.path.dirname(os.path.abspath(__file__))
bindir=rootdir+"/script"

html_header=""
html_footer=""
if os.path.isfile(rootdir+"/index.html"):
    fp=open(rootdir+"/index.html")
    txt=fp.read()
    fp.close()
    html_header=txt.split('<!-- CONTENT START -->')[0]
    html_footer=txt.split('<!-- CONTENT END -->')[-1]

def ExitWithError(msg,html_footer):
    print("ERROR!")
    print(msg)
    print("<p></p><a href=.>[Back]</a>")
    if len(html_footer):
        print(html_footer)
    else:
        print("</body> </html>")
    exit()

#### read cgi parameters ####

print("Content-type: text/html\n")
if len(html_header):
    print(html_header)
else:
    print('''<html>
<head>
<link rel="stylesheet" type="text/css" href="page.css" />
<title>BioLiP</title>
</head>
<body bgcolor="#F0FFF0">
<img src=images/BioLiP1.png ></br>
<p><a href=.>[Back to Home]</a></p>
''')

form = cgi.FieldStorage()
mode = 1
structure=form.getfirst("structure",'').strip()
if not structure:
    structure=form.getfirst("struct_file",'').strip()
if not structure:
    ExitWithError("empty input",html_footer)
if len(structure) and not "ATOM " in structure:
    if len(structure.split())!=1:
        ExitWithError("unknown input format",html_footer)
    elif ':' in structure:
        mode=3
    else:
        mode=2

outdir=rootdir
cmd="date +%N"
p=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
stdout,stderr=p.communicate()
outdir=rootdir+'/output/'+stdout.decode().strip()
if os.getenv("REMOTE_ADDR"):
    outdir+='.'+os.getenv("REMOTE_ADDR").replace(':','.')
outdir+=".fsearch"
if not os.path.isdir(outdir):
    os.makedirs(outdir)

#### download pdb if necessary ####
if mode==1:
    if "\nloop_" in structure:
        fp=open(outdir+"/input.cif",'w')
    else:
        fp=open(outdir+"/input.pdb",'w')
    fp.write(structure)
    fp.close()
elif mode==2:
    cmd="curl -s 'https://ftp.ebi.ac.uk/pub/databases/alphafold/'|grep -ohP '\\bv\\d+\\b'|sort -n|tail -1"
    p=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
    stdout,stderr=p.communicate()
    version=stdout.decode().strip()
    if len(version)==0:
        version="v4"

    http_link="https://alphafold.ebi.ac.uk/files/AF-"+structure+"-F1-model_"+version+".pdb"
    cmd="curl -s '"+http_link+"' -o "+outdir+"/input.pdb"
    os.system(cmd)
    if not os.path.isfile(outdir+"/input.pdb"):
        ExitWithError("Cannot download "+http_link,html_footer)
elif mode==3:
    pdbid,asym_id=structure.split(':')
    http_link="https://files.rcsb.org/download/"+pdbid+".cif.gz"
    cmd="curl '"+http_link+"' -o "+outdir+"/input.cif.gz"
    os.system(cmd)
    if not os.path.isfile(outdir+"/input.cif.gz"):
        ExitWithError("Cannot download "+http_link,html_footer)
    cmd="cd %s; %s/cif2chain input.cif.gz input.pdb %s"%(outdir,bindir,asym_id)
    os.system(cmd)
    if not os.path.isfile(outdir+"/input.pdb"):
        ExitWithError("Cannot extract chain %s from <a href=%s/input.cif.gz>%s</a>"%(asym_id,outdir,pdbid),html_footer)


#### run foldseek ####
infile="input.pdb"
if mode==1 and os.path.isfile(outdir+"/input.cif") and \
           not os.path.isfile(outdir+"/input.pdb"):
    infile="input.cif"
cmd="cd %s; %s/foldseek easy-search %s %s/foldseek/receptor_DB aln.m8 tmpFolder --threads 1 --format-output target,alnlen,nident,qlen,tlen,evalue -e 1 > /dev/null"%(outdir,bindir,infile,rootdir)
os.system(cmd)
filename=outdir+"/aln.m8"
if not os.path.isfile(filename):
    ExitWithError(cmd,html_footer)

m8_list=[]
m8_dict=dict()
fp=open(filename,'r')
for line in fp.read().splitlines():
    items=line.split('\t')
    target=items[0].split('.')[0]
    if target in m8_dict:
        continue
    evalue=items[-1]
    if len(m8_list)>10 and float(evalue)>0.001:
        break
    m8_list.append(target)
    m8_dict[target]=items[1:]
fp.close()

if len(m8_list)==0:
    ExitWithError("Cannot find foldseek hit",html_footer)

fp=open(outdir+"/list",'w')
fp.write(''.join([target+'\n' for target in m8_list]))
fp.close()

#### run TMalign ####

cmd="cd %s; %s/xyz_sfetch %s/foldseek/receptor_nr.xyz list | gzip - > xyz.gz"%(outdir,bindir,rootdir)
os.system(cmd)

cmd="cd %s; %s/USalign %s xyz.gz -infmt2 2 -fast -outfmt 2 |grep -v '^#'|cut -f2-|sed 's/^xyz.gz://g'"%(outdir,bindir,infile)
p=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
stdout,stderr=p.communicate()
lines_dict=dict()
for line in stdout.decode().splitlines():
    items=line.split('\t')
    lines_dict[items[0]]=items
fp.close()
lines=[lines_dict[target] for target in m8_list if target in lines_dict]

print('''
<table border="0" align=center width=100%>    
<tr BGCOLOR="#FF9900">
    <th width=5% ALIGN=center><strong> # </strong></th>
    <th width=10% ALIGN=center><strong> Hit<br>(length)</strong></th>
    <th width=10% ALIGN=center><strong> Aligned<br>length </strong></th>
    <th width=10% ALIGN=center><strong> TM-score (Identity)<br>normalized by query</strong> </th>           
    <th width=10% ALIGN=center><strong> TM-score (Identity)<br>normalized by hit</strong> </th>           
    <th width=10% ALIGN=center><strong> RMSD (Identity)<br>normalized by<br>aligned length</strong> </th>           
    <th width=10% ALIGN=center><strong> E-value </strong> </th>           
    <th width=35% ALIGN=center><strong> Homologs<br>to hit</strong> </th>           
</tr><tr ALIGN=center>
''')
hit2chain_dict=dict()
hit2clust_dict=dict()
if len(lines):
    cmd="zcat "+rootdir+"/data/protein.fasta.gz|grep '>'|cut -f1,2|sed 's/^>//g'"
    p=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
    stdout,stderr=p.communicate()
    for line in stdout.decode().splitlines():
        hit,chainID=line.split('\t')
        pdbid=hit[:-len(chainID)]
        hit2chain_dict[hit]=(pdbid,chainID)
    fp=gzip.open(rootdir+"/data/protein_nr.fasta.clust.gz",'rt')
    for line in fp.read().splitlines():
        rep,mem=line.split('\t')
        hit2clust_dict[rep]=mem.split(',')
    fp.close()

totalNum=0
for line in lines:
    items=line#.split('\t')
    sacc =items[0]
    TM1  =items[1]
    TM2  =items[2]
    RMSD =items[3]
    ID1  =items[4]
    ID2  =items[5]
    IDali=items[6]
    L1   =items[7]
    L2   =items[8]
    Lali =items[9]
    evalue="NA"
    if sacc in m8_dict:
        evalue=m8_dict[sacc][-1]
    totalNum+=1
    bgcolor=''
    if totalNum%2==0:
        bgcolor='BGCOLOR="#DEDEDE"'

    
    pdbid,chainID=hit2chain_dict[sacc]
    hit="<a href=qsearch.cgi?pdbid=%s&chain=%s target=_blank>%s:%s</a>"%(
        pdbid,chainID,pdbid,chainID)
    homolog_list=[]
    if sacc in hit2clust_dict:
        for mem in hit2clust_dict[sacc]:
            pdbid,chainID=hit2chain_dict[mem]
            homolog_list.append("<a href=qsearch.cgi?pdbid=%s&chain=%s target=_blank>%s:%s</a>"%(
                pdbid,chainID,pdbid,chainID))
    print('''
<tr %s ALIGN=center>
    <td>%d</td>
    <td>%s (%s)</td>
    <td>%s</td>
    <td>%s (%s)</td>
    <td>%s (%s)</td>
    <td>%s (%s)</td>
    <td>%s</td>
    <td>%s</td>
</tr>
'''%(bgcolor,
    totalNum,
    hit,L2,
    Lali,
    TM1,ID1,
    TM2,ID2,
    RMSD,IDali,
    evalue,
    ', '.join(homolog_list)))
print("</table><p></p><a href=.>[Back]</a>")
if len(html_footer):
    print(html_footer)
else:
    print("</body> </html>")
