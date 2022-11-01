#!/usr/bin/python3
import cgi
import cgitb; cgitb.enable()  # for troubleshooting
import os
from datetime import datetime
import subprocess

rootdir=os.path.dirname(os.path.abspath(__file__))

html_header=""
html_footer=""
if os.path.isfile(rootdir+"/index.html"):
    fp=open(rootdir+"/index.html")
    txt=fp.read()
    fp.close()
    html_header=txt.split('<!-- CONTENT REFRESH -->')[0]
    html_footer=txt.split('<!-- CONTENT END -->')[-1]

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

mtime=os.path.getmtime(os.path.join(rootdir,"download","lig_frequency.txt"))
mtime=datetime.fromtimestamp(mtime)
mtime="%d-%d-%d"%(mtime.year,mtime.month,mtime.day)

cmd="zcat %s/data/metal.tsv.gz|cut -f1"%rootdir
p=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
stdout,stderr=p.communicate()
stdout=stdout.decode()
metal_list=stdout.splitlines()
metal_set=set(metal_list)

cmd="zcat %s/data/lig_all.tsv.gz|tail -n +2|cut -f1,2,4,9-12"%rootdir
p=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
stdout,stderr=p.communicate()
stdout=stdout.decode()
lines=stdout.splitlines()

rnas      =0
dnas      =0
peptides  =0
metals    =0
regulars  =0
affinities=0
moads     =0
pdbbinds  =0
bindingdbs=0
manuals   =0
chain_list=[]
for line in lines:
    items=line.split('\t')
    chain_list.append(items[0]+items[1])

    ccd=items[2]
    if ccd=="peptide":
        peptides+=1
    elif ccd=="rna":
        rnas+=1
    elif ccd=="dna":
        dnas+=1
    elif ccd in metal_set:
        metals+=1
    else:
        regulars+=1

    has_manual   =len(items[3])>0
    has_moad     =len(items[4])>0
    has_pdbbind  =len(items[5])>0
    has_bindingdb=len(items[6])>0
    affinities+=has_manual or has_moad or has_pdbbind or has_bindingdb
    manuals   +=has_manual;
    moads     +=has_moad
    pdbbinds  +=has_pdbbind
    bindingdbs+=has_bindingdb

print('''
<p>
<h1><u>BioLiP in numbers</u></h1>
</p>

BioLiP is updated weekly and the current version (%s) contains:
<li>Number of protein receptors: <a href=browse.cgi>%d</a></li>
<li>Number of entries: <a href=browse.cgi>%d</a></li>
<li>Number of entries for regular ligands: <a href=qsearch.cgi?lig3=regular>%d</a></li>
<li>Number of entries for metal ligands: <a href=qsearch.cgi?lig3=metal>%d</a></li>
<li>Number of entries for peptide ligands: <a href=qsearch.cgi?lig3=peptide>%d</a></li>
<li>Number of entries for DNA ligands: <a href=qsearch.cgi?lig3=dna>%d</a></li>
<li>Number of entries for RNA ligands: <a href=qsearch.cgi?lig3=rna>%d</a></li>
<li>Number of entries with binding affinity data: <a href=qsearch.cgi?baff=baff>%d</a>
(<a href=qsearch.cgi?baff=moad>%d</a> from Binding MOAD, <a href=qsearch.cgi?baff=pdbbind>%d</a> from PDBbind-CN, <a href=qsearch.cgi?baff=bindingdb>%d</a> from BindingDB, and <a href=qsearch.cgi?baff=manual>%d</a> from manual survey of the original literature)
'''%(mtime, len(set(chain_list)), regulars+metals+peptides+dnas+rnas,
    regulars,metals,peptides,dnas,rnas,
    affinities,moads,pdbbinds,bindingdbs,manuals
))

if len(html_footer):
    print(html_footer)
else:
    print("</body> </html>")
