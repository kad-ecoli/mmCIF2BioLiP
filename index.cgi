#!/usr/bin/python3
import cgi
import cgitb; cgitb.enable()  # for troubleshooting
import os
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

filename=rootdir+"/data/index.txt"
if os.path.isfile(filename):
    fp=open(filename,'r')
    print(fp.read())
    fp.close()

cmd="ls -rt output/|grep -P '[a-z0-9]+_[A-Za-z0-9]+_[FPC]\.svg'|cut -f1,2 -d_|uniq|tail -1"
p=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
stdout,stderr=p.communicate()
stdout=stdout.decode().strip()
if len(stdout):
    pdbid,chain=stdout.split('_')[:2]
    cmd="ls output/%s_*_*_*.pdb.gz|wc -l"%(pdbid)
    p=subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
    stdout,stderr=p.communicate()
    bs=''
    if int(stdout.decode()):
        bs="&bs=BS01"
    print('''
<p>
<h1><span title="PDB $pdbid Chain $chain Binding Site BS01"><a href=pdb.cgi?pdb=$pdbid&chain=$chain$bs target=_blank>View a random BioLiP entry</a></span></h1>
</p>
'''.replace("$pdbid",pdbid).replace("$chain",chain).replace("$bs",bs))


if len(html_footer):
    print(html_footer)
else:
    print("</body> </html>")
